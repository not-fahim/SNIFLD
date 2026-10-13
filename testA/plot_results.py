"""Plot SNIF output. Usage: python plot_results.py testA.in.out"""
import sys
import re
from pathlib import Path
import matplotlib.pyplot as plt

filename = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('testA.in.out')
outdir = filename.parent / (filename.stem + '_plots')
outdir.mkdir(exist_ok=True)

# Each section has a header containing its name and number of data rows.
names = ['CELL_AVG_PHI', 'PHI_HAT', 'CELL_EDGE_PHI', 'CELL_EDGE_Jx',
         'CELL_AVG_PSI', 'SCHEME_RESIDUAL', 'SCHEME_RESIDUAL_HAT',
         'ITERATIVE_RESIDUAL', 'ITERATIVE_RESIDUAL_HAT',
         'ERROR_L_INF', 'RHO_L']
data = {name: [] for name in names}
lines = filename.read_text().splitlines()

for i, line in enumerate(lines):
    # Handles both "PHI_HAT40" and "CELL_AVG_PHI 40".
    match = re.fullmatch(r'\s*(' + '|'.join(names) + r')\s*(\d+)\s*', line)
    if not match:
        continue
    name, count = match.group(1), int(match.group(2))
    for row in lines[i + 1:]:
        parts = row.split()
        try:
            numbers = [float(v.replace('D', 'E')) for v in parts]
        except ValueError:
            continue  # Ignore column headings.
        needed = 4 if name == 'CELL_AVG_PSI' else (2 if name in ('ERROR_L_INF', 'RHO_L') else 3)
        if len(numbers) == needed:
            data[name].append(numbers)
        if len(data[name]) == count:
            break
    if len(data[name]) != count:
        print(f'Warning: {name}: expected {count} rows, found {len(data[name])}')


def save_plot(name, title, ylabel, xlabel='Position x (cm)', log=False):
    """Plot a three-column section: index, position, value."""
    if not data[name]:
        print(f'Skipping missing section: {name}')
        return
    x = [row[1] for row in data[name]]
    y = [row[2] for row in data[name]]
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, '-o', markersize=3)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if log:
        plt.yscale('symlog', linthresh=1e-14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / (name + '.png'), dpi=200)
    plt.close()


# Spatial solution plots
save_plot('CELL_AVG_PHI', 'Cell-average scalar flux', r'$\phi_i$')
save_plot('PHI_HAT', 'First spatial moment of scalar flux', r'$\hat{\phi}_i$')
save_plot('CELL_EDGE_PHI', 'Cell-edge scalar flux', r'$\phi_{i+1/2}$')
save_plot('CELL_EDGE_Jx', 'Cell-edge current', r'$J_{i+1/2}$')

# Angular flux: one line for each direction mu.
if data['CELL_AVG_PSI']:
    plt.figure(figsize=(9, 5))
    angles = sorted(set(row[1] for row in data['CELL_AVG_PSI']))
    for mu in angles:
        rows = [row for row in data['CELL_AVG_PSI'] if row[1] == mu]
        rows.sort(key=lambda row: row[2])
        plt.plot([r[2] for r in rows], [r[3] for r in rows],
                 label=f'mu = {mu:+.4f}')
    plt.xlabel('Position x (cm)')
    plt.ylabel(r'Cell-average angular flux $\bar{\psi}_{m,i}$')
    plt.title('Cell-average angular flux by direction')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(outdir / 'CELL_AVG_PSI.png', dpi=200)
    plt.close()

# Residuals (one plot for each and one combined plot).
residuals = ['SCHEME_RESIDUAL', 'SCHEME_RESIDUAL_HAT',
             'ITERATIVE_RESIDUAL', 'ITERATIVE_RESIDUAL_HAT']
for name in residuals:
    save_plot(name, name.replace('_', ' ').title(), 'Residual')

plt.figure(figsize=(9, 5))
for name in residuals:
    if data[name]:
        plt.plot([r[1] for r in data[name]], [r[2] for r in data[name]],
                 label=name)
plt.xlabel('Position x (cm)')
plt.ylabel('Residual')
plt.title('Four residuals')
plt.grid(True, alpha=0.3)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(outdir / 'ALL_RESIDUALS.png', dpi=200)
plt.close()

# Iteration history: columns are iteration number and value.
for name, ylabel in [('ERROR_L_INF', r'$\|\phi^{(l)}-\phi^{(l-1)}\|_\infty$'),
                     ('RHO_L', r'$\rho_\phi^{(l)}$')]:
    if data[name]:
        plt.figure(figsize=(8, 5))
        plt.plot([r[0] for r in data[name]], [r[1] for r in data[name]], '-o', markersize=3)
        if name == 'ERROR_L_INF':
            plt.yscale('log')
        plt.xlabel('Source iteration l')
        plt.ylabel(ylabel)
        plt.title(name.replace('_', ' '))
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(outdir / (name + '.png'), dpi=200)
        plt.close()

print('Input:', filename)
print('Number of source iterations:', len(data['ERROR_L_INF']))
print('Plots saved to:', outdir)
