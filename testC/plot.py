import pandas as pd
import matplotlib.pyplot as plt

filename = 'testC.in.out'

# Data storage
data = {
    'CELL_AVG_PHI': [],
    'CELL_EDGE_PHI': [],
    'CELL_AVG_Jx': [],
    'CELL_EDGE_Jx': [],
    'CELL_EDGE_PSI': [],
    'CELL_AVG_PSI': []
}

current_section = None
remaining_rows = 0

# 1. READ FILE
with open(filename, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if not parts:
            continue

        # Check if line starts with a section name (e.g. "CELL_AVG_PHI 40")
        if parts[0] in data:
            current_section = parts[0]
            # Get expected row count from the 2nd number on the header line
            remaining_rows = int(parts[1]) if len(parts) > 1 else 0
            continue

        # Read data only while remaining_rows > 0
        if current_section and remaining_rows > 0:
            try:
                nums = [float(p) for p in parts]

                # 3 columns: [index, x, value]
                if current_section in ['CELL_AVG_PHI', 'CELL_EDGE_PHI', 'CELL_AVG_Jx', 'CELL_EDGE_Jx'] and len(nums) == 3:
                    data[current_section].append({'x': nums[1], 'val': nums[2]})
                    remaining_rows -= 1

                # 4 columns: [index, mu, x, psi]
                elif current_section in ['CELL_EDGE_PSI', 'CELL_AVG_PSI'] and len(nums) == 4:
                    data[current_section].append({'mu': nums[1], 'x': nums[2], 'psi': nums[3]})
                    remaining_rows -= 1

            except ValueError:
                pass  # Skip non-numeric header lines (e.g. "i x_mid CELL_AVG_PHI")

            # Stop storing as soon as expected count is reached
            if remaining_rows == 0:
                current_section = None

# Convert lists to DataFrames
df_avg_phi = pd.DataFrame(data['CELL_AVG_PHI'])
df_edge_phi = pd.DataFrame(data['CELL_EDGE_PHI'])
df_avg_jx = pd.DataFrame(data['CELL_AVG_Jx'])
df_edge_jx = pd.DataFrame(data['CELL_EDGE_Jx'])
df_psi_edge = pd.DataFrame(data['CELL_EDGE_PSI'])
df_psi_avg = pd.DataFrame(data['CELL_AVG_PSI'])


# 2. PLOT SCALAR FLUX & CURRENT
plt.figure(figsize=(8, 6))

# Top plot: Scalar Flux
plt.subplot(2, 1, 1)
if not df_avg_phi.empty:
    plt.plot(df_avg_phi['x'], df_avg_phi['val'], 'o-', label=r'Cell Avg $\Phi$')
if not df_edge_phi.empty:
    plt.plot(df_edge_phi['x'], df_edge_phi['val'], 's--', label=r'Cell Edge $\Phi$')
plt.ylabel(r'Scalar Flux ($\Phi$)')
plt.title('Scalar Flux')
plt.grid(True)
plt.legend()

# Bottom plot: Current
plt.subplot(2, 1, 2)
if not df_avg_jx.empty:
    plt.plot(df_avg_jx['x'], df_avg_jx['val'], 'o-', label=r'Cell Avg $J_x$')
if not df_edge_jx.empty:
    plt.plot(df_edge_jx['x'], df_edge_jx['val'], 's--', label=r'Cell Edge $J_x$')
plt.xlabel('x')
plt.ylabel(r'Current ($J_x$)')
plt.title('Current')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig('testA_phi_J.png')


# 3. PLOT ANGULAR FLUX (1 curve per mu)
plt.figure(figsize=(8, 5))

if not df_psi_edge.empty:
    for mu, group in df_psi_edge.groupby('mu'):
        style = '--' if mu < 0 else '-'
        plt.plot(group['x'], group['psi'], style, label=rf'$\mu$ = {mu:+.4f}')

plt.xlabel('x')
plt.ylabel(r'Angular Flux ($\psi$)')
plt.title(r'Cell Edge Angular Flux ($\psi$) per Angle')
plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('testA_psi_edge.png')

# 3. PLOT ANGULAR FLUX (1 curve per mu)
plt.figure(figsize=(8, 5))

if not df_psi_avg.empty:
    for mu, group in df_psi_avg.groupby('mu'):
        style = '--' if mu < 0 else '-'
        plt.plot(group['x'], group['psi'], style, label=rf'$\mu$ = {mu:+.4f}')

plt.xlabel('x')
plt.ylabel(r'Angular Flux ($\psi$)')
plt.title(r'Cell Avg. Angular Flux ($\psi$) per Angle')
plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('testA_psi_avg.png')
