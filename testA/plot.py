import pandas as pd
import matplotlib.pyplot as plt

filename = 'testA.in.out'

# Data storage
dateset_names = {
    'CELL_AVG_PHI',
    'PHI_HAT',
    'CELL_EDGE_PHI',
    'CELL_EDGE_Jx',
    'CELL_AVG_PSI',
    'SCHEME_RESIDUAL',
    'SCHEME_RESIDUAL_HAT',
    'ITERATIVE_RESIDUAL', 
    'ITERATIVE_RESIDUAL_HAT',
    'ERROR_L_INF',
    'RHO_L'
}

data = {}
for name in dataset_names:
    data['name'] = []	

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
                if current_section in ['CELL_AVG_PHI', 'CELL_EDGE_PHI', 'CELL_EDGE_Jx', 'SCHEME_RESIDUAL', 'SCHEME_RESIDUAL_HAT', 'ITERATIVE_RESIDUAL', 'ITERATIVE_RESIDUAL_HAT'] and len(nums) == 3:
                    data[current_section].append({'x': nums[1], 'val': nums[2]})
                    remaining_rows -= 1

                # 4 columns: [index, mu, x, psi]
                elif current_section in ['CELL_AVG_PSI'] and len(nums) == 4:
                    data[current_section].append({'mu': nums[1], 'x': nums[2], 'psi': nums[3]})
                    remaining_rows -= 1

                # 2 columns: [index, error]
                elif current_section in ['ERROR_L_INF', 'RHO_L'] and len(nums) == 2:
                    data[current_section].append({'x': nums[0], 'val': numes[1]})
                    remaining_rows -= 1

            except ValueError:
                pass  # Skip non-numeric header lines (e.g. "i x_mid CELL_AVG_PHI")
            # Stop storing as soon as expected count is reached
            if remaining_rows == 0:
                current_section = None

# Convert lists to DataFrames
    'CELL_AVG_PHI',
    'PHI_HAT',
    'CELL_EDGE_PHI',
    'CELL_EDGE_Jx',
    'CELL_AVG_PSI',
    'SCHEME_RESIDUAL',
    'SCHEME_RESIDUAL_HAT',
    'ITERATIVE_RESIDUAL', 
    'ITERATIVE_RESIDUAL_HAT',
    'ERROR_L_INF',
    'RHO_L'
df_avg_phi = pd.DataFrame(data['CELL_AVG_PHI'])
df_phi_hat = pd.DataFrame(data['PHI_HAT'])
df_phi_edge= pd.DataFrame(data['CELL_EDGE_PHI'])
df_edge_jx = pd.DataFrame(data['CELL_EDGE_Jx'])
df_psi_avg = pd.DataFrame(data['CELL_AVG_PSI'])
df_scheme_res = pd.DataFrame(data['SCHEME_RESIDUAL',)
df_scheme_res_hat = pd.DataFrame(data['SCHEME_RESIDUAL_HAT',)
df_iter_res = pd.DataFrame(data['ITERATIVE_RESIDUAL',)
df_iter_res_hat = pd.DataFrame(data['ITERATIVE_RESIDUAL_HAT',)
df_error_l_inf = pd.DataFrame(data['ERROR_L_INF',)
df_rho_l = pd.DataFrame(data['RHO_L',)



#cell avg. scalar Flux
if not df_avg_phi.empty:
    plt.plot(df_avg_phi['x'], df_avg_phi['val'], 'o-', label=r'Cell Avg $\Phi$')
plt.ylabel(r'Scalar Flux ($\Phi$)')
plt.title('Cell Average Scalar Flux')
plt.grid(True)
plt.legend()
plt.savefig('PHI_AVG.png')

# Scalar Flux first spatial moment
if not df_phi_hat.empty:
    plt.plot(df_phi_hat['x'], df_phi_hat['val'], 'o-', label=r'Scalar Flux spatial moment $\hat \Phi$')
plt.ylabel(r'($\hat \Phi$)')
plt.title('First Spatial Moment of Scalar Flux')
plt.grid(True)
plt.legend()
plt.savefig('PHI_HAT.png')

# Scalar Flux edge value
if not df_phi_edge.empty:
    plt.plot(df_phi_edge['x'], df_phi_edge['val'], 'o-', label=r'Cell Edge $\Phi$')
plt.ylabel(r'Scalar Flux ($\Phi$)')
plt.title('Cell Edge Scalar Flux')
plt.grid(True)
plt.legend()
plt.savefig('PHI_EDGE.png')

# Scalar Flux edge value
if not df_edge_jx.empty:
    plt.plot(df_edge_jx['x'], df_edge_jx['val'], 'o-', label=r'Cell Edge $\Phi$')
plt.ylabel(r'Current ($J_x$)')
plt.title('Cell Edge Scalar Flux')
plt.grid(True)
plt.legend()
plt.savefig('PHI_EDGE.png')

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
