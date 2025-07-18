
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py


def plot_field(data, z_position):
    #data['Abs_E'] = np.sqrt(data['Ex'] ** 2 + data['Ey'] ** 2 + data['Ez'] ** 2)
    target_z = z_position
    closest_z = data['z'].iloc[(data['z'] - target_z).abs().argmin()]
    filtered_data =  data[data['z'] == closest_z ]

    x = filtered_data['x'].values
    y = filtered_data['y'].values
    z = filtered_data['Ey'].values
    #z = np.sqrt(filtered_data['Ex'] ** 2 + filtered_data['Ey'] ** 2 + filtered_data['Ez'] ** 2)
    levels = 30
    plt.tricontourf(x, y, z, levels=levels, cmap='plasma', vmin=-90, vmax=0)

    plt.colorbar(label='|E|')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Contour plot at z = {z_position}')
    plt.show()

def plot_current(filename):
    df = pd.read_csv(filename)
    info = filename.split('/')[-1]
    plt.plot(df['Time'], df['Current'], label=info)
    plt.yscale('symlog')
    #plt.ylim(-1e-5, 1e-6)
    
    
filenames = [
    'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_charge_10k_off_tom_pix.csv',
    'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_charge_10k_on_tom_pix.csv',
]
plt.figure()
for filename in filenames:
    plot_current(filename)
    
plt.xlabel('Time (ns)')
plt.ylabel('Induced Current (A)')

plt.legend()
plt.title('Induced Current Over Time')
plt.show()

#file_path = './3d/E-Field [Es].h5'
#z_value = -16  # Set your desired z index or value

# with h5py.File(file_path, 'r') as f:
#     # Load mesh grids
#     x = f['Mesh line x'][:]
#     y = f['Mesh line y'][:]
#     z = f['Mesh line z'][:]
#     # Find the index closest to the desired z value
#     z_idx = np.argmin(np.abs(z - z_value))
#     # Load the E-field at this z slice
#     E_field = f['E-Field'][:, :, z_idx]
#     # If E_field is structured, extract components
#     if hasattr(E_field, 'dtype') and E_field.dtype.fields:
#         Ex = E_field['x']
#         Ey = E_field['y']
#         Ez = E_field['z']
#         E_abs = np.sqrt(Ex**2 + Ey**2 + Ez**2)
#     else:
#         E_abs = np.abs(E_field)  # fallback

# # Plot the absolute field in the x-y plane at the chosen z
# plt.figure()
# print(x.shape, y.shape, E_abs.shape)
# print(z)
# plt.pcolormesh(y, np.arange(E_abs.shape[0]), E_abs, shading='auto',cmap='plasma', vmin=0, vmax=1000)
# plt.xlabel('y')
# plt.ylabel('x')

# plt.title(f'|E| at z={z[z_idx]:.2f}')
# plt.colorbar(label='|E|')
# plt.show()
    

