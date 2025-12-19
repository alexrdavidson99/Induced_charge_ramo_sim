
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py
import functools
from pathlib import Path
import re
from scipy.optimize import curve_fit
from scipy.signal import chirp, find_peaks, peak_widths
from scipy import integrate
from matplotlib.patches import Polygon

# import mplhep
# mplhep.style.use(mplhep.style.LHCb2)
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.special import erf
#Gaussian function
def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) 

def gaussian_convolved_tophat(x, amp, mean, stddev, width=0.55):
    return amp * (erf((x - mean + width / 2) / (np.sqrt(2) * stddev)) - erf((x - mean - width / 2) / (np.sqrt(2) * stddev))) / 2


def load_and_group_sum(DATA_DIR,shift, round_decimals=4):
    
    print(f"Data directory: {DATA_DIR}")
    # Get files for each pattern
    files = list(DATA_DIR.glob(f"*off*{shift}*.csv"))
    print(f"Files found: {files}")
    df = pd.concat([pd.read_csv(file) for file in files], ignore_index=True)
    df['Time'] = df['Time'].round(round_decimals)
    df.groupby('Time', as_index=False)['Current'].sum()
    #df.to_csv(DATA_DIR / f'summed_induced_current_from_tom_p_center__off_in_x.csv', index=False)
    return df.groupby('Time', as_index=False)['Current'].sum()

def extract_co_value(co,s):
    match = re.search(fr"{co}_(-?\d+)", s)
    if match:
        return float(match.group(1))
    match = re.search(r"z_neg_(\d+)", s)
    if match:
        return -float(match.group(1))
    return 0.0



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
    
    
# filenames = [
#     'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_charge_10k_off_tom_pix.csv',
#     'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_charge_10k_on_tom_pix.csv',
# ]
# plt.figure()
# for filename in filenames:
#     plot_current(filename)
    
# plt.xlabel('Time (ns)')
# plt.ylabel('Induced Current (A)')

# plt.legend()
# plt.title('Induced Current Over Time')
# plt.show()

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

#load charge tom_ p csv data for all electrons from 1 to 1309561
# concat all csv files into one dataframe

Path_to_induced_charge = 'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/second_run_1mill_1500v/'
#DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/elastic-0.135phat/Fast-int-600k-1500V')

plt.figure(figsize=(10, 6))
DATA_DIR = Path(Path_to_induced_charge)
print(f"Data directory: {DATA_DIR}")
files = DATA_DIR.glob(f"*z_620-1*.csv")
for file in files:
    print(f"Processing file: {file}")
    df = pd.read_csv(file)
    # Plot the current data
    plt.plot(df['Time'], df['Current'], label=file.name)

Path_to_induced_charge = 'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/second_run_1mill_1500v_new_field/'


DATA_DIR = Path(Path_to_induced_charge)
print(f"Data directory: {DATA_DIR}")
files = DATA_DIR.glob(f"*z_610*.csv")
print(f"Data directory: {files}")
for file in files:
    print(f"Processing file: {file}")
    df = pd.read_csv(file)
    # Plot the current data
    plt.plot(df['Time'], df['Current'], label=file.name)
#plt.legend()

plt.show()
# Path_to_induced_charge = 'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/'


# DATA_DIR = Path(Path_to_induced_charge)
# print(f"Data directory: {DATA_DIR}")
# files = DATA_DIR.glob(f"induced_charge_tom_P*off*z*.csv")
# for file in files:
#     print(f"Processing file: {file}")
#     df = pd.read_csv(file)
#     # Plot the current data
#     plt.plot(df['Time'], df['Current'], label=file.name)


# plt.figure(figsize=(10, 6))
# DATA_DIR = Path(Path_to_induced_charge)
# print(f"Data directory: {DATA_DIR}")

# DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/TORCH_127/')

# # list_of_positions = ['z_70']
# list_of_positions = [f'z_{(z)}' for z in range(-630, 670, 50)]


# peaks_array = []
# peaks_positions = []

# plt.figure(figsize=(10, 6))
# for position in list_of_positions:
#     grouped = load_and_group_sum(DATA_DIR_v,position)
#     print(f"Processing group sum : {position }")
#     grouped['Current'].min()
#     max_current = grouped['Current'].min()
#     peaks_array.append(max_current)
#     print(f"Max current for {position}: {max_current}")
#     z_co = extract_co_value("z",position)
#     peaks_positions.append(z_co)
#     #plt.plot(z_co, max_current, 'o', label=position)
#     plt.plot(grouped['Time'], grouped['Current']*1e3, label= position)



# list_of_positions = ['z_70']



# #files_set2 = list(DATA_DIR.glob("induced_charge_tom_P*off*z*.csv"))
# #z_50'
#DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/fast_1500V/')
#DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/elastic-0.135phat/Fast-int-600k-1500V')

# #----> pd-25 DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/TORCH_127_long_range/')
DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/second_run_1mill_1500V' )
# # #DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/') # no file = 1000V
# # #DATA_DIR_v = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/')
# # #list_of_positions = ['x', 'z_49.9','z_60', 'z_100','z_175','z_550', 'z_225', 'z_275', 'z_300', 'z_400', 
# # #                     'z_350', 'z_450', 'z_500', 'z_1100', 'z_1650', 'z_neg_60','z_neg_1040', 'z_neg_1590', 'z_neg_550','z_neg_410', 'z_neg_275','z_neg_100','z_neg_150','z_neg_200']
# #                      #'z_neg_100']
# # #list_of_positions = ['z_60','z_550', 'z_400','z_225', 'z_neg_550', 'z_neg_275', 'z_neg_100']
# # #list_of_positions = ['z_60','z_350','z_550','z_1100']
# # #list_of_positions = ['z_60','85','z_490','z_neg_490']
# # #list_of_positions = ['z_neg_490','z_neg_50', 'z_neg_150','z_neg_250','z_0','z_60','z_75','z_85','z_150','z_160','z_250','z_325','z_375','z_375','z_400','z_490','z_610']
# # #list_of_positions = ['z_neg_490','z_60','z_610']
list_of_positions = [f'z_{(z)}' for z in range(-480, 670, 50)]
# #list_of_positions = [f'x_{(x)}' for x in range(-3030, 2470, 250)]
# #list_of_positions = [f'z_{(z)}' for z in range(-630, 670, 50)]
# #list_of_positions = [f'z_{(z)}' for z in range(-530, 570, 50)]
# # #list_of_positions = [f'x_{(z)}' for z in range(-3030, 2580, 250)]
# # #list_of_positions = ['x','z_60','z_550', 'z_1100', 'z_1650', 'z_neg_60']
# #                      #'z_neg_100']
# #list_of_positions = ['z_120']
# #list_of_positions = [f'z_{(z)}' for z in range(120, 220, 50)]


peaks_array = []
peaks_positions = []
plt.figure(figsize=(10, 6))
for position in list_of_positions:
    grouped = load_and_group_sum(DATA_DIR_v,position)
    print(f"Processing group sum : {position }")
    grouped['Current'].min()
    max_current = grouped['Current'].min()
    peaks_array.append(max_current)
    print(f"Max current for {position}: {max_current}")
    z_co = extract_co_value("z",position)
    z_co = np.array(z_co)
    peaks_positions.append(z_co)
    plt.plot(z_co*0.001, max_current, 'o', label=position)
    #plt.plot(grouped['Time'], grouped['Current'], label= position)
    grouped['Time'] = grouped['Time']* 1e-9  # Convert time to ns
    x = grouped['Time']
    y = grouped['Current']
    a = 0       # lower limit
    b = 1.56e-9       # upper limit

        # Method 1: Using the trapezoidal rule
    area_trapz = np.trapz(y, x)

    # Method 2: Using Simpson's rule (more accurate)
    #area_simps = integrate.simps(y, x)

    # Method 3: Using definite integral directly
    #area_quad, _ = integrate.quad(f, a, b)

    # Print results
    print(f"Area using trapezoidal rule   : {area_trapz:.2e}")
    #print(f"Area using Simpson's rule    : {area_simps:.2e}")
    print(f"gain in electrons = {area_trapz/-1.6e-19:.2f} electrons")
#     #print(f"Area using scipy.integrate.quad : {area_quad:.6f}")
    grouped = pd.concat([pd.DataFrame({'Time': [0.0], 'Current': [0.0]}),grouped], ignore_index=True)
    
    grouped.to_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_{position}_sumed.csv', index=False,header=False)

# # list_of_positions = [f'z_{(z)}' for z in range(0, 120, 10)]


# # for position in list_of_positions:
# #     grouped = load_and_group_sum(DATA_DIR_v,position)
# #     print(f"Processing group sum : {position }")
# #     grouped['Current'].min()
# #     max_current = grouped['Current'].min()
# #     peaks_array.append(max_current)
# #     print(f"Max current for {position}: {max_current}")
# #     z_co = extract_co_value("z",position)
# #     peaks_positions.append(z_co)
# #     plt.plot(z_co, max_current, 'o', label=position)
# #     plt.plot(grouped['Time'], grouped['Current'], label= position)




# plt.ylabel('Peak current (mA)')
# plt.xlabel('x position (mm)')
# plt.title(' Simulated Induced Current at Different Positions')
# #plt.yscale('log')
# #plt.legend()
# #plt.xlim(-1500, 1500)

# #peaks v  peak  posotions  datafame 
# peaks_df = pd.DataFrame({'Position': peaks_positions, 'Peak Current (A)': peaks_array})
# DATA_DIR_P = Path('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/')
# peaks_df.to_csv(DATA_DIR_P / 'peaks_positions_1000V.csv', index=False)

popt, pcov = curve_fit(gaussian_convolved_tophat, peaks_positions, peaks_array, p0=[0.0008, 60, 100])
peaks_positions_array = np.array(peaks_positions)
peaks_array = np.array(peaks_array)
plt.figure(figsize=(10, 6))
x_fit = np.linspace(min(peaks_positions), max(peaks_positions), 1000)
#plt.plot(x_fit, gaussian_convolved_tophat(x_fit, *popt)*1e3, 'r--', label=f'Gaussian Fit FWHM = {2 * np.sqrt(2 * np.log(2)) * popt[2]:.0f}um')
plt.scatter(peaks_positions_array*1e-3, peaks_array*1e3,  label='Data from simulation')
print(peaks_positions_array)
index = np.where((peaks_positions_array == 620) | (peaks_positions_array == -480))[0]
print (f"Indices of highlighted points: {index}")
if len(index) > 0:

    plt.scatter(peaks_positions_array[index]*1e-3, peaks_array[index]*1e3, color='red', label='Highlighted Points')
print(f"FWHM: {2 * np.sqrt(2 * np.log(2)) * popt[2]}")  # Full Width at Half Maximum
#plt.title('Peak Current vs Position')



plt.xlabel('z position (mm)')

plt.ylabel('peak current (mA)')

plt.legend()
plt.xlim(-0.5, 0.65)
plt.ylim(-0.7, 0)


plt.figure(figsize=(10, 6))
summed_current_df = pd.read_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_z_620_sumed.csv', names=['Time', 'Current'])
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='z_610')
summed_current_df = pd.read_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_z_-480_sumed.csv', names=['Time', 'Current'])
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='z_-480')
summed_current_df = pd.read_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_z_70_sumed.csv', names=['Time', 'Current'])
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='z_60')
plt.ylabel('Induced Current (mA)')
plt.xlabel('Time (ns)')

plt.legend()
plt.xlim(0, 1.5)
# plt.title('Different Patterns of Induced Current')



# # # signal = -grouped1['Current'].values
# # # time = grouped1['Time'].values

# # # # Normalize the signal between 0 and 1
# # # signal_norm = (signal - np.min(signal)) / (np.max(signal) - np.min(signal))

# # # # Find indices where signal crosses 10% and 90%
# # # idx_10 = np.where(signal_norm >= 0.1)[0][0]
# # # idx_90 = np.where(signal_norm >= 0.9)[0][0]

# # # t_10 = time[idx_10]
# # # t_90 = time[idx_90]
# # # rise_time = t_90 - t_10

# # # print(f"10% rise time: {t_10}")
# # # print(f"90% rise time: {t_90}")
# # # print(f"Rise time (10% to 90%): {rise_time}")
# # # plt.vlines([t_10,t_90], ymin=0, ymax=-0.0008, colors='red', linestyles='dashed', label='Rise Time')
# # # plt.legend()

# # # V1500 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/1500V/summed_induced_current_from_tom_p_center__off_in_x_test.csv")
# # # V1000 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current/summed_induced_current_from_tom_p_center__off_in_x_test.csv")
# # # plt.figure(figsize=(10, 6))
# # # plt.plot(V1500['Time'], V1500['Current']*1e3, label='1500V', color='blue')
# # # plt.plot(V1000['Time'], V1000['Current']*1e3, label='1000V', color='orange')
# # # plt.xlabel('Time (ns)')
# # # plt.ylabel('Induced Current (mA)')
# # # plt.title('Induced Current at 1500V and 1000V')

# # V1500 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/peaks_positions_1500V.csv")
# # # V1000 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/peaks_positions_1000V.csv")
# # # plt.figure(figsize=(10, 6))
# # plt.scatter(V1500['Position'], V1500['Peak Current (A)'], label='1500V', color='blue')
# # # plt.scatter(V1000['Position'], V1000['Peak Current (A)'], label='1000V', color='orange')
# # # plt.xlabel('Position (um)')
# # # plt.ylabel('Induced Current (mA)')
# # # plt.title('Induced Current at 1500V and 1000V')

# plt.figure(figsize=(10, 6))
# V1500 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/x_z_data_from_python_1500v_Loffler_2022-400k-elastic.csv")
# # V1000 = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/x_z_data_from_python_1.3_million.csv")
# # # plt.figure(figsize=(10, 6)) 

# #plt.scatter(V1500['x_position'], V1500['z_position'], label='1500V', color='blue')
# # #plt.scatter(V1000['x_position'], V1000['z_position'], label='1000V', color='orange')
# plt.hist2d(V1500['x_position'], V1500['z_position'], bins=100, range=[[-2000, 2000], [-300, 300]], cmap='plasma')
# #plt.hist(V1500['x_position'], bins=100, alpha=0.5, label='1000V', color='orange')
# #plt.hist(V1500['z_position'], bins=100, alpha=0.5, label='1500V', color='blue')
# # #plt.hist(V1000['x_position'], bins=500, alpha=0.5, label='1000V', )
# # plt.hist(V1000['z_position'], bins=500, alpha=0.5, label='1500V',)
# box_coords = [(-1650, -275), (-1650, 275), (1650, 275), (1650, -275)]
# polygon = Polygon(box_coords, closed=True, edgecolor='g', linestyle='--', fill=None)
# plt.gca().add_patch(polygon)

# # # plt.xlabel('Position (um)')
# # # plt.ylabel('Position (um)')
# # # plt.title('Hit map at 1500V and 1000V')

plt.show()

