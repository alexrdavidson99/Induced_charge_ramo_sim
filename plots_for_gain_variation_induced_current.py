
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

import mplhep
mplhep.style.use(mplhep.style.LHCb2)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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


def load_and_group_sum_new(DATA_DIR,shift, round_decimals=4):
    
    print(f"Data directory: {DATA_DIR}")

    # Get files for each pattern
    files = list(DATA_DIR.glob(f"*{shift}*.csv"))
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
    
gains = [
    2889476,
    40232,
    85952,
    294858,
    383788,
    487319,
    1031399,
    1074941,
    1179721,
    1261111,
    1972558
]

colors = [
    "#0d0887",
    "#3a049a",
    "#6300a7",
    "#8b0aa5",
    "#b12a90",
    "#cc4778",
    "#e16462",
    "#f2844b",
    "#fca636",
    "#f0f921",
    "#ffff33"   # highest gain
]

colors = [
    "#e6194b",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#46f0f0",
    "#f032e6",
    "#bcf60c",
    "#fabebe",
    "#008080"
]

currents = []
times = []


plt.figure(figsize=(10*1.5, 6*1.5))

for gain in gains:
    Path_to_induced_charge = f'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/run_{gain}_electrons_1500v_new_field_with_4_pix'
    index = gains.index(gain)


    gain_1million = gain*1e-6
    DATA_DIR = Path(Path_to_induced_charge)
    print(f"Data directory: {DATA_DIR}")


    list_of_positions = [f'z-1035'] # [f'z_{(z)}' for z in range(-480, 670, 50)]
    for position in list_of_positions:


        

        grouped = load_and_group_sum_new(DATA_DIR,position)
        t = np.array(grouped['Time'])
        I = np.array(grouped['Current']) * 1e3

        times.append(t)
        currents.append(I)



        print(f"Processing group sum : {position }")
        plt.plot(grouped['Time'], grouped['Current']*1e3, color=colors[index], label= rf"Gain = {gain_1million:.2f}$\times10^6$ ")

n_gains = len(gains)
t_ref = times[0]

interp_currents = []

for t, I in zip(times, currents):
    I_interp = np.interp(t_ref, t, I)
    interp_currents.append(I_interp)

# compute average
avg_current = np.mean(interp_currents, axis=0)
avg_voltage = avg_current * 50  # Convert current to voltage using R=50 Ohms
# avg_data_datafame
avg_df = pd.DataFrame({
    'Time': t_ref,
    'Average Current (mA)': avg_current,
    'Average Voltage (mV)': avg_voltage
})
#avg_df.to_csv('average_current_voltage_from_10_runs.csv', index=False)
#plt.plot(t_ref, avg_current*50, color='black', linewidth=3, label="Average Current")



plt.legend(loc="lower right")
plt.ylim(-0.11, 0.045)
plt.xlabel('Time (ns)')
plt.ylabel('Induced current (mA)')
plt.grid()
plt.savefig('induced_current_comparison_gain_variation_next_to_n.pdf')


plt.figure(figsize=(10, 6))
for gain in gains:
    Path_to_induced_charge = f'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/run_{gain}_electrons_1500v_new_field_with_4_pix'



    
    DATA_DIR = Path(Path_to_induced_charge)
    print(f"Data directory: {DATA_DIR}")


    list_of_positions = [f'z-1035','z-485','z65','z615','z1165'] # [f'z_{(z)}' for z in range(-480, 670, 50)]
    for position in list_of_positions:
        

        grouped = load_and_group_sum_new(DATA_DIR,position)
        print(f"Processing group sum : {position }")
        plt.plot(grouped['Time'], grouped['Current']*1e3, label= position)

Path_to_induced_charge = 'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/run_1074941_electrons_1500v_new_field_with_4_pix/90_degree_clockwise_rotation'



plt.figure(figsize=(10, 6))
DATA_DIR = Path(Path_to_induced_charge)
print(f"Data directory: {DATA_DIR}")



colors=["#4363d8", "#26bbbb", "#10C5F2", "#46f0f0", "#3078C1"]

labels = ['Nex to Neighbour 1', 'Neighbour 1', 'Center', 'Neighbour 2', 'Next to Neighbour 2']

list_of_positions = [f'z-1035','z-485','z65','z615','z1165'] # [f'z_{(z)}' for z in range(-480, 670, 50)]
for position in list_of_positions:
    

    grouped = load_and_group_sum_new(DATA_DIR,position)
    print(f"Processing group sum : {position }")
    # add to grouped a 0 current at time 0 to make the plot start at 0
    grouped = pd.concat([pd.DataFrame({'Time': [0.0], 'Current': [0.0]}),grouped], ignore_index=True)
    plt.plot(grouped['Time'], grouped['Current']*1e3, label= labels[list_of_positions.index(position)], color=colors[list_of_positions.index(position)])
    grouped['Time'] = grouped['Time']* 1e-9  # Convert time to ns
    #save each grouped to csv
    grouped.to_csv(DATA_DIR / f'summed_induced_current_1500v_in_z_{position}_sumed.csv', index=False,header=False)


plt.xlabel('Time (ns)')
plt.ylabel('Induced current (mA)')
plt.grid()
#plt.xlim(0, 2)
plt.legend(loc="lower left", fontsize='22')
plt.savefig('induced_current_at_60_z_1.07million_e.png', dpi=300)

plt.show()



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
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='Neighbour 1', color="#5296C7")
summed_current_df = pd.read_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_z_-480_sumed.csv', names=['Time', 'Current'])
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='Neighbour 2', color='#10C5F2')
summed_current_df = pd.read_csv(DATA_DIR_v / f'summed_induced_current_127__center_120_1500v_in_z_=_65_x=_z_70_sumed.csv', names=['Time', 'Current'])
plt.plot(summed_current_df['Time']*1e9, summed_current_df['Current']*1e3, label='Center', color='#1E90FF')
plt.ylabel('Induced current (mA)')
plt.xlabel('Time (ns)')

#plt.show()

