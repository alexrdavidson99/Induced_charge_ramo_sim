import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import scipy
from matplotlib.pyplot import figure
from scipy.interpolate import griddata
import h5py
import time

start = time.time()


from main import step_position, step_velocity, solve_for_intercept_time


def find_closest_values(target_value, x, y, z):
    closest_value_x = np.argmin(np.abs(x - target_value[0]))
    closest_value_y = np.argmin(np.abs(y - target_value[1]))
    closest_value_z = np.argmin(np.abs(z - target_value[2]))



    return closest_value_x, closest_value_y, closest_value_z

c = scipy.constants.speed_of_light * 1e-3  # in um/ns
V = 1000  # electrods potential in V
d = 2500 # mcp anode gap in um
m = 511e3
E = V * (c ** 2) / (d * m)  # electric field acceleration in um/ns^2
E_field = V/(d*1e-6) 
print(f"E field = {E_field}")  
orientation = np.array([0., 1., 0])
a0 = E * orientation
# Load the E-Field data from the .h5 file
#file_path = 'C:/Users/lexda/Downloads/E-Field_higher_mesh[Es].h5'
#file_path = 'C:/Users/lexda/Downloads/E-Field [Es]_W_field_toms_anode_2.h5'
file_path = './3d/E-Field [Es].h5'
# start_position_data = pd.read_csv("C:/Users/lexda/Downloads/ascii_export_2_pore.csv", comment='#', skip_blank_lines=True, sep=';',
#                                       header=None, names=["Position [X]","Position [Y]","Position [Z]",
#                                                           "Position [ABS (XYZ)]","Time","Velocity [X]","Velocity [Y]",
#                                                           "Velocity [Z]"], index_col=False)

start_position_data = pd.read_csv("C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/ascii_export_of_15d_pore.csv", comment='#', skip_blank_lines=True, sep=';',
                                      header=None, names=["Position [X]","Position [Y]","Position [Z]",
                                                          "Position [ABS (XYZ)]","Time","Velocity [X]","Velocity [Y]",
                                                          "Velocity [Z]"], index_col=False)


# Load the E-Field data using the correct dataset names
with h5py.File(file_path, 'r') as f:
    x = f['Mesh line x'][:]
    y = f['Mesh line y'][:]
    z = f['Mesh line z'][:]

sum_data = []
position_start_data_x = []
position_start_data_z = []
position_end_data_x = []
position_end_data_z = []

electron_x = []
electron_y = []
electron_z = []

energies = []
alpha_list = []
vabs_list = []
shift_a_pad = 0
#print(len(start_position_data["Velocity [X]"]))
num_electrons = len(start_position_data)
middle = num_electrons // 2
indices = list(range(5)) + list(range(num_electrons - 5, num_electrons))
#indices = list(range(middle - 1, middle + 1))
indices = range(200000)  # Exclude the first and last 5 indices
for j in indices:

    index_electron = j
    target_coord = np.array([start_position_data["Position [X]"][index_electron]+ shift_a_pad,
                             start_position_data["Position [Y]"][index_electron],
                             start_position_data["Position [Z]"][index_electron] , ]) #
    position_start_data_x.append(target_coord[0])
    position_start_data_z.append(target_coord[2])


    v0 = np.array([start_position_data["Velocity [X]"][index_electron],
                   start_position_data["Velocity [Y]"][index_electron],
                   start_position_data["Velocity [Z]"][index_electron]]) # [m/s]
    
    vabs = np.sqrt(v0[0] ** 2 + v0[1] ** 2 + v0[2] ** 2)
    
    energy = ( vabs / scipy.constants.speed_of_light) ** 2 * 0.5 * m
    #print(f"Energy = {energy} eV")
    vabs_list.append(vabs)
    alpha = np.arctan(v0[1] / v0[0])
    alpha_list.append(alpha)

    energies.append(energy)
    #print(f"start_position = {target_coord}")
    closest_x, closest_y, closest_z = find_closest_values(target_coord, x, y, z)
    # print(f"Closest value to {target_coord[0]} in x: {x[closest_x]}")
    # print(f"Closest value to {target_coord[1]} in y: {y[closest_y]}")
    # print(f"Closest value to {target_coord[2]} in z: {z[closest_z]}")

    with h5py.File(file_path, 'r') as f:
        field_array = f['E-Field'][closest_z, closest_y, closest_x]
    field_x = field_array['x']
    field_y = field_array['y']
    field_z = field_array['z']

    # Create a numerical vector for the electric field
    field_vector = np.array([field_x, field_y, field_z])
    field_array = field_vector
    #print (f"field_array = {field_array*-1.6e-19}")

    current = -1.6e-19*np.dot(v0, field_array) #[m/s]  * [V/m] = [A]
    #print(f"Current = {current} A")

    x0 = target_coord

    v0 = v0*1e-3
    t = solve_for_intercept_time(x0, v0, a0, target_coord[1]+d) #d is distance to the anode in um

    #print(f"t = {t} time in ns")
    #print (f"start velocity = {v0} um/ns")
    #print (f"position_end = {step_position(x0, v0, a0, t)} um")
    position_end = step_position(x0, v0, a0, t)
    position_end_data_x.append(position_end[0])
    position_end_data_z.append(position_end[2])
    t_max = 5 # [ns]
    step = 1000 # number of steps to take in the time interval
    induced_current = []
    field_array_steps_y = []
    start_time_of_electron = np.array(start_position_data["Time"][index_electron]) * 1e9
    #print(f"start_time_of_electron = {start_time_of_electron}")

    for i in range (step):
        time_step = 0.001 * i
        if time_step >= t:
            #print(f"time step stopped at  = {time_step} ns")
            sum_data.append({"Time": time_step + start_time_of_electron, "Current": 0})
            #for k in range(10):
            #    sum_data.append({"Time": time_step + start_time_of_electron + (k/10)*2, "Current": 0})
            break

        xi = step_position(x0, v0, a0, time_step)
        vi = step_velocity(v0, a0, time_step)
        
        #xi = np.array([5.2, 3507, -30])
        #xi = np.array([-4985.1, 1515, -2010])
        #xi = np.array([-1379.8, 1515, -2010])
        closest_x, closest_y, closest_z = find_closest_values(xi, x, y, z)
        

        with h5py.File(file_path, 'r') as f:
            #field_array = f['E-Field'][closest_x, closest_y, closest_z]
            field_array = f['E-Field'][closest_z, closest_y, closest_x]
   
        # print(f"Closest value to {xi[0]} in x: {x[closest_x]}")
        # print(f"Closest value to {xi[1]} in y: {y[closest_y]}")
        # print(f"Closest value to {xi[2]} in z: {z[closest_z]}")
        field_array = np.array(field_array)
        field_x = field_array['x']
        field_y = field_array['y']
        
        field_z = field_array['z']

        # Create a numerical vector for the electric field
        field_vector = np.array([field_x, field_y, field_z])
        field_array = field_vector
        field_array_steps_y.append(field_array[1])
        #print (f"field_array = {field_array}")

        #energy = ((np.sqrt(vi[0] ** 2 + vi[1] ** 2 + vi[2] ** 2)) / c) ** 2 * 0.5 * m
        #print(f"Energy = {energy} eV")
        current = 1.6e-19 * np.dot(vi*1e3, field_array)
        # print(f"field_array = {field_array}")
        # print(f"velocity = {vi*1e3} m/s")
        induced_current.append(current)
        #print(f"Current = {current} A")


        sum_data.append({"Time": time_step+start_time_of_electron, "Current": current})
        electron_x.append(xi[0])
        electron_y.append(xi[1])
        electron_z.append(xi[2])

        #print(f"Energy = {energy} eV")

        #print(f"x = {xi}")

df = pd.DataFrame(sum_data)
zero_index = df[df['Current'] == 0].index[0]
plt.scatter(df["Time"][:zero_index+1], df["Current"][:zero_index+1])
plt.scatter(df["Time"][zero_index+1:], df["Current"][zero_index+1:])
#plt.plot(df["Time"][:-1], field_array_steps_y)
plt.xlabel("Time (ns)")
plt.ylabel("Current (A)")
plt.title("Induced Current")
result_sort = df.sort_values(["Time"], ascending=[True])
last_time = result_sort['Time'].iloc[-1]

#Generate additional zero current entries
additional_entries = [{'Time': last_time + 0.001 * i, 'Current': 0} for i in range(1, 501)]

# Append to result_sort
additional_df = pd.DataFrame(additional_entries)
result_sort = pd.concat([result_sort, additional_df], ignore_index=True)

df['Time'] = df['Time'].round(4) 
grouped_sum = df.groupby('Time', as_index=False)['Current'].sum()

print(grouped_sum[grouped_sum.isna().any(axis=1)])
# Optional: sort by Time if not already sorted
#grouped_sum = grouped_sum.sort_values(by='Time')

print("Rows with zero Current:")
print(grouped_sum[grouped_sum['Current'] == 0])

# Plotting
plt.figure()
plt.plot(grouped_sum['Time'], grouped_sum['Current'])
plt.xlabel("Time (ns)")
plt.ylabel("Current (A)")
plt.title("Induced Charge")
grouped_sum.to_csv('C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_charge.csv', index=False)

print("Rows where Time = 0.06:")
print(grouped_sum[grouped_sum['Time'] == 0.06])

# Sort the result_sort by Time
# result_sort = result_sort.sort_values(by='Time').reset_index(drop=True)

# rolling_sum = result_sort['Current'].rolling(window=250).sum()
# plt.figure()
# plt.plot(result_sort['Time'], rolling_sum)
# plt.xlabel("Time (ns)")
# plt.ylabel("Current (A)")
# plt.title("Induced Charge")
# plt.xlim(0, 10)
plt.figure()
plt.scatter(position_end_data_x, position_end_data_z, color='b',label='python',alpha=0.5)
plt.figure()
plt.hist2d(position_end_data_x, position_end_data_z, cmap='plasma', bins=50)
plt.figure()
box_coords = [(-143, -143), (-143, 143), (143, 143), (143, -143)]
polygon = Polygon(box_coords, closed=True, edgecolor='g', linestyle='--', fill=None)
plt.gca().add_patch(polygon)
data_from_csv = pd.read_csv('C:/Users/lexda/Downloads/compare_CST_with_python_CST_phase_space.txt',
                            names=["X","Z"], sep="\t+", comment='#')
plt.scatter(data_from_csv["X"], data_from_csv["Z"], color='r',label='CST',alpha=0.5)

plt.legend()
plt.xlabel('x[um]')
plt.ylabel('z[um]')
plt.title('Comparison between CST and python')
plt.figure()
plt.hist(data_from_csv["X"], bins=100)
plt.hist(data_from_csv["Z"], bins=100)
plt.xlabel('x[um]')
plt.figure()
plt.hist2d(alpha_list, vabs_list, cmap='plasma', bins=50)
plt.colorbar()
plt.xlabel('alpha[rad]')
plt.ylabel('vabs[m/s]')
figure()
alpha_list = [abs(alpha) for alpha in alpha_list]

#print (alpha_list)
alpha_list = np.array(alpha_list)
plt.hist(90-alpha_list*(180/np.pi), bins=100)
plt.xlabel('Energy[eV]')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(electron_x, electron_y, electron_z, label='Electron Trajectory')
# ax.set_xlabel('X [um]')
# ax.set_ylabel('Y [um]')
# ax.set_zlabel('Z [um]')
# ax.set_title('3D Trajectory of the Electron')
# ax.legend()

end = time.time()
print(f"total time {end - start} seconds")


plt.show()
