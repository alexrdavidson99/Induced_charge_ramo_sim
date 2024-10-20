import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import scipy
from matplotlib.pyplot import figure
from scipy.interpolate import griddata
import h5py


from main import step_position, step_velocity, solve_for_intercept_time




def find_closest_values(target_value, x, y, z):
    closest_value_x = np.argmin(np.abs(x - target_value[0]))
    closest_value_y = np.argmin(np.abs(y - target_value[1]))
    closest_value_z = np.argmin(np.abs(z - target_value[2]))



    return closest_value_x, closest_value_y, closest_value_z

c = scipy.constants.speed_of_light * 1e-3  # in um/ns
V =  3000  # electrods potential in V
d = 2850 - 459.3 # mcp anode gap in um
m = 511e3
E = V * (c ** 2) / (d * m)  # electric field acceleration in mm/ns^2
orientation = np.array([0., 1., 0])
a0 = E * orientation
# Load the E-Field data from the .h5 file
file_path = 'C:/Users/lexda/Downloads/E-Field_higher_mesh[Es].h5'
start_position_data = pd.read_csv("C:/Users/lexda/Downloads/ascii_export_close_to_pore.csv", comment='#', skip_blank_lines=True, sep=';',
                                      header=None, names=["Position [X]","Position [Y]","Position [Z]",
                                                          "Position [ABS (XYZ)]","Time","Velocity [X]","Velocity [Y]",
                                                          "Velocity [Z]"], index_col=False)

# Load the E-Field data using the correct dataset names
with h5py.File(file_path, 'r') as f:
    x = f['Mesh line x'][:]
    y = f['Mesh line y'][:]
    z = f['Mesh line z'][:]

print(z)


sum_data = []
position_start_data_x = []
position_start_data_z = []
position_end_data_x = []
position_end_data_z = []#

energies = []
alpha_list = []
vabs_list = []
shift_a_pad = 143*2+128   # x = 144 z= 229.1
# add 200 to shift on to the pad.
print(len(start_position_data["Velocity [X]"]))
for j in range(1):

    index_electron = j
    target_coord = np.array([start_position_data["Position [X]"][index_electron]-shift_a_pad,
                             start_position_data["Position [Y]"][index_electron]-shift_a_pad,
                             start_position_data["Position [Z]"][index_electron]]) #
    position_start_data_x.append(target_coord[0])
    position_start_data_z.append(target_coord[2])


    v0 = np.array([start_position_data["Velocity [X]"][index_electron],
                   start_position_data["Velocity [Y]"][index_electron],
                   start_position_data["Velocity [Z]"][index_electron]]) # [m/s]

    vabs = np.sqrt(v0[0] ** 2 + v0[1] ** 2 + v0[2] ** 2)
    energy = ( vabs / scipy.constants.speed_of_light) ** 2 * 0.5 * m
    vabs_list.append(vabs)
    alpha = np.arctan(v0[1] / v0[0])
    alpha_list.append(alpha)

    energies.append(energy)

    closest_x, closest_y, closest_z = find_closest_values(target_coord, x, y, z)
    print(f"Closest value to {target_coord[0]} in x: {x[closest_x]}")
    print(f"Closest value to {target_coord[1]} in y: {y[closest_y]}")
    print(f"Closest value to {target_coord[2]} in z: {z[closest_z]}")

    with h5py.File(file_path, 'r') as f:
        field_array = f['E-Field'][closest_x, closest_y, closest_z]
    field_x = field_array['x']
    field_y = field_array['y']
    field_z = field_array['z']

    # Create a numerical vector for the electric field
    field_vector = np.array([field_x, field_y, field_z])
    field_array = field_vector
    #print (f"field_array = {field_array*-1.6e-19}")

    current = -1.6e-19*np.dot(v0, field_array) #[m/s]  * [V/m] = [A]
    print(f"Current = {current} A")

    x0 = target_coord

    v0 = v0*1e-3
    t = solve_for_intercept_time(x0, v0, a0, 2850)#2850

    print(f"t = {t} time in ns")
    print (f"start velocity = {v0} um/ns")
    print (f"position_end = {step_position(x0, v0, a0, t)} um")
    position_end = step_position(x0, v0, a0, t)
    position_end_data_x.append(position_end[0])
    position_end_data_z.append(position_end[2])
    t_max = 5 # [ns]
    step = 5000
    induced_current = []
    field_array_steps_y = []
    start_time_of_electron = np.array(start_position_data["Time"][index_electron]) * 1e9
    print(f"start_time_of_electron = {start_time_of_electron}")

    for i in range (step):
        time_step = (t_max / step) * i
        if time_step >= t:
            print(f"time step stopped at  = {time_step} ns")
            sum_data.append({"Time": time_step + start_time_of_electron, "Current": 0})
            #for k in range(10):
            #    sum_data.append({"Time": time_step + start_time_of_electron + (k/10)*2, "Current": 0})
            break

        xi = step_position(x0, v0, a0, time_step)
        vi = step_velocity(v0, a0, time_step)

        closest_x, closest_y, closest_z = find_closest_values(xi, x, y, z)


        with h5py.File(file_path, 'r') as f:
            field_array = f['E-Field'][closest_x, closest_y, closest_z]

        print(f"Closest value to {xi[0]} in x: {x[closest_x]}")
        print(f"Closest value to {xi[1]} in y: {y[closest_y]}")
        print(f"Closest value to {xi[2]} in z: {z[closest_z]}")
        field_array = np.array(field_array)
        field_x = field_array['x']
        field_y = field_array['y']
        field_z = field_array['z']

        # Create a numerical vector for the electric field
        field_vector = np.array([field_x, field_y, field_z])
        field_array = field_vector
        field_array_steps_y.append(field_array[1])
        print (f"field_array = {field_array}")

        energy = ((np.sqrt(vi[0] ** 2 + vi[1] ** 2 + vi[2] ** 2)) / c) ** 2 * 0.5 * m
        current = 1.6e-19 * np.dot(vi*1e3, field_array)
        induced_current.append(current)
        print(f"Current = {current} A")


        sum_data.append({"Time": time_step+start_time_of_electron, "Current": current})

        print(f"Energy = {energy} eV")

        print(f"x = {xi}")

df = pd.DataFrame(sum_data)
zero_index = df[df['Current'] == 0].index[0]
plt.plot(df["Time"][:zero_index+1], df["Current"][:zero_index+1])
plt.plot(df["Time"][zero_index+1:], df["Current"][zero_index+1:])
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

# Sort the result_sort by Time
result_sort = result_sort.sort_values(by='Time').reset_index(drop=True)

rolling_sum = result_sort['Current'].rolling(window=250).sum()
plt.figure()
plt.plot(result_sort['Time'], rolling_sum)
plt.xlabel("Time (ns)")
plt.ylabel("Current (A)")
plt.title("Induced Charge")
plt.xlim(0, 10)
plt.figure()
plt.scatter(position_end_data_x, position_end_data_z, color='b',label='python',alpha=0.5)

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
print("test")
print (alpha_list)
alpha_list = np.array(alpha_list)
plt.hist(90-alpha_list*(180/np.pi), bins=100)
plt.xlabel('Energy[eV]')


plt.show()
