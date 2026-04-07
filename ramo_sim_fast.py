import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import h5py
import time
from scipy import constants
import mplhep
mplhep.style.use(mplhep.style.LHCb2)
import os



from main import step_position, step_velocity, solve_for_intercept_time

start = time.time()

# ------------------------------
# Helper: Fast nearest index lookup
# ------------------------------
def find_closest_values(target_value, x, y, z):
    closest_value_x = np.searchsorted(x, target_value[0])
    closest_value_y = np.searchsorted(y, target_value[1])
    closest_value_z = np.searchsorted(z, target_value[2])
    return closest_value_x, closest_value_y, closest_value_z

def find_closest_indices(arr, values):
    """
    Vectorized nearest-neighbor index lookup for non-uniform grids.
    arr: 1D sorted coordinate array
    values: 1D array of positions along the same axis
    Returns array of indices of closest points in arr for each value.
    """
    idx = np.searchsorted(arr, values)
    idx = np.clip(idx, 1, len(arr) - 1)
    left = arr[idx - 1]
    right = arr[idx]
    idx -= values - left < right - values
    return idx

# ------------------------------
# Constants
# ------------------------------
c = constants.speed_of_light * 1e-3  # in um/ns
V = 1500  # electrode potential in V
d = 2510  # mcp anode gap in um
m = 511e3
E = V * (c ** 2) / (d * m)  # electric field acceleration in um/ns^2
E_field = V / (d * 1e-6)
print(f"E field = {E_field}")
orientation = np.array([0., 1., 0])
a0 = E * orientation

# ------------------------------
# Load E-Field data
# ------------------------------
file_path = './3d/E-Field [Es]_unifrom_filed.h5'
with h5py.File(file_path, 'r') as f:
    x = f['Mesh line x'][:]
    y = f['Mesh line y'][:]
    z = f['Mesh line z'][:]
    E_raw = f['E-Field'][:]

    # Split into Ex, Ey, Ez for faster access
    Ex, Ey, Ez = E_raw['x'], E_raw['y'], E_raw['z']

# Create interpolators for smooth, fast E-field lookups


# ------------------------------
# Load start positions
#C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/ascii_export_of_15d_pore.csv
#C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/ascii_export_ramo_600k_elestic.csv
#C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/ascii_export_Torch_layout_127_pores_1.4_million.csv - first used for pd-25

#ascii_export_sey_Loffler_2022_487319.csv

# ------------------------------ DEFULT!!! ascii_export_sey_Loffler_2022_1.007m ----------------------
start_position_data = pd.read_csv(
    "C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/ascii_export_sey_Loffler_2022_40232.csv",
    comment='#', skip_blank_lines=True, sep=';',
    header=None,
    names=[
        "Position [X]", "Position [Y]", "Position [Z]",
        "Position [ABS (XYZ)]", "Time",
        "Velocity [X]", "Velocity [Y]", "Velocity [Z]"
    ],
    index_col=False
)

# ------------------------------
# Simulation setup
# ------------------------------
num_electrons = len(start_position_data)
print(f"Number of electrons = {num_electrons}")

shift_a_pads = [-1035,-485,65,615,1165]
#shift_a_pads = [-480,70,620]
#shift_a_pads = range(-1035, 1165, 50) #[-490,610] #range(-530, 580, 50) x ----> range(-3030, 2580, 250)
#shift_a_pads_z = [65] #65 is where loffler 2022 is maximum 
#shift_a_pads_x = range(-530, 2580, 250)  # 670 430 z=fine x=coarse
batch_size = 10000
total_range = num_electrons

# Preallocate results
position_end_data_x = []
position_end_data_z = []

# Store all currents and times globally
all_times = []
all_currents = []

# ------------------------------
# Main simulation loop
# ------------------------------
for shift_a_pad in shift_a_pads:
    for batch_start in range(0, total_range, batch_size):
        sum_times = []
        sum_currents = []

        for j in range(batch_start, min(batch_start + batch_size, total_range)):
            index_electron = j
            x0 = np.array([
                start_position_data["Position [X]"][index_electron],
                start_position_data["Position [Y]"][index_electron],
                start_position_data["Position [Z]"][index_electron] + shift_a_pad
            ])

            v0 = np.array([
                start_position_data["Velocity [X]"][index_electron],
                start_position_data["Velocity [Y]"][index_electron],
                start_position_data["Velocity [Z]"][index_electron]
            ])  # [m/s]
            energy = 0.5 * m * np.sum(v0)/(c**2)
              # in eV
            v0_um = v0 * 1e-3  # convert to um/ns
            target_coord = x0[1] + d # target y position at anode
            # Solve for intercept time
            t = solve_for_intercept_time(x0, v0_um, a0, target_coord)
            #print(f"Electron {j}, intercept time: {t:.4f} ns")

            # End position
            position_end = step_position(x0, v0_um, a0, t)
            position_end_data_x.append(position_end[0])
            position_end_data_z.append(position_end[2])

            # Vectorize kinematics
            #steps = 1000
            #time_steps = np.linspace(0, t, steps)
            # Fixed timestep
            dt = 0.001  # ns
            steps = int(np.ceil(t / dt)) + 1 # total number of steps until intercept

            # Precompute time steps (stop at t)
            time_steps = np.arange(0, steps * dt, dt)
            time_steps = time_steps[time_steps <= t]


            xi = x0 + v0_um * time_steps[:, None] + 0.5 * a0 * time_steps[:, None] ** 2
            vi = v0_um + a0 * time_steps[:, None]

            # Interpolate E-field at all steps (vectorized)
            closest_x, closest_y, closest_z = find_closest_values(xi, x, y, z)
            closest_x = find_closest_indices(x, xi[:, 0])
            closest_y = find_closest_indices(y, xi[:, 1])
            closest_z = find_closest_indices(z, xi[:, 2])

            Ex_here = Ex[closest_z, closest_y, closest_x]
            Ey_here = Ey[closest_z, closest_y, closest_x]
            Ez_here = Ez[closest_z, closest_y, closest_x]
            E_here = np.column_stack((Ex_here, Ey_here, Ez_here))

            # Induced current: q * v · E
            q = 1.6e-19
            currents = q * np.einsum("ij,ij->i", vi * 1e3, E_here)

            # Stop current after arrival time
            currents[time_steps >= t] = 0

            # Store results for this electron
            start_time_of_electron = np.array(start_position_data["Time"][index_electron]) * 1e9
            sum_times.extend(time_steps + start_time_of_electron)
            sum_currents.extend(currents)

        # Store batch results
        

        # ------------------------------
        # Combine results using NumPy instead of Pandas groupby
        # ------------------------------
        all_times = np.round(np.array(sum_times), 4)
        all_currents = np.array(sum_currents)

        # Aggregate currents per time using NumPy
        unique_times, inverse_idx = np.unique(all_times, return_inverse=True)
        grouped_currents = np.zeros_like(unique_times)
        np.add.at(grouped_currents, inverse_idx, all_currents)

        # Save results
        df_grouped = pd.DataFrame({"Time": unique_times, "Current": grouped_currents})
        # make folder if not exists
        
        # output_folder = f"C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/second_run_{num_electrons}_electrons_1500v_new_field_with_4_pix"
        # os.makedirs(output_folder, exist_ok=True)

        # df_grouped.to_csv(
        #     f'{output_folder}/fast_induced_charge_Loffler_2022_1500v_fix_{j-batch_size}_{j}_off_center_z_{shift_a_pad}-{num_electrons}-elastic-0.25_TORCH_17_pores-new-field-one-electron.csv', 
        #     index=False    
        # )

        output_folder = f"C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/induced_current_loffler/run_{num_electrons}_electrons_1500v_new_field_with_4_pix"
        os.makedirs(output_folder, exist_ok=True)

        filename = f"induced_charge_{j-batch_size}_{j}_z{shift_a_pad}_{num_electrons}.csv"

        df_grouped.to_csv(
            os.path.join(output_folder, filename),
            index=False
        )

        print(f"saved batch")

# ------------------------------
# Plot once at the end
# ------------------------------
plt.figure(figsize=(10, 6))
#add a 0 current after the final df_grouped current value to make the plot go back to 0 at the end
df_grouped = df_grouped._append({"Time": df_grouped["Time"].iloc[-1] + 0.001, "Current": 0}, ignore_index=True)


plt.plot(df_grouped["Time"], df_grouped["Current"]*1e3, color="darkblue")
plt.xlabel("Time (ns)")
plt.ylabel("Induced current (mA)")
#plt.title("Induced Charge")
plt.grid()
#plt.tight_layout()

plt.savefig(f'induce_current_plot_1500v_fix_1e_off_center_z_{shift_a_pad}_from-{num_electrons}-elastic-0.25.pdf')
plt.show()
end = time.time()
print(f"Total time: {end - start:.2f} seconds")
# ------------------------------

x_z_data = {"x_position": position_end_data_x,
           "z_position": position_end_data_z }
df_x_z = pd.DataFrame(x_z_data)
df_x_z.to_csv(f'C:/Users/lexda/PycharmProjects/Induced_charge_ramo_sim/x_z_data_from_python_1500v_Loffler_2022-{num_electrons}-elastic-1electron.csv', index=False)