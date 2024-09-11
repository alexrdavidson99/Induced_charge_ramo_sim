import numpy as np
import pandas as pd
import scipy.constants
from plots import plot_field
import matplotlib.pyplot as plt

def solve_for_intercept_time(x0, v0, acc, target_distance):
    """
    Solve for the intercept time when the particle reaches the target distance in the y direction.
    """
    # Define the polynomial coefficients for the equation of motion in the y direction
    coeffs = [
        0.5 * acc[1],  # t^2 term
        v0[1],  # t term
        x0[1] - target_distance  # constant term
    ]

    # Solve the quadratic equation
    roots = np.roots(coeffs)

    # Filter out complex roots and negative times
    real_roots = roots[np.isreal(roots) & (roots >= 0)]

    if len(real_roots) == 0:
        raise ValueError("No valid intercept time found.")

    # Return the smallest positive real root
    return np.min(real_roots)


def step_position(x0, v0, acc, time):
    """
    Position after time step
    """
    x0 = np.array(x0)
    v0 = np.array(v0)
    acc = np.array(acc)
    return x0 + v0 * time + 0.5 * acc * time ** 2

def step_velocity( v0, acc, time ):
    '''
    Velocity after time step
    '''
    return v0 + acc*time

def find_target(target_coord, nom_E_field):
    # Calculate the Euclidean distance between the target coordinate and all points in nom_E_field
    distances = np.sqrt(
        (nom_E_field['x'] - target_coord[0]) ** 2 +
        (nom_E_field['y'] - target_coord[1]) ** 2 +
        (nom_E_field['z'] - target_coord[2]) ** 2
    )

    # Find the index of the minimum distance
    min_index = distances.idxmin()
    print(f"{nom_E_field['x'][min_index]}")

    nom_E_field_pos_x = nom_E_field['Ex'][min_index]
    nom_E_field_pos_y = nom_E_field['Ey'][min_index]
    nom_E_field_pos_z = nom_E_field['Ez'][min_index]

    field_array = np.array([nom_E_field_pos_x, nom_E_field_pos_y, nom_E_field_pos_z])

    return field_array

if __name__ == "__main__":

    # constants
    c = scipy.constants.speed_of_light * 1e-3  # in mm/ns
    V = 1000  # electrods potential in V
    d = 2900  # pore depth in mm
    m = 511e3
    E = V * (c ** 2) / (d * m)  # electric field acceleration in mm/ns^2
    orientation = np.array([0., 1., 0])
    a0 = E * orientation

    # data C:/Users/lexda/Downloads/
    nom_E_field = pd.read_csv("C:/Users/lexda/Downloads/E-Field [Es].txt", names=["x", "y", "z", "Ex", "Ey", "Ez"], skiprows=3, sep='\s+')
    start_position_data = pd.read_csv("ascii_export(2).csv", comment='#', skip_blank_lines=True, sep=';',
                                      header=None, names=["Position [X]","Position [Y]","Position [Z]",
                                                          "Position [ABS (XYZ)]","Time","Velocity [X]","Velocity [Y]",
                                                          "Velocity [Z]"], index_col=False)
    # (start_position_data)
    #3.198215199407748e-14
    # -6.123234e-14 cut
    sum_data = []
    plot_field(nom_E_field,0)
    for j in range(len(start_position_data["Velocity [X]"])-812):

        index_electron = j
        target_coord = np.array([start_position_data["Position [X]"][index_electron],
                                 start_position_data["Position [Y]"][index_electron],
                                 start_position_data["Position [Z]"][index_electron]]) # [um]

        v0 = np.array([start_position_data["Velocity [X]"][index_electron],
                       start_position_data["Velocity [Y]"][index_electron],
                       start_position_data["Velocity [Z]"][index_electron]]) # [m/s]

        field_array = find_target(target_coord, nom_E_field)  # [V/m]
        print (f"field_array = {field_array}")

        current = -1.6e-19*np.dot(v0, field_array) #[m/s]  * [V/m] = [A]
        print(f"Current = {current} A")

        x0 = target_coord

        v0 = v0*1e-3
        t = solve_for_intercept_time(x0, v0, a0, d)

        print(f"t = {t} time in ns")
        print (f"start velocity = {v0} um/ns")
        print (f"position = {step_position(x0, v0, a0, t)} um")

        t_max = 0.8 # [ns]
        step = 500
        induced_current = []
        field_array_steps_y = []
        start_time_of_electron = np.array(start_position_data["Time"][index_electron]) * 1e9
        print(f"start_time_of_electron = {start_time_of_electron}")

        for i in range (step):
            time_step = (t_max / step) * i
            if time_step >= t:
                print(f"time step stopped at  = {time_step} ns")
                sum_data.append({"Time": time_step + start_time_of_electron, "Current": 0})
                break

            xi = step_position(x0, v0, a0, time_step)
            vi = step_velocity(v0, a0, time_step)

            field_array = find_target(xi, nom_E_field)
            field_array_steps_y.append(field_array[1])
            print (f"field_array = {field_array}")
            energy = ((np.sqrt(vi[0] ** 2 + vi[1] ** 2 + vi[2] ** 2)) / c) ** 2 * 0.5 * m
            current = 1.6e-19 * np.dot(vi*1e3, field_array)
            induced_current.append(current)
            print(f"Current = {current} A")


            sum_data.append({"Time": time_step+start_time_of_electron, "Current": current})

            print(f"t = {t} time in ns")

            print(f"Energy = {energy} eV")

            print(f"x = {xi}")

    df = pd.DataFrame(sum_data)
    zero_index = df[df['Current'] == 0].index[0]
    plt.plot(df["Time"][:zero_index+1], df["Current"][:zero_index+1])
    plt.plot(df["Time"][zero_index+1:], df["Current"][zero_index+1:])
    result_sort = df.sort_values(["Time"], ascending=[True])
    rolling_sum = result_sort['Current'].rolling(window=500).sum()
    plt.figure()
    plt.plot(result_sort['Time'], rolling_sum)
    plt.xlabel("Time (ns)")
    plt.ylabel("Current (A)")
    plt.title("Induced Charge")
    print (rolling_sum)
    print(result_sort )
    plt.show()
    pass