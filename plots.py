
import matplotlib.pyplot as plt
import numpy as np


def plot_field(data, z_position):
    #data['Abs_E'] = np.sqrt(data['Ex'] ** 2 + data['Ey'] ** 2 + data['Ez'] ** 2)
    filtered_data =  data[data['z'] == z_position]

    x = filtered_data['x'].values
    y = filtered_data['y'].values
    z = filtered_data['Ex'].values
    levels = 30
    plt.tricontourf(x, y, z, levels=levels, cmap='plasma')

    plt.colorbar(label='|E|')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Contour plot at z = {z_position}')
    plt.show()