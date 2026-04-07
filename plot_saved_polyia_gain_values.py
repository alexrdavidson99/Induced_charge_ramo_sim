import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import scipy
from matplotlib.pyplot import figure
from scipy.interpolate import griddata
import h5py
import time
# 3D plot
from mpl_toolkits.mplot3d import Axes3D

from scipy.optimize import curve_fit

from matplotlib.patches import Polygon
from scipy.special import erf
import mplhep
mplhep.style.use(mplhep.style.LHCb2)


#Gaussian function
def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) 

def gaussian_convolved_tophat(x, amp, mean, stddev, width=0.55): #3.31
    return amp * (erf((x - mean + width / 2) / (np.sqrt(2) * stddev)) - erf((x - mean - width / 2) / (np.sqrt(2) * stddev))) / 2


gain_data = pd.read_csv("C:/Users/lexda/Downloads/plot_saved_polyia_gain_values_fft_1500_back_cut_pos_bad_fit.txt")
#gain_data = pd.read_csv("C:/Users/lexda/Downloads/plot_saved_polyia_gain_values_nol_con_330_back.csv")
gain_data = pd.read_csv("C:/Users/lexda/Downloads/gain_data_polya_90_hour_values.csv")
#gain_data = pd.read_csv("C:/Users/lexda/Downloads/plot_saved_polyia_gain_values_nol_con_330_back_cut_out_gain_drop_at_92.1.txt")
colors = {
    "F5": "#e03ba8",
    "F6": "#3aa3a5",
    "F7": "#f79c23",
    "F8": "#2e86f7"
}


pixel_label = {
    "F5": 4,
    "F6": 3,
    "F7": 2,
    "F8": 1,
}

plt.figure(figsize=(13, 7))

mask_F5 = (gain_data["Pixel"] == "F5")
#if mask_F5.any():
    #idx_to_drop = gain_data[mask_F5].index[-1]  # second-to-last index
    #gain_data = gain_data.drop(idx_to_drop)
#     idx_to_drop = gain_data[mask_F5].index[7]  # second-to-last index
    
#     idx_to_drop = gain_data[mask_F5].index[-3]  # second-to-last index
#     gain_data = gain_data.drop(idx_to_drop)

mask_F7 = (gain_data["Pixel"] == "F5")
if mask_F7.any():
    pos_mask = np.isclose(gain_data.loc[mask_F7, "Position_x"], 81.8, atol=1e-3)
    if pos_mask.any():
        idx_to_drop = gain_data.loc[mask_F7].index[pos_mask]
        gain_data = gain_data.drop(idx_to_drop)
        print(f"Dropped {len(idx_to_drop)} row(s) for F7 at Position_x ~ 97.3")
    else:
        print("No F7 rows with Position_x ~ 97.3 found")
# For F8: drop last point
mask_F8 = (gain_data["Pixel"] == "F8")
if mask_F8.sum() >= 1:
    idx_to_drop = gain_data[mask_F8].index[-1]  # last index
    gain_data = gain_data.drop(idx_to_drop)

target_pixel = next((p for p, lab in pixel_label.items() if lab == 1), None)
if target_pixel is not None:
    mask_target = (gain_data["Pixel"] == "F6")
    n = mask_target.sum()
    if n >= 3:
        # zero-based: index 7 is the 8th entry for that pixel group
        idx_to_drop = gain_data[mask_target].index[2]
        gain_data = gain_data.drop(idx_to_drop)


# Group by Pixel and plot each separately
for pixel, df_pixel in gain_data.groupby("Pixel"):

    xdata = df_pixel["Position_x"].values
    ydata = df_pixel["Gain"].values

    # Initial guess: amp=max(y), mean=position of max(y), stddev=0.1
    max_index = np.argmax(ydata)
    p0 = [ydata[max_index], xdata[max_index], 0.1]

    
    popt, pcov = curve_fit(gaussian_convolved_tophat, xdata, ydata, p0=p0)

        # Generate fitted curve
    fit_x = np.linspace(min(xdata)-10, max(xdata)+10, 10000)
    fit_y = gaussian_convolved_tophat(fit_x, *popt)

    print(f"Pixel {pixel}: Fit parameters: amp={popt[0]:.3f}, mean={popt[1]:.3f}, stddev={popt[2]:.3f}")
    print(f"fwhm = {2.355 * popt[2]:.3f} mm")

    plt.errorbar(df_pixel["Position_x"], df_pixel["Gain"],
                yerr=df_pixel["error"],
                color=colors[pixel],
                marker='o',
                markersize=6,              # smaller points
                markeredgecolor='k',
                markeredgewidth=0.5,
                markerfacecolor=colors[pixel],
                ecolor='k',                # darker error bars
                elinewidth=1.5,
                capsize=3,
                linestyle='None',
                label=f"Channel {pixel_label[pixel]} data")
    if pixel != "F7":
        plt.plot(fit_x, fit_y, color=colors[pixel], linestyle='-' , label=f"Channel {pixel_label[pixel]} fit")

    half_max = np.max(fit_y)/2
    print(np.max(fit_y))
    indices = np.where(fit_y >= half_max)[0]
    fwhm_true = fit_x[indices[-1]]-fit_x[indices[0]]
    print(fit_x[indices[-1]])
    print(fit_x[indices[0]])
    
    print(f"true FWHM = {fwhm_true} mm ")



# pixel_boundary_left = 83-0.552/2
# pixel_boundary_left_1 = pixel_boundary_left - 0.552
# pixel_boundary_left_2 = pixel_boundary_left - 2*0.552
# plexl_boundary_left_3 = pixel_boundary_left - 3*0.552
# pixel_boundary_right = 83+0.552/2

pixel_boundary_left = 93.6-3.31/2
pixel_boundary_right = 93.6+3.31/2
plt.vlines([pixel_boundary_left, pixel_boundary_right], ymin=0, ymax=7e5, color='gray', linestyle='--', label='Pixel boundary') 
# plt.vlines([pixel_boundary_left_1, pixel_boundary_left_2], ymin=0, ymax=5e5, color='gray', linestyle='--')
# plt.vlines([plexl_boundary_left_3], ymin=0, ymax=5e5, color='gray', linestyle='--')

plt.xlabel("X Coordinate (mm)",fontsize=32)
plt.xlim(80, 84.7)
#plt.xlim(90, 97)
plt.ylabel(r"Gain ($e^{-}$)", fontsize=32)
#plt.title("Gain as a function of position")
plt.legend(fontsize=16, loc='upper right')
plt.savefig("polyia_gain_vs_position_1500_V_back_long_pix_fit.png", dpi=300)





plt.show()

