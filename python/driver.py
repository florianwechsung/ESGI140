import pyesgi140 as e
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from visualization import make_runner_video
from samples import ranked_samples
from scipy import interpolate

tmax = 75*60
dt = 1.0
rho_start = 4.

def run_simulation(sizes, generate_histogram=False, generate_video=False):
    vmaxs = ranked_samples(sum(sizes))
    tstarts = [0] * sizes[0] + [240] * sizes[1] + [516] * sizes[2]
    # 

    width_data = np.loadtxt("./../data/road_width.csv", delimiter=";")
    width = interpolate.interp1d(width_data[:,0], width_data[:, 1], kind="linear", bounds_error=False, fill_value=(width_data[0, 1], width_data[-1, 1]))
    width_callable = lambda x: width(x)

    t1 = datetime.now()
    settings = e.SpeedSettings(rho_1=0.3, rho_2=2., v_min=0.1)
    def slope(x):
        if x<= 0: return 0.
        elif x<1000: return 1.
        elif x<6000: return 0.
        elif x<7500: return -1.
        elif x<8500: return -2.
        elif x<9000: return 0.
        else: return 2

    res_t, res_x, res_v, res_rho = e.solve_particle_model(vmaxs, tstarts, width_callable, slope, settings, tmax=tmax, dt=dt, rho_start=rho_start)
    t2 = datetime.now()
    print("Time to solve: ", t2-t1)
    wasted_time = e.calculate_average_wasted_time(res_x, vmaxs, res_t, 10000);
    print("Average wasted time per runner: ", wasted_time, " when using wave sizes ", sizes)
    times_5k = e.get_finish_times(res_x, vmaxs, res_t, 5000)
    print("Fastest runner for 5K: ", min(times_5k)/60)
    times_10k = e.get_finish_times(res_x, vmaxs, res_t, 10000)
    print("Fastest runner for 10K: ", min(times_10k)/60)

    if generate_histogram:
        plt.figure()
        bins = list(range(int(min(times_5k)), int(max(times_5k)), 60))
        plt.hist(times_5k, bins=bins, color="blue")
        plt.hist(times_5k[0:sizes[0]], bins=bins, color="red", alpha=0.6)
        plt.hist(times_5k[sum(sizes[0:1]):sum(sizes[0:2])], bins=bins, color="orange", alpha=0.6)
        plt.hist(times_5k[sum(sizes[0:2]):sum(sizes[0:3])], bins=bins, color="green", alpha=0.6)
        plt.xlabel("Time (in s)")
        plt.ylabel("Number of runners")
        plt.title("5k finish times")
        plt.savefig("5k_times.pdf")
        plt.figure()

        bins = list(range(int(min(times_10k)), int(max(times_10k)), 60))
        plt.hist(times_10k, bins=bins, color="blue")
        plt.hist(times_10k[0:sizes[0]], bins=bins, color="red", alpha=0.6)
        plt.hist(times_10k[sum(sizes[0:1]):sum(sizes[0:2])], bins=bins, color="orange", alpha=0.6)
        plt.hist(times_10k[sum(sizes[0:2]):sum(sizes[0:3])], bins=bins, color="green", alpha=0.6)
        plt.xlabel("Time (in s)")
        plt.ylabel("Number of runners")
        plt.title("10k finish times")
        plt.savefig("10k_times.pdf")

    if generate_video:
        res_x = np.asarray(res_x)
        res_rho = np.asarray(res_rho)
        res_v = np.asarray(res_v)
        make_runner_video(res_t, res_x, res_v, res_rho, width, slope, tmax//10, xmax=10000)

    return wasted_time


# these were the wave sizes in the actual race
sizes = [1336, 2976, 1712] 
run_simulation(sizes, generate_video=True, generate_histogram=False)

# Uncomment line below if you want to run the simulation for many wave sizes
import sys; sys.exit(0)
sizes1 = range(1000, 3001, 100)
sizes2 = range(1000, 3001, 100)
wasted_times = np.zeros((len(sizes1), len(sizes2)))
for i, size1 in enumerate(sizes1):
    for j, size2 in enumerate(sizes2):
        sizes = [size1, size2, 6024-size1-size2]
        wasted_time = run_simulation(sizes)
        wasted_times[i, j] = wasted_time

wasted_times = np.asarray(wasted_times)
plt.figure()
plt.contourf(sizes1, sizes2, wasted_times, 20, cmap='cool')
plt.xlabel("Size of wave 1")
plt.ylabel("Size of wave 2")
plt.title("Lost time for different wave sizes")
plt.colorbar()
plt.savefig("contour.pdf")

