import pyesgi140 as e
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from visualization import make_runner_video
from samples import ranked_samples

n = 8000
vmaxs = ranked_samples(n)
# fractions = [0.3, 0.3, 0.4]
# times = [0, 150, 300]

tstarts = [0] * (n//3) + [300] * (n//3) + [600] * (n-n//3-n//3)

t1 = datetime.now()
constriction_at = 1000
width = e.Barrier(3000, 100, 200, 7.5, 3.0)

res_t, res_x, res_v, res_rho = e.solve_particle_model(vmaxs, tstarts, width)
t2 = datetime.now()
print("Time to solve: ", t2-t1)
res_x = np.asarray(res_x)
res_rho = np.asarray(res_rho)
res_v = np.asarray(res_v)
# np.save(
make_runner_video(res_t, res_x, res_v, res_rho, width)
import sys
sys.exit(1)
plt.figure()
for i in range(n):
    plt.plot(res_t, res_x[:, i])
plt.xlabel("Time")
plt.ylabel("Distance")
# plt.xlim((0,100))
# plt.ylim((0,500))
plt.savefig("start.pdf")
plt.show()
