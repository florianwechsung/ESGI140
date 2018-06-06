import pyesgi140 as e
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from visualization import make_runner_video

n = 1000

t1 = datetime.now()
constriction_at = 1000
def width(x):
    if x<constriction_at-100 or x>constriction_at + 100:
        return 2.
    else:
        return 1.
res_t, res_x, res_rho = e.solve_particle_model(n, width)
t2 = datetime.now()
print("Time to solve: ", t2-t1)
print(len(res_t))
print(res_t)
res_x = np.asarray(res_x)
res_rho = np.asarray(res_rho)
make_runner_video(res_t, res_x, res_rho, width)
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
