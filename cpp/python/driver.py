import pyesgi140 as e
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

n = 500

t1 = datetime.now()
constriction_at = 5000
def width(x):
    if x<constriction_at-100 or x>constriction_at + 100:
        return 7.5
    else:
        return 0.0075
res_t, res_x = e.solve_particle_model(n, width)
t2 = datetime.now()
print("Time to solve: ", t2-t1)
print(len(res_t))
print(res_t)
res_x = np.asarray(res_x).T
# plt.figure()
# ax = sns.tsplot(data=res_x, time=res_t, err_style="boot_traces")
# plt.show()
plt.figure()
for i in range(n):
    plt.plot(res_t, res_x.T[:, i])
plt.xlabel("Time")
plt.ylabel("Distance")
# plt.xlim((0,100))
# plt.ylim((0,500))
plt.savefig("start.pdf")
plt.show()
