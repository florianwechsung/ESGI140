import pyesgi140 as e
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

n = 1000

t1 = datetime.now()
res_t, res_x = e.solve_particle_model(n)
t2 = datetime.now()
print("Time to solve: ", t2-t1)
print(len(res_t))
res_x = np.asarray(res_x).T
# plt.figure()
# ax = sns.tsplot(data=res_x, time=res_t, err_style="boot_traces")
# plt.show()
plt.figure()
for i in range(n):
    plt.plot(res_t, res_x.T[:, i])
plt.show()
