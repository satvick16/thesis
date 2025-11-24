import numpy as np
import matplotlib.pyplot as plt

rdrv = 40
rs = 22.5
rh1 = 190
# rrep = 315
rrep = 40
z0 = 50
rterm = z0
k = (rdrv + rs) / rterm

# Use numpy.arange for floating point steps
rch_options = np.arange(0, 50, 0.5)
rh2 = [rh1 * (1 + k * rterm * ((1 / (rch + rterm)) + (1 / rh1))) - rrep for rch in rch_options]

plt.plot(rch_options, rh2)
plt.xlabel('R_ch (Ohms)')
plt.ylabel('R_h2 (Ohms)')
plt.title('R_h2 vs R_ch')
plt.grid(True)
plt.show()
