rdrv = 38.8636363636
rs = 22.5
rh1 = 190
rrep = 273
z0 = 50
rterm = z0
k = (rdrv + rs) / rterm
rch = 0.9

rh2 = rh1 * (1 + k * rterm * ((1 / (rch + rterm)) + (1 / rh1))) - rrep

print(rh2)
