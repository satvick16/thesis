import numpy as np

def parallel(r1, r2):
    return (r1 * r2) / (r1 + r2)

def find_rh_pairs(rch, rterm, k, rrep, rh1_range):
    results = []
    for rh1 in rh1_range:
        if rh1 >= 2 * rs:
            rhs = rh1 * (1 + k * rterm * ((1 / (rch + rterm)) + (1 / rh1)))
            rh2 = rhs - rrep
            if rh2 > 0:
                results.append((rh1, rh2))
    return results

rch = 9.5
rterm = 100
rdrv = 50
rs = rdrv
k = (rs + rdrv) / rterm
rrep = rdrv * 5

rh1_values = np.linspace(1, 200, 10000)
pairs = find_rh_pairs(rch, rterm, k, rrep, rh1_values)

for rh1, rh2 in pairs:
    rtia = 200
    rout = rdrv + rs
    rhac = rh1 + parallel(rrep + rh2, rtia)

    rtx = parallel(rout, rhac)
    print(rtx)

    if (rtx > 0.999 * rterm) and (rtx < 1.001 * rterm):
        print(f"R_h1 = {rh1:.2f} Ω, R_h2 = {rh2:.2f} Ω, Rtx = {rtx:.2f} Ω")
