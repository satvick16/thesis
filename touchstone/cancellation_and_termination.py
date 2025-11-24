import numpy as np

def parallel(r1, r2):
    return (r1 * r2) / (r1 + r2)

def find_rh_pairs(rch, rterm, k, rrep, rh1_range, rs):
    results = []
    for rh1 in rh1_range:
        if rh1 < 2 * rs:  # enforce rh1 at least double rs
            continue
        rhs = rh1 * (1 + k * rterm * ((1 / (rch + rterm)) + (1 / rh1)))
        rh2 = rhs - rrep
        if rh2 > 0:
            results.append((rh1, rh2))
    return results

rch = 9.5
rterm = 100
rh1_values = np.linspace(10, 1000, 1000)
rtia = 200

# Search ranges for rdrv and rrep
rdrv_values = np.linspace(10, 1000, 50)   # adjust as needed
rrep_values = np.linspace(10, 1000, 50)   # adjust as needed

valid_configs = []

for rdrv in rdrv_values:
    rs = rdrv
    rout = rdrv + rs
    k = (rs + rdrv) / rterm

    for rrep in rrep_values:
        if rrep <= rdrv:  # enforce rrep > rdrv
            continue

        pairs = find_rh_pairs(rch, rterm, k, rrep, rh1_values, rs)
        for rh1, rh2 in pairs:
            rhac = rh1 + parallel(rrep + rh2, rtia)
            rtx = parallel(rout, rhac)

            if 0.999 * rterm < rtx < 1.001 * rterm:
                valid_configs.append((rdrv, rrep, rh1, rh2, rtx))
                # break  # optionally stop at first valid pair per (rdrv, rrep)

# Print all valid configurations
for cfg in valid_configs:
    rdrv, rrep, rh1, rh2, rtx = cfg
    print(f"rdrv = {rdrv:.2f} Ω, rrep = {rrep:.2f} Ω, rh1 = {rh1:.2f} Ω, rh2 = {rh2:.2f} Ω, Rtx = {rtx:.2f} Ω")
