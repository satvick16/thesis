def parallel(r1, r2):
    return (r1 * r2) / (r1 + r2)

rch = 9.5
rterm = 100
rtia = 200

for rdrv in range(10, 1000, 1):
    rrep = rdrv * 5
    rs = rdrv
    k = (rs + rdrv) / rterm

    pairs = []
    for rh1 in range(1, 1000, 1):
        if rh1 >= 2 * rs:
            rhs = rh1 * (1 + k * rterm * ((1 / (rch + rterm)) + (1 / rh1)))
            rh2 = rhs - rrep
            if rh2 > 0:
                pairs.append((rh1, rh2))

    for rh1, rh2 in pairs:
        rout = rdrv + rs
        rhac = rh1 + parallel(rrep + rh2, rtia)

        rtx = parallel(rout, rhac)

        if (rtx > 0.99 * rterm) and (rtx < 1.01 * rterm):
            # if rdrv == 55:
            print(f"rdrv = {rdrv} Ω, rrep = {rrep} Ω, R_h1 = {rh1:.2f} Ω, R_h2 = {rh2:.2f} Ω, Rtx = {rtx:.2f} Ω")
