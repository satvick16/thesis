lambda_ = 0.1
mu = 0.04359
cox = 0.0197
L = 45e-9
vgs = 0.75
vth = 0.47
fsat = 0.8

Rdrv = 55
Rrep = 275

percent = 0.01

def main_driver():
    Rterm = 2 * Rdrv

    # top 2: ro || 1/gme = 2Rterm
    print("Main Driver Top 2")
    for W in range(1, 100000, 1):
        W = W * 1e-9
        gme = fsat * mu * cox * (W / L) * (vgs - vth)
        ro = 1 / (lambda_ * mu * cox * (W / L) * (vgs - vth) * (vgs - vth))
        # check if within 5% of target
        threshold = percent * 2 * Rterm
        if abs((ro * (1 / gme)) / (ro + (1 / gme)) - 2 * Rterm) < threshold:
            print(f"\tWidth: {W*1e9} nm")
            # print((ro * (1 / gme)) / (ro + (1 / gme)))
            # print(2 * Rterm)

    # bottom 2: ro = 2Rterm
    print("Main Driver Bottom 2")
    for W in range(1, 100000, 1):
        W = W * 1e-9
        ro = 1 / (lambda_ * mu * cox * (W / L) * (vgs - vth) * (vgs - vth))
        # check if within 5% of target
        threshold = percent * 2 * Rterm
        if abs(ro - 2 * Rterm) < threshold:
            print(f"\tWidth: {W*1e9} nm")
            # print(ro)
            # print(2 * Rterm)


def replica_driver():
    Rterm = 2 * Rrep

    # top 2: ro || 1/gme = 2Rterm
    print("Replica Driver Top 2")
    for W in range(1, 100000, 1):
        W = W * 1e-9
        gme = fsat * mu * cox * (W / L) * (vgs - vth)
        ro = 1 / (lambda_ * mu * cox * (W / L) * (vgs - vth) * (vgs - vth))
        # check if within 5% of target
        threshold = percent * 2 * Rterm
        if abs((ro * (1 / gme)) / (ro + (1 / gme)) - 2 * Rterm) < threshold:
            print(f"\tWidth: {W*1e9} nm")
            # print((ro * (1 / gme)) / (ro + (1 / gme)))
            # print(2 * Rterm)

    # bottom 2: ro = 2Rterm
    print("Replica Driver Bottom 2")
    for W in range(1, 100000, 1):
        W = W * 1e-9
        ro = 1 / (lambda_ * mu * cox * (W / L) * (vgs - vth) * (vgs - vth))
        # check if within 5% of target
        threshold = percent * 2 * Rterm
        if abs(ro - 2 * Rterm) < threshold:
            print(f"\tWidth: {W*1e9} nm")
            # print(ro)
            # print(2 * Rterm)
    

if __name__ == "__main__":
    main_driver()
    replica_driver()
