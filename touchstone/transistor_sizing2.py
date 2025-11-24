# Improved width calculator for desired on-resistance
# - Units: mobility in m^2/(V·s), Cox in F/m^2, L in m, voltages in V, R in ohms.
# - Two estimates returned:
#    1) linear-region R_on estimate (valid when VDS is small compared to VGS-VT)
#    2) saturation-equivalent estimate (useful if transistor enters saturation; gives width
#       such that Isat ~= Vdd/R)

import math

# device / bias parameters (as provided)
cox_n = 0.01973       # F/m^2 (NMOS)
cox_p = 0.018665      # F/m^2 (PMOS)

mu_n = 0.04359        # m^2/(V·s) (electron mobility) -> 435.9 cm^2/Vs
mu_p = 0.00432        # m^2/(V·s) (hole mobility)     -> 43.2 cm^2/Vs

Vtn = 0.471           # V (NMOS threshold)
Vtp = 0.423           # V (PMOS threshold magnitude; treat as positive magnitude)
Vgs = 0.75            # V (gate drive)
Vdd = 0.75            # V (supply)

L = 45e-9             # m (channel length)

# desired resistances
R_driver = 100.0      # ohm
R_replica = 200.0     # ohm

def widths_for_R(R_target, mu, cox, Vth, Vgs, L):
    Vov = Vgs - Vth  # overdrive
    if Vov <= 0:
        raise ValueError(f"Vgs ({Vgs}) <= Vth ({Vth}); transistor is off or marginal.")
    # 1) linear-region width (small Vds) from R_on ≈ 1 / (mu * Cox * (W/L) * Vov)
    #    => W_linear = L / (R * mu * Cox * Vov)
    W_linear = L / (R_target * mu * cox * Vov)

    # 2) saturation-equivalent width:
    #    Isat ≈ 0.5 * mu * Cox * (W/L) * Vov^2. If we want Vdd / Isat = R_target:
    #    => W_sat = (2 * Vdd * L) / (R_target * mu * Cox * Vov^2)
    W_sat = (2.0 * Vdd * L) / (R_target * mu * cox * Vov * Vov)

    # sanity-check reconstructed resistances
    R_from_linear = 1.0 / (mu * cox * (W_linear / L) * Vov)  # should ~= R_target
    Isat = 0.5 * mu * cox * (W_sat / L) * Vov * Vov
    R_from_sat_equiv = Vdd / Isat if Isat > 0 else math.inf

    return {
        "W_linear_m": W_linear,
        "W_linear_um": W_linear * 1e6,
        "R_check_linear_ohm": R_from_linear,
        "W_sat_m": W_sat,
        "W_sat_um": W_sat * 1e6,
        "R_check_sat_ohm": R_from_sat_equiv,
        "Isat_at_W_sat_A": Isat,
    }

# compute for driver (100 ohm) and replica (200 ohm) both NMOS and PMOS
for name, R in [("driver", R_driver), ("replica", R_replica)]:
    print(f"\n== {name} target R = {R} ohm ==")
    # NMOS
    try:
        n = widths_for_R(R, mu_n, cox_n, Vtn, Vgs, L)
        print(f"NMOS linear-width: {n['W_linear_m']:.3e} m  ({n['W_linear_um']:.3f} µm)")
        print(f"NMOS sat-equivalent width: {n['W_sat_m']:.3e} m  ({n['W_sat_um']:.3f} µm)")
        print(f"  check R (linear-model) = {n['R_check_linear_ohm']:.3f} ohm")
        print(f"  check R (sat-equivalent) = {n['R_check_sat_ohm']:.3f} ohm, Isat = {n['Isat_at_W_sat_A']:.3e} A")
    except ValueError as e:
        print("NMOS:", e)

    # PMOS (note: using magnitude of threshold and same Vgs magnitude; treat Vov = Vgs - |Vtp|)
    try:
        p = widths_for_R(R, mu_p, cox_p, Vtp, Vgs, L)
        print(f"PMOS linear-width: {p['W_linear_m']:.3e} m  ({p['W_linear_um']:.3f} µm)")
        print(f"PMOS sat-equivalent width: {p['W_sat_m']:.3e} m  ({p['W_sat_um']:.3f} µm)")
        print(f"  check R (linear-model) = {p['R_check_linear_ohm']:.3f} ohm")
        print(f"  check R (sat-equivalent) = {p['R_check_sat_ohm']:.3f} ohm, Isat = {p['Isat_at_W_sat_A']:.3e} A")
    except ValueError as e:
        print("PMOS:", e)
