import skrf as rf
import numpy as np

# === USER SETTINGS ===
s4p_path = "C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p"   # Path to your file
R0 = 50.0                       # Characteristic impedance per line
zs = R0                         # Source termination (Ω)
zl = R0                         # Load termination (Ω)
is_differential = False         # Set True if file is a 4-port differential channel

# === LOAD NETWORK ===
ntw = rf.Network(s4p_path)
print(f"Loaded {ntw.nports}-port network with freq range {ntw.f[0]/1e9:.3f}–{ntw.f[-1]/1e9:.3f} GHz")

# === IF DIFFERENTIAL, CONVERT TO MIXED-MODE 2-PORT ===
if is_differential and ntw.nports == 4:
    mm = ntw.copy()
    mm.renumber([0,1,2,3], [0,2,1,3])   # enforce correct ordering if needed
    mm = mm.s2mm()                      # convert to mixed-mode
    ntw = mm.subnetwork([0,1])          # keep differential-differential subnetwork
    print("Converted to differential mixed-mode 2-port.")

# === EXTRACT LOW-FREQ S-PARAMS ===
S = ntw.s[0]   # take lowest frequency point (index 0)
S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
S21_mag = abs(S21)
print(f"S21(0) = {S21_mag:.5f} magnitude, {np.angle(S21, deg=True):.2f}° phase")

# === COMPUTE REFLECTION COEFFICIENTS ===
z0 = 2 * R0    # matches your MATLAB code's "z0 = 2 * simSettings.z0"
Gamma1 = (2*zs - z0) / (2*zs + z0)
Gamma2 = (2*zl - z0) / (2*zl + z0)

# === COMPUTE H_ch USING SAME EQUATION AS MATLAB ===
DeltaS = S11*S22 - S12*S21
numerator = S21 * (1 - Gamma1) * (1 + Gamma2)
denominator = 1 - S11*Gamma1 - S22*Gamma2 + Gamma1*Gamma2*DeltaS
Hch0 = numerator / denominator
Hch_mag = abs(Hch0)
print(f"H_ch(0) = {Hch_mag:.5f}")

# === COMPUTE Rch ===
# For matched terminations (Gamma1=Gamma2=0), simplifies to 2*R0*(1/|Hch|-1)
Rch = R0 * (1/Hch_mag - 2) if (Gamma1 or Gamma2) else 2*R0*(1/Hch_mag - 1)
print(f"Estimated R_ch = {Rch:.3f} Ω")

