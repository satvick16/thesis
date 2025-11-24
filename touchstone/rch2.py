import numpy as np
import skrf as rf

def extract_channel_resistance(s4p_path):
    # Load the S4P file
    ntwk = rf.Network(s4p_path)
    z0 = ntwk.z0[0,0]  # reference impedance, usually 50Ω

    # Convert to Z-parameters
    z_params = ntwk.z
    freqs = ntwk.f

    # Pick DC (or the lowest frequency)
    Z = z_params[0]

    # Extract relevant elements
    Z11, Z12, Z21, Z22 = Z[0,0], Z[0,1], Z[1,0], Z[1,1]

    # Compute channel equivalent series resistance
    Rch = np.real(Z11 + Z22 - Z12 - Z21)

    print(f"Reference Impedance (Z0): {z0} Ω")
    print(f"Lowest frequency point: {freqs[0]/1e9:.3f} GHz")
    print(f"Channel Resistance Rch ≈ {Rch:.4f} Ω")

    return Rch

# Example usage
if __name__ == "__main__":
    path = "C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p"  # <-- replace with your file path
    extract_channel_resistance(path)
