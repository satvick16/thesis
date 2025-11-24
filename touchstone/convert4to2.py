import skrf as rf
import numpy as np
import sys

def convert_differential_to_single_ended(input_file, output_file):
    """
    Convert a differential 4-port S-parameter file (S4P)
    into a 2-port differential (S2P) Touchstone file.

    Expected port mapping:
        Port 1 -> TX+
        Port 2 -> TX-
        Port 3 -> RX+
        Port 4 -> RX-
    """

    # Load the 4-port S-parameter file
    s4p = rf.Network(input_file)
    if s4p.number_of_ports != 4:
        raise ValueError("Input file must have 4 ports (differential pair system).")

    # Convert single-ended (SE) -> mixed-mode (MM)
    # New API in scikit-rf >= 0.25
    mm = rf.network.se2gmm(s4p)

    # Extract the differential-to-differential (Sdd) 2x2 block
    # Mixed-mode order: [diff1, diff2, comm1, comm2]
    Sdd11 = mm.s[:, 0, 0]
    Sdd12 = mm.s[:, 0, 1]
    Sdd21 = mm.s[:, 1, 0]
    Sdd22 = mm.s[:, 1, 1]

    # Construct a 2-port differential S-matrix
    s2p_matrix = np.zeros((len(Sdd11), 2, 2), dtype=complex)
    s2p_matrix[:, 0, 0] = Sdd11
    s2p_matrix[:, 0, 1] = Sdd12
    s2p_matrix[:, 1, 0] = Sdd21
    s2p_matrix[:, 1, 1] = Sdd22

    # Create the new 2-port network
    s2p = rf.Network(frequency=s4p.frequency, s=s2p_matrix, z0=s4p.z0[0, 0])

    # Write the resulting S2P file
    s2p.write_touchstone(output_file)
    print(f"✅ Successfully converted {input_file} → {output_file}.s2p")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert4to2.py input.s4p output_basename")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_differential_to_single_ended(input_file, output_file)
