import random

def prbs7(length):
    """
    Generate a PRBS-7 sequence of specified length (in bits).
    """
    sequence = []
    for i in range(length):
        sequence.append(random.choice([0, 1]))

    return sequence

def write_ltspice_pwl(filename, sequence, timestep_ps=20.0):
    """
    Write PRBS sequence as PWL format for LTspice.
    Each step is timestep_ps long.
    """
    with open(filename, "w") as f:
        prev_bit = sequence[0]
        time_ps = 0.0
        for bit in sequence:
            if bit:
                bit = 0.75
            if time_ps > 0:
                f.write(f"{time_ps-1:.3f}p {prev_bit}\n")
            f.write(f"{time_ps:.3f}p {bit}\n")
            time_ps += timestep_ps
            prev_bit = bit

        # Repeat last point to hold final value
        if sequence[-1]:
            f.write(f"{time_ps:.3f}p 0.75\n")
        else:
            f.write(f"{time_ps:.3f}p 0\n")

# Parameters
num_bits = 127  # Full PRBS-7 period
timestep = 20.0  # in ps
output_file = "prbs7.txt"

# Generate and export
seq = prbs7(num_bits)
write_ltspice_pwl(output_file, seq, timestep)

print(f"PRBS-7 sequence written to {output_file}")
