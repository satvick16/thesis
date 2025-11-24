import skrf as rf
import numpy as np

# Input and output files
input_file = 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p'
output_file = 'channel.s4p'

# Desired reference impedance
new_z0 = 40  # Ohms

# Load the 4-port S4P
ntwk = rf.Network(input_file)

# Original impedance (assume uniform)
z0_old = ntwk.z0  # could be scalar or array of length 4
if np.isscalar(z0_old):
    z0_old = np.array([z0_old]*ntwk.number_of_ports)

# Make new z0 array
if np.isscalar(new_z0):
    z0_new = np.array([new_z0]*ntwk.number_of_ports)
else:
    z0_new = np.array(new_z0)

# Renormalize network
ntwk.renormalize(z0_new)

# Save new S4P
ntwk.write_touchstone(output_file)

print(f"Renormalized S4P saved as {output_file} with Z0 = {new_z0} Ω")
