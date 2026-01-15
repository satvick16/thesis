import skrf as rf
import matplotlib.pyplot as plt

# Load the s2p file
ntw = rf.Network("channel_two_port.s2p")

# Frequency axis in GHz
freq = ntw.f / 1e9

# Plot |S-parameters| in dB
plt.figure()
plt.plot(freq, ntw.s_db[:, 0, 0], label="S11")
plt.plot(freq, ntw.s_db[:, 1, 0], label="S21")
plt.plot(freq, ntw.s_db[:, 0, 1], label="S12")
plt.plot(freq, ntw.s_db[:, 1, 1], label="S22")

plt.xlabel("Frequency (GHz)")
plt.ylabel("Magnitude (dB)")
plt.title("S-Parameters")
plt.legend()
plt.grid(True)
plt.show()
