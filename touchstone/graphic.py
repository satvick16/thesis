import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

# -----------------------------
# Parameters
# -----------------------------
num_symbols = 500
sps = 64                  # samples per symbol
snr_db = 20
np.random.seed(1)

# -----------------------------
# Modulation levels (normalized power)
# -----------------------------
levels = {
    "NRZ": np.array([-1, 1]),
    "PAM3": np.array([-1, 0, 1]),
    "PAM4": np.array([-3, -1, 1, 3])
}

def normalize_power(x):
    return x / np.sqrt(np.mean(x**2))

# -----------------------------
# Channel (simple low-pass)
# -----------------------------
def lowpass(signal, cutoff=0.25):
    b, a = butter(5, cutoff)
    return lfilter(b, a, signal)

# -----------------------------
# Signal generation
# -----------------------------
def generate_signal(levels):
    symbols = np.random.choice(levels, num_symbols)
    symbols = normalize_power(symbols)

    # Rectangular pulse shaping
    tx = np.zeros(num_symbols * sps)
    tx[::sps] = symbols

    # Channel filtering
    tx = lowpass(tx)

    # Noise
    signal_power = np.mean(tx**2)
    noise_power = signal_power / (10**(snr_db / 10))
    noise = np.random.normal(0, np.sqrt(noise_power), len(tx))

    return tx + noise

signals = {k: generate_signal(v) for k, v in levels.items()}

# -----------------------------
# Time-domain plots
# -----------------------------
plt.figure(figsize=(12, 7))

for i, (name, sig) in enumerate(signals.items(), 1):
    plt.subplot(3, 1, i)
    plt.plot(sig[:12 * sps])
    plt.title(f"{name} Time-Domain Waveform")
    plt.ylabel("Amplitude")
    plt.grid(True)

plt.xlabel("Samples")
plt.tight_layout()
plt.show()

# -----------------------------
# Eye diagram function
# -----------------------------
def eye_diagram(signal, sps, title):
    span = 2 * sps
    for i in range(sps, len(signal) - span, sps):
        plt.plot(signal[i:i+span], color="blue", alpha=0.12)
    plt.title(title)
    plt.xlabel("Time")
    plt.grid(True)

# -----------------------------
# Eye diagrams
# -----------------------------
plt.figure(figsize=(12, 4))

for i, (name, sig) in enumerate(signals.items(), 1):
    plt.subplot(1, 3, i)
    eye_diagram(sig, sps, f"{name} Eye Diagram")

plt.tight_layout()
plt.show()
