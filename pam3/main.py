import numpy as np

# === PAM-3 Data Settings ===
bitrate = 50e9  # 50 Gbps
symbol_time = 1 / bitrate  # 20 ps
num_symbols = 320  # more symbols for better test
levels = [0, 1, 2]  # PAM-3 levels: maps to 0V, 0.375V, 0.75V
data = np.random.choice(levels, size=num_symbols)

# === Time Axis ===
time_step = 1e-12  # 1 ps resolution
t = np.arange(0, num_symbols * symbol_time, time_step)

# === sel0 / sel1 One-Hot Control Signals ===
sel0 = np.zeros_like(t)
sel1 = np.zeros_like(t)

for i, symbol in enumerate(data):
    start = int(i * symbol_time / time_step)
    end = int((i + 1) * symbol_time / time_step)
    if symbol == 0:
        sel0[start:end] = 1
    elif symbol == 1:
        sel1[start:end] = 1
    # symbol == 2 → both stay 0 → output = 0.75V

# === Save to File ===
def save_pwl(filename, t, v):
    with open(filename, 'w') as f:
        for ti, vi in zip(t, v):
            f.write(f"{ti:.12e}\t{vi:.1f}\n")

save_pwl("sel0.txt", t, sel0)
save_pwl("sel1.txt", t, sel1)
print("Saved: sel0.txt and sel1.txt")
