# generate 2 random bit sequences called H and L
# the permitted combinations of values at any given time step are:
# H=0, L=0
# H=0, L=1
# H=1, L=0
# they cannot both be 1 at the same time

import random

def generate_bit_sequences(length):
    H = ""
    L = ""
    notH = ""
    notL = ""
    for _ in range(length):
        choice = random.choice([(0,0), (0,1), (1,0)])
        H += str(choice[0])
        L += str(choice[1])
        notH += str(1 - choice[0])
        notL += str(1 - choice[1])
    return H, L, notH, notL

H, L, notH, notL = generate_bit_sequences(20)
print("H:", H)
print("L:", L)
print("notH:", notH)
print("notL:", notL)
