# generate 4 sequences of 30 bits each: H, L, !H, !L
# permitted values are: (H = 0, L = 0), (H = 0, L = 1), (H = 1, L = 0)
# not permitted: (H = 1, L = 1)

import random
def generate_bit_sequences(length=30):
    H = ""
    L = ""
    not_H = ""
    not_L = ""

    for _ in range(length):
        bit_H = random.randint(0, 1)
        if bit_H == 1:
            bit_L = 0
        else:
            bit_L = random.randint(0, 1)

        H += str(bit_H)
        L += str(bit_L)
        not_H += str(1 - bit_H)
        not_L += str(1 - bit_L)

    return H, L, not_H, not_L
H, L, not_H, not_L = generate_bit_sequences()
print("H:     ", H)
print("L:     ", L)
print("!H:    ", not_H)
print("!L:    ", not_L)
print("Length:", len(H))
print("Length:", len(L))
print("Length:", len(not_H))
print("Length:", len(not_L))
