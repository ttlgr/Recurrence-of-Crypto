from Crypto.Util.number import getPrime, inverse, bytes_to_long
from secret import flag
from random import randint
from math import floor

def hint(p, length, upper):
    gift = {}
    pstr = bin(p)[2:].zfill(upper + 1)
    for i in range(length):
        idx = randint(1, upper)
        if idx not in gift:
            gift[idx] = pstr[idx]
    return gift

p = getPrime(512)
q = getPrime(512)
n = p * q
e = getPrime(32)
d = inverse(e, (p - 1) * (q - 1))
c = pow(bytes_to_long(flag), e, n)
hint1 = hint(p, floor(0.42 * p.bit_length()), 511)
hint2 = hint(q, floor(0.42 * q.bit_length()), 511)
hint3 = hint(d, floor(0.42 * d.bit_length()), 1023)
print("n =", n)
print("e =", e)
print("c =", c)
print("hint1 =", hint1)
print("hint2 =", hint2)
print("hint3 =", hint3)