# Calculate BHP using Cullender-Smith method
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
import numpy as np

# %%
gSG = 0.71
d = 2.259
rel_roughness = 0.0006
L = 10000
theta = 0
phf = 800
Thf = 150
Twf = 200
qsc = 2


# %%
def zFactor(p, t, gSG):
    ppc = 677 + 15 * gSG - 37.5 * gSG ** 2
    tpc = 168 + 325 * gSG - 12.5 * gSG ** 2

    ppr = p / ppc
    tpr = t / tpc

    F = 0.3106 - 0.49 * tpr + 0.1824 * (tpr ** 2)
    E = 9 * (tpr - 1)
    D = 10 ** F
    C = 0.132 - 0.32 * math.log10(tpr)
    B = (0.62 - 0.23 * tpr) * ppr + ((0.066 / (tpr - 0.86)) - 0.037) * (
            ppr ** 2) + ((0.32 * (ppr ** 6)) / (10 ** E))
    A = 1.39 * (tpr - 0.92) ** 0.5 - 0.36 * tpr - 0.101

    z = A + ((1 - A) / (math.exp(B))) + C * (ppr ** D)
    return z


# %%
f = (1 / (1.74 - 2 * math.log10(2 * rel_roughness))) ** 2

# %%
depth = [x for x in range(0, L + 1, int(L / 2))]
T_hf_R = Thf + ((Twf - Thf) / L) * depth[0] + 460
T_mf_R = Thf + ((Twf - Thf) / L) * depth[1] + 460
T_wf_R = Thf + ((Twf - Thf) / L) * depth[-1] + 460
Z_hf = zFactor(phf, T_hf_R, gSG)
pzt_hf = phf / (Z_hf * T_hf_R)
I_hf = pzt_hf / (0.001 * math.cos(math.radians(theta)) * (pzt_hf ** 2) + (0.6666 * f * (qsc ** 2)) / (d ** 5))


# %%

def pCalc(piz_guess,t,I_prev,p_prev):
    p = piz_guess[0]
    I = piz_guess[1]
    z = piz_guess[2]

    ppc = 677 + 15 * gSG - 37.5 * gSG ** 2
    tpc = 168 + 325 * gSG - 12.5 * gSG ** 2

    tpr = t / tpc

    f0 = z - ((1.39 * (tpr - 0.92) ** 0.5 - 0.36 * tpr - 0.101) + (
                (1 - (1.39 * (tpr - 0.92) ** 0.5 - 0.36 * tpr - 0.101)) /
                (math.exp((0.62 - 0.23 * tpr) * (p / ppc) + ((0.066 / (tpr - 0.86)) - 0.037) * ((p / ppc) ** 2) + (
                            (0.32 * ((p / ppc) ** 6)) / (10 ** (9 * (tpr - 1))))))
                + (0.132 - 0.32 * math.log10(tpr)) * (
                            (p / ppc) ** (10 ** (0.3106 - 0.49 * tpr + 0.1824 * (tpr ** 2))))))
    f1 = I - ((p / (z * t)) / (
                0.001 * math.cos(math.radians(theta)) * ((p / (z * t)) ** 2) + (0.6666 * f * (qsc ** 2)) / (
                    d ** 5)))
    f2 = p - p_prev - (18.75 * gSG * L / (I_prev + I))
    return np.array([f0, f1, f2])


piz_guess0 = np.array([phf, I_hf, Z_hf])
piz = optimize.fsolve(pCalc, piz_guess0,(T_mf_R,I_hf,phf))
pmf = piz[0]
I_mf = piz[1]
Z_mf = piz[2]
print(piz, pCalc(piz,T_mf_R,I_mf,phf))

piz_guess1 = np.array([pmf, I_mf, Z_mf])
piz1 = optimize.fsolve(pCalc, piz_guess1,(T_wf_R,I_mf,pmf))
pwf = piz1[0]
I_wf = piz1[1]
Z_wf = piz1[2]
print(piz1, pCalc(piz1,T_wf_R,I_wf,pmf))


# %%




