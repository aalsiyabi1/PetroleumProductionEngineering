# Calculate BHP using Cullender-Smith method
import math
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize

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

#%%
depth = [x for x in range(0, L + 1, int(L / 2))]
T_hf_R = Thf + ((Twf - Thf) / L) * depth[0] + 460
T_mf_R = Thf + ((Twf - Thf) / L) * depth[1] + 460
T_wf_R = Thf + ((Twf - Thf) / L) * depth[-1] + 460
Z_hf = zFactor(phf,T_hf_R,gSG)
pzt_hf = phf/(Z_hf*T_hf_R)
I_hf = pzt_hf / (0.001 *math.cos(math.radians(theta)) * (pzt_hf**2) + (0.6666*f*(qsc**2))/(d**5))

# %%



#
# p_guess = phf+50
# z_guess = zFactor(p_guess,T_mf_R,gSG)
#
# obj_fun = lambda Y: (Y - phf - (18.75*gSG*L/(I_hf+((Y/(z_guess*T_mf_R)) /
#                                                    (0.001 *math.cos(math.radians(theta)) * ((Y/(z_guess*T_mf_R))**2) +
#                                                     (0.6666*f*(qsc**2))/(d**5))))))
# Y_soln = optimize.fsolve(obj_fun, 0)
# fY_soln = obj_fun(Y_soln)
# # print(Y_soln,fY_soln)
#
# pmf = Y_soln
#




#%%
'''
fY = lambda Y: (P[1] - P[0] - (18.75*gSG*L/(I[1]+I[0]))
Y_soln = optimize.fsolve(fY, 0)
fY_soln = fY(Y_soln)
print(Y_soln,fY_soln)
'''

