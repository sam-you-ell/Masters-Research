import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
# Fvalues = pd.read_excel('FermionCircuit20.xlsx')
# time = np.array(Fvalues['Time'])
# Ent1 = np.array(Fvalues['Ent1'])
# Ent2 = np.array(Fvalues['Ent2'])
# Ent3 = np.array(Fvalues['Ent3'])
# Ent4 = np.array(Fvalues['Ent4'])
# Ent5 = np.array(Fvalues['Ent5'])

# U2 = pd.read_excel('U2onFockState120.xlsx')
# timefourt = np.array(U2['Time'])
# Su2 = np.array(U2['Ave'])

# U3 = pd.read_excel('U3onFockState120.xlsx')
# Time = np.array(U3['Time'])
# Su3 = np.array(U3['Ave'])


Fermi20 = pd.read_excel('U2onFockState20.xlsx')
time = np.array(Fermi20['Time'])
Sf20 = np.array(Fermi20['Ave'])

Fermi40 = pd.read_excel('U2onFockState40.xlsx')
Time = np.array(Fermi40['Time'])
Sf40 = np.array(Fermi40['Ave'])

Fermi60 = pd.read_excel('U2onFockState60.xlsx')
Sf60 = np.array(Fermi60['Ave'])

Fermi80 = pd.read_excel('U2onFockState80.xlsx')
Time120 = np.array(Fermi80['Time'])
Sf80 = np.array(Fermi80['Ave'])

# Fock120 = pd.read_excel('U2onFockState120.xlsx')
# TimeFock120 = np.array(Fock120['Time'])
# SFock120 = np.array(Fock120['Ave'])

# # FockU3120 = pd.read_excel('U3onFockState120.xlsx')
# # TimeFock120 = np.array(FockU3120['Time'])
# # SFockU3120 = np.array(FockU3120['Ave'])

# Fermi120 = pd.read_excel('U2onFermiSea120.xlsx')
# TimeFerm120 = np.array(Fermi120['Time'])
# fermS120 = np.array(Fermi120['Ave'])


# U3Fermi120 = pd.read_excel('U3onFermiSea80.xlsx')
# Time1 = np.array(U3Fermi120['Time'])
# fermS1 = np.array(U3Fermi120['Ave'])

# U2Fermi80 = pd.read_excel('U2onFermiSea80.xlsx')
# Time2 = np.array(U2Fermi80['Time'])
# fermS2 = np.array(U2Fermi80['Ave'])


# plt.plot(time[:4000], Averaged_Entropy[:4000], label="L=20")
# plt.plot(timefourt[:4000], Ave40[:4000], label="L=40")
# plt.plot(timesixt[:4000], Ave60[:4000], label="L=60")

# bigvals = pd.read_excel('FermionCircuit100.xlsx')
# time = np.array(bigvals['Time'])
# fermS = np.array(bigvals['Ent1'])

Stabvals = pd.read_excel('StabilizerCircuits.xlsx')
stabS = np.array(Stabvals['Ave'])


plt.plot(time / 20**2, Sf20, label='L = 20')
plt.plot(Time / 40**2, Sf40, label='L = 40')
plt.plot(Time / 60**2, Sf60, label='L = 60')
plt.plot(Time / 80**2, Sf80, label='L = 80')
# plt.plot(Time120 / 120**2, Sf120, label='L = 80')

# plt.plot(TimeFock120, SFock120, label=r'Free-Fermion')
# plt.plot(TimeFerm120, fermS120,  color='blue', label=r'Super Clifford')

plt.xlabel(r"  $t/L^2$", rotation=0, loc='center')
plt.ylabel(r" $S_{A}/S_{\infty}$       ", rotation=0)

plt.xlim((0.0, 0.3))
# plt.ylim((None, 0.7))
# for pos in ['right', 'top']:
#     plt.gca().spines[pos].set_visible(False)

# plt.xticks(np.arange(min(TimeFerm120), max(TimeFerm120)+10, 5000))
plt.legend(loc='lower right')
plt.title(
    "Entanglement Entropy vs Time for Variable System Size")

plt.savefig("varysystem.pdf")

plt.show()


# Svalues = pd.read_excel('StabilizerCircuits.xlsx')
# time = np.array(Svalues['Time'])
# Stab_Entropy1 = np.array(Svalues['Ent1'])
# Stab_Entropy2 = np.array(Svalues['Ent2'])
# Stab_Entropy3 = np.array(Svalues['Ent3'])
# Stab_Entropy4 = np.array(Svalues['Ent4'])
# Stab_Entropy5 = np.array(Svalues['Ent5'])
# Stab_Entropy6 = np.array(Svalues['Ent6'])
# Stab_Entropy7 = np.array(Svalues['Ent7'])
# Stab_Entropy8 = np.array(Svalues['Ent8'])
# Stab_Entropy9 = np.array(Svalues['Ent9'])
# Stab_Entropy10 = np.array(Svalues['Ent10'])

# Averaged_Entropy_Stabilizers = (
#     Stab_Entropy1 + Stab_Entropy2 + Stab_Entropy3 + Stab_Entropy4 + Stab_Entropy5 + Stab_Entropy6 + Stab_Entropy7 + Stab_Entropy8 + Stab_Entropy9 + Stab_Entropy10) / 10
# Averaged_Entropy_Stabilizers = Averaged_Entropy_Stabilizers * \
#     math.log(2**60)/math.log(2)
