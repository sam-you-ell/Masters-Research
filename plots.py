import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# Fvalues = pd.read_excel('FermionCircuit20.xlsx')
# time = np.array(Fvalues['Time'])
# Ent1 = np.array(Fvalues['Ent1'])
# Ent2 = np.array(Fvalues['Ent2'])
# Ent3 = np.array(Fvalues['Ent3'])
# Ent4 = np.array(Fvalues['Ent4'])
# Ent5 = np.array(Fvalues['Ent5'])

# fourtvals = pd.read_excel('FermionCircuit.xlsx')
# timefourt = np.array(fourtvals['Time'])
# Ave40 = np.array(fourtvals['Ave40'])

# sixtvals = pd.read_excel('FermionCircuit60.xlsx')
# timesixt = np.array(sixtvals['Time'])
# Ave60 = np.array(sixtvals['Ave60'])

# Averaged_Entropy = (Ent1 + Ent2 + Ent3 + Ent4 + Ent5) / 5

# plt.plot(time[:4000], Averaged_Entropy[:4000], label="L=20")
# plt.plot(timefourt[:4000], Ave40[:4000], label="L=40")
# plt.plot(timesixt[:4000], Ave60[:4000], label="L=60")

bigvals = pd.read_excel('FermionCircuit100.xlsx')
time = np.array(bigvals['Time'])
fermS = np.array(bigvals['Ent1'])

Stabvals = pd.read_excel('StabilizerCircuits100.xlsx')
stabS = np.array(Stabvals['AveStab'])


plt.plot(time, fermS, label='Free Fermion')
plt.plot(time, stabS, label='Clifford')

plt.xlabel(r"  t", rotation=0, loc='center')
plt.ylabel(r" $S_{A}/S_{\infty}$       ", rotation=0)

for pos in ['right', 'top']:
    plt.gca().spines[pos].set_visible(False)

plt.xticks(np.arange(min(time), max(time)+10, 5000))
plt.legend()
# plt.title(
#     "Entanglement Entropy under Free Fermion Dynamics")


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
