import numpy as np
import matplotlib.pyplot as plt

# Molecular weights (kg/kmol)
M_CH4 = 16
M_O2 = 32
M_N2 = 28

# Air composition (molar ratio)
x_O2_air = 0.21
x_N2_air = 0.79

# Equivalence ratio range
phi = np.arange(0.5, 2.01, 0.01)

# Initialize arrays for mole fractions and mass fractions
chi_CH4 = np.zeros(len(phi))
chi_O2 = np.zeros(len(phi))
chi_N2 = np.zeros(len(phi))
mass_fraction_CH4 = []
mass_fraction_O2 = []
mass_fraction_N2 = []

# Calculate mole and mass fractions for each phi
for i in range(len(phi)):
    if phi[i] <= 1:  # Lean mixture (excess oxygen)
        n_CH4 = 1
        n_O2 = 2 / phi[i]
        n_N2 = n_O2 * (x_N2_air / x_O2_air)
    else:  # Rich mixture (excess fuel)
        n_CH4 = phi[i]
        n_O2 = 2
        n_N2 = n_O2 * (x_N2_air / x_O2_air)
    
    # Total moles of reactants
    total_moles = n_CH4 + n_O2 + n_N2

    # Mole fractions
    chi_CH4[i] = n_CH4 / total_moles
    chi_O2[i] = n_O2 / total_moles
    chi_N2[i] = n_N2 / total_moles

    # Mass fractions
    m_CH4 = chi_CH4[i] * M_CH4
    m_O2 = chi_O2[i] * M_O2
    m_N2 = chi_N2[i] * M_N2
    total_mass = m_CH4 + m_O2 + m_N2

    mass_fraction_CH4.append(m_CH4 / total_mass)
    mass_fraction_O2.append(m_O2 / total_mass)
    mass_fraction_N2.append(m_N2 / total_mass)

# Convert lists to numpy arrays for plotting
mass_fraction_CH4 = np.array(mass_fraction_CH4)
mass_fraction_O2 = np.array(mass_fraction_O2)
mass_fraction_N2 = np.array(mass_fraction_N2)

# Plotting molar fractions
plt.figure(figsize=(12, 8))
plt.plot(phi, chi_CH4, label='CH4 (molar fraction)', color='blue', linewidth=2)
plt.plot(phi, chi_O2, label='O2 (molar fraction)', color='red', linewidth=2)
plt.plot(phi, chi_N2, label='N2 (molar fraction)', color='green', linewidth=2)
plt.title('Reactants Molar Fraction Profiles', fontsize=16)
plt.xlabel('Equivalence Ratio (ϕ)', fontsize=14)
plt.ylabel('Molar Fraction', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()

# Plotting mass fractions
plt.figure(figsize=(12, 8))
plt.plot(phi, mass_fraction_CH4, label='CH4 (mass fraction)', color='blue', linewidth=2)
plt.plot(phi, mass_fraction_O2, label='O2 (mass fraction)', color='red', linewidth=2)
plt.plot(phi, mass_fraction_N2, label='N2 (mass fraction)', color='green', linewidth=2)
plt.title('Reactants Mass Fraction Profiles', fontsize=16)
plt.xlabel('Equivalence Ratio (ϕ)', fontsize=14)
plt.ylabel('Mass Fraction', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()




