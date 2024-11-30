import numpy as np
import matplotlib.pyplot as plt

# Define the equivalence ratio range
phi = np.arange(0.5, 2.01, 0.01)

# Initialize mole fractions
chi_CO2 = np.zeros_like(phi)
chi_H2O = np.zeros_like(phi)
chi_O2 = np.zeros_like(phi)
chi_N2 = np.zeros_like(phi)

# Compute mole fractions based on phi
for i in range(len(phi)):
    if phi[i] < 1:  # Lean mixture
        chi_CO2[i] = 1 / (10.52 - 2 * phi[i])
        chi_H2O[i] = chi_CO2[i]
        chi_O2[i] = (2 - 2 * phi[i]) / (10.52 - 2 * phi[i])
        chi_N2[i] = 7.52 / (10.52 - 2 * phi[i])
    else:  # Rich mixture
        chi_CO2[i] = 1 / (7.52 + phi[i])
        chi_H2O[i] = chi_CO2[i]
        chi_O2[i] = 0
        chi_N2[i] = 7.52 / (7.52 + phi[i])

# Plot mole fractions
plt.figure(figsize=(10, 8))
plt.plot(phi, chi_H2O, label='H2O (molar fract)', color='blue')
plt.plot(phi, chi_N2, label='N2 (molar fract)', color='green')
plt.plot(phi, chi_O2, label='O2 (molar fract)', color='orange')  # O2 added for lean mixtures

# Add title, labels, and legend
plt.title('Combustion Products Molar Fraction Profiles', fontsize=14)
plt.xlabel('Equivalence Ratio (Phi)', fontsize=12)
plt.ylabel('Molar Fraction', fontsize=12)
plt.legend(fontsize=10)
plt.grid()

# Show plot
plt.tight_layout()
plt.show()
