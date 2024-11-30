import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 8.314  # Universal gas constant in J/(mol·K)
T_ref = 298.15  # Reference temperature in Kelvin
T_comb = 1500  # Approximate temperature for combustion in Kelvin (example)

# NASA polynomials coefficients
coeffs_CH4 = [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01]
coeffs_O2 = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]
coeffs_N2 = [2.95257626E+00, 1.39690000E-03, -4.92631603E-07, 7.86010367E-11, -4.60755321E-15, -9.23948688E+02, 5.87188762E+00]
coeffs_CO2 = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
coeffs_H2O = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]

# Molecular weights (kg/kmol)
MW_CH4 = 16
MW_CO2 = 44
MW_H2O = 18

# Equivalence ratio range
phi_range = np.arange(0.5, 2.01, 0.01)

# Specific enthalpy calculation function
def specific_enthalpy(T, coeffs):
    h = R * T * (
        coeffs[0]
        + coeffs[1] * T / 2
        + coeffs[2] * T**2 / 3
        + coeffs[3] * T**3 / 4
        + coeffs[4] * T**4 / 5
        + coeffs[5] / T
    )
    return h

# Calculate specific enthalpies at reference and combustion temperatures
H_CH4_ref = specific_enthalpy(T_ref, coeffs_CH4)
H_O2_ref = specific_enthalpy(T_ref, coeffs_O2)
H_N2_ref = specific_enthalpy(T_ref, coeffs_N2)
H_CO2_ref = specific_enthalpy(T_ref, coeffs_CO2)
H_H2O_ref = specific_enthalpy(T_ref, coeffs_H2O)

H_CH4_comb = specific_enthalpy(T_comb, coeffs_CH4)
H_O2_comb = specific_enthalpy(T_comb, coeffs_O2)
H_N2_comb = specific_enthalpy(T_comb, coeffs_N2)
H_CO2_comb = specific_enthalpy(T_comb, coeffs_CO2)
H_H2O_comb = specific_enthalpy(T_comb, coeffs_H2O)

# Initialize arrays for combustion enthalpy per kmol and per kg of fuel
enthalpy_per_kmol = np.zeros(len(phi_range))
enthalpy_per_kg = np.zeros(len(phi_range))

# Calculate enthalpy change for each equivalence ratio
for i, phi in enumerate(phi_range):
    # Reactant moles
    n_CH4_react = 1  # Assume 1 kmol of CH4
    n_O2_react = (2 / phi) * n_CH4_react
    n_N2_react = (79 / 21) * n_O2_react

    # Product moles
    if phi <= 1:  # Lean mixture
        n_CO2_prod = 1 * n_CH4_react
        n_H2O_prod = 2 * n_CH4_react
        n_O2_prod = n_O2_react - 2 * n_CH4_react
        n_N2_prod = n_N2_react
        n_CH4_prod = 0
    else:  # Rich mixture
        n_CO2_prod = 1 * n_CH4_react
        n_H2O_prod = 2 * n_CH4_react
        n_O2_prod = 0
        n_N2_prod = n_N2_react
        n_CH4_prod = phi - 1

    # Enthalpy of reactants and products
    H_react = (
        n_CH4_react * H_CH4_comb
        + n_O2_react * H_O2_comb
        + n_N2_react * H_N2_comb
    )
    H_prod = (
        n_CO2_prod * H_CO2_comb
        + n_H2O_prod * H_H2O_comb
        + n_O2_prod * H_O2_comb
        + n_N2_prod * H_N2_comb
        + n_CH4_prod * H_CH4_comb
    )
    
    # Enthalpy change (J/kmol CH4)
    delta_H = H_prod - H_react

    # Store results
    enthalpy_per_kmol[i] = delta_H  # KJ/kmol of CH4
    enthalpy_per_kg[i] = delta_H / MW_CH4  # KJ/kg of CH4

# Convert enthalpy to MJ for better readability
enthalpy_per_kmol_KJ = enthalpy_per_kmol / 1e3
enthalpy_per_kg_KJ = enthalpy_per_kg / 1e3

# Plot results
plt.figure(figsize=(12, 6))
plt.plot(phi_range, enthalpy_per_kmol_KJ, label='Enthalpy of CH4 (KJ/Kmol)', color='blue')
plt.plot(phi_range, enthalpy_per_kg_KJ, label='Enthalpy of CH4 (KJ/Kg)', color='red')
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')

# Add labels, title, and legend
plt.title('Combustion Enthalpy vs Equivalence Ratio', fontsize=16)
plt.xlabel('Equivalence Ratio (Phi)', fontsize=14)
plt.ylabel('Combustion Enthalpy (MJ)', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
