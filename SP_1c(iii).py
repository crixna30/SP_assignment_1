import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Constants
R = 8.314  # Universal gas constant [J/mol-K]

# Temperature range and equivalence ratio range
phi_range = np.arange(0.5, 2.01, 0.01)

# Initial temperature
T0 = 298.15  # [K]

# coefficients for NASA polynomial (valid for T >= 1000 K)
coeffs_CH4 = [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01]
coeffs_O2 = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]
coeffs_N2 = [2.95257626E+00, 1.39690000E-03, -4.92631603E-07, 7.86010367E-11, -4.60755321E-15, -9.23948688E+02, 5.87188762E+00]
coeffs_CO2 = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
coeffs_H2O = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]




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


def thermal_balance(T, H_reactants, phi, coeffs_CH4, coeffs_CO2, coeffs_H2O, coeffs_O2, coeffs_N2):
    
    if phi <= 1:
        n_CO2 = 1
        n_H2O = 2
        n_O2 = 2 / phi - 2
        n_N2 = (79 / 21) * (2 / phi)
        n_CH4 = 0
    else:
        n_CO2 = 1
        n_H2O = 2
        n_O2 = 0
        n_N2 = (79 / 21) * 2
        n_CH4 = phi - 1

    h_CO2 = specific_enthalpy(T, coeffs_CO2)
    h_H2O = specific_enthalpy(T, coeffs_H2O)
    h_O2 = specific_enthalpy(T, coeffs_O2)
    h_N2 = specific_enthalpy(T, coeffs_N2)
    h_CH4 = specific_enthalpy(T, coeffs_CH4)
    H_products = n_CO2 * h_CO2 + n_H2O * h_H2O + n_O2 * h_O2 + n_N2 * h_N2 + n_CH4 * h_CH4

    return H_reactants - H_products


# Preallocate result storage
T_ad = np.zeros(len(phi_range))

# Loop over equivalence ratios
for i, phi in enumerate(phi_range):
    n_CH4 = 1  # 1 mol of CH4
    n_O2 = (2 / phi) * n_CH4
    n_N2 = (79 / 21) * n_O2

    h_CH4 = specific_enthalpy(T0, coeffs_CH4)
    h_O2 = specific_enthalpy(T0, coeffs_O2)
    h_N2 = specific_enthalpy(T0, coeffs_N2)
    H_reactants = n_CH4 * h_CH4 + n_O2 * h_O2 + n_N2 * h_N2

    T_guess = 2000  # Initial guess
    T_ad[i] = fsolve(thermal_balance, T_guess, args=(H_reactants, phi, coeffs_CH4, coeffs_CO2, coeffs_H2O, coeffs_O2, coeffs_N2))[0]

fig, ax = plt.subplots(figsize=(12, 6))

# Plot the adiabatic flame temperature as a red line
ax.plot(phi_range, T_ad, color='red', label='Adiabatic Flame Temperature (K)')
ax.set_xlabel('Equivalence Ratio (Phi)', fontsize=14)
ax.set_ylabel('Adiabatic Flame Temperature (K)', fontsize=14, color='red')
ax.set_title('Adiabatic Flame Temperature vs Equivalence Ratio', fontsize=16)

# Grid and legend
ax.grid(True, linestyle='--', alpha=0.7)
ax.legend(loc='upper right', fontsize=12)

# Final layout adjustments
fig.tight_layout()
plt.show()



