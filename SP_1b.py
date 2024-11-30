import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 8.314  # Universal gas constant (J/mol·K)

# Temperature (in K)
T = 1500  # Example temperature, can be adjusted

# NASA polynomial coefficients (valid for T >= 1000 K)
coeffs_CH4 = [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01]
coeffs_O2 = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]
coeffs_N2 = [2.95257626E+00, 1.39690000E-03, -4.92631603E-07, 7.86010367E-11, -4.60755321E-15, -9.23948688E+02, 5.87188762E+00]
coeffs_CO2 = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00]
coeffs_H2O = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]

# Molecular weights (kg/kmol)
M_CH4 = 16
M_O2 = 32
M_N2 = 28
M_CO2 = 44
M_H2O = 18

# Function to calculate specific enthalpy (kJ/kmol) using NASA polynomials
def absolute_enthalpy(T, coeffs):
    h = R * T * (
        coeffs[0]
        + coeffs[1] * T / 2
        + coeffs[2] * T**2 / 3
        + coeffs[3] * T**3 / 4
        + coeffs[4] * T**4 / 5
        + coeffs[5] / T
    ) / 1000  # Convert from J/kmol to kJ/kmol
    return h

# Calculate enthalpies of individual species at the given temperature
H_CH4 = absolute_enthalpy(T, coeffs_CH4)
H_O2 = absolute_enthalpy(T, coeffs_O2)
H_N2 = absolute_enthalpy(T, coeffs_N2)
H_CO2 = absolute_enthalpy(T, coeffs_CO2)
H_H2O = absolute_enthalpy(T, coeffs_H2O)

# Equivalence ratio range
phi = np.arange(0.50, 2.01, 0.01)

# Initialize arrays for mole fractions and enthalpy profiles
chi_CO2 = np.zeros(len(phi))
chi_H2O = np.zeros(len(phi))
chi_O2 = np.zeros(len(phi))
chi_N2 = np.zeros(len(phi))
H_molar = []  # Enthalpy per kmol
H_mass = []   # Enthalpy per kg

# Calculate mole fractions and enthalpy for each phi
for i in range(len(phi)):
    if phi[i] < 1:  # Lean mixture
        denominator = 10.52 - 2 * phi[i]
        if denominator <= 0:
            chi_CO2[i] = chi_H2O[i] = chi_O2[i] = chi_N2[i] = np.nan
        else:
            chi_CO2[i] = 1 / denominator
            chi_H2O[i] = chi_CO2[i]
            chi_O2[i] = (2 - 2 * phi[i]) / denominator
            chi_N2[i] = 7.52 / denominator
    else:  # Rich mixture
        denominator = 7.52 + phi[i]
        if denominator <= 0:
            chi_CO2[i] = chi_H2O[i] = chi_O2[i] = chi_N2[i] = np.nan
        else:
            chi_CO2[i] = 1 / denominator
            chi_H2O[i] = chi_CO2[i]
            chi_O2[i] = 0
            chi_N2[i] = 7.52 / denominator

    # Molar enthalpy (kJ/kmol)
    H_mixture_molar = (
        chi_CO2[i] * H_CO2
        + chi_H2O[i] * H_H2O
        + chi_O2[i] * H_O2
        + chi_N2[i] * H_N2
    )
    H_molar.append(H_mixture_molar)

    # Average molecular weight of the mixture (kg/kmol)
    M_mixture = (
        chi_CO2[i] * M_CO2
        + chi_H2O[i] * M_H2O
        + chi_O2[i] * M_O2
        + chi_N2[i] * M_N2
    )

    # Mass enthalpy (kJ/kg)
    if M_mixture == 0:
        H_mixture_mass = np.nan
    else:
        H_mixture_mass = H_mixture_molar / M_mixture
    H_mass.append(H_mixture_mass)

# Remove NaN values for plotting
phi_valid = phi[~np.isnan(H_molar)]
H_molar_valid = np.array(H_molar)[~np.isnan(H_molar)]
H_mass_valid = np.array(H_mass)[~np.isnan(H_mass)]

# Plot the results
plt.figure(figsize=(12, 8))
plt.plot(phi_valid, H_molar_valid, label='Enthalpy per kmol of fuel (kJ/kmol)', color='blue', linewidth=2)
plt.plot(phi_valid, H_mass_valid, label='Enthalpy per kg of fuel (kJ/kg)', color='red', linewidth=2)
plt.title(f'Absolute Enthalpy of Fresh Reacting Mixture at T = {T} K', fontsize=16)
plt.xlabel('Equivalence Ratio (ϕ)', fontsize=14)
plt.ylabel('Enthalpy', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.tight_layout()
plt.show()
