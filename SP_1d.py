import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Constants
R = 8.314  # Universal gas constant [J/mol-K]
T0 = 298.15  # Initial temperature [K]
phi_range = np.arange(0.5, 2.01, 0.01)  # Equivalence ratio range
n_species = 7  # CO2, H2O, CO, H2, O2, N2, CH4

# Initialize arrays
mole_fractions = np.zeros((len(phi_range), n_species))
T_ad_dissociation = np.zeros(len(phi_range))

# NASA Polynomial Coefficients (Valid for T >= 1000 K)
coeffs_CO2 = [3.85746029, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806]
coeffs_H2O = [3.03399249, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010]
coeffs_CO = [3.57953347, 3.05176840E-04, -1.69469055E-07, 4.55849056E-11, -4.55225040E-15, -1.43440860E+04, 3.50840928]
coeffs_H2 = [2.34433112, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238]
coeffs_O2 = [3.28253784, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129]
coeffs_N2 = [2.95257626, 1.39690000E-03, -4.92631603E-07, 7.86010367E-11, -4.60755321E-15, -9.23948688E+02, 5.87188762]
coeffs_CH4 = [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180]

# Specific enthalpy function
def specific_enthalpy(T, coeffs):
    """Calculate specific enthalpy using NASA polynomials."""
    h = R * T * (
        coeffs[0]
        + coeffs[1] * T / 2
        + coeffs[2] * T**2 / 3
        + coeffs[3] * T**3 / 4
        + coeffs[4] * T**4 / 5
        + coeffs[5] / T
    )
    return h

# Function to calculate equilibrium mole fractions
def equilibrium_mole_fractions(T, phi):
    """Calculate equilibrium mole fractions at temperature T and equivalence ratio phi."""
    if phi <= 1:
        # Lean combustion
        x_CO2 = 1 - 0.1 * np.exp(-3000 / T)
        x_H2O = 2 - 0.2 * np.exp(-2500 / T)
        x_O2 = 2 / phi - 2 + 0.1 * np.exp(-3000 / T) + 0.2 * np.exp(-2500 / T)
        x_N2 = (79 / 21) * (2 / phi)
        x_CO = 0.1 * np.exp(-3000 / T)
        x_H2 = 0.2 * np.exp(-2500 / T)
        x_CH4 = 0
    else:
        # Rich combustion
        x_CO2 = 1 - 0.1 * np.exp(-3000 / T)
        x_H2O = 2 - 0.2 * np.exp(-2500 / T)
        x_O2 = 0
        x_N2 = (79 / 21) * 2
        x_CO = 0.1 * np.exp(-3000 / T)
        x_H2 = 0.2 * np.exp(-2500 / T)
        x_CH4 = phi - 1
    total = x_CO2 + x_H2O + x_O2 + x_N2 + x_CO + x_H2 + x_CH4
    return np.array([x_CO2, x_H2O, x_CO, x_H2, x_O2, x_N2, x_CH4]) / total

def equilibrium_mole_fractions(T, phi):
    """Calculate equilibrium mole fractions at temperature T and equivalence ratio phi."""
    T = np.atleast_1d(T)[0]  # Ensure T is scalar
    
    if phi <= 1:
        # Lean combustion
        x_CO2 = 1 - 0.1 * np.exp(-3000 / T)
        x_H2O = 2 - 0.2 * np.exp(-2500 / T)
        x_O2 = 2 / phi - 2 + 0.1 * np.exp(-3000 / T) + 0.2 * np.exp(-2500 / T)
        x_N2 = (79 / 21) * (2 / phi)
        x_CO = 0.1 * np.exp(-3000 / T)
        x_H2 = 0.2 * np.exp(-2500 / T)
        x_CH4 = 0
    else:
        # Rich combustion
        x_CO2 = 1 - 0.1 * np.exp(-3000 / T)
        x_H2O = 2 - 0.2 * np.exp(-2500 / T)
        x_O2 = 0
        x_N2 = (79 / 21) * 2
        x_CO = 0.1 * np.exp(-3000 / T)
        x_H2 = 0.2 * np.exp(-2500 / T)
        x_CH4 = phi - 1
    total = x_CO2 + x_H2O + x_O2 + x_N2 + x_CO + x_H2 + x_CH4
    return np.array([x_CO2, x_H2O, x_CO, x_H2, x_O2, x_N2, x_CH4]) / total


# Main computation
for i, phi in enumerate(phi_range):
    # Initial guess for adiabatic flame temperature
    T_guess = 2000

    # Compute H_reactants
    n_CH4 = 1  # 1 mole of CH4
    n_O2 = (2 / phi) * n_CH4
    n_N2 = (79 / 21) * n_O2
    h_CH4 = specific_enthalpy(T0, coeffs_CH4)
    h_O2 = specific_enthalpy(T0, coeffs_O2)
    h_N2 = specific_enthalpy(T0, coeffs_N2)
    H_reactants = n_CH4 * h_CH4 + n_O2 * h_O2 + n_N2 * h_N2

    # Iterative solver
    def adiabatic_with_dissociation(T):
        x = equilibrium_mole_fractions(T, phi)
        x_CO2, x_H2O, x_CO, x_H2, x_O2, x_N2, x_CH4 = x

        # Compute enthalpy of products
        h_CO2 = specific_enthalpy(T, coeffs_CO2)
        h_H2O = specific_enthalpy(T, coeffs_H2O)
        h_CO = specific_enthalpy(T, coeffs_CO)
        h_H2 = specific_enthalpy(T, coeffs_H2)
        h_O2 = specific_enthalpy(T, coeffs_O2)
        h_N2 = specific_enthalpy(T, coeffs_N2)
        h_CH4 = specific_enthalpy(T, coeffs_CH4)

        H_products = (
            x_CO2 * h_CO2
            + x_H2O * h_H2O
            + x_CO * h_CO
            + x_H2 * h_H2
            + x_O2 * h_O2
            + x_N2 * h_N2
            + x_CH4 * h_CH4
        )
        return H_reactants - H_products

    T_ad_dissociation[i] = fsolve(adiabatic_with_dissociation, T_guess)[0]

# Plot results
plt.figure(figsize=(10, 8))
plt.plot(phi_range, T_ad_dissociation, color="red", label="Adiabatic Flame Temperature with Dissociation (K)")
plt.xlabel("Equivalence Ratio (ϕ)", fontsize=14)
plt.ylabel("Adiabatic Flame Temperature (K)", fontsize=14)
plt.title("Adiabatic Flame Temperature with Dissociation vs Equivalence Ratio", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()

