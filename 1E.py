import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Constants
P = 101325  # Atmospheric pressure in Pa
T_initial = 298.15  # Initial temperature in Kelvin
T_standard = 298.15  # Standard conditions temperature in Kelvin
phi_range = np.arange(0.5, 2.01, 0.01)  # Equivalence ratio range
species_list = ["CH4", "O2", "N2", "CO2", "H2O", "H2", "CO", "OH"]  # Relevant species
fuel = "CH4"  # Fuel
oxidizer = ["O2", "N2"]  # Oxidizer
molar_mass_fuel = 16.04  # Molar mass of CH4 in g/mol

# Prepare results storage
adiabatic_flame_temp = np.zeros(len(phi_range))
product_mole_fractions = np.zeros((len(phi_range), len(species_list)))
combustion_enthalpy_kmol = np.zeros(len(phi_range))  # Combustion enthalpy per kmol of CH4
combustion_enthalpy_kg = np.zeros(len(phi_range))  # Combustion enthalpy per kg of CH4

# Loop over equivalence ratios
for i, phi in enumerate(phi_range):
    # Define stoichiometry
    oxidizer_ratio = 2 / phi  # Adjust O2 for equivalence ratio
    nitrogen_ratio = oxidizer_ratio * 3.76  # Air composition: 21% O2, 79% N2 by volume
    reactants = f"{fuel}:{phi}, {oxidizer[0]}:{oxidizer_ratio}, {oxidizer[1]}:{nitrogen_ratio}"

    # Create a gas object for GRI-Mech 3.0
    gas = ct.Solution("gri30.yaml")
    gas.TPX = T_initial, P, reactants

    # Calculate enthalpy of reactants at initial conditions
    h_reactants = gas.enthalpy_mass

    # Equilibrate to find adiabatic flame temperature and products
    gas.equilibrate("HP")  # Adiabatic conditions: constant enthalpy and pressure
    adiabatic_flame_temp[i] = gas.T  # Flame temperature

    # Mole fractions of products
    for j, species in enumerate(species_list):
        product_mole_fractions[i, j] = gas[species].X[0]  # Mole fraction of each species

    # Reset to standard conditions
    gas.TPX = T_standard, P, reactants

    # Enthalpy of reactants at standard conditions
    h_reactants_std = gas.enthalpy_mass

    # Equilibrate at standard conditions
    gas.equilibrate("TP")  # Equilibrium at constant temperature and pressure
    h_products_std = gas.enthalpy_mass

    # Combustion enthalpy at standard conditions
    combustion_enthalpy_kmol[i] = (h_reactants_std - h_products_std) * molar_mass_fuel  # J/kmol
    combustion_enthalpy_kg[i] = combustion_enthalpy_kmol[i] / molar_mass_fuel * 1e3  # J/kg

# Plot 1: Adiabatic Flame Temperature vs Equivalence Ratio
plt.figure(figsize=(10, 8))
plt.plot(phi_range, adiabatic_flame_temp, label="Adiabatic Flame Temperature (K)", color="blue", linewidth=1.5)
plt.title("Adiabatic Flame Temperature vs Equivalence Ratio")
plt.xlabel("Equivalence Ratio (Phi)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.grid()
plt.show()

# Plot 2: Product Mole Fractions vs Equivalence Ratio
plt.figure(figsize=(10, 8))
for j, species in enumerate(species_list):
    plt.plot(phi_range, product_mole_fractions[:, j], label=species, linewidth=1.5)
plt.title("Product Mole Fractions vs Equivalence Ratio")
plt.xlabel("Equivalence Ratio (Phi)")
plt.ylabel("Mole Fraction")
plt.legend(loc="best")
plt.grid()
plt.show()

# Plot 3: Combustion Enthalpy vs Equivalence Ratio
plt.figure(figsize=(10, 8))
plt.plot(phi_range, combustion_enthalpy_kmol, label="Combustion Enthalpy (J/kmol of CH4)", color="blue", linewidth=1.5)
plt.plot(phi_range, combustion_enthalpy_kg, label="Combustion Enthalpy (J/kg of CH4)", color="red", linewidth=1.5)
plt.title("Combustion Enthalpy vs Equivalence Ratio")
plt.xlabel("Equivalence Ratio (Phi)")
plt.ylabel("Combustion Enthalpy")
plt.legend()
plt.grid()
plt.show()





