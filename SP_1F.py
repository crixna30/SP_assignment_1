import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Constants
P = 101325  # Atmospheric pressure in Pa
T_initial = 298.15  # Initial temperature in Kelvin
T_standard = 298.15  # Standard conditions in Kelvin
phi_range = np.arange(0.5, 2.01, 0.01)  # Equivalence ratio range
fuel = "CH4"  # Fuel
oxidizer = ["O2", "N2"]  # Oxidizer
molar_mass_fuel = 16.04  # Molar mass of CH4 in g/mol
species_list = ["CO2", "H2O", "O2", "CO", "H2", "OH", "O", "H", "N2", "N", "NO", "CH4"]  # Expanded species list

# Prepare results storage
adiabatic_flame_temp = np.zeros(len(phi_range))
product_mole_fractions = np.zeros((len(phi_range), len(species_list)))
combustion_enthalpy_kmol = np.zeros(len(phi_range))  # Combustion enthalpy per kmol of CH4
combustion_enthalpy_kg = np.zeros(len(phi_range))  # Combustion enthalpy per kg of CH4

# Loop over equivalence ratios
for i, phi in enumerate(phi_range):
    # Define reactant composition
    oxidizer_ratio = 2 / phi  # O2 amount adjusted for equivalence ratio
    nitrogen_ratio = oxidizer_ratio * 3.76  # Air composition: 21% O2, 79% N2
    reactants = f"{fuel}:{phi}, O2:{oxidizer_ratio}, N2:{nitrogen_ratio}"

    # Create gas object for GRI-Mech 3.0 mechanism
    gas = ct.Solution("gri30.yaml")
    gas.TPX = T_initial, P, reactants

    # Calculate enthalpy of reactants at initial conditions
    h_reactants = gas.enthalpy_mass

    # Equilibrate gas mixture to find adiabatic flame temperature and products
    gas.equilibrate("HP")  # Adiabatic flame: constant enthalpy and pressure
    adiabatic_flame_temp[i] = gas.T  # Store flame temperature

    # Get mole fractions of all specified products
    for j, species in enumerate(species_list):
        product_mole_fractions[i, j] = gas[species].X[0]

    # Reset to standard conditions
    gas.TPX = T_standard, P, reactants

    # Enthalpy of reactants at standard conditions
    h_reactants_std = gas.enthalpy_mass

    # Equilibrate at standard conditions
    gas.equilibrate("TP")  # Equilibrium at constant T and P
    h_products_std = gas.enthalpy_mass

    # Combustion enthalpy at standard conditions
    combustion_enthalpy_kmol[i] = (h_reactants_std - h_products_std) * molar_mass_fuel  # J/kmol
    combustion_enthalpy_kg[i] = combustion_enthalpy_kmol[i] / molar_mass_fuel * 1e3  # J/kg

# Plot 1: Adiabatic Flame Temperature vs Equivalence Ratio
plt.figure(figsize=(10, 8))
plt.plot(phi_range, adiabatic_flame_temp, color="blue", label="Adiabatic Flame Temperature (K)", linewidth=2)
plt.xlabel("Equivalence Ratio (Phi)", fontsize=14)
plt.ylabel("Adiabatic Flame Temperature (K)", fontsize=14)
plt.title("Adiabatic Flame Temperature vs Equivalence Ratio", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()

# Plot 2: Product Mole Fractions vs Equivalence Ratio
plt.figure(figsize=(10, 8))
for j, species in enumerate(species_list):
    plt.plot(phi_range, product_mole_fractions[:, j], linewidth=2, label=species)
plt.xlabel("Equivalence Ratio (Phi)", fontsize=14)
plt.ylabel("Mole Fraction", fontsize=14)
plt.title("Product Mole Fractions vs Equivalence Ratio", fontsize=16)
plt.legend(loc="upper right", fontsize=12, ncol=2)
plt.grid(True)
plt.show()

# Plot 3: Combustion Enthalpy vs Equivalence Ratio
plt.figure(figsize=(10, 8))
plt.plot(phi_range, combustion_enthalpy_kmol, label="Combustion Enthalpy (J/kmol of CH4)", color="blue", linewidth=2)
plt.plot(phi_range, combustion_enthalpy_kg, label="Combustion Enthalpy (J/kg of CH4)", color="red", linewidth=2)
plt.xlabel("Equivalence Ratio (Phi)", fontsize=14)
plt.ylabel("Combustion Enthalpy", fontsize=14)
plt.title("Combustion Enthalpy vs Equivalence Ratio", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()






