import matplotlib.pyplot as plt
import numpy as np

# Sample data for each temperature
current_density = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

# Experimental data
voltage_exp_750 = [0.9, 1.1, 1.25, 1.35, 1.45, 1.5, 1.55, 1.58, 1.6]
voltage_exp_850 = [0.93, 1.05, 1.2, 1.3, 1.35, 1.4, 1.45, 1.48, 1.5]
voltage_exp_950 = [0.9, 0.95, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35]

# Model data
voltage_mod_750 = [0.9, 1.0, 1.15, 1.25, 1.35, 1.4, 1.45, 1.48, 1.5]
voltage_mod_850 = [0.9, 1.02, 1.18, 1.28, 1.34, 1.38, 1.42, 1.45, 1.47]
voltage_mod_950 = [0.9, 0.98, 1.04, 1.08, 1.12, 1.18, 1.22, 1.25, 1.28]

current_density_850 = np.linspace(0., 1, 11)
voltage_exp_850 = [0.94, 0.97, 0.99, 1.01, 1.03, 1.05, 1.07, 1.085, 1.1, 1.12, 1.14]

import pickle

# Load data from temp_data.pkl
with open('temp_data.pkl', 'rb') as file:
       data = pickle.load(file)

current_density_850_mod = data['current_density']
voltage_mod_850 = data['V_cell']

# Plotting
plt.figure(figsize=(8, 6))

# Scatter plot for experimental data
# plt.scatter(current_density, voltage_exp_750, color="green", label="Exp 750°C", edgecolor="black", s=50)
# plt.scatter(current_density, voltage_exp_850, color="orange", label="Exp 850°C", edgecolor="black", s=50)
plt.scatter(current_density_850, voltage_exp_850, color="blue", label="Experimental 850°C", edgecolor="black", s=75, zorder=4)

# Line plot for model data
# plt.plot(current_density, voltage_mod_750, color="green", linestyle="--", label="Mod 750°C")
# plt.plot(current_density, voltage_mod_850, color="orange", linestyle="--", label="Mod 850°C")
plt.plot(current_density_850_mod, voltage_mod_850, color="blue", linestyle="--", label="Model 850°C", lw=3)

# Labels and legend
plt.xlabel("Current density [A/cm²]", fontsize=14)
plt.ylabel("Voltage [V]", fontsize=14)
plt.legend(fontsize=14)
plt.grid(True)
plt.ylim(0.9, 1.2)
plt.xlim(0, 1.)

# Display the plot
plt.savefig('validated_data.png', dpi=300)
