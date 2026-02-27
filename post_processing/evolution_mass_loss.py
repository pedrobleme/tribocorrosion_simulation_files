import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# Configuration
# --------------------------------------------------

# Material density [g/mm^3]
density = 0.00785  

# Initial volume [mm^3]
initial_volume = 60.18880  

# Files and corrosion rates
cases = {
    r"$\dot d_c = 0.13$ mm/year": "enh_1.csv",
    r"$\dot d_c = 0.471$ mm/year": "enh_2.csv",
    r"$\dot d_c = 0.85$ mm/year": "enh_3.csv"
}

# Academic matplotlib style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.8
})

fig, ax = plt.subplots(figsize=(7, 5))

# --------------------------------------------------
# Plotting
# --------------------------------------------------

for label, file_path in cases.items():
    data = pd.read_csv(file_path)

    # Time column
    time = data["Time"]

    # SDV2 column (corrosion damage variable)
    sdv2 = data.iloc[:, 2]

    # Corroded volume [mm^3]
    corroded_volume = sdv2 * initial_volume

    # Mass loss [kg]
    mass_loss = corroded_volume * density / 1000.0

    ax.plot(
        time,
        mass_loss,
        linestyle='-',
        label=label
    )

    print(f"Final mass loss for {label}: {mass_loss.iloc[-1]:.4f} kg")

# --------------------------------------------------
# Axes, grid and labels
# --------------------------------------------------

ax.set_xlabel("Time (s)")
ax.set_ylabel("Cumulative mass loss (kg)")

ax.set_ylim(bottom=0.0)

ax.grid(
    True,
    which='major',
    linestyle='--',
    linewidth=0.5,
    alpha=0.5
)

ax.legend(loc="upper left", frameon=True)

plt.tight_layout()
plt.show()