import numpy as np
import matplotlib.pyplot as plt

# =========================
# INPUT DATA
# =========================
# x and y experimental / numerical values
x = np.array([
    2.9946e-12,
    8.9839e-13,
    2.9946e-13,
    8.9839e-14
])

y = np.array([
    0.000124160647886,
    3.72461239575e-05,
    1.24156095889e-05,
    3.72472869004e-06
])

# Target y value (e.g. scaled experimental wear volume)
y_target = 7.1296e-6


# =========================
# LINEAR FIT
# =========================
a, b = np.polyfit(x, y, 1)
x_target = (y_target - b) / a

# =========================
# PLOT
# =========================
x_fit = np.linspace(0.0, 1.05 * x.max(), 300)
y_fit = a * x_fit + b

plt.figure(figsize=(7,5))

# Numerical results (thin markers)
plt.plot(
    x, y,
    linestyle='',
    marker='o',
    markersize=5,
    color='#1f4e79',
    label='Numerical results'
)

# Linear fit (clean solid line)
plt.plot(
    x_fit, y_fit,
    linewidth=2.0,
    color='#1f4e79',
    label='Linear fit'
)

# Experimental scaled volume
plt.axhline(
    y_target,
    linestyle='--',
    linewidth=1.2,
    color='gray',
    label=rf'$V_{{exp,scaled}} = {y_target:.2e}$'
)

# Calibrated point
plt.plot(
    x_target, y_target,
    marker='s',
    markersize=6,
    color='#2e8b57',
    label=rf'$k_{{cal}} = {x_target:.2e}$'
)

# Equation text close to the line
x_eq = 0.5 * x.max()
y_eq = 0.00005

equation_text = rf'$V = {a:.2e}\,k {b:+.2e}$'
plt.text(
    x_eq,
    y_eq * 0.95,
    equation_text,
    fontsize=11,
    color='black'
)

# Labels and layout
plt.xlabel(r'Wear coefficient $k$ [mm$^3$/Nmm]')
plt.ylabel(r'Wear volume loss $V$ [mm$^3$]')
plt.legend(frameon=True)
plt.grid(True, linewidth=0.6, alpha=0.7)
plt.tight_layout()
plt.show()

# =========================
# OUTPUT
# =========================
print('Linear fit:')
print(f'  V = {a:.4e} * k + {b:.4e}')
print()
print('Calibrated wear coefficient:')
print(f'  k_cal = {x_target:.4e} mm^3/(N·mm)')