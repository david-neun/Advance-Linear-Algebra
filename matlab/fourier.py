import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# Define f(x) and parameters
f = lambda x: x - 1
L = 1  # interval [-1, 1]
N = 50  # compute coefficients for n = 0, 1, ..., 50

# Compute a0
a0, _ = integrate.quad(f, -L, L)
a0 /= (2 * L)

# Compute a_n and b_n for n = 1, ..., N
a = np.zeros(N + 1)
b = np.zeros(N + 1)
a[0] = a0

for n in range(1, N + 1):
    integrand_a = lambda x, n=n: f(x) * np.cos(n * np.pi * x / L)
    integrand_b = lambda x, n=n: f(x) * np.sin(n * np.pi * x / L)
    a[n], _ = integrate.quad(integrand_a, -L, L)
    a[n] /= L
    b[n], _ = integrate.quad(integrand_b, -L, L)
    b[n] /= L

# Print coefficients
print("Fourier Coefficients for f(x) = x - 1 on [-1, 1]")
print("=" * 50)
print(f"{'n':>4s}  {'a_n':>20s}  {'b_n':>20s}")
print("-" * 50)
print(f"{'0':>4s}  {a[0]:>20.15f}  {'N/A':>20s}")
for n in range(1, N + 1):
    print(f"{n:>4d}  {a[n]:>20.15f}  {b[n]:>20.15f}")

# Fourier approximation function
def fourier_approx(x, a, b, L, N):
    result = a[0]  # a0
    for n in range(1, N + 1):
        result += a[n] * np.cos(n * np.pi * x / L) + b[n] * np.sin(n * np.pi * x / L)
    return result

# --- Plot on [-1, 1] and [-3, 3] ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot on [-1, 1]
x1 = np.linspace(-1, 1, 1000)
y1_exact = f(x1)
y1_approx = fourier_approx(x1, a, b, L, N)

axes[0].plot(x1, y1_exact, 'b-', linewidth=2, label=r'$f(x) = x - 1$')
axes[0].plot(x1, y1_approx, 'r--', linewidth=1.5, label=f'Fourier approx ($N={N}$)')
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[0].set_title(r'Fourier Approximation on $[-1,\,1]$')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Plot on [-3, 3]
x2 = np.linspace(-3, 3, 3000)
y2_exact = f(x2)
y2_approx = fourier_approx(x2, a, b, L, N)

axes[1].plot(x2, y2_exact, 'b-', linewidth=2, label=r'$f(x) = x - 1$')
axes[1].plot(x2, y2_approx, 'r--', linewidth=1.5, label=f'Fourier approx ($N={N}$)')
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')
axes[1].set_title(r'Fourier Approximation on $[-3,\,3]$')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fourier_plot.png', dpi=150)
print("\nPlot saved to fourier_plot.png")
