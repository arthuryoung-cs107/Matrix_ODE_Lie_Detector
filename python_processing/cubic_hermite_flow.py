import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------
# SYMBOLIC CONSTRUCTION
# ----------------------------------
# Define symbolic variables
h, x = sp.symbols('h x')
u_mh, up_mh, u_ph, up_ph = sp.symbols('u(-h) u\'(-h) u(h) u\'(h)')

# Define the inverse matrix A_inv
A_inv = sp.Matrix([
    [1/2, h/4, 1/2, -h/4],
    [-3/(4*h), -1/4, 3/(4*h), -1/4],
    [0, -1/(4*h), 0, 1/(4*h)],
    [1/(4*h**3), 1/(4*h**2), -1/(4*h**3), 1/(4*h**2)]
])

# Define the vector y
y = sp.Matrix([u_mh, up_mh, u_ph, up_ph])

# Compute the coefficient vector [a0, a1, a2, a3]
coeffs = A_inv * y
coeffs_simplified = [sp.simplify(c) for c in coeffs]

# ----------------------------------
# NUMERICAL EVALUATION
# ----------------------------------
# Define h and the target function
h_val = 1.0
u_func = lambda x: np.sin(x)
up_func = lambda x: np.cos(x)

# Sample u and u' at -h and h
umh = u_func(-h_val)
upmh = up_func(-h_val)
uph = u_func(h_val)
upph = up_func(h_val)

# Substitute numerical values
subs_dict = {h: h_val, u_mh: umh, up_mh: upmh, u_ph: uph, up_ph: upph}
numerical_coeffs = [float(c.subs(subs_dict).evalf()) for c in coeffs_simplified]
a0, a1, a2, a3 = numerical_coeffs

# Define interpolating polynomial and its derivative
p = lambda x: a0 + a1 * x + a2 * x**2 + a3 * x**3
dp = lambda x: a1 + 2 * a2 * x + 3 * a3 * x**2

# ----------------------------------
# PLOTTING INTERPOLANT
# ----------------------------------
x_vals = np.linspace(-1.5*h_val, 1.5*h_val, 400)
true_vals = u_func(x_vals)
interp_vals = p(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, true_vals, label='Original function $u(x) = \sin(x)$', linewidth=2)
plt.plot(x_vals, interp_vals, '--', label='Cubic Hermite Interpolant', linewidth=2)
plt.scatter([-h_val, h_val], [umh, uph], color='red', zorder=5, label='Interpolation Points')
plt.title("Cubic Hermite Interpolation using values and derivatives at $x = \pm h$")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# ----------------------------------
# PLOTTING DERIVATIVE
# ----------------------------------
dtrue_vals = up_func(x_vals)
dinterp_vals = dp(x_vals)

plt.figure(figsize=(10, 6))
plt.plot(x_vals, dtrue_vals, label="True derivative $u'(x) = \cos(x)$", linewidth=2)
plt.plot(x_vals, dinterp_vals, '--', label="Derivative of Hermite Interpolant", linewidth=2)
plt.scatter([-h_val, h_val], [upmh, upph], color='red', zorder=5, label='Matched Derivatives')
plt.title("Derivative of the Cubic Hermite Interpolant vs True Derivative")
plt.xlabel("x")
plt.ylabel("u'(x)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
