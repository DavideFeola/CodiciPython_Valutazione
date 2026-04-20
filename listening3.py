import numpy as np
import matplotlib.pyplot as plt

# Parametri
G_dB = 120
G = 10**(G_dB / 10)
alpha = 0.188

LM_values = np.arange(2, 11, 1)
pout_targets = [0.25, 0.10, 0.05, 0.01]
labels = ['25%', '10%', '5%', '1%']

def K_from_margin(M_dB):
    return (np.log(10) / 10) * M_dB

def pout_circ(sigma, G, K, alpha=0.188):
    x = K / (2 * sigma**2 * alpha * G)
    return (1 + x) * np.exp(-x)

def bisection_circ(G, M_dB, pout_target, alpha=0.188, tol=1e-15):
    K = K_from_margin(M_dB)
    a, b = 1e-9, 1e-3
    for _ in range(100):
        mid = (a + b) / 2
        if pout_circ(mid, G, K, alpha) > pout_target:
            b = mid
        else:
            a = mid
        if (b - a) < tol:
            break
    return (a + b) / 2

# Tabella
print(f"{'LM [dB]':<10}", end='')
for p in labels:
    print(f"{p:>12}", end='')
print()

for M in LM_values:
    print(f"{M:<10.2f}", end='')
    for p in pout_targets:
        s = bisection_circ(G, M, p) * 1e6
        print(f"{s:>12.3f}", end='')
    print()

plt.figure(figsize=(10, 6))
for p, lbl in zip(pout_targets, labels):
    sigma_vals = [bisection_circ(G, M, p) * 1e6 for M in LM_values]
    plt.plot(LM_values, sigma_vals, marker='o', label=f'$P_{{out}}$ = {lbl}')

plt.xlabel('Margine di collegamento LM [dB]')
plt.ylabel(r'$\sigma_{\theta,\max}$ [$\mu$rad]')
plt.legend(frameon=False)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
