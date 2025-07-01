import numpy as np
import matplotlib.pyplot as plt
import os

def plot_wavefunction_evolution(directory, interval=100):
    files = sorted(f for f in os.listdir(directory) if f.endswith('.csv'))
    for i, file in enumerate(files):
        data = np.loadtxt(os.path.join(directory, file), delimiter=",")
        x = data[:, 0]
        psi_real = data[:, 1]
        psi_imag = data[:, 2]
        psi_abs = np.sqrt(psi_real**2 + psi_imag**2)

        plt.figure(figsize=(8, 4))
        plt.plot(x, psi_abs, label=f'Time step {i * interval}')
        plt.xlabel('x')
        plt.ylabel('|ψ(x)|')
        plt.title('Wavefunction Amplitude')
        plt.grid(True)
        plt.legend()
        plt.savefig(f"{directory}/frame_{i:04d}.png")
        plt.close()

# Example usage
plot_wavefunction_evolution("results/", interval=100)
