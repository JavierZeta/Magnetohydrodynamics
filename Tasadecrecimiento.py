import numpy as np
import matplotlib.pyplot as plt

def plot_mri_growth_rate():
    # Parámetros físicos (normalizados)
    Omega = 1.0        # Velocidad angular local
    q = 1.5            # Perfil Kepleriano (q = -d ln Omega / d ln R)
    va = 1.0           # Velocidad de Alfvén (B0 / sqrt(mu0 * rho))
    
    # Frecuencia epicíclica al cuadrado: kappa^2 = 2 * (2 - q) * Omega^2
    kappa_sq = 2 * (2 - q) * (Omega**2)
    
    # Definir el rango de números de onda kz (normalizados por va/Omega)
    # El rango interesante suele ser de 0 hasta la estabilización
    kz_va = np.linspace(0, 2.0, 500)
    
    # Relación de dispersión: omega^4 - B*omega^2 + C = 0
    # B = kappa^2 + 2*(kz*va)^2
    # C = (kz*va)^2 * ((kz*va)^2 - 2*q*Omega^2)
    
    B = kappa_sq + 2 * (kz_va**2)
    C = (kz_va**2) * ((kz_va**2) - 2 * q * (Omega**2))
    
    # Resolvemos para omega^2 usando la fórmula cuadrática: 
    # omega^2 = [B +- sqrt(B^2 - 4C)] / 2
    discriminant = B**2 - 4*C
    omega_sq = (B - np.sqrt(discriminant)) / 2
    
    # La tasa de crecimiento gamma es la parte imaginaria de omega.
    # Si omega^2 < 0, entonces gamma = sqrt(-omega^2)
    gamma = np.where(omega_sq < 0, np.sqrt(np.abs(omega_sq)), 0)
    
    # Crear la gráfica
    plt.figure(figsize=(8, 5))
    plt.plot(kz_va, gamma, label=r'Modo inestable ($\gamma/\Omega$)', color='blue', lw=2)
    
    # Añadir elementos informativos
    plt.title('Tasa de Crecimiento de la MRI (Modelo Lineal)', fontsize=14)
    plt.xlabel(r'$k_z v_A / \Omega$ (Número de onda normalizado)', fontsize=12)
    plt.ylabel(r'$\gamma / \Omega$ (Tasa de crecimiento)', fontsize=12)
    
    # Marcar el máximo teórico (aprox 0.75 para kepleriano)
    max_gamma = np.max(gamma)
    kz_max = kz_va[np.argmax(gamma)]
    plt.scatter(kz_max, max_gamma, color='red')
    plt.annotate(f'Máximo: {max_gamma:.2f}', xy=(kz_max, max_gamma), 
                 xytext=(kz_max+0.1, max_gamma), arrowprops=dict(arrowstyle='->'))

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.axhline(0, color='black', lw=1)
    
    plt.show()

if __name__ == "__main__":
    plot_mri_growth_rate()