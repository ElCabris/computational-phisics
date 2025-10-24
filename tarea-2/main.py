import numpy as np
import matplotlib.pyplot as plt

# =============================================
# PARÁMETROS F\'ISICOS Y NUMÉRICOS
# =============================================

v0 = 20.0      # velocidad inicial [m/s] - típico lanzamiento de proyectil
theta = np.pi/4  # ángulo de lanzamiento (45 grados = π/4 radianes)
g = 9.81       # aceleración gravitacional [m/s²] - valor estándar
h = 0.1        # paso de integración [s] - determina la precisión
t_max = 3.0    # tiempo máximo de simulación [s] - para evitar loops infinitos

# Explicación: Estos parámetros controlan tanto la física del problema como la
# precisión del método numérico. Un h más pequeño da mayor precisión pero mayor costo computacional.

# =============================================
# CONDICIONES INICIALES
# =============================================

# Descomposición de la velocidad inicial en componentes x e y
# vx0 = v0 * cos(theta): componente horizontal (coseno = adyacente/hipotenusa)
# vy0 = v0 * sin(theta): componente vertical (seno = opuesto/hipotenusa)
vx0 = v0 * np.cos(theta)
vy0 = v0 * np.sin(theta)

# Listas para almacenar la evolución temporal de las variables
# Iniciamos en el origen (0,0) con las velocidades calculadas
x, y = [0.0], [0.0]          # Posiciones iniciales
vx, vy = [vx0], [vy0]        # Velocidades iniciales
t_values = [0.0]             # Tiempos (para tracking temporal)

print(f"Condiciones iniciales:")
print(f"Velocidad inicial: {v0:.2f} m/s")
print(f"Ángulo: {np.degrees(theta):.1f}°")
print(f"vx0: {vx0:.2f} m/s, vy0: {vy0:.2f} m/s")

# ======================================
# MÉTODO DE EULER - INTEGRACIÓN NUMÉRICA
# ======================================

t = 0.0
step = 0
print("\nIniciando integración numérica...")

# El loop continúa mientras el proyectil esté arriba del suelo (y >= 0)
while y[-1] >= 0 and t < t_max:
    
    # --- PASO 1: Calcular nuevas posiciones usando velocidades actuales ---
    # x_new = x_actual + Δt * vx_actual (movimiento horizontal uniforme)
    # y_new = y_actual + Δt * vy_actual (movimiento vertical acelerado)
    x_new = x[-1] + h * vx[-1]
    y_new = y[-1] + h * vy[-1]
    
    # --- PASO 2: Calcular nuevas velocidades ---
    # vx permanece constante (no hay fuerza horizontal)
    # vy disminuye debido a la gravedad: vy_new = vy_actual - g*Δt
    vx_new = vx[-1]  
    vy_new = vy[-1] - h * g
    
    # --- PASO 3: Almacenar los nuevos valores ---
    x.append(x_new)
    y.append(y_new)
    vx.append(vx_new)
    vy.append(vy_new)
    t += h
    t_values.append(t)
    step += 1
    
    # Información de progreso cada 10 pasos
    if step % 10 == 0:
        print(f"Paso {step}: t = {t:.2f}s, x = {x_new:.2f}m, y = {y_new:.2f}m")

print(f"Simulación completada: {step} pasos, tiempo final = {t:.2f}s")
print(f"Alcance numérico: {x[-1]:.2f} m")

# ==================================
# SOLUCIÓN ANALÍTICA EXACTA (PARA COMPARACIÓN)
# ==================================

# La solución exacta para movimiento parabólico sin rozamiento:
# x(t) = v0*cos(θ)*t
# y(t) = v0*sin(θ)*t - 0.5*g*t²
# Eliminando el tiempo: y(x) = x*tan(θ) - (g*x²)/(2*v0²*cos²(θ))

print("\nCalculando solución analítica...")

# Crear array de posiciones x para evaluar la trayectoria analítica
x_analytic = np.linspace(0, max(x), 200)

# Calcular las y correspondientes usando la ecuación de la trayectoria
# y(x) = x*tan(θ) - (g*x²)/(2*v0²*cos²(θ))
tan_theta = np.tan(theta)
cos2_theta = np.cos(theta)**2
denominador = 2 * v0**2 * cos2_theta

y_analytic = x_analytic * tan_theta - (g * x_analytic**2) / denominador

# También calculamos el alcance teórico exacto
alcance_teorico = (v0**2 * np.sin(2*theta)) / g
altura_max_teorica = (v0**2 * np.sin(theta)**2) / (2*g)

print(f"Alcance teórico: {alcance_teorico:.2f} m")
print(f"Altura máxima teórica: {altura_max_teorica:.2f} m")

# ===========================
# VISUALIZACIÓN Y COMPARACIÓN
# ===========================

plt.figure(figsize=(10, 6))

# Graficar solución numérica (método de Euler)
plt.plot(x, y, 'bo-', markersize=3, linewidth=1, label=f'Euler (h={h}s)')

# Graficar solución analítica exacta
plt.plot(x_analytic, y_analytic, 'r-', linewidth=2, label='Solución analítica')

# Configuración del gráfico
plt.xlabel("Distancia horizontal, x [m]", fontsize=12)
plt.ylabel("Altura, y [m]", fontsize=12)
plt.title("Trayectoria de Proyectil: Método de Euler vs Solución Analítica", fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.axis('equal')  # Misma escala en ambos ejes para ver la forma real

# Añadir información adicional en el gráfico
plt.text(0.05, 0.95, f'v₀ = {v0} m/s, θ = {np.degrees(theta):.0f}°', 
         transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.show()

# ===================
# ANÁLISIS DEL ERROR
# ==================

print("\n" + "="*50)
print("ANÁLISIS DE PRECISIÓN")
print("="*50)

# Calcular error en el alcance
error_alcance = abs(x[-1] - alcance_teorico)
error_relativo = (error_alcance / alcance_teorico) * 100

print(f"Alcance numérico: {x[-1]:.4f} m")
print(f"Alcance teórico:  {alcance_teorico:.4f} m")
print(f"Error absoluto:   {error_alcance:.4f} m")
print(f"Error relativo:   {error_relativo:.2f}%")

# Sugerencia para mejorar la precisión
if error_relativo > 1:
    print("\nRecomendación: Reducir el paso h para mayor precisión")
    print(f"Pruebe con h = 0.01 o h = 0.001")
else:
    print("\n¡Buena precisión! El método funciona bien para este paso h")
