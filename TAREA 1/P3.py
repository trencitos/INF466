import random
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Definición de Constantes y Funciones ---

# Alfabeto de 6 bases para la especie alienígena
ALPHABET = ['A', 'C', 'G', 'T', 'B', 'D']
SEQ_LENGTH = 200

# A. Generar una secuencia aleatoria
def generar_secuencia_aleatoria(length=SEQ_LENGTH, alphabet=ALPHABET):
    """Genera una secuencia aleatoria de la longitud y alfabeto dados."""
    return [random.choice(alphabet) for _ in range(length)]

# B. Función de Mutación
def aplicar_mutacion(sequence, alphabet=ALPHABET):
    """Aplica una mutación (inserción, borrado o reemplazo) equiprobable."""
    
    # 1. Elegir el tipo de mutación equiprobablemente
    mutation_type = random.choice(['insercion', 'borrado', 'reemplazo'])
    
    # 2. Elegir la posición al azar
    # En una secuencia de N elementos, hay N+1 lugares para inserción (incluyendo inicio/fin)
    # y N lugares para borrado/reemplazo.
    N = len(sequence)
    
    if mutation_type == 'insercion':
        # Posición entre 0 y N (N es el final)
        pos = random.randint(0, N)
        # Elegir la base a insertar equiprobablemente
        base_a_insertar = random.choice(alphabet)
        
        # Aplicar la inserción
        nueva_secuencia = sequence[:pos] + [base_a_insertar] + sequence[pos:]
        return nueva_secuencia
    
    elif mutation_type == 'borrado':
        if N == 0:
            return [] # No se puede borrar de una secuencia vacía
        # Posición entre 0 y N-1
        pos = random.randint(0, N - 1)
        
        # Aplicar el borrado
        nueva_secuencia = sequence[:pos] + sequence[pos+1:]
        return nueva_secuencia
        
    elif mutation_type == 'reemplazo':
        if N == 0:
            return []
        # Posición entre 0 y N-1
        pos = random.randint(0, N - 1)
        
        # Elegir la nueva base (cualquiera de las otras 5, equiprobablemente)
        base_actual = sequence[pos]
        otras_bases = [b for b in alphabet if b != base_actual]
        base_nueva = random.choice(otras_bases)
        
        # Aplicar el reemplazo
        nueva_secuencia = sequence[:] # Copia
        nueva_secuencia[pos] = base_nueva
        return nueva_secuencia
    
    return sequence # En caso de error, devuelve la secuencia original

# C. Distancia de Levenshtein (Implementando Needleman-Wunsch)
# Parámetros para Levenshtein (similar a Needleman-Wunsch):
# Match: 0, Mismatch: 1, Gap Penalty: 1
def distancia_levenshtein(s1, s2):
    """
    Calcula la Distancia de Levenshtein (mínimas operaciones de edición)
    utilizando el algoritmo de Programación Dinámica (similar a Needleman-Wunsch
    con costos de edición).
    """
    m, n = len(s1), len(s2)
    
    # Inicialización de la matriz (m+1 filas, n+1 columnas)
    D = np.zeros((m + 1, n + 1), dtype=int)

    # Llenar la primera fila y columna (costo de gaps iniciales)
    for i in range(m + 1):
        D[i, 0] = i # i borrados
    for j in range(n + 1):
        D[0, j] = j # j inserciones

    # Llenar el resto de la matriz
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Costo de sustitución: 0 si coinciden, 1 si no
            costo_sustitucion = 0 if s1[i - 1] == s2[j - 1] else 1
            
            # El valor es el mínimo de tres opciones:
            # 1. Borrado (D[i-1, j] + 1)
            # 2. Inserción (D[i, j-1] + 1)
            # 3. Sustitución (D[i-1, j-1] + costo_sustitucion)
            D[i, j] = min(
                D[i - 1, j] + 1,
                D[i, j - 1] + 1,
                D[i - 1, j - 1] + costo_sustitucion
            )

    # La distancia de Levenshtein es el valor en la esquina inferior derecha
    return D[m, n]

# --- 2. Experimentos y Análisis ---

MAX_MUTATIONS = 300
N_ITERATIONS_C = 10000

# Preparación para graficar
M_values = range(MAX_MUTATIONS + 1) # Valores de M de 0 a 300

# a) Muta una secuencia, compara con la original
def experimento_a():
    """Generar secuencia, aplicar M mutaciones y medir D(S_final, S_inicial)"""
    S_inicial = generar_secuencia_aleatoria()
    S_actual = S_inicial[:]
    D_values = []
    
    for M in M_values:
        if M > 0:
            S_actual = aplicar_mutacion(S_actual)
        
        D = distancia_levenshtein(S_actual, S_inicial)
        D_values.append(D)
        
    return D_values

# b) Genera 2 copias, muta independientemente y mide D'(S_a, S_b)
def experimento_b():
    """Generar 2 copias, aplicar M mutaciones INDEPENDIENTES y medir D'(S_a, S_b)"""
    S_0 = generar_secuencia_aleatoria()
    S_a = S_0[:]
    S_b = S_0[:]
    Dp_values = []
    
    for M in M_values:
        D_prime = distancia_levenshtein(S_a, S_b)
        Dp_values.append(D_prime)
        
        # Aplicar UNA mutación a cada copia de forma INDEPENDIENTE
        if M < MAX_MUTATIONS:
            S_a = aplicar_mutacion(S_a)
            S_b = aplicar_mutacion(S_b)
            
    return Dp_values

# c) Distribución de Distancias Aleatorias
def experimento_c():
    """Genera 10.000 pares aleatorios y mide la distribución de D."""
    D_aleatorias = []
    for _ in range(N_ITERATIONS_C):
        S1 = generar_secuencia_aleatoria()
        S2 = generar_secuencia_aleatoria()
        D_aleatorias.append(distancia_levenshtein(S1, S2))
        
    D_aleatorias = np.array(D_aleatorias)
    media = np.mean(D_aleatorias)
    std_dev = np.std(D_aleatorias)
    
    # Graficar el histograma
    plt.figure(figsize=(8, 5))
    plt.hist(D_aleatorias, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(media, color='red', linestyle='dashed', linewidth=1, label=f'Media: {media:.2f}')
    plt.title(f'c) Histograma de Distancia de Levenshtein (N={N_ITERATIONS_C} Pares Aleatorios)\nmedia={media:.2f}, std={std_dev:.2f}')
    plt.xlabel('Distancia de Levenshtein')
    plt.ylabel('Frecuencia')
    plt.legend()
    plt.grid(axis='y', alpha=0.5)
    plt.show()
    
    return media, std_dev

# --- Código Principal de Ejecución ---

if __name__ == "__main__":
    print(f"--- Pregunta 3: Parámetros del Experimento ---")
    print(f"Alfabeto de Bases: {ALPHABET}")
    print(f"Longitud de Secuencia: {SEQ_LENGTH} bases")
    print(f"Máximo de Mutaciones (M): {MAX_MUTATIONS}")
    print("\n" + "="*50 + "\n")

    # a) y b) Ejecutar simulaciones y graficar
    D_a = experimento_a()
    Dp_b = experimento_b()

    plt.figure(figsize=(10, 6))
    plt.plot(M_values, D_a, label='a) $D(S_{final}, S_{inicial})$: Mutación vs Original', color='blue')
    plt.plot(M_values, Dp_b, label="b) $D'(S_a, S_b)$: Divergencia entre 2 secuencias", color='orange')
    plt.title('a) y b) Distancia de Levenshtein vs. Número de Mutaciones (M)')
    plt.xlabel('Número de Mutaciones (M)')
    plt.ylabel('Distancia de Levenshtein (D o D\')')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

    # c) Ejecutar cálculo de distribución aleatoria
    print(f"--- c) Análisis de Distribución Aleatoria ({N_ITERATIONS_C} pares) ---")
    media_aleatoria, std_dev_aleatoria = experimento_c()
    print(f"Media ($\mu$) de Distancia Aleatoria: {media_aleatoria:.2f}")
    print(f"Desviación Estándar ($\sigma$) de Distancia Aleatoria: {std_dev_aleatoria:.2f}")
    
    # d) Análisis de la Twilight Zone
    print("\n" + "="*50 + "\n")
    print("--- d) Determinación de la Twilight Zone ---")
    
    # Buscamos el punto donde la divergencia (Dp_b) alcanza la media aleatoria (media_aleatoria)
    try:
        # Se busca el primer índice donde Dp_b (distancia entre linajes) es >= media_aleatoria
        # Se convierte a un arreglo numpy para usar funciones de filtrado.
        Dp_b_np = np.array(Dp_b)
        
        # Filtra los M que cumplen la condición
        M_values = np.array(range(MAX_MUTATIONS + 1))
        M_twilight = M_values[Dp_b_np >= media_aleatoria][0]
        
        print(f"El parentesco se vuelve indetectable o entra en la 'twilight zone' cuando la distancia (D')")
        print(f"entre los linajes mutantes se acerca a la distancia promedio de secuencias aleatorias ({media_aleatoria:.2f}).")
        print(f"Esto ocurre aproximadamente a partir de M = {M_twilight} mutaciones.")
        
    except IndexError:
        print(f"Nota: La divergencia D'(M) no alcanzó la media de distancia aleatoria ({media_aleatoria:.2f}) dentro de {MAX_MUTATIONS} mutaciones.")
        print("Esto puede pasar si la tasa de cambio es baja respecto a la longitud de la secuencia.")