from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy
import random  # Importa random para usar en mutaciones
import matplotlib.pyplot as plt  # Importar matplotlib para visualización

# Inicializa la población
poblacion = []
path = "C:\\secuenciasBFOA\\multiFasta.fasta"  # Asegúrate de usar \\ en las rutas
numeroDeBacterias = 20  # Aumentamos el número de bacterias
numRandomBacteria = 2   # Insertamos más bacterias aleatorias
iteraciones = 15        # Reducimos ligeramente el número de iteraciones
tumbo = 2               # Incrementamos el tumbo para más variabilidad 
nado = 4                # Incrementamos el nado para más exploración
chemio = chemiotaxis()
veryBest = bacteria(path)
tempBacteria = bacteria(path)
original = bacteria(path)
globalNFE = 0

# Listas para almacenar el progreso
fitness_history = []  # Para almacenar el mejor fitness en cada iteración
nfe_history = []      # Para almacenar el NFE global en cada iteración

# Parámetros iniciales
dAttr = 0.15  # Aumentamos ligeramente la atracción
wAttr = 0.3   # Aumentamos la fuerza de atracción
hRep = 0.05   # Aumentamos la repulsión
wRep = 15     # Aumentamos la fuerza de repulsión
mutation_rate = 0.1  # Tasa de mutación inicial
umbral_convergencia = 0.01  # Umbral para verificar convergencia
tasa_diversidad = 0.1  # Probabilidad de introducir diversidad
tasa_nueva = 0.2  # Proporción de nuevas bacterias a introducir
min_ciclo = 2  # Ciclo mínimo de repulsión

# Función para clonar la mejor bacteria
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

# Valida que las secuencias originales no tengan gaps
def validaSecuencias(path, veryBest):
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

# Ajusta los parámetros dinámicamente según el progreso del fitness
def ajustaParametros(prevFitness, currentFitness):
    global dAttr, wAttr, hRep, wRep
    mejoraFitness = currentFitness - prevFitness
    if mejoraFitness < umbral_convergencia:  # Umbral bajo de mejora para ajustar
        dAttr *= 1.05
        wAttr *= 1.05
        hRep = dAttr
        wRep *= 1.05
    else:
        dAttr *= 0.95
        wAttr *= 0.95
        hRep = dAttr
        wRep *= 0.95

# Genera la población inicial
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))

# Función para mutar una bacteria
def mutar_bacteria(bacteria):
    for i in range(len(bacteria.matrix.seqs)):
        if random.random() < mutation_rate:  # Aplica mutación según la tasa de mutación
            seq = list(bacteria.matrix.seqs[i])
            position = random.randint(0, len(seq) - 1)
            seq[position] = random.choice(['A', 'T', 'C', 'G'])  # Asumiendo secuencias de ADN
            bacteria.matrix.seqs[i] = ''.join(seq)

# Función para introducir diversidad en la población
def introducir_diversidad(bacterias, tasa_diversidad):
    if random.random() < tasa_diversidad:
        num_nuevas_bacterias = int(len(bacterias) * tasa_nueva)
        nuevas_bacterias = [bacteria(path) for _ in range(num_nuevas_bacterias)]
        bacterias.extend(nuevas_bacterias)

# Algoritmo principal
prevBestFitness = float('-inf')
for iteracion in range(iteraciones):
    for b in poblacion:
        b.tumboNado(tumbo)  # Movimiento aleatorio
        b.autoEvalua()      # Evaluación del fitness
        
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE  # Acumula NFE global

    best = max(poblacion, key=lambda x: x.fitness)
    
    if best.fitness > prevBestFitness:
        clonaBest(veryBest, best)
        ajustaParametros(prevBestFitness, best.fitness)
        prevBestFitness = best.fitness
    else:
        print("Convergencia alcanzada. Deteniendo ajustes de parámetros.")
        break  # Salir del bucle si la convergencia se alcanza

    print("Iteración:", iteracion, "Fitness:", veryBest.fitness, "NFE:", globalNFE)

    # Almacena el progreso
    fitness_history.append(veryBest.fitness)
    nfe_history.append(globalNFE)

    # Aplicar mutaciones a la población
    for b in poblacion:
        mutar_bacteria(b)

    # Introducir diversidad en la población
    introducir_diversidad(poblacion, tasa_diversidad)

    # Funciones de eliminación y clonación
    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
    print("Población:", len(poblacion))

# Muestra el genoma de la mejor bacteria y valida las secuencias
veryBest.showGenome()
validaSecuencias(path, veryBest)

# Visualización de los resultados
plt.figure(figsize=(12, 6))

# Graficar fitness
plt.subplot(2, 1, 1)
plt.plot(fitness_history, marker='o', color='blue', label='Fitness')
plt.title('Progreso del Fitness a lo largo de las Iteraciones')
plt.xlabel('Iteraciones')
plt.ylabel('Fitness')
plt.grid()
plt.legend()

# Graficar NFE
plt.subplot(2, 1, 2)
plt.plot(nfe_history, marker='s', color='red', label='NFE')
plt.title('Progreso de NFE a lo largo de las Iteraciones')
plt.xlabel('Iteraciones')
plt.ylabel('NFE')
plt.grid()
plt.legend()

plt.tight_layout()  # Ajustar el layout
plt.show()  # Mostrar los gráficos
