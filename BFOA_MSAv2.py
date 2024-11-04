from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy
import random
import matplotlib.pyplot as plt

# Inicializa la población
poblacion = []
path = "C:\\secuenciasBFOA\\multiFasta.fasta"
numeroDeBacterias = 20
iteraciones = 15
tumbo = 2
chemio = chemiotaxis()
veryBest = bacteria(path)
original = bacteria(path)
globalNFE = 0

# Listas para almacenar el progreso
fitness_history = []
nfe_history = []

# Parámetros iniciales
dAttr = 0.15
wAttr = 0.3
hRep = 0.05
wRep = 15
mutation_rate = 0.2  # Aumenta la tasa de mutación
umbral_convergencia = 0.01

# Clonar la mejor bacteria
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

# Ajusta los parámetros dinámicamente según el progreso del fitness
def ajustaParametros(prevFitness, currentFitness):
    global dAttr, wAttr, hRep, wRep
    mejoraFitness = currentFitness - prevFitness
    if mejoraFitness < umbral_convergencia:
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
        if random.random() < mutation_rate:
            seq = list(bacteria.matrix.seqs[i])
            position = random.randint(0, len(seq) - 1)
            seq[position] = random.choice(['A', 'T', 'C', 'G'])  # Suponiendo ADN, ajusta según sea necesario
            bacteria.matrix.seqs[i] = ''.join(seq)

# Algoritmo principal
prevBestFitness = float('-inf')
for iteracion in range(iteraciones):
    for b in poblacion:
        b.tumboNado(tumbo)
        b.autoEvalua()

    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE
    best = max(poblacion, key=lambda x: x.fitness)

    if best.fitness > prevBestFitness:
        clonaBest(veryBest, best)
        ajustaParametros(prevBestFitness, best.fitness)
        prevBestFitness = best.fitness
    else:
        print("Convergencia alcanzada. Deteniendo ajustes de parámetros.")
        break

    print("Iteración:", iteracion, "Fitness:", veryBest.fitness, "NFE:", globalNFE)
    
    # Almacena el progreso
    fitness_history.append(veryBest.fitness)
    nfe_history.append(globalNFE)
    
    # Aplicar mutaciones a la población
    for b in poblacion:
        mutar_bacteria(b)

    # Funciones de eliminación y clonación
    chemio.eliminarClonar(path, poblacion)

# Muestra el genoma de la mejor bacteria
veryBest.showGenome()

# Visualización de los resultados
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(fitness_history, marker='o', color='blue', label='Fitness')
plt.title('Progreso del Fitness a lo largo de las Iteraciones')
plt.xlabel('Iteraciones')
plt.ylabel('Fitness')
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(nfe_history, marker='s', color='red', label='NFE')
plt.title('Progreso de NFE a lo largo de las Iteraciones')
plt.xlabel('Iteraciones')
plt.ylabel('NFE')
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
