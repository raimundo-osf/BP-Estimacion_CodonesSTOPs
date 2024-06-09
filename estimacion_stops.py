import time
from collections import Counter
import pandas as pd

# Función para cargar el genoma desde un archivo fasta
def cargar_genoma(archivo_fasta):
    with open(archivo_fasta, 'r') as f:
        secuencia = ''.join(line.strip() for line in f if not line.startswith('>'))
    return secuencia

# Función para contar los codones de STOP usando una ventana deslizante
def contar_codones_stop(secuencia):
    codones_stop = ["TAG", "TGA", "TAA"]
    contador = Counter()
    for i in range(len(secuencia) - 2):  # Ventana deslizante de longitud 3
        codon = secuencia[i:i+3]
        if codon in codones_stop:
            contador[codon] += 1
    return contador

# Función para calcular la frecuencia de nucleótidos
def calcular_frecuencias(secuencia):
    total = len(secuencia)
    frecuencias = Counter(secuencia)
    return {nt: frecuencias[nt]/total for nt in "ATCG"}

# Estimaciones
def estimar_stop_1(largo_genoma, largo_patron):
    return (largo_genoma - largo_patron + 1) / 4**largo_patron

def estimar_stop_2(largo_genoma, largo_patron, frecuencias):
    P = frecuencias['A'] * frecuencias['T'] * frecuencias['C'] * frecuencias['G']
    return (largo_genoma - largo_patron + 1) * P

# Función principal
def main():
    archivo_fasta = "/home/raimundoosf/Desktop/MMSB/GCF_000006945.2_ASM694v2_genomic.fna"  # Cambia esto por la ruta a tu archivo FASTA
    secuencia = cargar_genoma(archivo_fasta)
    largo_genoma = len(secuencia)
    frecuencias = calcular_frecuencias(secuencia)
    print("frecuencia relativa de nucleotidos:", frecuencias)
    num_stops_real = contar_codones_stop(secuencia)
    print("codones cant:", num_stops_real)
    
    # Estimaciones
    t0 = time.time()
    estimacion1 = estimar_stop_1(largo_genoma, 3)
    t1 = time.time()
    error_relativo1 = abs(sum(num_stops_real.values()) - estimacion1) / sum(num_stops_real.values())
    tiempo1 = t1 - t0
    
    t2 = time.time()
    estimacion2 = estimar_stop_2(largo_genoma, 3, frecuencias)
    t3 = time.time()
    error_relativo2 = abs(sum(num_stops_real.values()) - estimacion2) / sum(num_stops_real.values())
    tiempo2 = t3 - t2

    # Resultados
    resultados = []

    for codon in ["TAG", "TGA", "TAA"]:
        resultados.append([
            "Especie Desconocida", frecuencias['G'] + frecuencias['C'], codon, num_stops_real[codon], 
            estimacion1, error_relativo1, tiempo1, 
            estimacion2, error_relativo2, tiempo2
        ])

    # Crear el DataFrame
    columnas = ["Nombre de la especie", "%G+C", "Codón", "Número de STOPs real", 
                "Estimación 1", "Error relativo 1", "Tiempo de cálculo 1", 
                "Estimación 2", "Error relativo 2", "Tiempo de cálculo 2"]
    df = pd.DataFrame(resultados, columns=columnas)

    # Guardar el DataFrame en un archivo Excel
    df.to_excel("resultados_stop_codons.xlsx", index=False)
    print("Resultados guardados en 'resultados_stop_codons.xlsx'")

if __name__ == "__main__":
    main()
