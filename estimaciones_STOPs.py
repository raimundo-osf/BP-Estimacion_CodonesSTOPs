from collections import Counter

def contar_codones_stop(secuencia):
    codones_stop = ["TAG", "TGA", "TAA"]
    contador = Counter()
    for i in range(len(secuencia) - 2):  # Ventana deslizante de longitud 3
        codon = secuencia[i:i+3]
        if codon in codones_stop:
            contador[codon] += 1
    return contador

def conteo_nts(secuencia):
    nts = Counter(secuencia)
    return {nt: nts[nt] for nt in "ATCG"}

def conteo_dupletes(secuencia):
    dupletes = Counter()
    total = len(secuencia) - 1 
    for i in range(total): # Ventana deslizante de longitud 2
        duplete = secuencia[i:i+2] 
        dupletes[duplete] += 1
    return {duplete: count for duplete, count in dupletes.items()}

def conteo_tripletes(secuencia):
    tripletes = Counter()
    total = len(secuencia) - 2 
    for i in range(total): # Ventana deslizante de longitud 3
        triplete = secuencia[i:i+3]
        tripletes[triplete] += 1
    return {triplete: count for triplete, count in tripletes.items()}

# Método de Estimación Sin información alguna
def metodo_sin_informacion(largo_genoma):
    return (largo_genoma - 2) / 4**3

# Método de Estimación con Frecuencias Independientes
def metodo_frecuencias_independientes(largo_genoma, codon, cantidad_nts):
    P = (cantidad_nts[codon[0]] * cantidad_nts[codon[1]] * cantidad_nts[codon[2]]) / largo_genoma ** 3
    return (largo_genoma - 2) * P 

# Método de Estimación con Dependencia Markoviana (Estimaciones 3 y 4)
def metodo_markoviano(largo_genoma, codon, cantidad_nts, cantidad_dupletes):
    P = cantidad_nts[codon[0]] / largo_genoma
    P *= cantidad_dupletes[codon[0:2]] / cantidad_nts[codon[0]]
    P *= cantidad_dupletes[codon[1:3]] / cantidad_nts[codon[1]]
    return (largo_genoma - 2) * P

# Método de Estimación con Dependencia Total (Estimaciones 5 y 6)
def metodo_dependencia_total(largo_genoma, codon, cantidad_nts, cantidad_dupletes, cantidad_tripletes):
    P = cantidad_nts[codon[0]] / largo_genoma
    P *= cantidad_dupletes[codon[0:2]] / cantidad_nts[codon[0]]
    P *= cantidad_tripletes[codon[0:3]] / cantidad_dupletes[codon[0:2]]
    return (largo_genoma - 2) * P

