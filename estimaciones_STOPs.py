from collections import Counter

def contar_codones_stop(secuencia):
    codones_stop = ["TAG", "TGA", "TAA", "GAT", "AGT", "AAT"] #incluye codones inversos para estimaciones 4 y 6
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

# Estimacion 1 | Sin información alguna
def estimacion_1(largo_genoma):
    return (largo_genoma - 2) / 4**3

# Estimacion 2 | Agregando la información de la frecuencia de As, Ts, Cs y Gs
def estimacion_2(largo_genoma, codon, cantidad_nts):
    P = (cantidad_nts[codon[0]] * cantidad_nts[codon[1]] * cantidad_nts[codon[2]]) / largo_genoma ** 3
    return (largo_genoma - 2) * P 

# Estimacion 3 | Agregando la información del nucleótido anterior (Markov) | Leyendo del primero hacia el ultimo
def estimacion_3(largo_genoma, codon, cantidad_nts, cantidad_dupletes):
    P = cantidad_nts[codon[0]] / largo_genoma
    P *= cantidad_dupletes[codon[0:2]] / cantidad_nts[codon[0]]
    P *= cantidad_dupletes[codon[1:3]] / cantidad_nts[codon[1]]
    return (largo_genoma - 2) * P

# Estimacion 4 | Agregando la información del nucleótido anterior (Markov) | Leyendo del ultimo hacia el primero
def estimacion_4(largo_genoma, codon, cantidad_nts, cantidad_dupletes):
    codon_inverso = codon[::-1]
    P = cantidad_nts[codon_inverso[0]] / largo_genoma
    P *= cantidad_dupletes[codon_inverso[0:2]] / cantidad_nts[codon_inverso[0]]
    P *= cantidad_dupletes[codon_inverso[1:3]] / cantidad_nts[codon_inverso[1]]
    return (largo_genoma - 2) * P

# Estimacion 5 | Agregando la información de todos los nucleótidos anteriores | Leyendo del primero hacia el ultimo
def estimacion_5(largo_genoma, codon, cantidad_nts, cantidad_dupletes, cantidad_tripletes):
    P = cantidad_nts[codon[0]] / largo_genoma
    P *= cantidad_dupletes[codon[0:2]] / cantidad_nts[codon[0]]
    P *= cantidad_tripletes[codon[0:3]] / cantidad_dupletes[codon[0:2]]
    return (largo_genoma - 2) * P

# Estimacion 6 | Agregando la información de todos los nucleótidos anteriores | Leyendo del ultimo hacia el primero
def estimacion_6(largo_genoma, codon, cantidad_nts, cantidad_dupletes, cantidad_tripletes):
    codon_inverso = codon[::-1]
    P = cantidad_nts[codon_inverso[0]] / largo_genoma
    P *= cantidad_dupletes[codon_inverso[0:2]] / cantidad_nts[codon_inverso[0]]
    P *= cantidad_tripletes[codon_inverso[0:3]] / cantidad_dupletes[codon_inverso[0:2]]
    return (largo_genoma - 2) * P
