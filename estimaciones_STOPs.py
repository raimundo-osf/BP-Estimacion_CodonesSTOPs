from collections import Counter

def contar_codones_stop(secuencia):
    codones_stop = ["TAG", "TGA", "TAA"]
    contador = Counter()
    for i in range(len(secuencia) - 2):  # Ventana deslizante de longitud 3
        codon = secuencia[i:i+3]
        if codon in codones_stop:
            contador[codon] += 1
    return contador

def frecuencia_nts(secuencia):
    total = len(secuencia)
    probabilidad_nts = Counter(secuencia)
    return {nt: probabilidad_nts[nt]/total for nt in "ATCG"}

def frecuencia_dupletes(secuencia):
    probabilidad_dupletes = Counter()
    total = len(secuencia) - 1 
    for i in range(total): # Ventana deslizante de longitud 2
        duplete = secuencia[i:i+2] 
        probabilidad_dupletes[duplete] += 1
    return {duplete: count/total for duplete, count in probabilidad_dupletes.items()}

def frecuencia_tripletes(secuencia):
    probabilidad_tripletes = Counter()
    total = len(secuencia) - 2 
    for i in range(total): # Ventana deslizante de longitud 3
        triplete = secuencia[i:i+3]
        probabilidad_tripletes[triplete] += 1
    return {triplete: count/total for triplete, count in probabilidad_tripletes.items()}

# Estimacion 1 | Sin información alguna
def estimacion_1(largo_genoma, largo_codon):
    return (largo_genoma - largo_codon + 1) / 4**largo_codon

# Estimacion 2 | Agregando la información de la frecuencia de As, Ts, Cs y Gs
def estimacion_2(largo_genoma, codon, probabilidad_nts):
    P = probabilidad_nts[codon[0]] * probabilidad_nts[codon[1]] * probabilidad_nts[codon[2]]
    return (largo_genoma - len(codon) + 1) * P

# Estimacion 3 | Agregando la información del nucleótido anterior (Markov) | Leyendo del primero hacia el ultimo
def estimacion_3(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes):
    P = probabilidad_nts[codon[0]]
    P *= probabilidad_dupletes[codon[0:2]] / probabilidad_nts[codon[0]]
    P *= probabilidad_dupletes[codon[1:3]] / probabilidad_nts[codon[1]]
    return (largo_genoma - len(codon) + 1) * P

# Estimacion 4 | Agregando la información del nucleótido anterior (Markov) | Leyendo del ultimo hacia el primero
def estimacion_4(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes):
    P = probabilidad_nts[codon[2]]
    P *= probabilidad_dupletes[codon[1:3]] / probabilidad_nts[codon[2]]
    P *= probabilidad_dupletes[codon[0:2]] / probabilidad_nts[codon[1]]
    return (largo_genoma - len(codon) + 1) * P

# Estimacion 5 | Agregando la información de todos los nucleótidos anteriores | Leyendo del primero hacia el ultimo
def estimacion_5(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes, probabilidad_tripletes):
    P = probabilidad_nts[codon[0]]
    P *= probabilidad_dupletes[codon[0:2]] / probabilidad_nts[codon[0]]
    P *= probabilidad_tripletes[codon[0:3]] / probabilidad_dupletes[codon[0:2]]
    return (largo_genoma - len(codon) + 1) * P

# Estimacion 6 | Agregando la información de todos los nucleótidos anteriores | Leyendo del ultimo hacia el primero
def estimacion_6(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes, probabilidad_tripletes):
    P = probabilidad_nts[codon[2]]
    P *= probabilidad_dupletes[codon[1:3]] / probabilidad_nts[codon[2]]
    P *= probabilidad_tripletes[codon[0:3]] / probabilidad_dupletes[codon[1:3]]
    return (largo_genoma - len(codon) + 1) * P
