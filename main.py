import time
import pandas as pd
from estimaciones_STOPs import (
    contar_codones_stop, 
    frecuencia_nts, 
    frecuencia_dupletes, 
    frecuencia_tripletes, 
    estimacion_1, 
    estimacion_2, 
    estimacion_3, 
    estimacion_4, 
    estimacion_5, 
    estimacion_6
)

def cargar_multifasta(archivo_fasta):
    secuencias = {}
    contador_especies = {}
    with open(archivo_fasta, 'r') as f:
        nombre_especie = ""
        secuencia = ""
        for line in f:
            if line.startswith('>'):
                if secuencia:
                    # Si el nombre de la especie ya existe, agregar un sufijo numérico
                    if nombre_especie in secuencias:
                        if nombre_especie not in contador_especies:
                            contador_especies[nombre_especie] = 1
                        contador_especies[nombre_especie] += 1
                        nombre_especie_con_sufijo = f"{nombre_especie}_{contador_especies[nombre_especie]}"
                        secuencias[nombre_especie_con_sufijo] = secuencia
                    else:
                        secuencias[nombre_especie] = secuencia
                nombre_especie = " ".join(line.strip().split()[1:3])
                secuencia = ""
            else:
                secuencia += line.strip()
        if secuencia:
            if nombre_especie in secuencias:
                if nombre_especie not in contador_especies:
                    contador_especies[nombre_especie] = 1
                contador_especies[nombre_especie] += 1
                nombre_especie_con_sufijo = f"{nombre_especie}_{contador_especies[nombre_especie]}"
                secuencias[nombre_especie_con_sufijo] = secuencia
            else:
                secuencias[nombre_especie] = secuencia
    return secuencias

def main():
    archivo_multifasta = "/home/raimundoosf/Desktop/MMSB/BounsPoint1/multifasta_genomas.fasta" 
    secuencias = cargar_multifasta(archivo_multifasta)
    print(len(secuencias))
    resultados = []

    for nombre_especie, secuencia in secuencias.items():
        largo_genoma = len(secuencia)
        probabilidad_nts = frecuencia_nts(secuencia)
        probabilidad_dupletes = frecuencia_dupletes(secuencia)
        probabilidad_tripletes = frecuencia_tripletes(secuencia)
        num_stops_real = contar_codones_stop(secuencia)

        # Estimación 1 (misma para cada codón)
        t0 = time.time()
        estimacion1 = estimacion_1(largo_genoma, 3)
        t1 = time.time()
        tiempo1 = t1 - t0
        
        # Estimación 2 para cada codón STOP
        estimaciones2 = {}
        tiempos2 = {}
        errores_relativos2 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t2 = time.time()
            estimacion2 = estimacion_2(largo_genoma, codon, probabilidad_nts)
            t3 = time.time()
            estimaciones2[codon] = estimacion2
            tiempos2[codon] = t3 - t2
            errores_relativos2[codon] = abs(num_stops_real[codon] - estimacion2) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Estimación 3 para cada codón STOP
        estimaciones3 = {}
        tiempos3 = {}
        errores_relativos3 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t4 = time.time()
            estimacion3 = estimacion_3(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes)
            t5 = time.time()
            estimaciones3[codon] = estimacion3
            tiempos3[codon] = t5 - t4
            errores_relativos3[codon] = abs(num_stops_real[codon] - estimacion3) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Estimación 4 para cada codón STOP
        estimaciones4 = {}
        tiempos4 = {}
        errores_relativos4 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t6 = time.time()
            estimacion4 = estimacion_4(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes)
            t7 = time.time()
            estimaciones4[codon] = estimacion4
            tiempos4[codon] = t7 - t6
            errores_relativos4[codon] = abs(num_stops_real[codon] - estimacion4) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
        
        # Estimación 5 para cada codón STOP
        estimaciones5 = {}
        tiempos5 = {}
        errores_relativos5 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t8 = time.time()
            estimacion5 = estimacion_5(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes, probabilidad_tripletes)
            t9 = time.time()
            estimaciones5[codon] = estimacion5
            tiempos5[codon] = t9 - t8
            errores_relativos5[codon] = abs(num_stops_real[codon] - estimacion5) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Estimación 6 para cada codón STOP
        estimaciones6 = {}
        tiempos6 = {}
        errores_relativos6 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t10 = time.time()
            estimacion6 = estimacion_6(largo_genoma, codon, probabilidad_nts, probabilidad_dupletes, probabilidad_tripletes)
            t11 = time.time()
            estimaciones6[codon] = estimacion6
            tiempos6[codon] = t11 - t10
            errores_relativos6[codon] = abs(num_stops_real[codon] - estimacion6) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
        
        # Resultados
        for codon in ["TAG", "TGA", "TAA"]:
            error_relativo1 = abs(num_stops_real[codon] - estimacion1) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
            resultados.append([
                nombre_especie, probabilidad_nts['G'] + probabilidad_nts['C'], codon, num_stops_real[codon], 
                estimacion1, error_relativo1, tiempo1, 
                estimaciones2[codon], errores_relativos2[codon], tiempos2[codon],
                estimaciones3[codon], errores_relativos3[codon], tiempos3[codon],
                estimaciones4[codon], errores_relativos4[codon], tiempos4[codon],
                estimaciones5[codon], errores_relativos5[codon], tiempos5[codon],
                estimaciones6[codon], errores_relativos6[codon], tiempos6[codon]
            ])
            
    # Crear el DataFrame
    columnas = ["Nombre de la especie", "%G+C", "Codón", "Número de STOPs real", 
                "Estimación 1", "Error relativo 1", "Tiempo de cálculo 1", 
                "Estimación 2", "Error relativo 2", "Tiempo de cálculo 2",
                "Estimación 3", "Error relativo 3", "Tiempo de cálculo 3",
                "Estimación 4", "Error relativo 4", "Tiempo de cálculo 4",
                "Estimación 5", "Error relativo 5", "Tiempo de cálculo 5",
                "Estimación 6", "Error relativo 6", "Tiempo de cálculo 6"]
    df = pd.DataFrame(resultados, columns=columnas)

    # Calcular promedios y desviaciones estándar
    promedios = df.iloc[:, 4:].mean()
    desviaciones = df.iloc[:, 4:].std()

    # Agregar promedios y desviaciones al DataFrame
    promedios_row = ["Promedios", "", "", ""] + list(promedios)
    desviaciones_row = ["Desviaciones estándar", "", "", ""] + list(desviaciones)
    df.loc[len(df)] = promedios_row
    df.loc[len(df)] = desviaciones_row

    # Guardar el DataFrame en un archivo Excel
    df.to_excel("resultados_stop_codons.xlsx", index=False)
    print("Resultados guardados en 'resultados_stop_codons.xlsx'")
    
if __name__ == "__main__":
    main()
