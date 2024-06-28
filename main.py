import time
import pandas as pd
from estimaciones_STOPs import (
    contar_codones_stop, 
    conteo_nts, 
    conteo_dupletes, 
    conteo_tripletes, 
    metodo_sin_informacion, 
    metodo_frecuencias_independientes, 
    metodo_markoviano, 
    metodo_dependencia_total
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

        ''' Conteo de nts, dupletes y tripletes | Leyendo secuencia de izquierda a derecha '''

        # Calcular cantidad de cada nucleótido
        t_nts_start = time.time()
        cantidad_nts = conteo_nts(secuencia)
        t_nts_end = time.time()
        t_calculo_conteo_nts = t_nts_end - t_nts_start

        # Calcular cantidad de cada duplete
        t_dupletes_start = time.time()
        cantidad_dupletes = conteo_dupletes(secuencia)
        t_dupletes_end = time.time()
        t_calculo_conteo_dupletes = t_dupletes_end - t_dupletes_start

        # Calcular cantidad de cada triplete
        t_tripletes_start = time.time()
        cantidad_tripletes = conteo_tripletes(secuencia)
        t_tripletes_end = time.time()
        t_calculo_conteo_tripletes = t_tripletes_end - t_tripletes_start

        ''' Conteo de nts, dupletes y tripletes | Leyendo secuencia de derecha a izquierda'''

        # Se invierte secuencia, para conteo de derecha a izquierda de secuencia original
        secuencia_inv = secuencia[::-1]

        # Proporcion de nts en secuencia invertida debe ser igual a la proporcion de nts en secuencia normal (lectura de izq. a der.), tal que no se realiza conteo de nts nuevamente.

        # En realidad, la frecuencia de los inversos (dupletes y tripletes) debe ser la msima (Ej: #TG = #GT), pero para evitar confusiones al realizar los algoritmos asociados a las estimaciones, se decidio calcular nuevamente la cantidad de dupletes y tripletes para la secuencia inversa (lectura de derecha a izquierda)

        # Calcular cantidad de cada duplete | Secuencia inv.
        t_dupletes_start = time.time()
        cantidad_dupletes_sec_inv = conteo_dupletes(secuencia_inv)
        t_dupletes_end = time.time()
        t_calculo_conteo_dupletes_sec_inv = t_dupletes_end - t_dupletes_start

        # Calcular cantidad de cada triplete | Secuencia inv.
        t_tripletes_start = time.time()
        cantidad_tripletes_sec_inv = conteo_tripletes(secuencia_inv)
        t_tripletes_end = time.time()
        t_calculo_conteo_tripletes_sec_inv = t_tripletes_end - t_tripletes_start

        # Calcular cantidad de real de codones
        # Nuevamente, cantidad de codones STOPs debe ser igual a la de sus inversos. Asi, solo se cuentan los codones STOPs "TAG", "TGA", "TAA" y no sus inversos ("GAT", "AGT", "AAT")
        num_stops_real = contar_codones_stop(secuencia)

        ''' Estimaciones '''

        # Método De Estimación 1 (misma para cada codón)
        t0 = time.time()
        estimacion1 = metodo_sin_informacion(largo_genoma)
        t1 = time.time()
        tiempo1 = t1 - t0
        
        # Método De Estimación 2 para cada codón STOP
        estimaciones2 = {}
        tiempos2 = {}
        errores_relativos2 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t2 = time.time()
            estimacion2 = metodo_frecuencias_independientes(largo_genoma, codon, cantidad_nts)
            t3 = time.time()
            estimaciones2[codon] = estimacion2
            tiempos2[codon] = t3 - t2 + t_calculo_conteo_nts
            errores_relativos2[codon] = abs(num_stops_real[codon] - estimacion2) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Método De Estimación 3 para cada codón STOP
        estimaciones3 = {}
        tiempos3 = {}
        errores_relativos3 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t4 = time.time()
            estimacion3 = metodo_markoviano(largo_genoma, codon, cantidad_nts, cantidad_dupletes)
            t5 = time.time()
            estimaciones3[codon] = estimacion3
            tiempos3[codon] = t5 - t4 + t_calculo_conteo_nts + t_calculo_conteo_dupletes
            errores_relativos3[codon] = abs(num_stops_real[codon] - estimacion3) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Método De Estimación 4 para cada codón STOP
        estimaciones4 = {}
        tiempos4 = {}
        errores_relativos4 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t6 = time.time()
            codon_inverso = codon[::-1]
            # Estimación se da sobre codon inverso -lectura desde ultimo nt hacia 1er nt-
            estimacion4 = metodo_markoviano(largo_genoma, codon_inverso, cantidad_nts, cantidad_dupletes_sec_inv)
            t7 = time.time()
            estimaciones4[codon] = estimacion4
            tiempos4[codon] = t7 - t6 + t_calculo_conteo_nts + t_calculo_conteo_dupletes_sec_inv
            errores_relativos4[codon] = abs(num_stops_real[codon] - estimacion4) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
        
        # Método De Estimación 5 para cada codón STOP
        estimaciones5 = {}
        tiempos5 = {}
        errores_relativos5 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t8 = time.time()
            estimacion5 = metodo_dependencia_total(largo_genoma, codon, cantidad_nts, cantidad_dupletes, cantidad_tripletes)
            t9 = time.time()
            estimaciones5[codon] = estimacion5
            tiempos5[codon] = t9 - t8 + t_calculo_conteo_nts + t_calculo_conteo_dupletes + t_calculo_conteo_tripletes
            errores_relativos5[codon] = abs(num_stops_real[codon] - estimacion5) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')

        # Método De Estimación 6 para cada codón STOP
        estimaciones6 = {}
        tiempos6 = {}
        errores_relativos6 = {}
        for codon in ["TAG", "TGA", "TAA"]:
            t10 = time.time()
            # Estimación se da sobre codon inverso -lectura desde ultimo nt hacia 1er nt-
            codon_inverso = codon[::-1]
            estimacion6 = metodo_dependencia_total(largo_genoma, codon_inverso, cantidad_nts, cantidad_dupletes_sec_inv, cantidad_tripletes_sec_inv)
            t11 = time.time()
            estimaciones6[codon] = estimacion6
            tiempos6[codon] = t11 - t10 + t_calculo_conteo_nts + t_calculo_conteo_dupletes_sec_inv + t_calculo_conteo_tripletes_sec_inv
            errores_relativos6[codon] = abs(num_stops_real[codon] - estimacion6) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
        
        # Resultados
        for codon in ["TAG", "TGA", "TAA"]:
            error_relativo1 = abs(num_stops_real[codon] - estimacion1) / num_stops_real[codon] if num_stops_real[codon] != 0 else float('inf')
            porcentaje_GC = (cantidad_nts['G'] + cantidad_nts['C']) / largo_genoma
            resultados.append([
                nombre_especie, porcentaje_GC, codon, num_stops_real[codon],
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
