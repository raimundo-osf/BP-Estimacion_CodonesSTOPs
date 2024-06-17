import random
import time
from Bio import Entrez

# Correo electrónico
Entrez.email = "raioliva21@gmail.com"

# Función para obtener IDs de genomas
def obtener_ids_genomas():
    search_term = "bacteria[orgn] AND complete genome[title]"
    handle = Entrez.esearch(db="nucleotide", retmax=10000, term=search_term)
    genome_id = Entrez.read(handle)['IdList']
    return genome_id

# Función para descargar genoma en formato FASTA
def descargar_genoma_fasta(genome_id):
    handle = Entrez.efetch(db='nucleotide', id=genome_id, rettype='fasta', retmode='text')
    fasta_data = handle.read()
    return fasta_data

# Obtener 10000 IDs de genomas de bacterias
genome_ids = obtener_ids_genomas()
print(len(genome_ids))

# Seleccionar 100 IDs al azar
genome_ids_aleatorios = random.sample(genome_ids, 100)

# Abrir el archivo multifasta
with open("multifasta_genomas.fasta", "w") as multifasta_file:
    # Descargar y escribir las secuencias en el archivo multifasta
    for i, genome_id in enumerate(genome_ids_aleatorios):
        try:
            fasta_data = descargar_genoma_fasta(genome_id)
            multifasta_file.write(fasta_data)
            print(f"Descargado y agregado al archivo multifasta: genoma ID {genome_id}")
        except Exception as e:
            print(f"Error al descargar el genoma ID {genome_id}: {e}")
        time.sleep(1)  # Para evitar problemas de tasa de solicitud

print("Descarga completa de 100 genomas en formato multifasta.")
