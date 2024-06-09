import random
import time
from Bio import Entrez

# Configura tu correo electrónico
Entrez.email = "tu_correo@ejemplo.com"

# Función para obtener IDs de genomas
def obtener_ids_genomas():
    handle = Entrez.esearch(db="assembly", term="bacteria", retmax=1000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Función para descargar genoma en formato FASTA
def descargar_genoma_fasta(genome_id):
    # Primero obtenemos los accesiones asociados al ID del genoma
    handle = Entrez.esummary(db="assembly", id=genome_id, report="full")
    record = Entrez.read(handle)
    handle.close()
    
    # Acceder al campo 'FtpPath_GenBank' que contiene el enlace FTP al genoma
    ftp_path = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
    if not ftp_path:
        raise ValueError(f"No se encontró el enlace FTP para el ID del genoma: {genome_id}")
    
    # Obtener el archivo FASTA del genoma
    fasta_url = ftp_path + "/" + ftp_path.split('/')[-1] + "_genomic.fna.gz"
    
    # Descargar y descomprimir el archivo FASTA
    import urllib.request
    import gzip
    
    fasta_file = "temp.fna.gz"
    urllib.request.urlretrieve(fasta_url, fasta_file)
    
    # Descomprimir el archivo y devolver el contenido
    with gzip.open(fasta_file, 'rt') as f_in:
        fasta_data = f_in.read()
    
    return fasta_data

# Obtener 1000 IDs de genomas de bacterias
genome_ids = obtener_ids_genomas()

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
