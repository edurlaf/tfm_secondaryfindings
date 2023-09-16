# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 22:33:14 2023

@author: kindi
"""
import subprocess
import json

# Definir las funciones para cada módulo y opción
def run_intervar(vcf_path, output_file, category, assembly):
    """
    Ejecuta el programa Intervar para anotar variantes genéticas.
    
    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        output_file (str): Ruta al archivo de salida de Intervar.
        category (str): Categoría de genes para la anotación.
        assembly (str): Ensamblaje genómico a utilizar.
    
    Raises:
        subprocess.CalledProcessError: Si hay un error al ejecutar Intervar.
    """
    try:
        # Ruta vcf interseccion y directorio de salida
        input_vcf = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
        
        # Construir el comando para ejecutar Intervar
        cmd = [
            "/home/sagarruxki/InterVar/Intervar.py",
            "-b", "hg19",# assembly en clinvar es 37, aqui 19. decidir .   además no sé si se puede 38 en intervar
            "-i", input_vcf,
            "--input_type", "VCF",
            "-o", output_file
        ]

        # Ejecutar el comando y capturar la salida
        # Cambiar el directorio de trabajo solo para el comando Intervar
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd="/home/sagarruxki/InterVar") as process:
            output, _ = process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar InterVar: {e.output}")
        
def map_review_status(review_status):
    """
    Mapea el estado de revisión (review status) a un nivel de evidencia (evidence level) equivalente.
    
    Args:
        review_status (str): Estado de revisión en texto.
    
    Returns:
        int: Nivel de evidencia equivalente en número.
    """
    # Mapear "Review status" a número de estrellas
    mapping = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting interpretations": 1,
        "criteria provided, single submitter": 1,
        "no assertion for the individual variant": 0,
        "no assertion criteria provided": 0,
        "no assertion provided": 0
    }
    return mapping.get(review_status.lower(), 0)  # Valor predeterminado es 0 si no se encuentra en el mapeo

def run_clinvar_filtering(vcf_path, evidence_level, category, clinvar_path):
    """
    Filtra variantes del archivo VCF basado en datos de CLINVAR y un nivel de evidencia.
    
    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        evidence_level (int): Nivel de evidencia deseado.
        category (str): Categoría de genes para la anotación.
        clinvar_path (str): Ruta al archivo de base de datos de CLINVAR.
    
    Raises:
        Exception: Si ocurre un error durante el filtrado de variantes.
    """
    try:
        # Cargar el catálogo de genes desde el archivo JSON
        catalog_genes_path = "/home/sagarruxki/" + category + "_risk_genes.json"
        with open(catalog_genes_path, "r") as json_file:
                catalog_genes = json.load(json_file)
                
        # Leer la base de datos de CLINVAR
        clinvar_db = {}  # Un diccionario para almacenar la información de CLINVAR
        
        with open(clinvar_path, "r") as db_file:
            header = db_file.readline().strip().split("\t")
            for line in open(clinvar_path):
                line = line.rstrip()
                if line == "":
                    continue
                fields = line.strip().split("\t")
                variant = fields[9] + ":" + fields[14] + ":" + fields[15] + ":" + fields[16]
                gene = fields[2]
                clinical_significance = fields[3]
                rs_id = fields[4]
                review_status = fields[12]
                clinvar_id = fields[5]
                clinvar_db[variant] = {
                    "Gene": gene,
                    "ClinicalSignificance": clinical_significance,
                    "rs": rs_id,
                    "ReviewStatus": review_status,
                    "ClinvarID": clinvar_id
                }
        
        # Filtrar las variantes según CLINVAR y el nivel de evidencia
        filtered_variants = []
        
        input_vcf = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
        
        with open(input_vcf, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("#"):
                    continue  # Saltar líneas de encabezado del VCF
                fields = line.strip().split("\t")
                variant = fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4]

                if variant in clinvar_db and clinvar_db[variant]["ClinicalSignificance"] in ["Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"]:
                    stars = map_review_status(clinvar_db[variant]["ReviewStatus"])
                    if stars >= evidence_level:                   
                        # Formatear la línea en el formato deseado
                        formatted_line = f"{variant}\t{clinvar_db[variant]['Gene']}\t{clinvar_db[variant]['ClinicalSignificance']}\t{clinvar_db[variant]['rs']}\t{clinvar_db[variant]['ReviewStatus']}\t{clinvar_db[variant]['ClinvarID']}\n"
                        filtered_variants.append(formatted_line)
        
        # Guardar las variantes filtradas en un nuevo archivo VCF
        output_vcf = vcf_path.split("normalized")[0] + category + "_clinvar.vcf"
        with open(output_vcf, "w") as output_file:
            for line in filtered_variants:
                output_file.write(line)
        
        print(f"Variantes filtradas guardadas en '{output_vcf}'")

    except Exception as e:
        print(f"Error al filtrar variantes: {e}")


def run_personal_risk_module(vcf_path, assembly, mode, evidence_level, category, clinvar_path):
    """
    Ejecuta el módulo de riesgo personal según el modo seleccionado.
    
    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        assembly (str): Ensamblaje genómico a utilizar.
        mode (str): Modo de ejecución ("basic" o "advanced").
        evidence_level (int): Nivel de evidencia deseado.
        category (str): Categoría de genes para la anotación.
    """
    if mode == "basic":
        run_acmg(vcf_path, output_dir, category)
    elif mode == "advanced":
        run_acmg(vcf_path)
        run_clinvar_filtering(vcf_path, evidence_level, clinvar_path)
    else:
        print("Modo no válido. Elija 'basic' o 'advanced'.")

