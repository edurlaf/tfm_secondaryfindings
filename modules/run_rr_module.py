# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 18:15:14 2023

@author: kindi
"""

import subprocess
import json
import pysam
import csv

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
            "./InterVar/Intervar.py",
            "-b", "hg19",# assembly en clinvar es 37, aqui 19. decidir .   además no sé si se puede 38 en intervar
            "-i", input_vcf,
            "--input_type", "VCF",
            "-o", output_file
        ]

        # Ejecutar el comando y capturar la salida
        # Cambiar el directorio de trabajo solo para el comando Intervar
        with subprocess.Popen(cmd, stderr=subprocess.STDOUT, text=True, cwd="./InterVar") as process:
            output, _ = process.communicate()

    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar InterVar: {e.output}")
        
def parse_intervar_output(vcf_path, category, mode):
    """
    Procesa el archivo de salida de InterVar y extrae los campos necesarios.

    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        category (str): Categoría de genes para la anotación.
        mode (str): Modo de ejecución ("basico" o "avanzado").

    Returns:
        list: Una lista de diccionarios con los campos extraídos.
    """
    
    #intervar_output_file = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
    intervar_output_file = vcf_path.split('.hard-filtered')[0] + '_' + category + '.hg19_multianno.txt.intervar'
    intervar_results = {}
    
    try:    
        with open(intervar_output_file, "r") as intervar_file:
            for line in intervar_file:
                # Ignorar líneas que no son variantes
                if not line.startswith("#"):
                    fields = line.strip().split("\t")
                    # Extraer los campos necesarios y formatearlos como se requiere
                    variant = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
                    ref_gene = fields[8]
                    avsnp = fields[9]
                    #print(fields[13].split(":")[1].split("PVS")[0])
                    classification = fields[13].split(": ")[1].split(" PVS")[0]
                    orpha = fields[32]
                    other_info = fields[-1]
                    
                    # Filtrar solo las variantes patogénicas y posiblemente patogénicas
                    if classification in ["Pathogenic", "Likely pathogenic"] or mode == 'advanced':
                    
                        # Crear un diccionario con los campos extraídos
                        intervar_results[variant] = {
                            "Gene": ref_gene,
                            "Rs": avsnp,
                            "Classification": classification,
                            "Orpha": orpha,
                            "GT": other_info
                        }
    except Exception as e:
        raise Exception(f"Error al procesar la salida de InterVar: {e}")

    return intervar_results
        
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

def run_clinvar_filtering(evidence_level, clinvar_path):
    """
    Filtra variantes de la base de datos de CLINVAR según un nivel de evidencia dado.

    Args:
        evidence_level (int): El nivel de evidencia deseado para filtrar las variantes.
        clinvar_path (str): Ruta al archivo de la base de datos de CLINVAR.

    Returns:
        dict: Un diccionario que contiene las variantes filtradas de CLINVAR y su información relacionada,
              organizadas por variantes.
    
    Raises:
        Exception: Si ocurre un error durante el filtrado de variantes de CLINVAR.
    """
    try:              
        # Leer la base de datos de CLINVAR
        clinvar_dct = {}  # Un diccionario para almacenar la información de CLINVAR
        
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
                stars = map_review_status(review_status)
                if stars >= int(evidence_level):  
                    clinvar_id = fields[5]
                    clinvar_dct[variant] = {
                        "Gene": gene,
                        "ClinicalSignificance": clinical_significance,
                        "rs": rs_id,
                        "ReviewStatus": '(' + str(stars) + ') ' + review_status,
                        "ClinvarID": clinvar_id
                    }
        
        return(clinvar_dct)

    except Exception as e:
        print(f"Error al filtrar variantes: {e}")
        
def combine_results(vcf_norm, category, intervar_results, clinvar_dct):
    """
    Combina los resultados de Intervar y ClinVar en una sola línea por variante.

    Args:
        vcf_norm (str): Ruta al archivo VCF normalizado.
        category (str): Categoría de genes para la anotación.
        intervar_results (dict): Resultados de Intervar.
        clinvar_dct (dict): Base de datos de ClinVar.

    Returns:
        dict: Un diccionario con los resultados combinados.
    """
    combined_results = {}
    
    # para compbinar los resultados de intervar y clinvar, aunque lo ideal sería recorrer
    # las listas, intervar cambia la anotación, eliminando el nt de referencia en las idnels
    # por lo que no es posible encontrar esa variante en clinvar (no siempre hay rs disponible)
    # así que lo único que se me ocurre es recorrer el vcf de interseccion
    


    # Archivo VCF de intersección
    vcf_path = vcf_norm.split('normalized')[0] + category + '_intersection.vcf'
    
    try:
        # Abrir el archivo VCF de intersección
        #with pysam.VariantFile(vcf_path, "r") as vcf_file:
        with open(vcf_path, "r") as vcf_file:
            #for variant in vcf_file:
            # Recorrer el archivo línea por línea
            for line in vcf_file:
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]  # Supongo una sola alternativa porque he splitteado los multiallelic sites
        
                # Construye la clave de búsqueda
                variant_key = f"{chrom}:{pos}:{ref}:{alt}"
        
                # Busca la variante en los resultados de Intervar
                if len(ref) > len(alt):
                    if len(alt) == 1:
                        #variant_int = chrom + ':' + str(int(pos) + 1) + ':' + ref[1:] + ':-'
                        variant_int = f"{chrom}:{int(pos) + 1}:{ref[1:]}:-"
                    else:
                        variant_int = f"{chrom}:{str(int(pos) + 1)}:{ref[1:]}:{alt[1:]}"
                else:
                    variant_int = f"{chrom}:{pos}:{ref}:{alt}"
                        
                intervar_info = intervar_results.get(variant_int)
        
                # Busca la variante en los resultados de ClinVar
                clinvar_info = clinvar_dct.get(variant_key)
    
                # Combina la información si es "Pathogenic" o "Likely pathogenic" en alguno de los dos
                if (intervar_info and intervar_info["Classification"] in ["Pathogenic", "Likely pathogenic"]) or (clinvar_info and clinvar_info["ClinicalSignificance"] in ["Pathogenic", "Likely pathogenic"]):
                    combined_results[variant_key] = {
                        "Gene": clinvar_info["Gene"],
                        "Genotype": intervar_info["GT"],
                        "rs": intervar_info["Rs"] if intervar_info["Rs"] != '.' else clinvar_info["rs"],
                        "IntervarClassification": intervar_info["Classification"],
                        "ClinvarClinicalSignificance": clinvar_info["ClinicalSignificance"],
                        "ReviewStatus": clinvar_info["ReviewStatus"],
                        "ClinvarID": clinvar_info["ClinvarID"],
                        "Orpha": intervar_info["Orpha"]
                }
    except Exception as e:
        raise Exception(f"Error al combinar resultados: {e}")
    
    return(combined_results)
    
def write_combined_results_to_tsv(combined_results, vcf_path, category):
    """
    Escribe los resultados combinados en un archivo TSV.

    Args:
        combined_results (dict): Diccionario con los resultados combinados.
        vcf_path (str): Ruta al archivo VCF normalizado.
        category (str): Categoría de genes para la anotación.
    """
    try:
        # Abrir el archivo TSV para escritura
        output_tsv = vcf_path.split('normalized')[0] + category + '_combinedresults.tsv'
        with open(output_tsv, "w", newline="") as tsv_file:
            fieldnames = ["Variant", "Gene", "Genotype", "rs", "IntervarClassification", "ClinvarClinicalSignificance", "ReviewStatus", "ClinvarID", "Orpha"]
            writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")
        
            # Escribir el encabezado del archivo TSV
            writer.writeheader()
        
            # Iterar sobre las claves (variantes) en el diccionario result
            for variant, info in combined_results.items():
                # Asegurarse de que el diccionario 'info' tenga todas las claves necesarias
                row = {
                    "Variant": variant,
                    "Gene": info.get("Gene", ""),
                    "Genotype": info.get("Genotype", ""),
                    "rs": info.get("rs", ""),
                    "IntervarClassification": info.get("IntervarClassification", ""),
                    "ClinvarClinicalSignificance": info.get("ClinvarClinicalSignificance", ""),
                    "ReviewStatus": info.get("ReviewStatus", ""),
                    "ClinvarID": info.get("ClinvarID", ""),
                    "Orpha": info.get("Orpha", "")
                }
    
                # Escribir la fila en el archivo TSV
                writer.writerow(row)    
     
    except Exception as e:
        raise Exception(f"Error al escribir resultados en archivo TSV: {e}")

def run_reproductive_risk_module(vcf_path, assembly, mode, evidence_level, category, clinvar_path): #sobra category, este modulo es especifico de rr
    """
    Ejecuta el módulo de riesgo personal según el modo seleccionado.
    
    Args:
        vcf_path (str): Ruta al archivo VCF de entrada.
        assembly (str): Ensamblaje genómico a utilizar.
        mode (str): Modo de ejecución ("basic" o "advanced").
        evidence_level (int): Nivel de evidencia deseado.
        category (str): Categoría de genes para la anotación.
        clinvar_path (str): Ruta al archivo de base de datos de CLINVAR.
    """
    if mode == "basic":
        run_intervar(vcf_path, output_dir, category, assembly)
        intervar_results = parse_intervar_output(output_dir, category, mode)
        return(intervar_results)
        #check_inheritance()
    elif mode == "advanced":
        run_intervar(vcf_path, output_dir, category, assembly)
        intervar_results = parse_intervar_output(vcf_path, category, mode)
        clinvar_dct = run_clinvar_filtering(evidence_level, clinvar_path)
        combined_results = combine_results(vcf_path, category, intervar_results, clinvar_dct)
        write_combined_results_to_tsv(combined_results, vcf_path, category)
        return(combined_results)
        #check_inheritance()
    else:
        print("Modo no válido. Elija 'basic' o 'advanced'.")
        return None