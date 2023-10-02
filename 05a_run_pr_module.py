# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 22:33:14 2023

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
        
def process_intervar_output(vcf_path, category):
    """
    Procesa el archivo de salida de InterVar y extrae los campos necesarios.

    Args:
        intervar_output_file (str): Ruta al archivo de salida de InterVar.

    Returns:
        list: Una lista de diccionarios con los campos extraídos.
    """
    
    #intervar_output_file = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
    intervar_output_file = "/home/sagarruxki/input_vcf_example/NM769.01_pr.hg19_multianno.txt - copia.intervar" #arreglar
    intervar_results = {}
    
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
                #if classification in ["Pathogenic", "Likely pathogenic"]:
                #if classification in ["Pathogenic", "Likely pathogenic", "Benign"]:   #sólo por 
                
                # Crear un diccionario con los campos extraídos
                intervar_results[variant] = {
                    "Gene": ref_gene,
                    "Rs": avsnp,
                    "Classification": classification,
                    "Orpha": orpha,
                    "GT": other_info
                }

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
        # catalog_genes_path = "/home/sagarruxki/" + category + "_risk_genes.json"
        # with open(catalog_genes_path, "r") as json_file:
        #         catalog_genes = json.load(json_file)
                
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
        
        # # Filtrar las variantes según CLINVAR y el nivel de evidencia
        # #filtered_variants = []
        # filtered_variants = {}
        
        # input_vcf = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
        
        # with open(input_vcf, "r") as vcf_file:
        #     for line in vcf_file:
        #         if line.startswith("#"):
        #             continue  # Saltar líneas de encabezado del VCF
        #         fields = line.strip().split("\t")
        #         variant = fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4]
        #         gt = fields[9].split(':')[0]

        #         #if variant in clinvar_db and clinvar_db[variant]["ClinicalSignificance"] in ["Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"]:
        #         if variant in clinvar_db:
        #             stars = map_review_status(clinvar_db[variant]["ReviewStatus"])
        #             if stars >= evidence_level:                   
        #                 # # Formatear la línea en el formato deseado
        #                 # formatted_line = f"{variant}\t{clinvar_db[variant]['Gene']}\t{gt}\t{clinvar_db[variant]['ClinicalSignificance']}\t{clinvar_db[variant]['rs']}\t{clinvar_db[variant]['ReviewStatus']}\t{clinvar_db[variant]['ClinvarID']}\n"
        #                 # filtered_variants.append(formatted_line)
                        
        #                 # Guardar la información de la variante en el diccionario
        #                 filtered_variants[variant] = {
        #                     "Gene": clinvar_db[variant]['Gene'],
        #                     "Genotype": gt,
        #                     "ClinicalSignificance": clinvar_db[variant]['ClinicalSignificance'],
        #                     "rs": clinvar_db[variant]['rs'],
        #                     "ReviewStatus": clinvar_db[variant]['ReviewStatus'],
        #                     "ClinvarID": clinvar_db[variant]['ClinvarID']
        #                 }
                                
        # # Guardar las variantes filtradas en un nuevo archivo VCF
        # output_vcf = vcf_path.split("normalized")[0] + category + "_clinvar.vcf"
        # with open(output_vcf, "w") as output_file:
        #     # for line in filtered_variants:
        #     #     output_file.write(line)
        #     for variant, info in filtered_variants.items():
        #         formatted_line = f"{variant}\t{info['Gene']}\t{info['Genotype']}\t{info['ClinicalSignificance']}\t{info['rs']}\t{info['ReviewStatus']}\t{info['ClinvarID']}\n"
        #         output_file.write(formatted_line)
        
        # print(f"Variantes filtradas guardadas en '{output_vcf}'")
        
        # return filtered_variants
        return(clinvar_dct)

    except Exception as e:
        print(f"Error al filtrar variantes: {e}")
        
# def load_clinvar_database(clinvar_path):
#     """
#     Carga la base de datos de ClinVar desde un archivo de texto en un diccionario.

#     Args:
#         clinvar_path (str): Ruta al archivo de base de datos de CLINVAR.

#     Returns:
#         dict: Un diccionario con la información de CLINVAR.
#     """
#     clinvar_db = {}

#     with open(clinvar_path, "r") as db_file:
#         header = db_file.readline().strip().split("\t")
#         for line in db_file:
#             line = line.rstrip()
#             if line == "":
#                 continue
#             fields = line.strip().split("\t")
#             variant = fields[9] + ":" + fields[14] + ":" + fields[15] + ":" + fields[16]
#             gene = fields[2]
#             clinical_significance = fields[3]
#             rs_id = fields[4]
#             review_status = fields[12]
#             clinvar_id = fields[5]
#             clinvar_db[variant] = {
#                 "Gene": gene,
#                 "ClinicalSignificance": clinical_significance,
#                 "rs": rs_id,
#                 "ReviewStatus": review_status,
#                 "ClinvarID": clinvar_id
#             }

#     return clinvar_db

# # Cargar la base de datos de ClinVar
# clinvar_db = load_clinvar_database(clinvar_path)
        
def combine_results(vcf_norm, category, intervar_results, clinvar_dct):
    """
    Combina los resultados de Intervar y ClinVar en una sola línea por variante.

    Args:
        intervar_results (list): Resultados de Intervar.
        clinvar_db (dict): Base de datos de ClinVar.

    Returns:
        list: Una lista de diccionarios con la información combinada.
    """
    combined_results = {}
    
    # para compbinar los resultados de intervar y clinvar, aunque lo ideal sería recorrer
    # las listas, intervar cambia la anotación, eliminando el nt de referencia en las idnels
    # por lo que no es posible encontrar esa variante en clinvar (no siempre hay rs disponible)
    # así que lo único que se me ocurre es recorrer el vcf de interseccion
    


    # Archivo VCF de intersección
    vcf_path = vcf_norm.split('normalized')[0] + category + '_intersection.vcf'
    
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
                    variant_int = chrom + ':' + str(int(pos) + 1) + ':' + ref[1:] + ':' + alt[1:]
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
    
def write_combined_results_to_tsv(combined_results, output_tsv):
    """
    Escribe los resultados combinados en un archivo TSV.

    Args:
        combined_results (list): Lista de diccionarios con los resultados combinados.
        output_tsv (str): Ruta al archivo de salida TSV.
    """
    
    # Abrir el archivo TSV para escritura
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
    
    
    
    with open(output_tsv, "w", newline="") as tsv_file:
        fieldnames = [
            "Category",
            "Variant",
            "Gene",
            "Genotype",
            "rs",
            "IntervarClassification",
            "ClinvarClinicalSignificance",
            "ReviewStatus",
            "ClinvarID",
            "Orpha"
        ]
        writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter="\t")
        
        # Escribir la fila de encabezado
        writer.writeheader()
        
        # Escribir los resultados combinados
        for result in combined_results:
            writer.writerow(result)

 
        


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
        run_intervar(vcf_path, output_dir, category)
        intervar_results = process_intervar_output(output_dir)
        return(intervar_results)
        #check_inheritance()
    elif mode == "advanced":
        run_intervar(vcf_path, output_file, category, assembly)
        intervar_results = process_intervar_output(output_dir)
        run_clinvar_filtering(vcf_path, evidence_level, clinvar_path)
        combine_results()
        write_combined_results_to_tsv(combined_results, output_tsv_path)
        return(combined_results)
        #check_inheritance()
    else:
        print("Modo no válido. Elija 'basic' o 'advanced'.")
        
    

