# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:21:45 2023

@author: kindi
"""
import subprocess
import json
import pysam
import csv

# Definir las funciones para cada módulo y opción
def annotate_fg_variants(dir_path, vcf_path, output_file, assembly):
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
    # Cargar archivo fg_json con variantes farmacogenéticas
    fg_json_path = dir_path + 'fg_risk_variants_' + assembly + '.json'
    with open(fg_json_path, 'r') as file:
        fg_json = json.load(file)

    annotated_variants = []

    # Abrir el archivo VCF
    vcf_int_path = vcf_path.split('.vcf')[0] + '_fg_intersection.vcf'
    
    with open(vcf_int_path, 'r') as vcf_file:
    #with pysam.VariantFile(vcf_path) as vcf_file:
        for line in vcf_file:
            fields = line.split('\t')
            variant_key = f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
            #variant_key = f"{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}"
            variant_info = fg_json['variants'].get(variant_key, None)

            if variant_info:
                annotated_variant = {
                    "Variant": variant_key,
                    "GT": fields[9].split(':')[0],
                    "Gene": variant_info["gene_symbol"],
                    "rs": variant_info["rs"]
                }
            else:
                annotated_variant = {
                    "Variant": variant_key,
                    "GT": fields[9].split(':')[0],
                    "Gene": "No encontrado",
                    "rs": "No encontrado",
                }

            annotated_variants.append(annotated_variant)
            
            file = './annotated_variants.txt'
            with open(file, 'w') as fout:
                json.dump(annotated_variants, fout)
            

    return annotated_variants

def assign_cyp2c9_diplotype(variants):
    # Gen
    gene = 'CYP2C9'
    
    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    
    variants_gene = {}
    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']
            #break  # Se encontró una variante en el gen objetivo, así que puedes salir del bucle
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        if 'rs1799853' in variants_gene.keys():
            if variants_gene['rs1799853'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                else:
                    if 'rs1057910' in variants_gene.keys():
                        diplotype = '*2/*3'
                    else:
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        elif 'rs1057910' in variants_gene.keys():
            if variants_gene['rs1057910'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3'
                else:
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    return(diplotype)                    
        
def assign_cyp2c19_diplotype(variants):
    # Gen
    gene = 'CYP2C19'
    
    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    
    variants_gene = {}
    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']
            #break  # Se encontró una variante en el gen objetivo, así que puedes salir del bucle
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False or (len(variants_gene) == 1 and 'rs3758581' in variants_gene.keys()):
        diplotype = '*1/*1'
        
    else:
        if 'rs12769205' in variants_gene.keys() and 'rs4244285' in variants_gene.keys():
            if variants_gene['rs12769205'] == '1/1'  and variants_gene['rs4244285'] == '1/1': # en realidad habría que comprobar si es len >2, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 2 or (len(variants_gene) == 3 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*2'
                else:
                    if 'rs4986893' in variants_gene.keys():
                        diplotype = '*2/*3'
                    elif 'rs12248560' in variants_gene.keys():
                        diplotype = '*2/*17'
                    else:
                        print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                        
        elif 'rs4986893' in variants_gene.keys():
            if variants_gene['rs4986893'] == '1/1':
                diplotype = '*3/*3'
            else:
                if len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                    diplotype = '*1/*3'
                elif 'rs12248560' in variants_gene.keys():
                    diplotype = '*3/*17'
                else:
                    print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # me falta el *17
        elif 'rs12248560' in variants_gene.keys():
            if variants_gene['rs12248560'] == '1/1':
                diplotype = '*17/*17'
            elif len(variants_gene) == 1 or (len(variants_gene) == 2 and 'rs3758581' in variants_gene.keys()):
                diplotype = '*1/*17'
            else:
                print('Se han encontrado variantes en CYP2C19 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                
    return(diplotype)

def assign_dpyd_diplotype(variants):
    # Gen
    gene = 'DPYD'
    
    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    
    variants_gene = {}
    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']
            #break  # Se encontró una variante en el gen objetivo, así que puedes salir del bucle
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2A:
        if 'rs3918290' in variants_gene.keys():
            if variants_gene['rs3918290'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2A/*2A'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2A'
                elif len(variants_gene) == 2:
                    if 'rs55886062' in variants_gene.keys():
                        diplotype = '*2A/*13'
                    elif 'rs67376798' in variants_gene.keys():
                        diplotype = '*2A/c.2846A>T'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*2a/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *13:
        elif 'rs55886062' in variants_gene.keys():
            if variants_gene['rs55886062'] == '1/1':
                diplotype = '*13/*13'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*13'
                elif len(variants_gene) == 2:
                    elif 'rs67376798' in variants_gene.keys():
                        diplotype = '*13/c.2846A>T'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = '*13/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs67376798' in variants_gene.keys():
            if variants_gene['rs67376798'] == '1/1':
                diplotype = 'c.2846A>T/c.2846A>T'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/c.2846A>T'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    if 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
                        diplotype = 'c.2846A>T/HapB3'
                    else:
                        print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        # Chequear c.2846A>T:
        elif 'rs75017182' in variants_gene.keys() and 'rs56038477' in variants_gene.keys():
            if variants_gene['rs75017182'] == '1/1' and variants_gene['rs56038477'] == '1/1':
                diplotype = 'HapB3/HapB3'
            else:
                if len(variants_gene) == 2:
                    diplotype = '*1/HapB3'
                else len(variants_gene) == 3:
                    print('Se han encontrado variantes en DPYD no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    return(diplotype)  
        
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