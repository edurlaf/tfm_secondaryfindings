# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:26:14 2023

@author: kindi
"""

import json
import pandas as pd


# # Cargar el catálogo de genes de riesgo personal desde el archivo JSON
# def load_json(json_file):
#     with open(json_file, 'r') as file:
#         gene_cat = json.load(file)
#     return gene_cat

def combine_variant_and_gene_info(variant_info, gene_info):
    """
    Combina la información de una variante y un gen.

    Args:
        variant_info (dict): Información de la variante.
        gene_info (dict): Información del gen.

    Returns:
        dict: Información combinada de la variante y el gen.
    """
    combined_info = {
        "Gene": variant_info["Gene"],
        "Genotype": variant_info["Genotype"],
        "rs": variant_info.get("rs", ""),
        "IntervarClassification": variant_info["IntervarClassification"],
        "ClinvarClinicalSignificance": variant_info.get("ClinvarClinicalSignificance", ""),
        "ReviewStatus": variant_info.get("ReviewStatus", ""),
        "ClinvarID": variant_info.get("ClinvarID", ""),
        "Orpha": variant_info.get("Orpha", ""),
        "Phenotype": gene_info["phenotype"],
        "ACMG_version": gene_info.get("ACMG_version", ""),  # Usar get para manejar la falta de 'ACMG_version'
        "OMIM_disorder": gene_info["OMIM_disorder"],
        "inheritance": gene_info["inheritance"],
        "variants_to_report": gene_info.get("variants_to_report", "")  # Usar get para manejar la falta de 'variants_to_report'
    }
    return combined_info

def check_inheritance(results, category, categories_path):
    """
    Comprueba la herencia de las variantes en función de la categoría de genes y genera un diccionario de variantes
    a informar según las reglas de herencia definidas en el archivo JSON de genes.

    Args:
        results (dict): Un diccionario de resultados de variantes.
        category (str): Categoría de genes para la anotación.

    Returns:
        dict: Un diccionario de variantes a informar siguiendo las reglas de herencia definidas.
    """    

    try:
        # Cargar el archivo JSON de categoría de genes    
        genes_cat_path = f"{categories_path}{category.upper()}/{category}_risk_genes.json"
        genes_cat = None
        with open(genes_cat_path, "r") as genes_cat_file:
            genes_cat = json.load(genes_cat_file)
        
        # Crear diccionario de variantes a informar 
        reported_variants = {}
    
        # Recorrer claves del diccionario de resultados combinados
        for variant_key, variant_info in results.items():
            variant_gene = variant_info["Gene"]
            # Obtener el modo de herencia del gen
            for gene in genes_cat['genes']:
                if gene['gene_symbol'] ==  variant_gene:
                    inher = gene["inheritance"]
                          
                    # Si herencia es dominante (AD), semidominante o ligado al X, o si la categoría es RR, se informa la variante
                    if inher in ['AD', 'SD', 'XL'] or category == 'rr':
                        # Combina la información de la variante y el gen
                        combined_info = combine_variant_and_gene_info(variant_info, gene)
                        # Agrega la información combinada a reported_variants
                        reported_variants[variant_key] = combined_info
                    
                    # Si la herencia es recesiva, comprobar genotipo y/o otras variantes en mismo gen
                    elif inher == 'AR':
                        # if gene = 'HFE':
                        #     "variants_to_report": "HFE p.C282Yl\n homozygotes only"
                        # Si la variante está en homocigosis, se reporta
                        if variant_info["Genotype"] == 'hom':
                            # Combina la información de la variante y el gen
                            combined_info = combine_variant_and_gene_info(variant_info, gene)
                            # Agrega la información combinada a reported_variants
                            reported_variants[variant_key] = combined_info
                            
                        # Si la variante está en heterocigosis, sólo se reportará si se encuentra otra variante en el mismo gen    
                        elif variant_info["Genotype"] == 'het':
                            # Verifica si hay otra variante en el mismo gen
                            other_variant_in_gene = False
                            for other_variant_key in results:
                                other_variant_info = results[other_variant_key]
                                if (
                                    other_variant_info["Gene"] == variant_gene
                                    and other_variant_key != variant_key
                                ):
                                    other_combined_info = combine_variant_and_gene_info(other_variant_info, gene)                                
                                    combined_info = combine_variant_and_gene_info(variant_info, gene)
    
                                    # Agrega la información de ambas variantes a reported_variants
                                    reported_variants[variant_key] = combined_info
                                    reported_variants[other_variant_key] = other_combined_info
        return(reported_variants)     
    except Exception as e:
        print(f"Error en check_inheritance: {str(e)}")
        return {}
 
def get_hpo_dct(categories_path):

    hpo_file = f"{categories_path}phenotype_to_genes.txt"
    hpo_data = {}
    # Abrir el archivo de texto en modo lectura
    with open('hpo_file', 'r') as file:
        # Salta la primera línea si contiene encabezados
        next(file)
        
        # Iterar a través de las líneas del archivo
        for line in file:
            # Dividir cada línea en columnas usando tabulaciones o espacios en blanco como delimitadores
            fields = line.strip().split('\t')
            
            # Extraer los valores de cada columna
            hpo_id, hpo_name, ncbi_gene_id, gene_symbol, disease_id = fields
            
            # Comprobar si el HPO ya está en el diccionario, y si no, crea un nuevo diccionario
            if hpo_id not in hpo_data:
                hpo_data[hpo_id] = {
                    "hpo_name": hpo_name,
                    "gene_info": []
                }
            
            # Agrega información del gen y la enfermedad al diccionario correspondiente
            hpo_data[hpo_id]["gene_info"].append({
                "ncbi_gene_id": ncbi_gene_id,
                "gene_symbol": gene_symbol,
                "disease_id": disease_id
            })
               
    return(hpo_data)

def check_diagnosis(user_hpos, reported_results):
    """
    Comprueba si los HPOs del usuario coinciden con los resultados reportados.

    Args:
        user_hpos (list): Lista de HPOs proporcionada por el usuario.
        reported_results (dict): Resultados de diagnóstico reportados.

    Returns:
        list: Lista de resultados coincidentes.
    """
    matching_results = []

    for hpo in user_hpos:
        for result in reported_results[hpo]["gene_info"]:
            # Comprueba si gene_symbol o disease_id coinciden con los del resultado
            if (result["gene_symbol"] in user_hpos or result["disease_id"] in user_hpos):
                matching_results.append({
                    "hpo_id": hpo,
                    "hpo_name": reported_results[hpo]["hpo_name"],
                    "ncbi_gene_id": result["ncbi_gene_id"],
                    "gene_symbol": result["gene_symbol"],
                    "disease_id": result["disease_id"]
                })

    return matching_results

def write_report(pr_results, rr_results, fg_results, haplot_results, categories_path, out_path, categories):
    """
    Escribe los resultados de las categorías PR, RR y FG en un archivo Excel.

    Args:
        pr_results (dict): Resultados de la categoría PR.
        rr_results (dict): Resultados de la categoría RR.
        fg_results (dict): Resultados de la categoría FG.
        out_path (str): Ruta al archivo de salida.
    """
    try:
        outfile = f"{out_path}final_results.xlsx"
        # Crear un objeto ExcelWriter para el archivo Excel
        with pd.ExcelWriter(outfile) as writer:
            for category in categories:
                if category == 'pr':
                    reported_results = check_inheritance(pr_results, category, categories_path)
                    #pr_final = check_diagnosis(pr_reported_results)  # pendiente de desarrollar, warning si los términos orpha se corresponden con los hpo del paciente
                    results_df =  pd.DataFrame.from_dict(reported_results, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)
                    
                elif category == 'rr':
                    reported_results = check_inheritance(rr_results, category, categories_path)
                    #rr_final = check_diagnosis(rr_reported_results)  # pendiente de desarrollar, warning si los términos orpha se corresponden con los hpo del paciente
                    results_df =  pd.DataFrame.from_dict(reported_results, orient='index')
                    results_df.to_excel(writer, sheet_name= category.upper() + ' results', index=True)
        
                else:
                    fg_df = pd.DataFrame(fg_results)
                    fg_df.to_excel(writer, sheet_name='FG results', index=False)
                    
                    haplot_df = pd.DataFrame(haplot_results)
                    haplot_df.to_excel(writer, sheet_name='FG Diplotype-Phenotype', index=False)
                
        print(f"Los resultados se han guardado en '{outfile}'.")
        
    except Exception as e:
        print(f"Error en write_report: {str(e)}")
