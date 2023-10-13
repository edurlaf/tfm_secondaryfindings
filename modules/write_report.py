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

def check_inheritance(results, category, dir_path):
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
        genes_cat_path = dir_path + category + '_risk_genes.json'
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
                
    
def write_report(pr_results, rr_results, fg_results, haplot_results, out_path, dir_path):
    """
    Escribe los resultados de las categorías PR, RR y FG en un archivo Excel.

    Args:
        pr_results (dict): Resultados de la categoría PR.
        rr_results (dict): Resultados de la categoría RR.
        fg_results (dict): Resultados de la categoría FG.
        out_path (str): Ruta al archivo de salida.
    """
    try:
        for category in categories:
            if category in ['pr', 'rr']:
                reported_results = check_inheritance(results, category, dir_path)
                #pr_final = check_diagnosis(pr_reported_results)  # pendiente de desarrollar, warning si los términos orpha se corresponden con los hpo del paciente
                results_df =  pd.DataFrame.from_dict(reported_results, orient='index')
                results_df.to_excel(out_path, sheet_name= category.upper() + ' results', index=True)
    
            else: #pendiente
                fg_df = pd.DataFrame(fg_results)
                fg_df.to_excel(out_path, sheet_name='FG Results', index=False)
                
                haplot_df = pd.DataFrame(haplot_results)
                haplot_df.to_excel(out_path, sheet_name='FG Diplotype-Phenotypes', index=False)
        print(f"Los resultados se han guardado en '{out_path}'.")
        
    except Exception as e:
        print(f"Error en write_report: {str(e)}")
