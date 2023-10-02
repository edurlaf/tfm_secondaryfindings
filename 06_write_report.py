# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:26:14 2023

@author: kindi
"""

import json
import pandas as pd


# Cargar el catálogo de genes de riesgo personal desde el archivo JSON
def load_json(json_file):
    with open(json_file, 'r') as file:
        gene_cat = json.load(file)
    return gene_cat

def check_inheritance(results, genes_cat):
    #for gene_cat in genes_cat['genes']:
    reported_variants = {}
    #for result in results:
    for variant_key, variant_info in results.items():
        variant_gene = variant_info["Gene"]
        # variant_info = results[result]
        # variant_gene = variant_info["Gene"]
        #print(variant_info["Gene"])

        for gene in genes_cat['genes']:
            #print(gene['gene_symbol'])
            if gene['gene_symbol'] ==  variant_gene:
                #print(variant_gene)
                #print(gene["inheritance"])
                inher = gene["inheritance"]
                
        
                # Si herencia es dominante (AD) ### todo esto está sin temrinar
                if inher == 'AD' or inher == 'SD' or inher == 'XL':
                    # Combina la información de la variante y el gen
                    combined_info = {
                        "Gene": variant_gene,
                        "Genotype": variant_info["Genotype"],
                        "rs": variant_info["rs"],
                        "IntervarClassification": variant_info["IntervarClassification"],
                        "ClinvarClinicalSignificance": variant_info["ClinvarClinicalSignificance"],
                        "ReviewStatus": variant_info["ReviewStatus"],
                        "ClinvarID": variant_info["ClinvarID"],
                        "Orpha": variant_info["Orpha"],
                        "Phenotype": gene["phenotype"],
                        "ACMG_version": gene["ACMG_version"],
                        "OMIM_disorder": gene["OMIM_disorder"],
                        "inheritance": gene["inheritance"],
                        "variants_to_report": gene["variants_to_report"]
                    }
                    #print(combined_info)
                    # Agrega la información combinada a reported_variants
                    reported_variants[variant_key] = combined_info
  

                elif inher == 'AR':
                    # if gene = 'HFE':
                    #     "variants_to_report": "HFE p.C282Yl\n homozygotes only"
                    # Si la variante está en homocigosis, se reporta
                    if variant_info["Genotype"] == 'hom':
                        # Combina la información de la variante y el gen
                        combined_info = {
                            "Gene": variant_gene,
                            "Genotype": variant_info["Genotype"],
                            "rs": variant_info["rs"],
                            "IntervarClassification": variant_info["IntervarClassification"],
                            "ClinvarClinicalSignificance": variant_info["ClinvarClinicalSignificance"],
                            "ReviewStatus": variant_info["ReviewStatus"],
                            "ClinvarID": variant_info["ClinvarID"],
                            "Orpha": variant_info["Orpha"],
                            "Phenotype": gene["phenotype"],
                            "ACMG_version": gene["ACMG_version"],
                            "OMIM_disorder": gene["OMIM_disorder"],
                            "inheritance": gene["inheritance"],
                            "variants_to_report": gene["variants_to_report"]
                        }
                        # Agrega la información combinada a reported_variants
                        reported_variants[variant_key] = combined_info
                        
                    # Si la variante está en heterocigosis, sólo se portará si se encuentra otra variante en el mismo gen    
                    elif variant_info["Genotype"] == 'het':
                        # Verifica si hay otra variante en el mismo gen
                        other_variant_in_gene = False
                        for other_variant_key in results:
                            other_variant_info = results[other_variant_key]
                            if (
                                other_variant_info["Gene"] == variant_gene
                                and other_variant_key != variant_key
                            ):
                                other_combined_info = {
                                    "Gene": variant_gene,
                                    "Genotype": other_variant_info["Genotype"],
                                    "rs": other_variant_info["rs"],
                                    "IntervarClassification": other_variant_info["IntervarClassification"],
                                    "ClinvarClinicalSignificance": other_variant_info["ClinvarClinicalSignificance"],
                                    "ReviewStatus": other_variant_info["ReviewStatus"],
                                    "ClinvarID": other_variant_info["ClinvarID"],
                                    "Orpha": other_variant_info["Orpha"],
                                    "Phenotype": gene["phenotype"],
                                    "ACMG_version": gene["ACMG_version"],
                                    "OMIM_disorder": gene["OMIM_disorder"],
                                    "inheritance": gene["inheritance"],
                                    "variants_to_report": gene["variants_to_report"]
                                }
                                
                                # Combina la información de la variante y el gen
                                combined_info = {
                                    "Gene": variant_gene,
                                    "Genotype": variant_info["Genotype"],
                                    "rs": variant_info["rs"],
                                    "IntervarClassification": variant_info["IntervarClassification"],
                                    "ClinvarClinicalSignificance": variant_info["ClinvarClinicalSignificance"],
                                    "ReviewStatus": variant_info["ReviewStatus"],
                                    "ClinvarID": variant_info["ClinvarID"],
                                    "Orpha": variant_info["Orpha"],
                                    "Phenotype": gene["phenotype"],
                                    "ACMG_version": gene["ACMG_version"],
                                    "OMIM_disorder": gene["OMIM_disorder"],
                                    "inheritance": gene["inheritance"],
                                    "variants_to_report": gene["variants_to_report"]
                                }
                                # Agrega la información de ambas variantes a reported_variants
                                reported_variants[variant_key] = combined_info
                                reported_variants[other_variant_key] = other_combined_info
    return(reported_variants)     

                
    
def write_report(pr_results, rr_results, fg_results, out_path):
    # Crear un escritor de Excel
    writer = pd.ExcelWriter(out_path, engine='xlsxwriter')
    
    for category in categories:
        if category == 'PR':
            pr_reported_results = check_inheritance(pr_results, gene_cat)
            pr_final = check_diagnosis(pr_reported_results)
            pr_df =  pd.DataFrame.from_dict(pr_reported_results, orient='index')
            pr_df = dict_to_df_results
            pr_df.to_excel(writer, sheet_name='PR Results', index=False)
            
        elif category == 'RR':
            rr_results = check_inheritance(rr_results)
            rr_final = check_diagnosis(rr_results)
            rr_df = dict_to_rr_results
            rr_df.to_excel(writer, sheet_name='RR Results', index=False)

        else:
            

    fg_df = pd.DataFrame(fg_results)
    fg_df.to_excel(writer, sheet_name='FG Results', index=False)


    # Guarda el archivo Excel
    writer.save()