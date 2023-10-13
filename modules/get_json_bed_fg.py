# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 22:45:05 2023

@author: kindi
"""
"""
falta haer csv con assembly 38, y su bed y json correspondientes. También falta ver si funciona bien.
"""
import os
import csv
import json
#from pybedtools import BedTool #DEBERÍA MIRARLO PARA ORDENAR LOS CROMOSOMAS
from natsort import natsorted


def generate_json_from_fg_csv(csv_file, assembly, dir_path):
    """
    Genera un archivo JSON a partir de un archivo CSV de datos farmacogenéticos.
    
    Args:
        csv_file (str): Ruta al archivo CSV de entrada.
        assembly (str): Versión de ensamblaje (por ejemplo, "37").
    
    Returns:
        None
    """
    try:
        genes_data = []
        
        # Create dictionary and gene list to store info
        variants_dct = {
        "category": "Hallazgos secundarios de riesgo farmacogenético",
        "variants": {}
        }
        
        variants_lst = []
    
        # Read CSV file and store it in the dictionary
        with open(csv_file, 'r') as file:
            csv_reader = csv.DictReader(file)
            for row in csv_reader:
                alts = row['Alternative'].split(',')
                if len(alts) > 1:
                    for i in alts:
                        variant = row['Chromosome'] + ':' + row['Position'] + ':' + row['Reference'] + ':' + i
                        variant_info = {
                            'gene_symbol': row['Gene'],
                            'rs': row['dbSNP'],
                            }
                        variants_dct["variants"][variant] = variant_info
                else:
                    variant = row['Chromosome'] + ':' + row['Position'] + ':' + row['Reference'] + ':' + row['Alternative']
                    variant_info = {
                        'gene_symbol': row['Gene'],
                        'rs': row['dbSNP'],
                        }
                    variants_dct["variants"][variant] = variant_info
                #variants_lst.append(variant)
        
        #return(variants_dct, variants_lst)
            # Write JSON file
            out_json = f'{dir_path}fg_risk_variants_{assembly}.json'
            with open(out_json, 'w') as json_file:
                json.dump(variants_dct, json_file, indent = 4)
    
    # try:
    #     with open(csv_file, 'r') as file:
    #         csv_reader = csv.DictReader(file)
    #         current_gene = None
    
    #         for row in csv_reader:
    #             if current_gene != row['Gene']:
    #                 if current_gene is not None:
    #                     genes_data.append(gene_info)
    #                 current_gene = row['Gene']
    #                 gene_info = {
    #                     "gene_symbol": current_gene,
    #                     "variants": []
    #                 }
                
    #             variant_info = {
    #                 "chromosome": row['Chromosome'],
    #                 "position": row['Position'],
    #                 "reference": row['Reference'],
    #                 "alternative": row['Alternative'],
    #                 "dbSNP": row['dbSNP']
    #             }
    #             gene_info["variants"].append(variant_info)
            
    #         if current_gene is not None:
    #             genes_data.append(gene_info)
        
    #     json_data = {
    #         "category": "Hallazgos secundarios de riesgo farmacogenetico",
    #         "genes": genes_data
    #     }
    
    #     # Write JSON file
    #     out_json = '/home/sagarruxki/farmacogenetico_risk_genes_' + assembly + '.json'
    #     with open(out_json, 'w') as json_file:
    #         json.dump(json_data, json_file, indent = 4)
        print(f"JSON file '{out_json}' generated successfully.")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def generate_bed_from_fg_csv(csv_file, assembly):
    """
    Genera un archivo BED a partir de un archivo CSV de datos farmacogenéticos.
    
    Args:
        csv_file (str): Ruta al archivo CSV de entrada.
        assembly (str): Versión de ensamblaje (por ejemplo, "37").
    
    Returns:
        None
    """
    bed_data = []
    
    try:
        with open(csv_file, 'r') as file:
            csv_reader = csv.DictReader(file)
            
            for row in csv_reader:
                chr_name = row['Chromosome']
                position = int(row['Position']) - 1  # BED uses 0-based coordinates
                ref_allele = row['Reference']
                alt_alleles = row['Alternative'].split(',')  # Split multiple alleles
                
                for alt_allele in alt_alleles:
                    length_diff = len(alt_allele) - len(ref_allele)
                    start_pos = position
                    end_pos = position + max(len(ref_allele), len(alt_allele))  # Take longer allele length
                    bed_entry = (chr_name, start_pos, end_pos, row['Gene'], ref_allele, alt_allele)
                    bed_data.append(bed_entry)
        
        # Sort bed_data by chromosome, start position, end position, and allele
        sorted_bed_data = natsorted(bed_data)
        
        bed_filename = f"farmacogenetico_variants_grch{assembly}.bed"
        
        with open(bed_filename, "w") as bed_file:
            for entry in sorted_bed_data:
                bed_file.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\t{entry[4]}\t{entry[5]}\n")
        
        print(f"BED file '{bed_filename}' generated and sorted.")
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def get_json_bed_fg(assembly, dir_path):
    """
    Función principal que procesa un archivo CSV de datos farmacogenéticos y genera archivos JSON y BED.

    Args:
        assembly (str): Versión de ensamblaje (por ejemplo, "37").

    Returns:
        None
    """
    try:
        csv_file = f"{dir_path}fg_risk_genes.csv"
        generate_json_from_fg_csv(csv_file, assembly, dir_path)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
    
    #assembly = "GRCh37"  # Replace with the appropriate assembly name
    
    

    
    
