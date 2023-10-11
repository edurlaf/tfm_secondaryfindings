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
def annotate_fg_variants(dir_path, vcf_path, assembly):
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

def check_gene_variants(variants, gene):
    # Flag para verificar si se encontró una variante en el gen objetivo
    found = False
    # Crear diccionario para almacenar genotipo de variantes del gene       #igual mejor guardar la variante también, para el caso de indels con mismo rs y diferente nº de repeticiones
    variants_gene = {}
    
    # Iterar a través de las variantes
    for variant in variants:
        if variant['Gene'] == gene:
            found = True
            variants_gene[variant['rs']] = variant['GT']
    
    return(found, variants_gene)

def assign_cyp2c9_diplotype(variants):
    # Gen
    gene = 'CYP2C9'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
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
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
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
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
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

def assign_nudt15_diplotype(variants):
    # Gen
    gene = 'NUDT15'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2A:
        if 'rs116855232' in variants_gene.keys():
            if len(variants_gene) == 1:
                if variants_gene['rs116855232'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                    diplotype = '*3/*3'
                else:
                    diplotype = '*1/*3'
            elif len(variants_gene) == 2:
                if 'rs746071566' in variants_gene.keys():
                    if variants_gene['rs746071566'] == '1/1' and variants_gene['rs116855232'] == '1/1':
                        diplotype = '*2/*2'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '1/1': #hbaría que chequear variante si es inserción o del
                        diplotype = '*2/*3'
                    elif variants_gene['rs746071566'] == '0/1' and variants_gene['rs116855232'] == '0/1':
                        diplotype = '*1/*2'
                    else:
                        print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
            else:
                print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
        else:
            print('Se han encontrado variantes en NUDT15 no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    return(diplotype)  
 
def assign_tpmt_diplotype(variants):
    # Gen
    gene = 'TPMT'
    
    # Comprobar variantes presentes en gen
    found, variants_gene = check_gene_variants(variants, gene)
            
    # Si no se encontró ninguna variante, diplotipo *1/*1
    if found == False:
        diplotype = '*1/*1'
        
    else:
        # Chequear *2:
        if 'rs1800462' in variants_gene.keys():
            if variants_gene['rs1800462'] == '1/1': # en realidad habría que comprobar si es len >1, porque puede que esté en fase con otro alelo
                diplotype = '*2/*2'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*2'
                elif len(variants_gene) == 2:
                    if 'rs1800460' in variants_gene.keys():
                        diplotype = '*2/*3B'
                    elif 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*2/*4'    
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*2/*3A'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *3B:
        elif 'rs1800460' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3B/*3B'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3B'
                elif len(variants_gene) == 2:
                    if 'rs1142345' in variants_gene.keys() and (variants_gene['rs1800460'] == '0/1') and (variants_gene['rs1142345'] == '0/1'):
                        diplotype = '*1/*3A o *3B/*3C'
                    elif 'rs1800584' in variants_gene.keys():
                        diplotype = '*3B/*4' 
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # Chequear *3C:
        elif 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*3C/*3C'
            else:
                if len(variants_gene) == 1:
                    diplotype = '*1/*3C'
                elif len(variants_gene) == 2:
                    if 'rs1800584' in variants_gene.keys():
                        diplotype = '*3C/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*3B'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')

        # Chequear *4:
        elif 'rs1800584' in variants_gene.keys():
            if variants_gene['rs1800584'] == '1/1' and len(variants_gene) == 1:
                diplotype = '*4/*4'
            else:
                if len(variants_gene) == 1 and variants_gene['rs1800584'] == '0/1':
                    diplotype = '*1/*4'
                elif len(variants_gene) == 2:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                elif len(variants_gene) == 3:
                    if 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
                        diplotype = '*3A/*4'
                    else:
                        print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                else:
                    print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
                    
        # Chequear *3A:  IGUAL REORGANIZO TODO TPMT Y ESTO AL INICIO
        elif 'rs1800460' in variants_gene.keys() and 'rs1142345' in variants_gene.keys():
            if variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3A'
            elif variants_gene['rs1800460'] == '1/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3B'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '1/1' and len(variants_gene) == 2:
                diplotype = '*3A/*3C'
            elif variants_gene['rs1800460'] == '0/1' and variants_gene['rs1142345'] == '0/1' and len(variants_gene) == 2:
                diplotype = '*1/*3A o *3A/*3C'
            else:
                print('Se han encontrado variantes en TPMT no consideradas en la asignación de haplotipos en esta herramienta. Revisar manualmente.')
    return(diplotype)  
       

def run_pharmacogenomic_risk_module(dir_path, vcf_path, assembly): #sobra category, este modulo es especifico de rr
    """
    Ejecuta el módulo de riesgo farmacogenético.
    
    Args:
        dir_path (str)
        vcf_path (str): Ruta al archivo VCF de entrada.
        assembly (str): Ensamblaje genómico a utilizar.
    """
    # Anotar variantes fg presentes en el vcf
    fg_results = annotate_fg_variants(dir_path, vcf_path, assembly)
    
    # Asignar los diplotipos de cada gen
    cyp2c19_diplotype = assign_cyp2c9_diplotype(fg_results)
    dpyd_diplotype = assign_dpyd_diplotype(fg_results)
    nudt15_diplotype = assign_nudt15_diplotype(fg_results)
    tpmt_diplotype = assign_tpmt_diplotype(fg_results)
    

