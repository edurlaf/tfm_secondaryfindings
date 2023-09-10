# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 22:15:15 2023

@author: kindi
"""
import os
import gzip
import subprocess

def normalize_vcf_with_bcftools(input_vcf_path):
    # norm requiere el archivo input comprimido. opto por split multiallelic (-m -). si quisiéramos hacer left aligned, necesitamos un fasta de referencia
    try:
        # Rutas de los archivos de entrada y salida
        output_vcf_path = input_vcf_path.split(".vcf")[0] + "_normalized.vcf"
        
        # Verificar si el archivo VCF está comprimido
        if not input_vcf_path.endswith(".gz"):
            # Comprimir el archivo si no está comprimido
            compressed_vcf_path = input_vcf_path + ".gz"
            with open(input_vcf_path, 'rb') as f_in:
                with gzip.open(compressed_vcf_path, 'wb') as f_out:
                    f_out.writelines(f_in)
            input_vcf_path = compressed_vcf_path
        
        # Verificar si el archivo VCF necesita ser indexado
        if not (os.path.exists(input_vcf_path + ".csi") or os.path.exists(input_vcf_path + ".tbi")):
            print('no está indexado')
            # Indexar el archivo VCF si no está indexado
            #index_command = f"bcftools index {input_vcf_path}"
            #os.system(index_command)
            index_command = ["bcftools", "index", input_vcf_path]
            subprocess.run(index_command, check=True)
      
        # Comando para normalizar con bcftools
        #bcftools_command = f"bcftools norm -Ov -o {output_vcf_path} {input_vcf_path}"
        bcftools_command = ["bcftools", "norm", "-O", "v", "-m", "-", "-o", output_vcf_path, input_vcf_path]
        
        # # Ejecutar el comando utilizando os.system
        # os.system(bcftools_command)
        # Ejecutar el comando utilizando subprocess
        subprocess.run(bcftools_command, check=True)

        print("Normalización con bcftools completada.")
    except Exception as e:
        print(f"Error durante la normalización con bcftools: {e}")
