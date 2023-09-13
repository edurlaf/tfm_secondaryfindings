# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 23:13:23 2023

@author: kindi
"""

from pybedtools import BedTool

def intersect_vcf_with_bed(vcf_path, category, assembly):
    try:
        # Rutas del archivo BED y del VCF de salida:
        if category == "PR":
            bed_path = "personal_genes_" + assembly + ".bed"
            
        elif category == "RR":
            bed_path = "reproductivo_genes_" + assembly + ".bed"
            
        elif category == "FG":
            bed_path = "farmacogenetico_variants_" + assembly + ".bed"
        
        else:
            print(category + " no es una de las categorías a elegir (PR, RR o FG)")
            return
        
        output_vcf_path = vcf_path.split("normalized")[0] + category + "_intersection.vcf"
            
        # Cargar el archivo VCF y el archivo BED utilizando pybedtools
        vcf = BedTool(vcf_path)
        bed = BedTool(bed_path)
        
        # Realizar la intersección utilizando pybedtools
        intersected_variants = vcf.intersect(bed, u=True)
        
        # Guardar las variantes intersectadas en un nuevo archivo VCF
        intersected_variants.saveas(output_vcf_path)
        
        print(f"Intersección completada. Variantes intersectadas guardadas en {output_vcf_path}")
    except Exception as e:
        print(f"Error durante la intersección: {e}")


# faltaría el header, pero realmente no sé si hace falta
