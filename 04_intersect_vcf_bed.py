# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 23:13:23 2023

@author: kindi
"""

from pybedtools import BedTool

def intersect_vcf_with_bed(vcf_path, category, assembly):
    """
    Intersecta un archivo VCF con un archivo BED y guarda las variantes intersectadas en un nuevo archivo VCF.
    
    Args:
        vcf_path (str): Ruta al archivo VCF que se va a intersectar.
        category (str): Categoría de genes o variantes a utilizar para la intersección (PR, RR o FG).
        assembly (str): Versión del ensamblaje genómico a utilizar (por ejemplo, "grch38").
    
    Returns:
        None
    
    Raises:
        Exception: Si ocurre un error durante la intersección o si la categoría no es válida.
    """
    try:
        # Definir un diccionario de categorías a rutas de archivos BED
        category_to_bed = {
            "PR": f"personal_genes_{assembly}.bed",
            "RR": f"reproductivo_genes_{assembly}.bed",
            "FG": f"farmacogenetico_variants_{assembly}.bed"
        }

        # Verifica si la categoría es válida
        if category not in category_to_bed:
            raise ValueError(f"{category} no es una categoría válida (PR, RR o FG)")
            return
        
        # Obtener la ruta del archivo BED correspondiente y del archivo de salida
        bed_path = category_to_bed[category]
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
