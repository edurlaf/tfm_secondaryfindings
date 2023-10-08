# -*- coding: utf-8 -*-

"""
Herramienta para el manejo automático de hallazgos secundarios.

Esta herramienta permite a los usuarios analizar archivos VCF para el manejo automático de hallazgos secundarios relacionados con riesgo personal, riesgo reproductivo y farmacogenético.

@Dependencies InterVar and AnnoVar 
@Usage python3 secondary_findings.py input_file.vcf -o <out-path-dir> --mode <Option: 'basic' or 'advanced'> --evidence <integer> --assembly <Option: '37' or '38'> 
@Arguments:
    -vcf (str): Ruta al archivo VCF de entrada.
    -outpath (str): Ruta al directorio donde se guardarán los resultados.
    -mode (str): Modo de análisis (básico o avanzado).
    -evidence (int): Nivel de evidencia de ClinVar para el modo avanzado (1-4).#comprobar que lo he puesto de memoria
    -assembly (str): Ensamblaje genómico a utilizar (GRCh37 o GRCh38).

@Author Edurne Urrutia Lafuente
@Date 2023/08/01
@email edurlaf@gmail.com
@github github.com/edurlaf
"""

#import sys
import os
import json
from myFunctions.functions import arguments, get_json_bed, get_json_bed_fg, get_clinvar, normalize_vcf, pr_module, rr_module, fg_module, write_report

def main():
    
    """
    Read config file
    """
    # Leer el archivo de configuración config.json
    with open("./config.json", "r") as config_file:
        config_data = json.load(config_file)
    
    # Obtener los valores del archivo de configuración
    dir_path = config_data["dir_path"]
       
    
    """
    Get the arguments
    """
    args = arguments()
    
    # Argumentos del usuario
    vcf_file = args.vcf_file
    out_path = args.outpath
    mode = args.mode
    evidence = args.evidence
    assembly = args.assembly.lower()
    
    # Obtener la preferencia del usuario para las categorías a analizar (PR, RR, FG)
    categories_usr = input("Elija las categorías a analizar (PR, RR, FG separados por comas): ")
    categories = [category.strip().lower() for category in categories_usr.split(",")]
 
    # habría que chequear el outpath?
    
    
    """
    Generate JSON and BED files 
    """
    # Comprobar si los archivos JSON existen
    # Catálogo de riesgo personal
    if not os.path.exists(dir_path + "pr_risk_genes.json"):
        print("Generando archivos JSON y BED para riesgo personal.") # mejor dentro de la función get_json
        get_json_bed("rp", assembly, dir_path)
        
    # Catálogo de riesgo reproductivo
    if not os.path.exists(dir_path + "rr_risk_genes.json"):
        print("Generando archivos JSON y BED para riesgo reproductivo.") # mejor dentro de la función get_json
        get_json_bed("rr", assembly, dir_path)
        
    # Catálogo de riesgo farmacogenético
    if not os.path.exists(dir_path + "fg_risk_genes.json"):
        print("Generando archivos JSON y BED para riesgo farmacogenético.") # mejor dentro de la función get_json
        get_json_bed_fg(assembly, dir_path)
        
    
    """
    Get ClinVar database
    """
    # Si el modo es avanzado, comprobar si se ha descargado la BD ClinVar
    if mode == 'advanced':
        clinvar_files = [file for file in os.listdir(dir_path) if file.startswith("clinvar_database_")]
        
        # Si hay archivos clinvar, seleccionar el más reciente
        if clinvar_files:
            clinvar_files.sort(reverse=True)
            last_clinvar = clinvar_files[0]

            # Obtener la versión del nombre del archivo
            last_version = last_clinvar.split('_')[3].split('.')[0]      
            #print(f"El archivo ClinVar más reciente encontrado es {archivo_mas_reciente}.")
            print(f"La versión actual del archivo ClinVar es {last_version}.")
        
            # Preguntar al usuario si desea actualizar el archivo
            answr = input("¿Deseas actualizarlo? (S/N): ")       
            if answr.lower() == "s":
                # Descargar el archivo actualizado
                get_clinvar(dir_path)
            else:
                print("No se actualizará el archivo ClinVar.")
        # Si no se encuentran archivos ClinVar, descargarlo y guardarlo
        else:
            print("No se encontraron archivos ClinVar en el directorio.")
            print("El archivo ClinVar se descargará.")
            get_clinvar(dir_path)


    """
    Comprobar dependencias?
    """
    # Comprobar si InterVar está en el path
    if not os.path.exists(dir_path + "InterVar/"):
        print("InterVar no está instalado. Por favor, instálalo para continuar.")

    
    """
    Normalizar VCF de entrada
    """
    norm_vcf = normalize_vcf(vcf_file, dir_path)
    
    """
    Realizar la intersección con los archivos BED
    """
    for category in categories:
        intersect_vcf_bed(norm_vcf, category, assembly)
    
    """
    Ejecutar los módulos que correspondan:
    """
    # Verificar y ejecutar los módulos elegidos por el usuario
    if "pr" in categories:
        # Ejecutar el módulo de riesgo personal (PR)
        print("Ejecutando módulo de riesgo personal...")
        pr_results = run_personal_risk_module(vcf_path, assembly, mode, evidence, category, clinvar_path)
    
    if "rr" in categories:
        # Ejecutar el módulo de riesgo reproductivo (RR)
        print("Ejecutando módulo de riesgo reproductivo...")
        rr_results = run_reproductive_risk_module(vcf_path, assembly, mode, evidence, category, clinvar_path)
        
    if "fg" in categories:
        # Ejecutar el módulo farmacogenético (FG)
        print("Ejecutando módulo farmacogenético...")
        fg_results = fg_module(assembly, dir_path)
    
    # Informar al usuario que los módulos han sido ejecutados
    print("Módulos de análisis completados.")
    
    
    """
    Generar el informe de salida
    """
    out_file = write_report(pr_results, rr_results, fg_results, out_path)
    print("Informe de resultados generado.\n ---Finalizado---")

    
        
if __name__ == "__main__":
    main()
