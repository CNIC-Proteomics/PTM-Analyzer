import zipfile
import os

# Lista de nombres de las carpetas que quieres comprimir
carpetas = ['bin', 'src', 'config', 'INSTALLATION.md', 'MODULES.md', 'README.md', 'TESTS.md', 'changelog.md', 'install_packages.R', 'python_requirements.txt']

exclude = ["\\test", "\\samples", "\\__pycache__"]

# Nombre del archivo ZIP de salida
nombre_zip = 'ReportAnalysis-v1.05.zip'

# Crear un archivo ZIP
with zipfile.ZipFile(nombre_zip, 'w', zipfile.ZIP_DEFLATED) as archivo_zip:
    for carpeta in carpetas:
        # Recorre cada carpeta y agrega su contenido al archivo ZIP
        for root, dirs, files in os.walk(carpeta):
            for archivo in files:
                if any([i in root for i in exclude]):
                    continue
                print(archivo)
                ruta_completa = os.path.join(root, archivo)
                # Determina la ruta relativa al archivo ZIP
                #ruta_zip = os.path.relpath(ruta_completa, carpeta)
                archivo_zip.write(ruta_completa, ruta_completa)

print(f'Archivo ZIP "{nombre_zip}" creado exitosamente.')
