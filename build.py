import zipfile
import os


# List of folder and file names you want to compress
folders = ['bin', 'src', 'config', 'INSTALLATION.md', 'MODULES.md', 'README.md', 'TESTS.md', 'changelog.md', 'install_packages.R', 'python_requirements.txt']

exclude = ["\\test", "\\samples", "\\__pycache__"]

# Name of the output ZIP file
zip_name = 'PTM-Analyzer-v1.08.zip'

# Create ZIP file
with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as archivo_zip:
    for item in folders:

        # Case 1: If it is a file => add directly
        if os.path.isfile(item):
            print(f"Adding: {item}")
            archivo_zip.write(item, item)
            continue

        # Case 2: If it is a directory => traverse with os.walk
        if os.path.isdir(item):
            # Traverse each directory and add its contents to the ZIP file
            for root, dirs, files in os.walk(item):
                for archivo in files:
                    if any([i in root for i in exclude]): # apply exclusions
                        continue

                    ruta_completa = os.path.join(root, archivo)
                    ruta_zip = ruta_completa  # you can adjust this if you want relative paths
                    print(f"Adding: {ruta_completa}")
                    archivo_zip.write(ruta_completa, ruta_zip)

print(f'ZIP file"{zip_name}" created successfully.')
