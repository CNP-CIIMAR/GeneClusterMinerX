####Autor: Leandro de Mattos Pereira, Junior Researcher 
######## CNP Laboratory, Pedro Leao - Team Leader - Date: June, 06, 2023

import os
import shutil

output_dir = "./output_antismash"
log_file = "log.txt"

# Verifica se o diretório de saída existe
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Verifica os arquivos já processados
processed_files = set()
for result_dir in os.listdir(output_dir):
    if result_dir.startswith("Result_") and os.path.isdir(os.path.join(output_dir, result_dir)):
        fasta_file = result_dir.replace("Result_", "")
        processed_files.add(fasta_file)

# Função para salvar as mensagens no arquivo de log
def save_log(message):
    with open(log_file, "a") as f:
        f.write(message + "\n")

# Processa o arquivo .fna que não foi completamente processado
for fasta in processed_files:
    result_dir = os.path.join(output_dir, "Result_" + fasta)
    index_html_path = os.path.join(result_dir, "index.html")

    if not os.path.exists(index_html_path):
        message = f"O arquivo {fasta} foi parcialmente processado anteriormente e o index.html está faltando. Refazendo o processamento."
        save_log(message)
        shutil.rmtree(result_dir)  # Remove o diretório

        # Processa novamente o arquivo .fna parcialmente processado
        command = f"antismash --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --cb-subclusters --cb-knownclusters --rre --genefinding-tool prodigal {fasta} --fullhmmer --output-dir {result_dir}"
        os.system(command)

# Processa os arquivos .fna restantes
for file in os.listdir("."):
    if file.endswith(".fasta"):
        fasta = file
        result_dir = os.path.join(output_dir, "Result_" + fasta)

        # Verifica se o arquivo já foi processado
        if fasta in processed_files:
            index_html_path = os.path.join(result_dir, "index.html")

            if os.path.exists(index_html_path):
                message = f"O arquivo {fasta} já foi processado anteriormente e o index.html existe."
                save_log(message)
                continue

        # Processa o arquivo .fna se ainda não foi processado completamente
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
            command = f"antismash --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --cb-subclusters --cb-knownclusters --rre --genefinding-tool prodigal {fasta} --fullhmmer --output-dir {result_dir}"
            os.system(command)
            message = f"Processado o arquivo {fasta}."
            save_log(message)
