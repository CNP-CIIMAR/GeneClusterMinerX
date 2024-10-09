import os
import shutil
import subprocess

output_dir = "./output_antismash"
log_file = "log.txt"

# Verifica se o diretório de saída existe, caso contrário, usa o existente
if not os.path.exists(output_dir):
    try:
        os.mkdir(output_dir)
        print(f"Criado o diretório: {output_dir}")
    except FileExistsError:
        print(f"O diretório {output_dir} já existe. Continuando...")
else:
    print(f"O diretório {output_dir} já existe. Usando o diretório existente.")

# Função para ajustar permissões com sudo (executado apenas uma vez)
def adjust_permissions(directory):
    command = f"sudo chmod -R u+rwX {directory}"
    subprocess.run(command, shell=True, check=True)

# Função para salvar as mensagens no arquivo de log
def save_log(message):
    with open(log_file, "a") as f:
        f.write(message + "\n")
        f.flush()  # Garante que o conteúdo seja gravado no arquivo imediatamente

# Função para rodar o antismash no ambiente Conda correto e capturar as saídas
def run_antismash(fasta, result_dir):
    save_log(f"Iniciando processamento do arquivo {fasta}.")
    command = f"bash -c 'source /home/mattoslmp/anaconda3/bin/activate antismash && antismash --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --cb-subclusters --cb-knownclusters --rre --genefinding-tool prodigal {fasta} -c 8 --fullhmmer --output-dir {result_dir}'"
    
    # Loga o comando que será executado
    save_log(f"Executando comando: {command}")
    
    try:
        # Usando subprocess.Popen para capturar stdout e stderr
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Leitura da saída em tempo real
        for stdout_line in iter(process.stdout.readline, ""):
            save_log(stdout_line.strip())  # Escreve cada linha de stdout no log
        for stderr_line in iter(process.stderr.readline, ""):
            save_log(stderr_line.strip())  # Escreve cada linha de stderr no log

        process.stdout.close()
        process.stderr.close()
        
        return_code = process.wait()
        if return_code != 0:
            save_log(f"Erro ao processar o arquivo {fasta}: Retornou código {return_code}")
        else:
            save_log(f"Processamento do arquivo {fasta} concluído com sucesso.")
    except subprocess.CalledProcessError as e:
        save_log(f"Erro ao processar o arquivo {fasta}: {e}")

# Definir a variável processed_files como um conjunto vazio
processed_files = set()

# Ajustar permissões apenas uma vez
adjust_permissions(output_dir)

# Verifica os arquivos já processados
for result_dir in os.listdir(output_dir):
    if result_dir.startswith("Result_") and os.path.isdir(os.path.join(output_dir, result_dir)):
        fasta_file = result_dir.replace("Result_", "")
        processed_files.add(fasta_file)

# Processa o arquivo .fna que não foi completamente processado
for fasta in processed_files:
    result_dir = os.path.join(output_dir, "Result_" + fasta)
    index_html_path = os.path.join(result_dir, "index.html")

    if not os.path.exists(index_html_path):
        message = f"O arquivo {fasta} foi parcialmente processado anteriormente e o index.html está faltando. Refazendo o processamento."
        save_log(message)
        
        # Verifica se o diretório existe antes de removê-lo
        if os.path.exists(result_dir):
            shutil.rmtree(result_dir)  # Remove o diretório

        # Processa novamente o arquivo .fna parcialmente processado
        run_antismash(fasta, result_dir)
        message = f"Reprocessado o arquivo {fasta}."
        save_log(message)

# Processa os arquivos .fna restantes
for file in os.listdir("."):
    if file.endswith(".fna"):
        fasta = file
        result_dir = os.path.join(output_dir, "Result_" + fasta)

        # Verifica se o arquivo já foi processado
        if fasta in processed_files:
            index_html_path = os.path.join(result_dir, "index.html")

            if os.path.exists(index_html_path):
                message = f"O arquivo {fasta} já foi processado anteriormente e o index.html existe."
                save_log(message)
                continue

        # Verifica se o diretório já existe antes de criá-lo
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        else:
            save_log(f"Diretório {result_dir} já existe. Continuando...")
        
        run_antismash(fasta, result_dir)
        message = f"Processado o arquivo {fasta}."
        save_log(message)

