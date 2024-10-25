import os
import shutil
import subprocess
import multiprocessing

output_dir = "./output_antismash"
log_file = "log.txt"

# Verifica se o diretório de saída existe, caso contrário, cria-o
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

# Cria um Lock para sincronizar o acesso ao arquivo de log
log_lock = multiprocessing.Lock()

# Função para salvar as mensagens no arquivo de log
def save_log(message):
    with log_lock:
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
            if stdout_line:
                save_log(stdout_line.strip())  # Escreve cada linha de stdout no log
        for stderr_line in iter(process.stderr.readline, ""):
            if stderr_line:
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

# Função que processa um único arquivo .fna
def process_file(fasta):
    result_dir = os.path.join(output_dir, "Result_" + fasta)
    index_html_path = os.path.join(result_dir, "index.html")

    # Verifica se o arquivo já foi processado completamente
    if os.path.exists(index_html_path):
        message = f"O arquivo {fasta} já foi processado anteriormente e o index.html existe."
        save_log(message)
        return

    # Remove diretórios parcialmente processados
    if os.path.exists(result_dir):
        shutil.rmtree(result_dir)
        save_log(f"Removido diretório parcialmente processado: {result_dir}")

    # Cria o diretório de resultado
    os.mkdir(result_dir)
    
    # Executa o antismash para o arquivo
    run_antismash(fasta, result_dir)
    message = f"Processado o arquivo {fasta}."
    save_log(message)

if __name__ == "__main__":
    # Ajustar permissões apenas uma vez
    adjust_permissions(output_dir)

    # Lista todos os arquivos .fna
    all_fna_files = [file for file in os.listdir(".") if file.endswith(".fna")]

    # Identifica os arquivos já processados completamente
    processed_files = set()
    for result_dir_name in os.listdir(output_dir):
        if result_dir_name.startswith("Result_") and os.path.isdir(os.path.join(output_dir, result_dir_name)):
            fasta_file = result_dir_name.replace("Result_", "")
            index_html_path = os.path.join(output_dir, result_dir_name, "index.html")
            if os.path.exists(index_html_path):
                processed_files.add(fasta_file)
            else:
                # Remove diretórios parcialmente processados
                shutil.rmtree(os.path.join(output_dir, result_dir_name))
                save_log(f"Removido diretório parcialmente processado: {result_dir_name}")

    # Determina os arquivos que precisam ser processados
    files_to_process = [f for f in all_fna_files if f not in processed_files]

    # Define o número de processos paralelos (ajustado para 80)
    num_processes = 80  # <-- Aqui você define para 80 processos

    # Cria um pool de processos para processar os arquivos em paralelo
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(process_file, files_to_process)
