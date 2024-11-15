# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import multiprocessing
import logging
import sys
from functools import partial
import argparse

# Configuração do logging
logging.basicConfig(
    filename='plantismash.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

output_dir = "./output_plantismash"
log_file = "plantismash.log"

# Função para criar diretório se não existir
def create_directory(directory):
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
            logging.info("Criado o diretório: {}".format(directory))
        except Exception as e:
            logging.error("Erro ao criar o diretório {}: {}".format(directory, e))
            sys.exit(1)
    else:
        logging.info("Usando o diretório existente: {}".format(directory))

create_directory(output_dir)

# Função para rodar o PlantiSMASH no ambiente Conda correto e capturar as saídas
def run_antismash(fasta, result_dir, conda_env, plantismash_script, cpu):
    lock_file = os.path.join(result_dir, ".lock")
    if os.path.exists(lock_file):
        logging.warning("Arquivo {} já está sendo processado por outro processo.".format(fasta))
        return

    # Cria um arquivo de lock
    try:
        with open(lock_file, 'w') as lf:
            lf.write("Locked by PID {}\n".format(os.getpid()))
    except Exception as e:
        logging.error("Não foi possível criar o arquivo de lock para {}: {}".format(fasta, e))
        return

    logging.info("Iniciando processamento do arquivo {}.".format(fasta))

    command = (
        "bash -c 'source {}/bin/activate plantismash && "
        "python2 {} --input-type nucl --taxon plants --cdhit --clusterblast --subclusterblast "
        "--knownclusterblast --smcogs --inclusive --borderpredict --full-hmmer --asf "
        "--cpu {} --outputfolder {} {}'".format(
            conda_env, plantismash_script, cpu, result_dir, fasta
        )
    )

    logging.info("Executando comando: {}".format(command))

    try:
        # Usando subprocess.Popen para capturar stdout e stderr
        process = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )

        # Captura a saída em tempo real
        for stdout_line in iter(process.stdout.readline, ""):
            if stdout_line:
                logging.info("[{}] {}".format(fasta, stdout_line.strip()))
        for stderr_line in iter(process.stderr.readline, ""):
            if stderr_line:
                logging.error("[{}] {}".format(fasta, stderr_line.strip()))

        process.stdout.close()
        process.stderr.close()

        return_code = process.wait()
        if return_code != 0:
            logging.error("Erro ao processar o arquivo {}: Retornou código {}".format(fasta, return_code))
        else:
            logging.info("Processamento do arquivo {} concluído com sucesso.".format(fasta))
    except subprocess.CalledProcessError as e:
        logging.error("Erro ao processar o arquivo {}: {}".format(fasta, e))
    except Exception as e:
        logging.error("Erro inesperado ao processar o arquivo {}: {}".format(fasta, e))
    finally:
        # Remove o arquivo de lock após o processamento
        try:
            os.remove(lock_file)
        except Exception as e:
            logging.warning("Não foi possível remover o arquivo de lock para {}: {}".format(fasta, e))

# Função que processa um único arquivo .fna
def process_file(fasta, output_dir, conda_env, plantismash_script, cpu):
    # Sanitiza o nome do diretório de resultado
    safe_fasta = os.path.splitext(fasta)[0].replace("/", "_").replace("\\", "_")
    result_dir = os.path.join(output_dir, "Result_" + safe_fasta)
    index_html_path = os.path.join(result_dir, "index.html")

    # Verifica se o arquivo já foi processado completamente
    if os.path.exists(index_html_path):
        logging.info("O arquivo {} já foi processado anteriormente e o index.html existe.".format(fasta))
        return

    # Remove diretórios parcialmente processados
    if os.path.exists(result_dir):
        try:
            shutil.rmtree(result_dir)
            logging.info("Removido diretório parcialmente processado: {}".format(result_dir))
        except Exception as e:
            logging.error("Erro ao remover diretório {}: {}".format(result_dir, e))
            return

    # Cria o diretório de resultado
    try:
        os.makedirs(result_dir)
    except Exception as e:
        logging.error("Erro ao criar diretório de resultado {}: {}".format(result_dir, e))
        return

    # Executa o PlantiSMASH para o arquivo
    run_antismash(fasta, result_dir, conda_env, plantismash_script, cpu)

# Função para identificar e preparar os arquivos a serem processados
def get_files_to_process(input_dir, output_dir):
    all_fna_files = [file for file in os.listdir(input_dir) if file.endswith(".fna")]
    processed_files = set()

    for result_dir_name in os.listdir(output_dir):
        if result_dir_name.startswith("Result_") and os.path.isdir(os.path.join(output_dir, result_dir_name)):
            fasta_file = result_dir_name.replace("Result_", "")
            index_html_path = os.path.join(output_dir, result_dir_name, "index.html")
            if os.path.exists(index_html_path):
                processed_files.add(fasta_file)
            else:
                # Remove diretórios parcialmente processados
                try:
                    shutil.rmtree(os.path.join(output_dir, result_dir_name))
                    logging.info("Removido diretório parcialmente processado: {}".format(result_dir_name))
                except Exception as e:
                    logging.error("Erro ao remover diretório {}: {}".format(result_dir_name, e))

    # Determina os arquivos que precisam ser processados
    files_to_process = [f for f in all_fna_files if f not in processed_files]

    # Remove possíveis duplicatas
    files_to_process = list(set(files_to_process))

    logging.info("Total de arquivos para processar: {}".format(len(files_to_process)))
    return files_to_process

# Função principal
def main():
    # Configuração dos argumentos de linha de comando
    parser = argparse.ArgumentParser(description="Script para processar arquivos .fna com PlantiSMASH.")
    parser.add_argument(
        '-p', '--num-processes',
        type=int,
        default=96,
        help="Número de processos paralelos (default: 96)"
    )
    parser.add_argument(
        '-c', '--cpus-per-process',
        type=int,
        default=98,
        help="Número de CPUs por processo (default: 98)"
    )
    parser.add_argument(
        '-i', '--input-dir',
        type=str,
        default=".",
        help="Diretório de entrada contendo arquivos .fna (default: diretório atual)"
    )
    parser.add_argument(
        '-e', '--conda-env',
        type=str,
        default="/home/mattoslmp/anaconda3",
        help="Caminho para o ambiente Conda (default: /home/mattoslmp/anaconda3)"
    )
    parser.add_argument(
        '-s', '--plantismash-script',
        type=str,
        default="/home/mattoslmp/plantismash-1.0/run_antismash.py",
        help="Caminho para o script run_antismash.py do PlantiSMASH (default: /home/mattoslmp/plantismash-1.0/run_antismash.py)"
    )

    args = parser.parse_args()

    # Validação dos argumentos
    if args.num_processes < 1:
        logging.error("O número de processos deve ser pelo menos 1.")
        sys.exit(1)
    if args.cpus_per_process < 1:
        logging.error("O número de CPUs por processo deve ser pelo menos 1.")
        sys.exit(1)
    total_requested_cpus = args.num_processes * args.cpus_per_process
    total_available_cpus = multiprocessing.cpu_count()
    if total_requested_cpus > total_available_cpus:
        logging.warning(
            "O número total de CPUs solicitadas ({}) excede o disponível ({}) no sistema.".format(
                total_requested_cpus, total_available_cpus
            )
        )
        logging.warning("Pode haver sobrecarga de CPU. Considere reduzir o número de processos ou CPUs por processo.")

    logging.info("Configurações recebidas:")
    logging.info("Número de processos paralelos: {}".format(args.num_processes))
    logging.info("Número de CPUs por processo: {}".format(args.cpus_per_process))
    logging.info("Diretório de entrada: {}".format(args.input_dir))
    logging.info("Ambiente Conda: {}".format(args.conda_env))
    logging.info("Script do PlantiSMASH: {}".format(args.plantismash_script))

    # Identifica os arquivos a serem processados
    files_to_process = get_files_to_process(args.input_dir, output_dir)

    if not files_to_process:
        logging.info("Nenhum arquivo para processar. Encerrando o script.")
        sys.exit(0)

    # Configura a função parcial com os argumentos fixos
    process_func = partial(
        process_file,
        output_dir=output_dir,
        conda_env=args.conda_env,
        plantismash_script=args.plantismash_script,
        cpu=args.cpus_per_process
    )

    # Cria um pool de processos
    pool = multiprocessing.Pool(processes=args.num_processes)
    try:
        pool.map(process_func, files_to_process)
    except KeyboardInterrupt:
        logging.warning("Execução interrompida pelo usuário.")
        pool.terminate()
        pool.join()
        sys.exit(1)
    except Exception as e:
        logging.error("Erro durante o processamento paralelo: {}".format(e))
        pool.terminate()
        pool.join()
        sys.exit(1)
    else:
        pool.close()
        pool.join()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.warning("Execução interrompida pelo usuário.")
        sys.exit(1)
    except Exception as e:
        logging.error("Ocorreu um erro inesperado: {}".format(e))
        sys.exit(1)
