#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Autor: Leandro de Mattos Pereira
Junior Researcher, CNP Laboratory
Pedro Leao - Team Leader
Data: Junho, 06, 2023 (Atualizado: Outubro, 03, 2024)
Descrição: Script para executar antiSMASH em múltiplos arquivos .fna ou .fasta dentro de um diretório,
           gerando um diretório de resultados para cada arquivo de entrada, com suporte a processamento paralelo.
"""

import os
import shutil
import argparse
import subprocess
from pathlib import Path
import logging
import sys
import multiprocessing
from functools import partial

# Diretório de saída e arquivo de log padrão
DEFAULT_LOG_FILE = "logs.txt"

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Script para executar antiSMASH em múltiplos arquivos .fna ou .fasta dentro de um diretório com opções personalizadas."
    )

    # Argumentos posicionais
    parser.add_argument(
        "input_dir",
        type=str,
        help="Diretório contendo arquivos .fna ou .fasta para análise."
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Diretório onde os resultados serão salvos."
    )

    # Opções de ajuda do antiSMASH
    parser.add_argument(
        "--antismash-help",
        action="store_true",
        help="Exibe a ajuda do antiSMASH e sai."
    )

    # Opções básicas de análise
    parser.add_argument(
        "-t", "--taxon",
        choices=["bacteria", "fungi"],
        default="bacteria",
        help="Classificação taxonômica da sequência de entrada (padrão: bacteria)."
    )
    parser.add_argument(
        "-c", "--cpus",
        type=int,
        default=4,  # Valor padrão ajustado para 4, pode ser alterado conforme a necessidade
        help="Número de CPUs a serem usadas em paralelo pelo antiSMASH (padrão: 4)."
    )
    parser.add_argument(
        "--databases",
        type=str,
        default="/home/mattoslmp/anaconda3/envs/antismash/lib/python3.9/site-packages/antismash/databases",
        help="Diretório raiz das bases de dados usadas pelo antiSMASH."
    )

    # Opções de saída adicionais
    parser.add_argument(
        "--output-basename",
        type=str,
        help="Nome base para arquivos de saída dentro do diretório de saída."
    )
    parser.add_argument(
        "--html-title",
        type=str,
        help="Título personalizado para a página de saída em HTML."
    )
    parser.add_argument(
        "--html-description",
        type=str,
        help="Descrição personalizada para adicionar à saída."
    )
    parser.add_argument(
        "--html-start-compact",
        action="store_true",
        help="Usa a visualização compacta por padrão na página de visão geral."
    )
    parser.add_argument(
        "--html-ncbi-context",
        dest="html_ncbi_context",
        action="store_true",
        help="Mostra links para o contexto genômico NCBI dos genes."
    )
    parser.add_argument(
        "--no-html-ncbi-context",
        dest="html_ncbi_context",
        action="store_false",
        help="Não mostra links para o contexto genômico NCBI dos genes."
    )
    parser.set_defaults(html_ncbi_context=False)

    # Análises adicionais
    parser.add_argument(
        "--all",
        action="store_true",
        help="Ativa todas as análises disponíveis no antiSMASH."
    )
    parser.add_argument("--fullhmmer", action="store_true", help="Executa análise HMMer em todo o genoma usando perfis Pfam.")
    parser.add_argument("--cassis", action="store_true", help="Predição baseada em motivos de regiões de clusters de genes de SM.")
    parser.add_argument("--clusterhmmer", action="store_true", help="Executa análise HMMer limitada a clusters usando perfis Pfam.")
    parser.add_argument("--tigrfam", action="store_true", help="Anota clusters usando perfis TIGRFam.")
    parser.add_argument("--asf", action="store_true", help="Executa análise de sítios ativos.")
    parser.add_argument("--cc-mibig", action="store_true", help="Compara clusters identificados com o banco de dados MIBiG.")
    parser.add_argument("--cb-general", action="store_true", help="Compara clusters identificados com uma base de dados de clusters preditos pelo antiSMASH.")
    parser.add_argument("--cb-subclusters", action="store_true", help="Compara clusters identificados com subclusters conhecidos que sintetizam precursores.")
    parser.add_argument("--cb-knownclusters", action="store_true", help="Compara clusters com clusters conhecidos da base de dados MIBiG.")
    parser.add_argument("--pfam2go", action="store_true", help="Mapeia Pfam para Gene Ontology.")
    parser.add_argument("--rre", action="store_true", help="Executa o RREFinder em modo de precisão em todos os clusters de RiPP.")
    parser.add_argument("--smcog-trees", action="store_true", help="Gera árvores filogenéticas de grupos ortólogos de clusters de metabólitos secundários.")
    parser.add_argument("--tfbs", action="store_true", help="Executa o localizador de sítios de ligação de fatores de transcrição (TFBS) em todos os clusters.")
    parser.add_argument("--tta-threshold", type=float, default=0.65, help="Menor conteúdo GC para anotar códons TTA (padrão: 0.65).")

    # Opções de previsão de genes
    parser.add_argument(
        "--genefinding-tool",
        choices=["glimmerhmm", "prodigal", "prodigal-m", "none", "error"],
        default="error",
        help="Ferramenta de previsão de genes a ser usada (padrão: error)."
    )
    parser.add_argument(
        "--genefinding-gff3",
        type=str,
        help="Especifica um arquivo GFF3 para extrair características."
    )

    # Opções de log
    parser.add_argument(
        "--log-file",
        type=str,
        default=DEFAULT_LOG_FILE,
        help=f"Arquivo de log para salvar mensagens (padrão: {DEFAULT_LOG_FILE})."
    )

    # Opção para definir o número de processos paralelos
    parser.add_argument(
        "--parallel-processes",
        type=int,
        default=4,
        help="Número de processos paralelos para executar o antiSMASH em diferentes arquivos (padrão: 4)."
    )

    args = parser.parse_args()

    # Se o usuário solicitar ajuda do antiSMASH, exiba e saia
    if args.antismash_help:
        help_command = ["antismash", "--help"]
        subprocess.run(help_command)
        sys.exit()

    # Se --all for especificado, habilita todas as análises adicionais
    if args.all:
        args.fullhmmer = True
        args.cassis = True
        args.clusterhmmer = True
        args.tigrfam = True
        args.asf = True
        args.cc_mibig = True
        args.cb_general = True
        args.cb_subclusters = True
        args.cb_knownclusters = True
        args.pfam2go = True
        args.rre = True
        args.smcog_trees = True
        args.tfbs = True

    return args

def setup_logging(log_file):
    """
    Configura o módulo de logging para registrar mensagens com timestamps.
    """
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Também adiciona um handler para exibir as mensagens no console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', '%Y-%m-%d %H:%M:%S')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def process_file(fasta_path, args, lock):
    """
    Processa um único arquivo fasta usando o antiSMASH.
    """
    try:
        fasta = fasta_path.name
        # Remove a extensão do arquivo para nomear o diretório de resultados
        fasta_stem = fasta_path.stem
        result_dir = args.output_dir / f"Result_{fasta_stem}"
        index_html_path = result_dir / "index.html"

        # Verifica se o diretório de resultados já existe
        if result_dir.exists():
            if index_html_path.exists():
                with lock:
                    logging.info(f"O arquivo '{fasta}' já foi processado anteriormente e o 'index.html' existe. Pulando.")
                return (fasta, "já processado")
            else:
                with lock:
                    logging.warning(f"O arquivo '{fasta}' foi parcialmente processado anteriormente e o 'index.html' está faltando. Refazendo o processamento.")
                try:
                    shutil.rmtree(result_dir)
                    with lock:
                        logging.info(f"Diretório '{result_dir}' removido para reprocessamento.")
                except Exception as e:
                    with lock:
                        logging.error(f"Erro ao remover o diretório '{result_dir}': {e}")
                    return (fasta, "falha")

        # Cria o diretório de resultados
        try:
            result_dir.mkdir(parents=True, exist_ok=True)
            with lock:
                logging.info(f"Diretório de resultados criado: '{result_dir}'.")
        except Exception as e:
            with lock:
                logging.error(f"Erro ao criar o diretório '{result_dir}': {e}")
            return (fasta, "falha")

        # Ajusta a taxonomia se a ferramenta de previsão de genes for glimmerhmm
        taxon = args.taxon  # Inicializa com o valor fornecido
        if args.genefinding_tool == "glimmerhmm":
            taxon = "fungi"
            with lock:
                logging.info("Ferramenta de previsão de genes 'glimmerhmm' selecionada. Ajustando taxonomia para 'fungi'.")

        # Monta o comando antiSMASH
        command = [
            "antismash",
            str(fasta_path),
            "--taxon", taxon,
            "--cpus", str(args.cpus),
            "--databases", args.databases,
            "--output-dir", str(result_dir)
        ]

        # Adiciona opções de saída
        if args.output_basename:
            command += ["--output-basename", args.output_basename]
        if args.html_title:
            command += ["--html-title", args.html_title]
        if args.html_description:
            command += ["--html-description", args.html_description]
        if args.html_start_compact:
            command.append("--html-start-compact")
        if args.html_ncbi_context:
            command.append("--html-ncbi-context")
        else:
            command.append("--no-html-ncbi-context")

        # Adiciona análises adicionais
        if args.fullhmmer:
            command.append("--fullhmmer")
        if args.cassis:
            # Verifica se a ferramenta de previsão de genes é 'prodigal'
            if args.genefinding_tool == "prodigal":
                with lock:
                    logging.warning("CASSIS desativado porque a ferramenta de previsão de genes é 'prodigal'.")
            else:
                command.append("--cassis")
        if args.clusterhmmer:
            command.append("--clusterhmmer")
        if args.tigrfam:
            command.append("--tigrfam")
        if args.asf:
            command.append("--asf")
        if args.cc_mibig:
            command.append("--cc-mibig")
        if args.cb_general:
            command.append("--cb-general")
        if args.cb_subclusters:
            command.append("--cb-subclusters")
        if args.cb_knownclusters:
            command.append("--cb-knownclusters")
        if args.pfam2go:
            command.append("--pfam2go")
        if args.rre:
            command.append("--rre")
        if args.smcog_trees:
            command.append("--smcog-trees")
        if args.tfbs:
            command.append("--tfbs")
        if args.tta_threshold:
            command += ["--tta-threshold", str(args.tta_threshold)]

        # Adiciona opções de previsão de genes
        command += ["--genefinding-tool", args.genefinding_tool]
        if args.genefinding_gff3:
            command += ["--genefinding-gff3", args.genefinding_gff3]

        # Executa o comando antiSMASH
        try:
            with lock:
                logging.info(f"Iniciando o processamento do arquivo '{fasta}' com antiSMASH.")
            subprocess.run(command, check=True)
            with lock:
                logging.info(f"Processado o arquivo '{fasta}' com sucesso.")
            return (fasta, "sucesso")
        except subprocess.CalledProcessError as e:
            with lock:
                logging.error(f"Erro ao processar o arquivo '{fasta}': {e}")
            return (fasta, "falha")
        except Exception as e:
            with lock:
                logging.error(f"Erro inesperado ao processar o arquivo '{fasta}': {e}")
            return (fasta, "falha")
    except Exception as e:
        with lock:
            logging.error(f"Erro inesperado no processamento do arquivo '{fasta_path}': {e}")
        return (fasta_path.name, "falha")

def main():
    args = parse_arguments()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    log_file = Path(args.log_file).resolve()
    args.output_dir = output_dir  # Atualiza args.output_dir para ser um Path

    # Configura o logging
    setup_logging(log_file)
    logging.info("Início do processamento.")

    # Verifica se o diretório de entrada existe
    if not input_dir.is_dir():
        logging.error(f"O diretório de entrada '{input_dir}' não existe ou não é um diretório.")
        sys.exit(1)

    # Cria o diretório de saída se não existir
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Diretório de saída configurado em '{output_dir}'.")
    except Exception as e:
        logging.error(f"Erro ao criar o diretório de saída '{output_dir}': {e}")
        sys.exit(1)

    # Lista todos os arquivos .fna e .fasta no diretório de entrada
    fna_files = list(input_dir.glob("*.fna")) + list(input_dir.glob("*.fasta"))
    if not fna_files:
        logging.warning(f"Nenhum arquivo .fna ou .fasta encontrado no diretório '{input_dir}'.")
        sys.exit(1)

    logging.info(f"Encontrados {len(fna_files)} arquivos para processamento.")

    # Inicializa contadores
    manager = multiprocessing.Manager()
    lock = manager.Lock()
    results = []

    # Prepara a função parcial para multiprocessing
    process_file_partial = partial(process_file, args=args, lock=lock)

    # Define o número de processos paralelos
    num_processes = args.parallel_processes

    # Cria um pool de processos
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(process_file_partial, fna_files)

    # Calcula o resumo
    success_count = sum(1 for _, status in results if status == "sucesso")
    failure_count = sum(1 for _, status in results if status == "falha")
    already_processed_count = sum(1 for _, status in results if status == "já processado")

    # Finaliza o processamento com um resumo
    logging.info("Processamento concluído.")
    logging.info(f"Total de genomas processados com sucesso: {success_count}")
    logging.info(f"Total de genomas que falharam no processamento: {failure_count}")
    logging.info(f"Total de genomas já processados anteriormente: {already_processed_count}")

    # Opcional: Exibe o resumo no console
    print("\nResumo do Processamento:")
    print(f"Total de genomas processados com sucesso: {success_count}")
    print(f"Total de genomas que falharam no processamento: {failure_count}")
    print(f"Total de genomas já processados anteriormente: {already_processed_count}")

if __name__ == "__main__":
    main()
