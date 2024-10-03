# GeneClusterMineX v 2.0 :rocket:
_____________________________________________________________________________________________________________________________________________________

 - **v 1.0.0: automation of processing of several Genomes fasta/fna files by antismash**



run_antismash.py

## Automated Secondary Metabolite Analysis with antiSMASH
**run_antismash.py is a Python script designed to streamline the execution of antiSMASH on multiple sequence files (.fna or .fasta) within a directory. The script automates the creation of output directories for each input file, manages detailed logs, and offers flexibility to customize analyses according to user needs.**

- Batch Processing: Executes antiSMASH on multiple sequence files simultaneously.
## Organized Results: Creates specific output directories for each input file, named as Result_input_filename.
Comprehensive Logging: Records detailed logs for each processing step, including timestamps and error messages.
## Flexible Analyses: Allows activating all available antiSMASH analysis tools or selecting specific analyses.
## GlimmerHMM Support: Integrates the GlimmerHMM gene prediction tool, automatically adjusting taxonomy to fungi.


## Requirements

- Operating System: Linux or Unix-based
- Python: Version 3.6 or higher
- antiSMASH: Properly installed and configured
- Appropriate Permissions: Read permissions for input files and write permissions for output and log directories

1. Clone the repository:

 ```git clone https://github.com/mattoslmp/CNP-Ciimar.git ```

2. Navigate to the cloned directory.

3. Create a Conda environment and install antiSMASH:

- conda create -n antiSMASH-env
- conda activate antiSMASH-env
- install dependencies:
- conda install hmmer2 hmmer diamond fasttree prodigal blast muscle glimmerhmm
- install python (versions 3.9, 3.10, and 3.11 tested, any version >= 3.9.0 should work):
- conda install -c conda-forge python
- wget https://dl.secondarymetabolites.org/releases/7.0.0/antismash-7.0.0.tar.gz
- tar -zxf antismash-7.0.0.tar.gz
- pip install ./antismash-7.0.0
- download-antismash-databases


## Usage
- Directory Structure
Input Directory (input_dir): Contains all .fna or .fasta files you wish to analyze.

```bash
./fna_inputs/
├── sample1.fna
├── sample2.fasta
└── ...

```

## Output Directory (output_dir): Where the results will be stored. Each input file will have its own results subdirectory.

---

## Output and Logs

After running the script, the output directory will be organized as follows:

```bash
./output_antismash/
├── log.txt
├── Result_sample1/
│   ├── index.html
│   ├── ... (other output files)
├── Result_sample2/
│   ├── index.html
│   ├── ... (other output files)
└── ...



```shell
 ./run_antismash.py <input_dir> <output_dir> [options]
```

## How run antismahs for all tools available:

```shell
nohup python ./run_antismash.py ./ ./output_antismash/ --all --cpus 8 --genefinding-tool glimmerhmm  > run_output2.log 3>&1 &
```


- <input_dir>: Path to the directory containing .fna or .fasta files.
- <output_dir>: Path to the directory where results will be saved.


Available Options
Positional Arguments:

### Available Options

#### Positional Arguments:

- **input_dir**: Input directory with `.fna` or `.fasta` files.
- **output_dir**: Output directory to store results.

#### General Options:

- `--antismash-help`: Displays antiSMASH help and exits.
- `-t`, `--taxon`: Taxonomic classification of the input sequence. Options: `bacteria` (default), `fungi`.
- `-c`, `--cpus`: Number of CPUs to use in parallel (default: 4).
- `--databases`: Root directory of databases used by antiSMASH.

#### Additional Output Options:

- `--output-basename`: Base name for output files within the output directory.
- `--html-title`: Custom title for the HTML output page.
- `--html-description`: Custom description to add to the output.
- `--html-start-compact`: Uses compact view by default on the overview page.
- `--html-ncbi-context`: Shows links to NCBI genomic context of genes.
- `--no-html-ncbi-context`: Does not show links to NCBI genomic context of genes.

#### Additional Analyses:

- `--all`: Activates all available antiSMASH analyses.
- `--fullhmmer`: Executes HMMer analysis across the entire genome using Pfam profiles.
- `--cassis`: Prediction based on motifs of SM gene cluster regions.
- `--clusterhmmer`: Executes HMMer analysis limited to clusters using Pfam profiles.
- `--tigrfam`: Annotates clusters using TIGRFam profiles.
- `--asf`: Executes active site analysis.
- `--cc-mibig`: Compares identified clusters with the MIBiG database.
- `--cb-general`: Compares identified clusters with a database of antiSMASH-predicted clusters.
- `--cb-subclusters`: Compares identified clusters with known subclusters synthesizing precursors.
- `--cb-knownclusters`: Compares clusters with known clusters from the MIBiG database.
- `--pfam2go`: Maps Pfam to Gene Ontology.
- `--rre`: Executes RREFinder in precision mode on all RiPP clusters.
- `--smcog-trees`: Generates phylogenetic trees of orthologous groups of secondary metabolite clusters.
- `--tfbs`: Executes transcription factor binding site locator on all clusters.
- `--tta-threshold`: Minimum GC content to annotate TTA codons (default: 0.65).

#### Gene Prediction Options:

- `--genefinding-tool`: Gene prediction tool to use. Options: `glimmerhmm`, `prodigal`, `prodigal-m`, `none`, `error` (default).
- `--genefinding-gff3`: Specifies a GFF3 file to extract features.

#### Logging Options:

- `--log-file`: Log file to save messages (default: `log.txt` in the output directory).


---

## Monitoring and Management

### Checking Running Processes

You can verify if the script is running using `ps` or `pgrep`. For example:

```bash
ps aux | grep run_antismash.py

Or:

```shell
pgrep -fl run_antismash.py
```

## Monitoring Logs in Real-Time
Internal Script Logs (output_log.txt):

```shell
tail -f ./output_antismash/output_log.txt
```

```shell
tail -f run_output.log
```

## Stopping the Process
If you need to terminate the background process, follow these steps:

## Identify the PID (Process ID):

```shell
pgrep -f run_antismash.py
```

##Kill the Process:

```shell
kill -9 PID
```

# Replace PID with the actual process number.

## Final Considerations

- File Extensions Compatibility: The script is configured to process files with .fna and .fasta extensions. Ensure your input files have one of these extensions to be recognized by the script.

- System Resources: Specifying a high number of CPUs with the --cpus option can speed up processing but ensure your system has sufficient resources to handle the load without overloading.

Detailed Logs: Utilize the log files to monitor progress and identify any issues during processing.

Future Updates: If antiSMASH introduces new analysis options in the future, you will need to update the script to include these new options, both in argparse and in the logic that activates options with --all.

File Extensions: If you have compressed files (like .fasta.gz), you will need to decompress them before processing, as the current script does not support compressed files.

## License
This project is licensed under the MIT License.

## Comments:

- Make sure you have activated the `antiSMASH-env` Conda environment before running the script, using the `conda activate antiSMASH-env` command.
- Make sure you have the necessary permissions to create directories and delete files/directories.

## Contribution:

Contributions and suggestions are welcome! Feel free to open an issue or submit a pull request with improvements to this script.

## Please cite:

-  antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures, and visualisation
   Kai Blin, Simon Shaw, Hannah E Augustijn, Zachary L Reitz, Friederike Biermann, Mohammad Alanjary, Artem Fetter, Barbara R Terlouw, William W. Metcalf, Eric JN 
   Helfrich, Gilles P van Wezel, Marnix. H Medema & Tilmann Weber. Nucleic Acids Research (2023) doi: 10.1093/nar/gkad344.
- If GeneClusterMineX was useful provide credits for the Authors: Leandro de Mattos Pereira and Pedro Leão.


  ## More about - Authors:

- [Junior Researcher, Leandro de Mattos Pereira](https://mattoslmp.github.io)

- [CNP team, Dr. Pedro Leão, Researcher Leader](https://leaolab.wixsite.com/leaolab)

