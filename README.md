# GeneClusterMineX v 1.0.0 :rocket:
_____________________________________________________________________________________________________________________________________________________

 - **v 1.0.0: automation of processing of several Genomes fasta/fna files by antismash**


This Python script automates the processing of **`.fna`** files using antiSMASH (**obs: if you file have sufix .fasta, change the python code**). 
First the code will create the **output_antismash directory** and later the program performed the antismash and store the Results of each run of Genome files (*.fna ou.fasta*) in subdirectory called subdirectory **`Result_*`name_genome**. 
But before each processing, the code will checks whether each **FASTA Genome file** has been previously processed based on the existence of the **`Result_*`** directory inside of directory: **output_antismash** and the `index.html` file inside of each subdirectory `Result_*`. If a file has not been processed or has been processed partially, the script will runs antiSMASH again to generate the results. **This .html file check was created because sometimes when we are running on a HPC server, the process can be interrupted for unknown reasons, it is a way to ensure that all genomes have been analyzed.**

## Prerequisites:

- Ubuntu operating system 20.04
- Python 3.11
- Conda (Anaconda or Miniconda)

## Requirements:

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

4. Edit the `antismash_run_fna.py` file and adjust the input directory of the `.fna` files (variable `input_directory`) according to where your files are located.

5. Run the Python script.

## Results:

- Results processed by antiSMASH will be stored in directories prefixed with `Result_` in the specified output directory followed by the name of the fasta genome file.
- The script will print messages indicating which files have already been processed previously and whether the `index.html` file exists.
- Previous result directories that do not have the `index.html` file will be deleted and processed again.

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

