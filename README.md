# ReAlign-P: A vertical iterative realignment method for protein multiple sequence alignment

ReAlign-P is a tool written in C++17 for realigning the protein multiple sequence alignment. It runs on Linux.

## 🔨Installation and Usage

### 1.1 Linux/WSL(Windows Subsystem for Linux ) - from Anaconda
1.Install WSL for Windows. Instructional video [1](https://www.youtube.com/watch?v=X-DHaQLrBi8&t=5s) or [2](http://lab.malab.cn/%7Etfr/1.mp4) (Copyright belongs to the original work).

2.Download and install Anaconda. Download Anaconda for different systems [here](https://www.anaconda.com/products/distribution#Downloads). Instructional video of anaconda installation [1](https://www.youtube.com/watch?v=AshsPB3KT-E) or [2](http://lab.malab.cn/%7Etfr/Install_anaconda_in_Linux.mp4) (Copyright belongs to the original work).

3.Install ReAlign-P.
```bash
#1 Create and activate a conda environment for ReAlign-P
conda create -n realign_p_env
conda activate realign_p_env

#2 Add channels to conda
conda config --add channels malab

#3 Install ReAlign-P
conda install -c malab realign_p

#4 Test ReAlign-P
realign_p -h
```

### 1.2 Linux/WSL(Windows Subsystem for Linux ) - from the source code

1.Download and Compile the source code. (Make sure your version of gcc >= 9.4.0)
```shell
#1 Download
git clone https://github.com/malabz/ReAlign-P.git

#2 Open the folder
cd ReAlign-P

#3 Compile
make

#4 Test ReAlign-P
./realign_p -h
```

2.Install the required alignment tools mafft and muscle5, we recommend using Conda to install.
```shell
# 1 Create and activate a conda environment for msa
conda create -n msa
conda activate msa
# 2 Install mafft and muscle
conda install -c bioconda mafft muscle=5.2
```

### 2 Usage
```
Usage: ./realign_p -i <input_file> [-o <output_file>] [-s <score_matrix>] [-l <length>] [-m <msa>]

Options:
  -i <input_file>    (required) Path to the input file containing sequence data.
  -o <output_file>   (optional) Path to the output file for storing results. Default is 'realign_p_result.fasta'.
  -s <score_matrix>   (optional) Score matrix to use, options are 'blosum62' or 'pam250'. Default is blosum62.
  -m <msa>           (optional) MSA tool to use, options are 'mafft' or 'muscle5'. Default is 'muscle5'.

Examples:
  ./realign_p -i data.fasta -o results.fasta -s pam250 -m mafft
  ./realign_p -i data.fasta -m muscle5

Note:
  - The '-i' option is required.
  - The '-m' option only supports 'mafft', and 'muscle5'.
```

## 🔬Test dataset and the use case
### 1. Information about the test dataset
Dataset|Nummber|Avg Number|Avg Length
:---:|:---:|:---:|:---:
BAliBASE|386|28.71|338.28
OXBench|395|8.33|138.58
PREFAB4|1682|45.19|233.51
SABRE|423|5.72|171.22

<!-- ### 2. The use case
```shell
# Download data
wget http://lab.malab.cn/soft/ReAlign-N/data/16s_like.tar.gz

# Unzip data
tar -zxvf 16s_like.tar.gz

# Get the folder path
cd 16s_like

# Run ReAlign-P

``` -->

## 📍Reminder
1. Currently ReAlign-P is **ONLY** available for Protein. 
2. Please ensure that the sequence ID entered into ReAlign is unique.
3. MAFFT and MUSCLE5 installation are required for the utilization of ReAlign-P. 

## 🖥️Environment
System|GCC version
:---:|:---:
Linux|GCC 9.4.0
WSL|GCC 9.4.0

## 🙏References & Acknowledgements
We would like to acknowledge the following msa tools that contributed to the development of ReAlign-P:

- **[MAFFT](https://mafft.cbrc.jp/alignment/server/index.html)**: This is a widely used multiple sequence alignment tool known for its high accuracy and scalability.

- **[MUSCLE5](https://www.drive5.com/muscle/)**: This is a novel algorithm which constructs an ensemble of high-accuracy alignment with diverse biases by perturbing a hidden Markov model and permuting its guide tree. 

## 🔖Citation

## 👋Contacts
The software tools are developed and maintained by 🧑‍🏫[ZOU's lab](http://lab.malab.cn/~zq/en/index.html).

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/ReAlign-P/issues) or email us at 👉[📩](zhai1xiao@gmail.com).

More tools and infomation can visit our [github](https://github.com/malabz).
