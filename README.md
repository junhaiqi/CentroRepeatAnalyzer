
## Overview of CentroRepeatAnalyzer
CentroRepeatAnalyzer is a pipeline based on [centroAnno](https://github.com/junhaiqi/centroAnno.git), which can directly analyze repeats and higher-order tandem repeat (HOR) structures of centromere regions from noisy sequencing data.
## Table of contents

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [Output](#output)
  * [Acknowledgments](#acknowledgments)
  * [License](#license)
  * [Cite](#cite)


## Requirements
CentroRepeatAnalyzer runs Linux and only requires a few python packages, such as numpy and pyfastx.


## Installation

```bash
git clone https://github.com/junhaiqi/CentroRepeatAnalyzer.git
pip install numpy pyfastx tqdm
```
We provide some required executable files, they are in "bin". If they cannot be used, users can download the required software packages ( [centroAnno](https://github.com/junhaiqi/centroAnno.git), [cd-hit](https://github.com/weizhongli/cdhit.git) ) and compile them to get the corresponding executable files and put them in the "bin".



## Usage

Basic command:
```bash
cd CentroRepeatAnalyzer
python main.py --h

usage: main.py [-h] [-k KMER_SIZE] [-l LEN_CUTOFF] [-c COUNT_CUTOFF] [-r REP_REDIO_CUTOFF] [-p] [-s CLUSTER_SIZE_CUTOFF] [-d CLUS_IDENTITY] [-e REP_IDENTITY] [-m MONOMER_TEMPLATE_FILE] [-y] [-t THREAD_NUM]
               fastX_file output_folder_path

positional arguments:
  fastX_file            It is the input filename in fasta format, required, can be in .gz format.
  output_folder_path    It is the path to the output folder, required.

optional arguments:
  -h, --help            show this help message and exit
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        It is the size of k-mer, which is used to filter repeat sequences and infer repeats (default = 9).
  -l LEN_CUTOFF, --len-cutoff LEN_CUTOFF
                        It is the length cutoff, reads longer than this value are considered repeat sequences (default = 5000).
  -c COUNT_CUTOFF, --count-cutoff COUNT_CUTOFF
                        It is the count cutoff, repeats that appear more than this value are considered high-quality repeats (default = 10).
  -r REP_REDIO_CUTOFF, --rep-redio-cutoff REP_REDIO_CUTOFF
                        It is the redio cutoff of non-unique k-mer, sequences greater than this value are considered potential repeat sequences (default = 0.5).
  -p, --close-hpc       It is the bool flag, when it appears, homopolymer compression technology will not be applied (default = False).
  -s CLUSTER_SIZE_CUTOFF, --cluster-size-cutoff CLUSTER_SIZE_CUTOFF
                        It is the size cutoff, the clusters which size larger this values are used to inferred repeats (default = 5).
  -d CLUS_IDENTITY, --clus-identity CLUS_IDENTITY
                        It is the sequence identity threshold for cdhit (default = 0.9).
  -e REP_IDENTITY, --rep-identity REP_IDENTITY
                        It is the identity threshold to filter reads with high-quality annotation results (default = 0.85).
  -m MONOMER_TEMPLATE_FILE, --monomer-template-file MONOMER_TEMPLATE_FILE
                        It is the fasta file which includes monomer templates (default = None), If it is given, the program will scan out all reads potentially related to these monomers and complete more
                        accurate annotations based on these monomers.
  -y, --only-infer      It is the bool flag, when it appears, it means that there is already annotation results, and the annotation results are stored in $output_folder_path (default = False).
  -t THREAD_NUM, --thread-num THREAD_NUM
                        It is the number of threads used to clustering (default = 4).
```

## Example
Analyze the tendem repeats and HORs on sequencing data without template information:
```bash
python main.py testData/subset_SRR11292120.fa test -p -k 9 -s 3
```

Analyze the tandem repeats and HORs on sequencing data with monomer template information for centromeric reads:

```bash
python main.py testData/subset_SRR11292120.fa test_mono -p -k 9 -s 3 -m testData/AlphaSat.fa
```

## Output
If the above example runs successfully, the folder 'test' and 'test_mono' will be the following files:

| File   | Description |
   |  :----:  | :----:  |
   | *_decomposedResult.csv  | Load the decompose results at monomer level |
   | *_horDecomposedResult.csv  | Load the decompose results at HOR level |
   | *_monomerTemplates.fa  | Load all monomers inferred by centroAnno |
   | *_HORs.fa  | Load all HORs inferred by centroAnno |
   | potential_repeat_reads.fa  | Load all potential reads which include repeats |
   | potential_repeats.fa  | Load all potential repeats|
   | rep_clusters  | Load all representative sequences produced by cd-hit clustering for "potential_repeats.fa" |
   | rep_clusters.clstr  | Load results produced by cd-hit clustering for "potential_repeats.fa" |
   | representative_repeats.fa  | Load all representative repeats |

'test_mono' also has some unique output:

| File   | Description |
   |  :----:  | :----:  |
   | monomers-related_reads.fa  | Load all reads which include some monomers from "testData/AlphaSat.fa" |
   | monomers-related_refine_annotations  | Loads all more detailed annotation results for reads related to the monomers from "testData/AlphaSat.fa"  |

## Acknowledgments
None.

## License 
MIT License.

## Cite
None.
