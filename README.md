# Chi-T

Pipeline creates tRNAs for use in genetic code expansion. tRNAs are designed to be active, orthogonal to the E. coli machinery, and recognised by their corresponding synthetase.

## Requirements

RNAFold from the ViennaRNA package must be installed and added to the PATH. The ViennaRNA package can be installed from 'here'. \n
Python 3 must be installed and added to the PATH. \n
Python packages are found in requirements.txt \n

## Installation

git clone .....
requirements.txt

## Usage
### tRNA database
The tRNA database to be used can be downloaded from ... in tab separated format. To format in a way Chi-T can use, use the script cleanup.py with the following command.
An alignment file is supplied in the distribution as d_align.xlsx.
" python cleanup.py <tRNA_file.tab> <alignment_file.xlsx> "

To use Chi-T call the script main.py with the corresponding arguments and options. Use python main.py -h to get a description of all possible inputs. 
An example use case is shown below

python main.py clean_trnas.csv my_synths.csv Pyrococcus_synth_112 Arg -o outputs/Pyr_112_output -l -a CTA TGA CGA -s
