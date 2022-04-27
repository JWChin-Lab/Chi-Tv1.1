# Chi-T

 /ˈkaɪtiː/ ('kai-tee')

Pipeline creates tRNAs for use in genetic code expansion. tRNAs are designed to be active, orthogonal to the E. coli machinery, and recognised by their corresponding synthetase.

## Requirements

RNAFold from the ViennaRNA package must be installed and added to the PATH. The ViennaRNA package can be installed from [here](https://www.tbi.univie.ac.at/RNA/#download).

Python 3 must be installed and added to the PATH.

## Installation

Clone the GitHub repository using

```git clone https://github.com/zyzzyva23/Chi-T.git```

Or download directly from [https://github.com/zyzzyva23/Chi-T](https://github.com/zyzzyva23/Chi-T)

Once installed, required Python packages can be installed with the command:

```pip install -r requirements.txt```

## Usage
### tRNA database
The tRNA database to be used can be downloaded from [here](http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi) in tab separated format. 

To format in a way Chi-T can use, use the script cleanup.py with the following command.

```> python cleanup.py <tRNA_file.tab> <alignment_file.xlsx> ```

A D-loop alignment file is supplied in the distribution as d_align.xlsx.

### Synthetase File

To supply Chi-T with synthetases, you must supply an excel file in the format:

Synth Name | Synth ID | tRNA ID | Genome ID

* Synth name - Arbitrary Name for your own purposes

*  Synth ID - a protein identifier e.g. a UNIPROT ID although currently not need

### tRNA Designs

To use Chi-T call the script main.py with the corresponding arguments and options. Use python main.py -h to get a description of all possible inputs. 
An example use case is shown below

```> python main.py clean_trnas.csv my_synths.xlsx Pyrococcus_synth_112 Arg -o outputs/Pyr_112_output -l -a CTA TGA CGA -s```
