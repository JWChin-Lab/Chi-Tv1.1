# Chi-T

 /ˈkaɪ ˈtiː/ ('kai-tee')

Pipeline creates tRNAs for use in genetic code expansion. tRNAs are designed to be active, orthogonal to the _E. coli_ machinery, and recognised by their corresponding synthetase.

## Requirements

RNAFold from the ViennaRNA package must be installed and added to the PATH. The ViennaRNA package can be installed from [here](https://www.tbi.univie.ac.at/RNA/#download).

Python 3 must be installed and added to the PATH. Package versions are compatible with Python 3.7 (and likely other Python 3.x which have not been tested).

Chi-T has been tested successfully on Windows 10, macOS, and Ubuntu 18.04. To run with the default tRNADB-CE dataset, we recommend a minimum of 8 GB RAM.

## Installation

Clone the GitHub repository using

```> git clone https://github.com/zyzzyva23/Chi-T.git```

Or download directly from [https://github.com/zyzzyva23/Chi-T](https://github.com/zyzzyva23/Chi-T)

Once installed, required Python packages can be installed with the command:

```> pip install -r requirements.txt```

You should note many of these package versions are now 'old', and so conflicts may arise if you attempt to mix the versions in the requirements with other newer versions (especially note the incompatibility with many versions of pandas with numpy 2.0.0 or later)

To test everything is working, navigate to the Chi-T directory and run the following command for a quick demo run on a reduced dataset:

```python main.py test_files/merged_test.csv test_files/Vir_pro_synth.xlsx Vir_pro Pro -o test_output/run_test -a CTA TGA CGA```

Total installation time should be on the order of minutes, and this demo should take a few seconds to run without errors, saving files of tRNAs at every filtering stage and outputting a log file. Due to stochasticity within Chi-T there is no defined expected output.

## Usage

A typical Chi-T run in Windows 10 on a 16-thread Ryzen 3700X CPU and 16 GB RAM took 1-2 hours. Usage without multiprocessing capabilities will take much longer.

### tRNA database
The tRNA database to be used can be downloaded from [here](http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi) in tab separated format. 

To format in a way Chi-T can use, use the script cleanup.py with the following command.

```> python cleanup.py <tRNA_file.tab> <alignment_file.xlsx> ```

A D-loop alignment file is supplied in the distribution as d_align.xlsx. More aligned D-loop sequences can be added manually, although Chi-T also automatically aligns D-loops not in the D-loop alignment file by anchoring to the XGGX motif.

### Synthetase File

To supply Chi-T with synthetases, you must supply an excel file in the format:

Synth Name | Synth ID | tRNA ID | Genome ID

* Synth name - Arbitrary Name for your own purposes

* Synth ID - A protein identifier e.g. a UNIPROT ID (although currently unused)

* tRNA ID - Identifier for the tRNA sequence corresponding to the tRNADB-CE entry

* Genome ID - A genome identifier for the organism (although currently unused)

Supplying synthetases to RS-ID requires this format also.

For an example of formatting, see the file test_files/Vir_pro_synth.xlsx

### Identitying synthetases to test using RS-ID

RS-ID can aid in synthetase selection by filtering for tRNAs with similar identity elements and identity parts. Supplying a synthetase, or several synthetases, restricts the search to sequences that pass similarity thresholds to the supplied synthetases' tRNAs. Running the script without supplying a synthetase returns a UMAP projection and clusters of all unique identity sequences in the isoacceptor class given.

To run, a usearch executable must be downloaded, and its file location given in the command line. This executable can be downloaded [here](https://www.drive5.com/usearch/download.html).

An example use case is shown below:

```> python synth_finder.py clean_trnas.csv Arg -o outputs/Pyr_112_output -sf my_synths.xlsx -sn Pyrococcus_synth_112 -d 0.2 -i 1```

#### Arguments

* positional arguments:
    * file
        * Clean tRNADB-CE file
    * amino_acid
        * Isoacceptor of choice

* Optional arguments:
    * -h, --help
        * Show help message and exit
    * -a, --anticodon
        * Change all anticodon sequences - has little effect on process
    * -o, --output_directory
        * Directory to save all files
    * -sf, --synth_file:
        * If specifying a synthetase, must be found in this file in excel format (formatting of file above)
    * -sn, --synth_name:
        * Name of user-defined synthetase - name must be found in synth_file
    * -u, --usearch:
        * File path of usearch executable
    * -d, --distance:
        * Upper threshold for USEARCH distance - default <= 0.2
    * -i, --identity:
        * Upper threshold for number of identity element nucleotide differences - default <= 1


### tRNA Designs

To use Chi-T call the script main.py with the corresponding arguments and options. Use python main.py -h to get a description of all possible inputs. 
An example use case is shown below:

```> python main.py clean_trnas.csv my_synths.xlsx Pyrococcus_synth_112 Arg -o outputs/Pyr_112_output -l -a CTA TGA CGA```

#### Arguments

* positional arguments:
    * file
        * Clean tRNADB-CE file
    * synth_file
        * File containing synthetase information (see above)
    * synth_name
        * Names of synthetases matching to entries in synth_file
        * Chi-T will iterate through each synthetase in order
        * Input is names of synthetases separated by a space e.g. Pro_04 Pro_05 Pro_08
    * amino_acid 
        * Amino Acid specified for tRNA generation (use three letter code e.g. Arg)

* optional arguments:
    * -h, --help: 
        * show help message and exit
    * -o, --output_directory: OUTPUT_DIRECTORY
        * Directory to store output files
    * -ip, --id_part_change: ID_PART_CHANGE
        * Identity parts that should be chimerified (except ID element)
        * Using this will select sequences containing all identity elements, but not necessarily the wild-type sequence
        * Input: names of parts separated by a space e.g. tRNA1-7_66-72* tRNA27-31_39-43*
    * -cp, --cluster_parts: CLUSTER_PARTS
        * Number of parts for each part type to cluster
        * Default is 200
        * i.e. 200 best scoring (lowest ID score) sequences from each part type to be clustered
        * Different isoacceptors have different numbers of sequences in the database, so this number requires tuning for different isoacceptor classes
    * -cm, --cluster_min: CLUSTER_MIN
        * Minimum number of remaining parts for clustering.
        * If the number of sequences in the part pool is lower than this number, clustering doesn't occur in order to keep numbers high enough.
        * Instead, 30 best scoring sequences are sampled randomly, and the most diverse sequences from this pool used for chimerification.
        * Separate from -cp and -cm arguments, if the number of sequences for a part is < 15, all sequences are used for chimera generation in every iteration.
    * -s, --subtle:          
        * If true, only parts <=2 mutations away from the native sequence are chosen (or a reference sequence if one is supplied)
    * -r, --reference: REFERENCE
        * If using subtle mode, this is the reference tRNA. Otherwise wild-type sequence used as reference.
    * -l, --length_filt: LENGTH_FILT
        * Filter chimeras if 79 nts or longer (or specified value).
        * If -l flag present without an integer following, then 79 nts used as a cutoff.
        * If integer supplied, this integer used as length cutoff.
    * -cf, --cervettini_filt: CERVETTINI_FILT
        * Parameters for Cervettini Filtering (employed in Cervettini et al., 2020) in order start_stringency, minimum_stringency, target number of chimeras, and step size
            * start_stringency: Initial ID score threshold to use. Default 0.5
            * minimum_stringency: Lowest ID score threshold to use before filtering gives up trying to reduce chimera number. Useful since ID score is not linear with orthogonality. Going too low can limit diversity by filtering chimeras that are alike (with similarly low scores). Default 0.2
            * target number of chimeras: As ID score threshold reduces, if the number of remaining chimeras goes below this number, then stop filtering. Default 2500000
            * step size: Size of ID score step to reduce each time. Default 0.05
        * Parameters must be in this order
            * To change parameters to 0.4, 0, 1000000, 0.1 then enter -cf 0.4 0 1000000 0.1
    * -a, --anticodons: ANTICODONS
        * Anticodons to iterate through in the order specified.
        * Can enter any number of anticodons >= 1.
        * First anticodon listed will be the one designs are finalised with e.g. CTA for TAG.
    * -f, --frequency: FREQUENCY
        * Average frequency from RNAFold across anticodons to use as threshold in final filtering step.
        * tRNAs with frequency lower than threshold removed.
        * Default 0.3.
    * -d, --diversity: DIVERSITY
        * Average diversity from RNAFold across anticodons to use as threshold in final filtering step.
        * tRNAs with diversity higher than threshold removed.
        * Default 5.0
    * -F, --final_frequency: FINAL_FREQUENCY
        * Average frequency threshold across anticodons
        * Default 0.35
    * -D, --final_diversity: FINAL_DIVERSITY
        * Average diversity threshold across anticodons
    * -n, --num_iterations: NUM_ITERATIONS
        * Number of times to iterate through Chi-T per synthetase.
        * After each iteration, used parts are removed from the initial part pool (except the wild-type sequence)
    * -m, --automatic: 
        * No user input required - without flag then user is required to approve multiple stages (mainly to prevent a silly amount of tRNAs crashing your computer).
    * -i, --initial:
        * If true, save initial chimeras to csv (only if needed for analysis as can be a large file).
    * -p, --pattern: PATTERN
        * Specify a csv file with synth name and regex string column.
        * Example file provided in tRNA_patterns.csv
    * -t, --num_tRNAs: NUM_TRNAS
        * Number of designs to output for each synthetase.
        * Default 4.

#### Example Usage with Arguments

```> python main.py tRNA_database.csv synthetase_file.xlsx metSynth1 metSynth2 Met –o output_folder –ip tRNA1-7_66-72* –cp 70 –cm 30 –cf 0.3 0 1500000 0.1 –-subtle –-reference reference_file.csv –-length_filt 78 –-anticodons CTA TGA CGA –-frequency 0.5 –-diversity 3 –n 3  –m –i –p pattern_file.csv –t 6```

### Creating Oligos

oligo_maker.py script included to produce csv file with designs.

```> python oligo_maker.py <tRNA design file> <output csv file> --forward_5 NNNNN --forward_3 NNNNN --reverse_5 NNNN --reverse_3 NNNN --anticodon CTA```

* tRNA design file: csv output file of main.py containing final designs
* output csv file: name of file oligos should be written to.
* --forward_5: sequence to add to 5' end of forward (sense) primer. Default GGCCGC
* --forward_3: sequence to add to 3' end of forward primer. Default CTGCA
* --reverse_5: sequence to add to 5' end of reverse primer. Default GATCTGCAG
* --reverse_3: sequence to add to 3' end of reverse primer. Default GC
* --anticodon: Designs to be written to output file will have this anticodon. Default CTA

### Checking a tRNA is in the database

tRNA_check.py script included to check if a tRNA ID is in the clean csv file of all tRNAs. Can be useful since many of the tRNA IDs specified in suppl. tables of Cervettini (2020) are not found in updated versions of tRNADB-CE.

```> python tRNA_check.py <tRNA DB file> <tRNA_id>```

* tRNA DB file: Clean tRNADB-CE csv file (processed through cleanup.py)
* tRNA_id: tRNADB-CE ID (or any) to check

If present, script will return the entry. Otherwise will return 'False'

### Adding a tRNA to the database file

tRNA_adder.py script used to add tRNAs to your cleaned tRNADB-CE csv file. Useful as tRNADB-CE does not currently contain eukaryotic tRNAs.

```> python tRNA_adder.py <tRNA_file> <alignment_file> --new_name updated_tRNAs.csv```

* tRNA_file: Clean tRNADB-CE csv file
* alignment_file: xlsx file containing D-loop alignments. Provided as d_align.xlsx
* -n, --new_name: Optional name of new csv file.

Script will ask for user input to build new entry manually.
