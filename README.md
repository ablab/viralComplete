# viralComplete: BLAST-based viral completeness verification

viralComplete is intended for completeness verification of novel viral contigs. It heavily relies on following assumptions:
1) Virus genome size is consistent across the viral family.
2) If a newly constructed viral contig is complete and belongs to a known family of viruses then 
its gene content should be similar to the gene content of a known virus.

We thus compute the “similarity” of a given contig (based on the Naive Bayesian Classifier)
to each known virus from the RefSeq database, and check whether the most similar 
known virus have length similar to the contig length. 


### Requirements

viralComplete is a Python script, thus, installation is not required. However, it has the following dependencies:

* Python 3,
* BLAST,
* BioPython,
* Prodigal (https://github.com/hyattpd/Prodigal, available via conda).

To work properly, viralComplete require Prodigal in your PATH environment variable.


## Usage 

    ./viralcomplete.py 
            -f Input fasta file
            -o output_directory 

            Optional arguments:
            -h, --help  Show the help message and exit
            -t          Number of threads
            -thr        Completeness threshold (0.0-1.0), default = 0.9



Output file: comma-separated table <input_file>_result_table.csv

Output format: contig name, prediction result, log probability, most probable virus name
 
Fasta files with prediction results can be found in the _Prediction_results_fasta_ folder

viralComplete considers a virus complete if its length exceeds a specified threshold (% of the length of the closest known reference, 90% by default). Working with novel viruses, you may decrease the detection threshold (provided as an optional argument), to reduce number of potenital false negatives.
