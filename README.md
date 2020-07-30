# blast_filter

It is a simple python (python3) script devised to filter and group sequences according to their similarity with reference sequences. The filtering is based on BLAST and LAST results. It also creates a tab delimited file summarizing the results.

# Usage

`python blast_filter.py <identity> <aln_len> <file.blast> <file.fasta> <file.db>`

where:
    * `identity` is the minimum identity score (above which a sequence will be filtered).
    * `aln_len` is the minimum length of the aligned segment (in base pairs).
    * `file.blast` is a file in standard blast format ('Hit Table (text)' on NCBI).
    * `file.fasta` is a fasta file containing all the query sequences.
    * `file.db` is a fasta file with standardized description lines (please, refer to the example file) containing the reference sequences.
