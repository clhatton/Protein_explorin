import pandas as pd
from snakemake import shell

# Read the CSV file with UniProt IDs
uniprot_ids = pd.read_csv("uniprot_ids_2.csv", header=None)[0].tolist()

rule all:
    input:
        expand("output/{uniprot_id}_docked.pdb", uniprot_id=uniprot_ids)

rule download_alphafold_pdb:
    output:
        "pdb/{uniprot_id}.pdb"
    params:
        uniprot_id="{uniprot_id}"
    shell:
        """
        curl -L -o {output} https://alphafold.ebi.ac.uk/files/AF-{params.uniprot_id}-F1-model_v4.pdb
        """

rule run_rosetta_docking:
    input:
        "c_term_tail.pdb" ,
        "pdb/{uniprot_id}.pdb"
    output:
        "output/{uniprot_id}_docked.pdb"
    params:
        uniprot_id="{uniprot_id}"
    shell:
        """
        /zata/zippy/hattonc/rosetta/source/bin/docking_protocol.default.linuxgccrelease -s {input[0]} -s {input[1]} \
        -nstruct 10 \
        -out:path:all output/ \
        -ignore_unrecognized_res \
        -missing_density_to_jump
        """
rule cleanup_pdb:
    input:
        "output/{uniprot_id}_docked.pdb"
    shell:
        "rm pdb/{wildcards.uniprot_id}.pdb"
