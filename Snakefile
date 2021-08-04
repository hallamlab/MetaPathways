## -*- python -*-

#### Snakemake file for MetaPathways

## This is currently used for orchestrating the staging of reference databases used by end-users.
## In the future this will be used to integrate new functionality into MetaPathways.

### Local Definitions:

## User-configurations are set in the config YAML file, not in this Snakemake file:
configfile: "snakemake_config.yaml"

## Database versions:
arb_release = "138.1"

## Mapping files to URL prefixes:
file2url_prefix = {}
file2url_prefix['SILVA_' + arb_release + '_LSUParc_tax_silva.fasta.gz'] = 'https://www.arb-silva.de/fileadmin/silva_databases/release_' + arb_release + '/Exports'
file2url_prefix['SILVA_' + arb_release + '_SSUParc_tax_silva.fasta.gz'] = 'https://www.arb-silva.de/fileadmin/silva_databases/release_' + arb_release + '/Exports'

hosted_functional_dbs = ["cazy",
                         "kegg",
                         "metacyc"]

thrird_party_functional_dbs = ["refseq"]


### DB Staging Rules:

rule stage_blast_full:
    input:
        expand(config["ref_db_dir"] + '/taxonomic/formatted/SILVA_{arb_release}_SSURef_tax_silva.ndb', arb_release=arb_release),
        expand(config["ref_db_dir"] + '/taxonomic/formatted/SILVA_{arb_release}_LSURef_tax_silva.ndb', arb_release=arb_release),
        config["ref_db_dir"] + '/functional/formatted/uniprot_sprot.pdb',
        config["ref_db_dir"] + '/functional/formatted/kegg-uniprot-2018-12-20.pdb',
        #config["ref_db_dir"] + '/functional/formatted/cazy-2020-06-01.pdb',
        config["ref_db_dir"] + '/functional/formatted/metacyc-2020-08-10.pdb',
        config["ref_db_dir"] + '/functional/formatted/refseq_protein.26.psq'
        
        
#rule stage_blast_lite:
#    input:

#rule stage_fast_full:
#    input:

#rule stage_fast_lite:
#    input:


rule download_only:
    input:
        config["ref_db_dir"] + '/taxonomic/SILVA_' + arb_release + '_SSURef_tax_silva.fasta.gz',
        config["ref_db_dir"] + '/functional/uniprot_sprot.fasta.gz',
        config["ref_db_dir"] + '/functional/kegg-uniprot-2018-12-20',
        config["ref_db_dir"] + '/functional/cazy-2020-06-01',
        config["ref_db_dir"] + '/functional/metacyc-2020-08-10',
        config["ref_db_dir"] + '/functional/formatted/refseq_protein.26.psq'

rule create_dirs_local_files:
    params:
        target_dir = config["ref_db_dir"]
    output:
        config["ref_db_dir"] + '/Dsignal'
    shell:
        """
        mkdir -p {params.target_dir}/functional/formatted
        mkdir -p {params.target_dir}/taxonomic/formatted
        cp tests/data/ref_data/Dsignal {params.target_dir}/
        cp tests/data/ref_data/TPCsignal {params.target_dir}/
        ## Stage Functional Categories Files:
        cp -R tests/data/ref_data/functional_categories {params.target_dir}/
        cp -R tests/data/ref_data/ncbi_tree {params.target_dir}/
        """
        
rule fetch_silva_db:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config["ref_db_dir"],
        arb_release = arb_release
    output:
        ssu_silva_file = expand(config["ref_db_dir"] + '/taxonomic/SILVA_{arb_release}_SSURef_tax_silva.fasta', arb_release=arb_release),
        lsu_silva_file = expand(config["ref_db_dir"] + '/taxonomic/SILVA_{arb_release}_LSURef_tax_silva.fasta', arb_release=arb_release)
    shell:
        """
        cd {params.target_dir}/taxonomic
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_{params.arb_release}/Exports/SILVA_{params.arb_release}_LSURef_tax_silva.fasta.gz
        wget https://www.arb-silva.de/fileadmin/silva_databases/release_{params.arb_release}/Exports/SILVA_{params.arb_release}_SSURef_tax_silva.fasta.gz
        gunzip {output.ssu_silva_file}
        gunzip {output.lsu_silva_file}
        """

rule make_silva_blast_db:
    input:
        ssu_silva_file = expand(config["ref_db_dir"] + '/taxonomic/SILVA_{arb_release}_SSURef_tax_silva.fasta', arb_release=arb_release),
        lsu_silva_file = expand(config["ref_db_dir"] + '/taxonomic/SILVA_{arb_release}_LSURef_tax_silva.fasta', arb_release=arb_release)
    params:
        target_dir = config["ref_db_dir"],
        arb_release = arb_release
    output:
        expand(config["ref_db_dir"] + '/taxonomic/formatted/SILVA_{arb_release}_SSURef_tax_silva.ndb', arb_release=arb_release),
        expand(config["ref_db_dir"] + '/taxonomic/formatted/SILVA_{arb_release}_LSURef_tax_silva.ndb', arb_release=arb_release)
    shell:
        """
        cd {params.target_dir}/taxonomic
        makeblastdb -in SILVA_{arb_release}_LSURef_tax_silva.fasta -dbtype nucl -parse_seqids -out formatted/SILVA_{arb_release}_LSURef_tax_silva
        makeblastdb -in SILVA_{arb_release}_SSURef_tax_silva.fasta -dbtype nucl -parse_seqids -out formatted/SILVA_{arb_release}_SSURef_tax_silva
        """
    
rule fetch_uniprot_swissprot_db:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/uniprot_sprot.fasta'
    shell:
        """
        cd {params.target_dir}/functional
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip uniprot_sprot.fasta.gz
        """

rule make_swissprot_blast_db:
    input:
        config["ref_db_dir"] + '/functional/uniprot_sprot.fasta'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/formatted/uniprot_sprot.pdb'
    shell:
        """
        cd {params.target_dir}/functional
        makeblastdb -in {input} -dbtype prot -parse_seqids -out formatted/uniprot_sprot
        """

        
rule fetch_kegg_uniprot_db:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/kegg-uniprot-2018-12-20'
    shell:
        """
        cd {params.target_dir}/functional
        wget https://ndownloader.figshare.com/files/27422531 -O kegg-uniprot-2018-12-20.gz
        gunzip kegg-uniprot-2018-12-20.gz
        """

rule make_kegg_blast_db:
    input:
        config["ref_db_dir"] + '/functional/kegg-uniprot-2018-12-20'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/formatted/kegg-uniprot-2018-12-20.pdb'
    shell:
        """
        cd {params.target_dir}/functional
        makeblastdb -in kegg-uniprot-2018-12-20 -dbtype prot -parse_seqids -out formatted/kegg-uniprot-2018-12-20
        """


        
rule fetch_cazy_db:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/cazy-2020-06-01.fasta'
    shell:
        """
        cd {params.target_dir}/functional
        wget https://ndownloader.figshare.com/files/27421229 -O cazy-2020-06-01.fasta.gz
        gunzip cazy-2020-06-01.fasta.gz
        """

rule make_cazy_blast_db:
    input:
        config["ref_db_dir"] + '/functional/cazy-2020-06-01.fasta'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/formatted/cazy-2020-06-01.pdb'
    shell:
        """
        cd {params.target_dir}/functional
        makeblastdb -in {input} -dbtype prot -parse_seqids -out formatted/cazy-2020-06-01
        """


        
rule fetch_metacyc_db:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/metacyc-2020-08-10.fasta'
    shell:
        """
        cd {params.target_dir}/functional
        wget https://ndownloader.figshare.com/files/27419069 -O metacyc-2020-08-10.fasta
        """

rule make_metacyc_blast_db:
    input:
        config["ref_db_dir"] + '/functional/metacyc-2020-08-10.fasta'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/formatted/metacyc-2020-08-10.pdb'
    shell:
        """
        cd {params.target_dir}/functional
        makeblastdb -in {input} -dbtype prot -parse_seqids -out formatted/metacyc-2020-08-10
        """


rule fetch_refseq_via_update_blastdb:
    input:
        config["ref_db_dir"] + '/Dsignal'
    params:
        target_dir = config['ref_db_dir']
    output:
        config["ref_db_dir"] + '/functional/formatted/refseq_protein.26.psq'
    shell:
        """
        cd {params.target_dir}/functional/formatted
        update_blastdb.pl --blastdb_version 5 --decompress refseq_protein
        """
