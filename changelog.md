# Changelog

## 3.5.2 — feature candidate

- Nextflow controllers for local and Slurm workflows, multi-sample manifests, resource budgets, database builds and licensed Pathway Tools container preparation.
- Read-mapping/checkpoint fixes and sequence-backed PGDB staging.
- SwissProt, UniRef and eggNOG taxonomy scoped to each annotation database and target; independent within-database LCA.
- Searchable reports with a sample-level landing page, wide abundance tables, source provenance, CSV exports and GitHub feedback links.
- Bundled three-sample CAMI reviewer inputs and documented local validation.

### Updates

**November 27, 2014**: [MetaPathways v2.5 released](https://github.com/hallamlab/metapathways2/releases/tag/v2.5) with upgrades to the pipeline:
    

* LAST homology searches with BLAST-equivalent output and E-values
* Reads per kilobase per million mapped (RPKM) coverage measure for Contig annotations calculated from raw reads (`.fastq`) or mapping files (`.SAM`) using [bwa](http://bio-bwa.sourceforge.net)
* Addition of the [CAZy sequence database](http://www.cazy.org) as a new compatible functional hierachy
* GUI Keyword-search from annotation subsetting and projection onto different functional hierarcies (KEGG, COG, SEED, MetaCyc, and now CAZy)

See [the release page](https://github.com/hallamlab/metapathways2/releases/tag/v2.5) and [the wiki](https://github.com/hallamlab/metapathways2/wiki) for more information.
