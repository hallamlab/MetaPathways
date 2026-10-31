# MetaPathways containers

Run metagenomic annotation, read abundance, genome splitting and reports with the versioned `quay.io/hallamlab/metapathways:4.0.3` image. MAGSplitter, Camelot and the three-sample test dataset are included. Licensed Pathway Tools and production references are supplied separately.

Follow the [Docker and Apptainer guide](https://hallamlab-metapathways.readthedocs.io/en/latest/containers.html) for installation, the three-sample test, and persistent input/output storage. See the [GitHub quick start](https://github.com/hallamlab/MetaPathways#quick-start) and [report an issue](https://github.com/hallamlab/MetaPathways/issues).

For FAST searches on Lustre or similar shared storage, MetaPathways 4.0.3 adds `--fast_local_db` and `--func_search_time_limit`. See the [Slurm and scratch guide](https://hallamlab-metapathways.readthedocs.io/en/latest/resources.html#fast-databases-on-shared-storage).
