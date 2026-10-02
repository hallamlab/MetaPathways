# Data flow: inputs to connected results

The arrows show which data products feed each calculation. They do not require independent branches to run serially. MP preserves sample, contig, feature and entity identifiers so the final tables can be joined and explored.

```mermaid
flowchart TB
    A[Assembly FASTA] --> QC[Preprocessed contigs]
    QC --> F[Predicted CDS and RNA features]
    F --> AN[Functional and taxonomic annotations]
    DB[MPDB reference sequences and mappings] --> AN
    READS[Optional reads] --> MAP[Read mapping and feature counting]
    QC --> MAP
    F --> MAP
    MAP --> AB[Contig and feature abundance]
    AN --> PI[Community Pathway Tools inputs]
    PI --> SPLIT[Genome-specific inputs]
    GM[Optional contig-to-genome map] --> SPLIT
    QC --> SEQ[Sequence-backed PGDB staging]
    PI --> SEQ
    SPLIT --> SEQ
    SEQ --> PT[Optional licensed Pathway Tools inference]
    PT --> PW[Pathway, reaction and gene tables]
    AN --> REPORT[Report tables and relational index]
    AB --> REPORT
    GM --> REPORT
    PW --> REPORT
    REPORT --> PORTAL[Explorer searches, subsets and CSV exports]
```

- Without reads, MP does not compute read abundance.
- Without a genome map, community analysis can still proceed; genome splitting is omitted.
- With `--skip_ptools`, annotation, abundance and reports remain available, but PGDB inference is omitted.
- PGDB staging joins annotations back to actual contig sequences and feature coordinates. Both community and genome PGDBs need sequence-backed inputs.

Abundance is calculated from the reads and features, independently of pathway inference. A pathway–gene association can be joined to abundance through the corresponding feature IDs; these are different measurements, not interchangeable values. Each annotation's taxonomy belongs to that annotation's reference database. See the [result schema](results-schema.md) for table keys and safe joins.

Follow the [complete workflow](analysis.md) for commands, or the [stage reference](workflow.md) for individual products and dependencies.
