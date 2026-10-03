# Data flow: inputs to connected results

For the tool-by-tool diagrams and citations, see the [detailed workflow](detailed-workflow.md).

The arrows show which data products feed each calculation. They do not require independent branches to run serially. MP preserves sample, contig, feature and entity identifiers so the final tables can be joined and explored.

```mermaid
%%{init: {"theme":"base","fontFamily":"Times New Roman, Times, serif","themeVariables":{"fontFamily":"Times New Roman, Times, serif","fontSize":"16px","primaryColor":"#CCCCCC","primaryTextColor":"#111111","primaryBorderColor":"#666666","secondaryColor":"#DAE8FC","tertiaryColor":"#F5F5F5","lineColor":"#808080","edgeLabelBackground":"#FFFFFF","background":"#FFFFFF","defaultLinkColor":"#808080","textColor":"#111111"},"flowchart":{"htmlLabels":false,"curve":"linear"}}}%%
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
    classDef module fill:#CCCCCC,stroke:#111111,stroke-width:1.5px,color:#111111;
    classDef compute fill:#F5F5F5,stroke:#666666,stroke-width:2px,color:#111111;
    classDef input fill:#DAE8FC,stroke:#6C8EBF,stroke-width:2px,color:#111111;
    classDef output fill:#D5E8D4,stroke:#82B366,stroke-width:2px,color:#111111;
    classDef data fill:#FFFFFF,stroke:#666666,stroke-width:1.5px,color:#111111;
    class A,DB,READS,GM input;
    class QC,F,AN,PI,AB,PW data;
    class MAP,SPLIT,SEQ,PT compute;
    class REPORT module;
    class PORTAL output;
```

- Without reads, MP does not compute read abundance.
- Without a genome map, community analysis can still proceed; genome splitting is omitted.
- With `--skip_ptools`, annotation, abundance and reports remain available, but PGDB inference is omitted.
- PGDB staging joins annotations back to actual contig sequences and feature coordinates. Both community and genome PGDBs need sequence-backed inputs.

Abundance is calculated from the reads and features, independently of pathway inference. A pathway–gene association can be joined to abundance through the corresponding feature IDs; these are different measurements, not interchangeable values. Each annotation's taxonomy belongs to that annotation's reference database. See the [result schema](results-schema.md) for table keys and safe joins.

Follow the [complete workflow](analysis.md) for commands, or the [stage reference](workflow.md) for individual products and dependencies.
