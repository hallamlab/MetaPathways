# What MetaPathways does

For the tool-by-tool diagrams and citations, see the [detailed workflow](detailed-workflow.md).

The [main workflow figure](index.md) and the diagrams below share the appnote's gray modules, blue inputs, green outputs and serif typography. Diamonds denote compute steps, not decisions; the schema diagram retains its relationship notation.

MetaPathways connects assembly annotations, read abundance and optional pathway inference in a shared set of sample and genome results. It accepts one sample or a collection and runs locally or through Slurm.

```mermaid
%%{init: {"theme":"base","fontFamily":"Times New Roman, Times, serif","themeVariables":{"fontFamily":"Times New Roman, Times, serif","fontSize":"16px","primaryColor":"#CCCCCC","primaryTextColor":"#111111","primaryBorderColor":"#666666","secondaryColor":"#DAE8FC","tertiaryColor":"#F5F5F5","lineColor":"#333333","edgeLabelBackground":"#FFFFFF","background":"#FFFFFF"},"flowchart":{"htmlLabels":false,"curve":"linear"}}}%%
flowchart LR
    A[Assemblies] --> B[Functional and taxonomic annotation]
    R[Reads] --> C[Read abundance]
    A --> C
    B --> D[Community and genome pathways]
    G[Genome assignments] --> D
    P[Licensed Pathway Tools] --> D
    B --> E[Reports and explorer]
    C --> E
    D --> E
    E --> F[Filtered tables and CSV exports]
    classDef module fill:#CCCCCC,stroke:#111111,stroke-width:1.5px,color:#111111;
    classDef compute fill:#F5F5F5,stroke:#666666,stroke-width:2px,color:#111111;
    classDef input fill:#DAE8FC,stroke:#6C8EBF,stroke-width:2px,color:#111111;
    classDef output fill:#D5E8D4,stroke:#82B366,stroke-width:2px,color:#111111;
    classDef data fill:#FFFFFF,stroke:#666666,stroke-width:1.5px,color:#111111;
    class A,R,G,P input;
    class B,C,D,E module;
    class F output;
```

Assemblies are required. Reads add abundance measurements; genome assignments let MP split community annotations into genome-specific inputs. MP does not assemble reads or perform genome binning. Reports remain available when optional inputs are omitted.

**For PGDBs, complete the [Pathway Tools installation guide](pathway-tools.md) first.** Use `--skip_ptools` to run without pathway inference. The [data-flow diagram](data-flow.md) shows these dependencies in more detail.

Start with [installation and the included three-sample test](installation.md), then follow the [complete workflow](analysis.md).
