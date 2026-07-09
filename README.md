# ENCODE Atlas — BPNet-family model file index

## Tables

| file | family | models | what the models predict |
|---|---|--:|---|
| `bpnet_model_files.tsv` | BPNet | 2,339 | base-resolution TF binding (ChIP-seq/-nexus) |
| `chrombpnet_model_files.tsv` | ChromBPNet | 1,512 | base-resolution chromatin accessibility (ATAC/DNase) |
| `procapnet_model_files.tsv` | ProCapNet | 6 | transcription initiation (PRO-cap) |
| `reporternet_model_files.tsv` | ReporterNet | 8 | massively-parallel reporter activity |

Each row is a trained model (ENCODE experiment accession), each column an output object type, and each cell the ENCODE file accession (`ENCFF…`) to download for that pair.
