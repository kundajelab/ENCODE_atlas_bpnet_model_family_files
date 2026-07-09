# ENCODE Atlas — BPNet-family model file index

## Tables

| file | family | models | what the models predict |
|---|---|--:|---|
| `bpnet_model_files.tsv` | BPNet | 2,339 | base-resolution TF binding (ChIP-seq/-nexus) |
| `chrombpnet_model_files.tsv` | ChromBPNet | 1,512 | base-resolution chromatin accessibility (ATAC/DNase) |
| `procapnet_model_files.tsv` | ProCapNet | 6 | transcription initiation (PRO-cap) |
| `reporternet_model_files.tsv` | ReporterNet | 8 | massively-parallel reporter activity |

Each row is a trained model (ENCODE experiment accession), each column an output object type, and each cell the ENCODE file accession (`ENCFF…`) to download for that pair.

## Client — `encode_atlas.py`

A single-file, dependency-free (stdlib-only) Python client any agent or user can use to resolve and download the models and their downstream products, and to reach the broader ENCODE portal. It fetches the tables from this repo automatically (cached), so you only need the one file.

```python
import encode_atlas as ea

# what's released for a model, and where to download each object
ea.model_files("ENCSR032RGS")                    # family + metadata + all object types -> ENCFF + URL
ea.download("ENCSR032RGS", "model", "./models/") # the trained model tar; also: motifs, contrib_counts,
                                                 #   contrib_profile, signal_predicted, motif_hits, ...

# search OUR model index (the four family tables)
ea.search(family="chrombpnet", biosample="K562") # -> [{accession, family, assay, target, tissue, qc}]

# reach the BROADER ENCODE portal (raw + processed data for ANY experiment, not just modelled ones)
ea.portal_files("ENCSR000EOT")                   # every file: reads, alignments, peaks, signal, ...
ea.portal_search(assay_title="ATAC-seq", biosample="liver")
```

Object types: `model`, `motifs`, `motif_hits`, `contrib_counts`, `contrib_profile`,
`signal_predicted`, `signal_observed`, `signal_biascorrected`, `regions`, `regions_train_test`,
`metrics` (see `ea.object_types()`). Download a file directly at
`https://www.encodeproject.org/files/<ENCFF>/@@download/`.
