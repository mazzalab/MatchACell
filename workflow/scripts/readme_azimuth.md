# Modulo Azimuth (panhumanpy) nella pipeline di annotazione

## Scopo
`Azimuth.py` aggiunge Pan-human Azimuth (pacchetto `panhumanpy`) come uno dei tool della pipeline di annotazione. Interfaccia e output sono gli stessi di `CIA.py`, così i risultati dei vari tool si possono confrontare a valle.

Azimuth è un classificatore a rete neurale addestrato su un riferimento pan-tessuto umano. Per ogni cellula restituisce label gerarchiche a tre livelli:

| Colonna | Granularità | Esempio |
|---|---|---|
| `azimuth_broad` | molto grossolana | Immune cell |
| `azimuth_medium` | intermedia | T cell |
| `azimuth_fine` | fine | Treg cell |
| `final_level_labels` | il livello più profondo raggiunto con confidenza, diverso da cellula a cellula | — |

Il modello restituisce anche `final_level_confidence` e, tramite `map_to_cell_ontology`, i termini e gli ID di Cell Ontology (CL) per i tre livelli.

## Input richiesto
- Un file `.h5ad` con i **raw counts interi in `layers["counts"]`** (se il layer manca, si usa `.X`, purché contenga counts interi), su **tutti i geni** (niente subset HVG) e con i simboli genici come `var_names`.
- Le cellule devono essere già filtrate per QC.
- In `.X` può esserci qualunque cosa (di solito i dati normalizzati), così come UMAP e leiden: Azimuth non li usa. Lo script costruisce internamente un AnnData minimale dai counts e lascia intatto l'oggetto originale.
- Una UMAP in `obsm["X_umap"]`, usata per i plot.
- Un file verdict, da cui `extract_best_res` (in `functions_annot`) ricava la colonna leiden alla risoluzione scelta.

Errore da evitare: dare ad Azimuth un oggetto scalato e ridotto ai geni HVG, come `pbmc3k_processed`. Le predizioni non sono valide.

## Utilizzo
```bash
python Azimuth.py -i input.h5ad -v verdict.txt -o results/azimuth -t 70 \
    [-l azimuth_fine] [-c <soglia>] [-m v1]
```

| Argomento | Significato | Default |
|---|---|---|
| `-i` | h5ad di input | obbligatorio |
| `-v` | file verdict (risoluzione leiden); richiede `functions_annot` | — |
| `--cluster_col` | colonna `.obs` dei cluster, in alternativa a `-v` | — |
| `-o` | cartella di output | obbligatorio |
| `-t` | % minima di cellule di un cluster perché una label venga assegnata al cluster | obbligatorio |
| `-l` | livello usato come predizione: `azimuth_broad`, `azimuth_medium`, `azimuth_fine`, `final_level_labels` | `azimuth_fine` |
| `-c` | soglia facoltativa su `final_level_confidence`: le cellule sotto soglia diventano `Unassigned` | nessun filtro |
| `-m` | versione del modello Azimuth | default del pacchetto |

## Output (in `-o`)
- `azimuth_prediction_cell.csv`: per ogni cellula, tutte le colonne prodotte da Azimuth (livelli, confidenza, termini e ID CL) più la predizione finale `Azimuth prediction`.
- `umap_azimuth_clusters_cell.png`: la UMAP dell'oggetto di input colorata per predizione, con `Unassigned` in grigio.
- Gli output di `multi_umap`: confronto tra leiden e predizione.
- `cluster_to_celltype_mapping.xlsx`: crosstab leiden × predizione con la colonna `Azimuth_cluster_annotation`. Contiene le label che superano `-t`% nel cluster, separate da virgola, oppure `Unknown`.
- `Azimuth_annotated.h5ad`: l'oggetto di input con le colonne Azimuth aggiunte.

## Scelte fatte
- **Si usano solo le annotazioni.** Embedding (`azimuth_embed`) e UMAP di Azimuth sono passi facoltativi, di sola visualizzazione, e non vengono calcolati. I plot usano la UMAP già presente nell'oggetto, uguale per tutti i tool. La UMAP di Azimuth è calcolata sullo spazio latente della rete e quindi ha una disposizione diversa da quella scanpy.
- **`Unassigned` non dipende da una soglia nostra.** Raccoglie solo le cellule che Azimuth lascia senza label al livello scelto, cioè i valori `False`/NaN. Il filtro sulla confidenza è disattivato di default. La vignetta ufficiale usa `final_level_confidence > 0.5` come esempio di QC, ma non è una soglia validata. Un eventuale filtro va deciso nel passo di confronto, con lo stesso criterio per tutti i tool.
- **Nessuno smoothing per cluster.** Azimuth classifica le singole cellule. Il passaggio al cluster avviene solo nel file di mapping, con la stessa regola di soglia di CIA.

## Avvertenze
- **`final_level_labels` mescola livelli di granularità diversi** (per esempio "B cell" accanto a "Treg cell"), quindi non è adatto a confronti diretti con altri tool. Per confrontare i tool, `azimuth_medium` è di solito il livello più vicino ai cluster tipici. A livello fine compaiono label poco plausibili per il tessuto in esame, su poche cellule.
- **`final_level_confidence` si riferisce al livello finale raggiunto da ogni cellula.** Se usi `-c` insieme a `-l azimuth_medium`/`azimuth_fine`, la soglia non riguarda necessariamente il livello scelto.
- **Fissare `-m`** per avere risultati riproducibili: il modello di default cambia con le versioni principali del pacchetto.

## Ambiente
Ambiente conda `azimuth` (Linux, solo CPU). I pacchetti Python sono installati con pip.

| Pacchetto | Versione |
|---|---|
| python | 3.12.13 |
| panhumanpy | 1.0.0 (modello di default: `v1`) |
| tensorflow / keras | 2.17.0 / 3.15.1 |
| scanpy | 1.11.5 |
| anndata | 0.12.19 |
| numpy | 1.26.4 |
| pandas | 2.3.3 |
| scipy | 1.17.1 |
| scikit-learn | 1.6.0 |
| matplotlib | 3.11.2 |
| openpyxl | 3.1.5 (serve per il file `.xlsx`) |

.

Note:
- **Messaggi di TensorFlow all'avvio** (CUDA/cuDNN/TensorRT non trovati): sono solo informativi e indicano che il modello gira su CPU.
- **Pesi del modello**: vengono scaricati da Zenodo al primo utilizzo e salvati in `~/.cache/panhumanpy/`. Su nodi senza accesso a internet vanno scaricati prima.
- **Versione del modello**: per avere risultati riproducibili, passare esplicitamente `-m v1`, la versione di default di panhumanpy 1.0.0.

## Dataset di prova
Per i test si usa PBMC3k: counts da `sc.datasets.pbmc3k()` (in `layers["counts"]`), ristretti alle 2638 cellule di `sc.datasets.pbmc3k_processed()`, da cui si copiano `louvain` e la UMAP (`obsm["X_umap"]`). Senza file verdict e senza `functions_annot` lo script si lancia così:
```bash
python Azimuth.py -i data/test_input.h5ad --cluster_col louvain -o results/test -t 70 -m v1
```


## Riferimenti
- panhumanpy: https://satijalab.org/pan_human_azimuth/python/panhumanpy
