#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MatchACell | CELLama adapter (CLI).

SCOPO
  Legge un dataset non annotato (.h5ad o matrice 10x .h5) e una reference
  AnnData con etichette note; utilizza le funzioni ORIGINALI di CELLama per
  gli embedding e per trasferire le etichette (nearest-neighbor cosine).
  Salva un .h5ad annotato, una UMAP PNG, un CSV e metadati JSON.

PRIMO AVVIO (dalla shell, senza creare manualmente un ambiente):
  python run_cellama.py --input query.h5ad --reference ref.h5ad \\
      --reference-label cell_type --output risultati/

COME FUNZIONA L'AMBIENTE
  1. Con Python standard >=3.9, questo script crea automaticamente un venv
     isolato in ~/.cache/matchacell/cellama-venv (oppure --venv-dir PATH).
  2. Si riavvia dentro il venv, installa CELLama originale dal GitHub e le
     eventuali dipendenze mancanti; ai run successivi riusa l'ambiente.
  3. Non esegue 'conda activate' (impossibile modificare l'ambiente della
     shell genitore). Se si usa un ambiente preconfigurato, specificare
     --no-auto-venv. --no-auto-install impedisce le installazioni pip.

NOTA: l'installazione iniziale richiede Internet, git e permessi di scrittura.
      Il modello sentence-transformers potrebbe essere scaricato al primo uso.
      I pacchetti originali di CELLama possono avere vincoli di versione.
      Non è stato effettuato un test completo su dati reali del progetto.

CELLama originale: https://github.com/portrai-io/CELLama
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import venv
from datetime import datetime, timezone
from pathlib import Path

CELLAMA_GIT = 'git+https://github.com/portrai-io/CELLama.git'
DEFAULT_VENV = Path.home() / '.cache' / 'matchacell' / 'cellama-venv'
# Moduli indispensabili per questa pipeline; i requisiti del progetto
# originale vengono installati da setup.py quando pip installa CELLama.
NEEDED = {
    'numpy': 'numpy', 'pandas': 'pandas', 'scipy': 'scipy',
    'scanpy': 'scanpy', 'sklearn': 'scikit-learn',
    'sentence_transformers': 'sentence-transformers',
    'langchain': 'langchain==0.1.16', 'torch': 'torch',
    'xgboost': 'xgboost', 'faiss': 'faiss-cpu',
    'matplotlib': 'matplotlib', 'umap': 'umap-learn',
    'h5py': 'h5py',
}


def say(msg):
    print(msg, flush=True)


def pip_install(requirements, enabled):
    """Installazioni esplicite, senza eseguire comandi di shell arbitrari."""
    if not enabled:
        raise RuntimeError('Dipendenze mancanti: ' + ', '.join(requirements)
                           + '. Rimuovere --no-auto-install per installarle.')
    say('[INSTALL] ' + ', '.join(requirements))
    subprocess.run([sys.executable, '-m', 'pip', 'install', *requirements], check=True)
    importlib.invalidate_caches()


def exists(module):
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def bootstrap_venv(args):
    """Crea/riutilizza un venv e riavvia il processo al suo interno.

    Non modifica l'ambiente Python o Conda dell'utente. Un marker di ambiente
    evita cicli di riavvio. --no-auto-venv fa usare il Python corrente.
    """
    if not args.auto_venv or os.environ.get('MATCHACELL_IN_VENV') == '1':
        return
    folder = args.venv_dir.expanduser().resolve()
    vpy = folder / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
    if not vpy.is_file():
        say(f'[ENV] Creo ambiente virtuale: {folder}')
        folder.parent.mkdir(parents=True, exist_ok=True)
        venv.EnvBuilder(with_pip=True, clear=False).create(str(folder))
    else:
        say(f'[ENV] Ambiente virtuale già presente: {folder}')
    say('[ENV] Eseguo lo script nel Python isolato')
    child_env = dict(os.environ, MATCHACELL_IN_VENV='1')
    # Propaga esattamente gli argomenti originali; nessun 'activate' necessario.
    completed = subprocess.run([str(vpy), str(Path(__file__).resolve()), *sys.argv[1:]], env=child_env)
    sys.exit(completed.returncode)


def dependency_check(args):
    """Installa una volta CELLama originale e ripara moduli mancanti.

    Nota: trovare un modulo non garantisce al 100% la compatibilità di versione.
    Un test di import finale rileva almeno incompatibilità evidenti.
    """
    say('[0/5] Controllo dipendenze e CELLama')
    if sys.version_info < (3, 9):
        raise RuntimeError('Serve Python >= 3.9.')
    if sys.version_info >= (3, 12):
        say('[WARN] Per i pin storici di CELLama è preferibile Python 3.9–3.11.')

    if args.cellama_dir:
        source = args.cellama_dir.expanduser().resolve()
        if not (source / 'setup.py').is_file() or not (source / 'cellama.py').is_file():
            raise ValueError(f'Non è una cartella CELLama valida: {source}')
        cellama_spec = str(source)
    else:
        cellama_spec = CELLAMA_GIT

    # Installare il repo originale tramite pip install importa i suoi moduli,
    # incluso _nn_model.py, senza copiarli/modificarli nella nostra directory.
    if not exists('cellama'):
        pip_install([cellama_spec], args.auto_install)
    else:
        say('  [OK] CELLama già installato (codice originale)')

    missing = [spec for module, spec in NEEDED.items() if not exists(module)]
    if missing:
        pip_install(missing, args.auto_install)
    else:
        say('  [OK] Tutti i moduli richiesti risultano presenti')

    # L'import reale valida i requisiti più critici. Fallire con errore chiaro
    # è più sicuro che aggiornare silenziosamente decine di versioni.
    check = ('from cellama import lm_cell_integrated_embed, '
             'setup_knn_search, find_nearest_cell_types; '
             'import scanpy, matplotlib, faiss, sentence_transformers')
    say('  Verifica compatibilità degli import ...')
    subprocess.run([sys.executable, '-c', check], check=True)
    say('  [OK] Import completati')


def read_data(path, sc):
    """.h5ad = AnnData; .h5 = SOLO matrice 10x Genomics standard."""
    if not path.is_file():
        raise FileNotFoundError(f'File non trovato: {path}')
    if path.suffix.lower() == '.h5ad':
        return sc.read_h5ad(path)
    if path.suffix.lower() == '.h5':
        try:
            return sc.read_10x_h5(path, gex_only=True)
        except Exception as exc:
            raise ValueError(f'{path} non è un H5 10x leggibile; convertire in .h5ad') from exc
    raise ValueError('Formati accettati: .h5ad e .h5 (10x)')


def run_analysis(args):
    # Import pesanti solo DOPO creazione venv e verifica installazione.
    import matplotlib
    matplotlib.use('Agg')  # Funziona anche via SSH su server senza display.
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    # FUNZIONI ORIGINALI DEL REPOSITORY CELLama (non sono riscritte qui).
    from cellama import lm_cell_integrated_embed, setup_knn_search, find_nearest_cell_types

    args.output.mkdir(parents=True, exist_ok=True)
    say('[1/5] Caricamento dati')
    target = read_data(args.input, sc)
    reference = read_data(args.reference, sc)
    if args.reference_label not in reference.obs:
        raise KeyError(f'Colonna reference {args.reference_label!r} assente. '
                       f'Disponibili: {list(reference.obs.columns)}')
    if reference.obs[args.reference_label].isna().any():
        raise ValueError('La reference contiene tipi cellulari mancanti (NaN).')
    if target.n_obs < 4 or reference.n_obs < 1:
        raise ValueError('Servono almeno 4 cellule query e 1 cellula reference.')
    if not target.var_names.is_unique or not reference.var_names.is_unique:
        raise ValueError('Gene names duplicati: rendere univoci var_names.')
    common = reference.var_names.intersection(target.var_names, sort=False)
    if len(common) < args.min_common_genes:
        raise ValueError(f'Solo {len(common)} geni comuni; richiesti almeno {args.min_common_genes}. '
                         'Controllare gli identificatori genici (symbol vs Ensembl).')
    if target.X is None or reference.X is None:
        raise ValueError('Entrambi i dataset devono avere la matrice di espressione in .X')
    say(f'  Query: {target.n_obs:,} cellule; reference: {reference.n_obs:,}; geni comuni: {len(common):,}')

    say('[2/5] Embedding con CELLama originale')
    # Usiamo geni condivisi per rendere i testi comparabili nella reference/query.
    # CELLama seleziona i geni più espressi per cellula e genera un testo;
    # il language model converte poi ogni descrizione in un vettore numerico.
    # CELLama.get_top_genes() usa .toarray(): serve una matrice CSR.
    ref_for_model = reference[:, common].copy()
    query_for_model = target[:, common].copy()
    for a in (ref_for_model, query_for_model):
        a.X = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(np.asarray(a.X))
    ref_emb, query_emb = lm_cell_integrated_embed(
        ref_for_model, query_for_model, model_name=args.model,
        top_k=min(args.top_genes, len(common)), use_intersect=False,
        obs_features=None,
    )
    ref_matrix = np.asarray(ref_emb.X, dtype=np.float32)
    query_matrix = np.asarray(query_emb.X, dtype=np.float32)
    if not np.isfinite(ref_matrix).all() or not np.isfinite(query_matrix).all():
        raise ValueError('Embedding CELLama contengono valori non finiti.')
    target.obsm['X_cellama'] = query_matrix

    say('[3/5] Trasferimento del cell type (kNN originale)')
    # k=1 (nearest-neighbor cosine), come helper originale CELLama.
    # reset_index è necessario perché l'helper usa etichette via series[i].
    knn = setup_knn_search(ref_matrix, n_neighbors=1)
    labels = reference.obs[args.reference_label].astype(str).reset_index(drop=True)
    predicted, distances = find_nearest_cell_types(knn, query_matrix, labels)
    target.obs['cellama_cell_type'] = pd.Categorical(predicted)
    # Distanza cosine: minore = riferimento più simile; NON è probabilità.
    target.obs['cellama_nn_distance'] = np.asarray(distances, dtype=float)

    say('[4/5] Calcolo UMAP dagli embedding CELLama')
    # UMAP descrittiva: non prova di identità cellulare o accuratezza.
    sc.pp.neighbors(target, use_rep='X_cellama',
                    n_neighbors=min(args.umap_neighbors, target.n_obs - 1),
                    key_added='cellama_neighbors')
    sc.tl.umap(target, neighbors_key='cellama_neighbors', random_state=args.random_state)

    say('[5/5] Salvataggio output')
    annotated = args.output / 'cellama_annotated.h5ad'
    target.write_h5ad(annotated, compression='gzip')
    fig = sc.pl.umap(target, color='cellama_cell_type', show=False,
                     title='CELLama | predicted cell types', frameon=False, return_fig=True)
    fig.savefig(args.output / 'cellama_umap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    target.obs['cellama_cell_type'].value_counts().rename_axis('cell_type').reset_index(
        name='n_cells').to_csv(args.output / 'cell_type_counts.csv', index=False)
    info = {
        'input': str(args.input.resolve()), 'reference': str(args.reference.resolve()),
        'reference_label': args.reference_label, 'model': args.model,
        'top_genes': args.top_genes, 'cells_query': target.n_obs,
        'cells_reference': reference.n_obs, 'common_genes': len(common),
        'annotation_method': 'Original CELLama embedding + cosine nearest-neighbor k=1',
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
    }
    (args.output / 'run_info.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    say(f'FINITO! File: {annotated}\nUMAP: {args.output / "cellama_umap.png"}')


def main():
    parser = argparse.ArgumentParser(
        description='MatchACell: CELLama originale, annotazione reference-based e UMAP.')
    parser.add_argument('--input', type=Path, help='Dataset da annotare (.h5ad o 10x .h5)')
    parser.add_argument('--reference', type=Path, help='Dataset reference con tipi cellulari noti (.h5ad)')
    parser.add_argument('--reference-label', default='cell_type',
                        help='Nome colonna con i tipi cellulari nella reference (default cell_type)')
    parser.add_argument('--output', type=Path, help='Directory output')
    parser.add_argument('--model', default='all-MiniLM-L6-v2', help='Sentence-transformer CELLama')
    parser.add_argument('--top-genes', type=int, default=30, help='Geni più espressi per cellula')
    parser.add_argument('--umap-neighbors', type=int, default=15)
    parser.add_argument('--random-state', type=int, default=42)
    parser.add_argument('--min-common-genes', type=int, default=30)
    parser.add_argument('--cellama-dir', type=Path,
                        help='Opzionale: checkout locale IMMUTATO della repository originale')
    parser.add_argument('--venv-dir', type=Path, default=DEFAULT_VENV,
                        help=f'Dove creare/riutilizzare ambiente virtuale (default {DEFAULT_VENV})')
    parser.add_argument('--no-auto-venv', dest='auto_venv', action='store_false',
                        help='Usa ambiente Python corrente (consigliato su cluster con ambienti gestiti)')
    parser.add_argument('--no-auto-install', dest='auto_install', action='store_false',
                        help='Non installare pacchetti: fallisci se mancano')
    parser.add_argument('--check-deps', action='store_true',
                        help='Prepara/controlla ambiente senza eseguire annotazione')
    parser.set_defaults(auto_venv=True, auto_install=True)
    args = parser.parse_args()
    if not args.check_deps and any(v is None for v in (args.input, args.reference, args.output)):
        parser.error('Servono --input, --reference e --output (tranne con --check-deps)')
    if args.top_genes < 1 or args.umap_neighbors < 2 or args.min_common_genes < 1:
        parser.error('--top-genes >=1, --umap-neighbors >=2, --min-common-genes >=1')
    if sys.version_info < (3, 9):
        raise RuntimeError('Python >=3.9 necessario per lo script.')
    bootstrap_venv(args)
    dependency_check(args)
    if args.check_deps:
        say('Controllo dipendenze concluso.')
        return
    run_analysis(args)


if __name__ == '__main__':
    try:
        main()
    except (Exception, subprocess.CalledProcessError) as e:
        say(f'ERRORE: {e}')
        sys.exit(1)
