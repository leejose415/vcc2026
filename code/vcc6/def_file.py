"""The D/E/F competition file from the three per-context C+G runs (Jose 2026-10-28): each context's prediction without
its control cells, obs['context'] set (run_line.submission_file), concatenated into one h5ad (3 x 300 x 400 = 360,000
cells, raw integer counts, 18,533 genes), then vcc prep against the final phase's gene list and pert_counts.csv
(contexts D, E, F from the manifest beside it). Exit codes checked; the report goes to <out>/prep.json.

  ~/venv/bin/python code/vcc6/def_file.py --runs results/stage1/2026-10-28_DEF --out results/stage1/2026-10-28_DEF/final
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

V = Path.home() / 'vcc'
sys.path.insert(0, str(V / 'code/vcc6'))
FINAL = V / 'data_final/extracted'


def main(a):
    import anndata as ad
    import run_line as RL
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    parts = []
    for x in ('D', 'E', 'F'):
        p = out / f'submission_{x}.h5ad'
        n_kd, n_ctrl = RL.submission_file(Path(a.runs) / x / 'e2e' / 'salt_none' / 'prediction.h5ad', p, context=x)
        print(f'{x}: {n_kd} knockdown cells kept, {n_ctrl} control cells removed', flush=True)
        parts.append(ad.read_h5ad(p))
    full = ad.concat(parts, join='outer', merge='same')
    full.obs_names_make_unique()
    full.obs['context'] = full.obs['context'].astype('category')
    h5 = out / 'prediction_DEF.h5ad'
    full.write_h5ad(h5, compression='gzip')
    print(f'combined: {full.n_obs} cells x {full.n_vars} genes, contexts {full.obs.context.value_counts().to_dict()}', flush=True)
    del full, parts
    vcc = out / 'prediction_DEF.vcc'
    cmd = [str(Path(sys.executable).parent / 'vcc'), 'prep', str(h5), '-g', str(FINAL / 'gene_names.csv'),
           '--perts', str(FINAL / 'pert_counts.csv'), '-o', str(vcc), '--force']
    r = subprocess.run(cmd, capture_output=True, text=True)
    rep = dict(command=' '.join(cmd), returncode=r.returncode, stdout=r.stdout[-6000:], stderr=r.stderr[-6000:],
               vcc=str(vcc), vcc_bytes=vcc.stat().st_size if vcc.exists() else None, h5ad=str(h5), h5ad_bytes=h5.stat().st_size)
    (out / 'prep.json').write_text(json.dumps(rep, indent=1) + '\n')
    print(r.stdout[-3000:]); print(r.stderr[-3000:], file=sys.stderr)
    print(f'vcc prep exit {r.returncode}')
    sys.exit(r.returncode)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--runs', required=True); p.add_argument('--out', required=True)
    main(p.parse_args())
