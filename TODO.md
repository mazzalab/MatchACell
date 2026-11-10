# TODO

## Cluster stability: per-resolution CSV checkpoints, Excel at the end

**Where:** `workflow/scripts/matchacell_cluster_stability.py`, `run_all_stability()`

**Problem:** `Cluster_Stability_Summary.xlsx` is opened with `pd.ExcelWriter` before
the bootstrap loop starts. Each resolution's sheet is only flushed to disk when the
writer closes, after the last resolution. Until then the file is 0 bytes. On large
inputs the loop runs for many hours: about 2 h 10 min per resolution with 200
iterations on about 330k cells on CPU, so about 17 h for 8 resolutions. If the job
dies partway (crash, node failure, PBS walltime, kill), the resolutions that already
finished are lost as well.

**Proposed change:**
- [x] When each resolution finishes, write its scores right away to
      `stability/csv/stability_<resolution>.csv`, with the columns `Cluster` and
      `Jaccard_Stability`. Done: the CPU bootstrap scores every resolution on the
      same subsamples, so it rewrites every CSV with the running means after each
      batch of iterations (plus an `n_iter_done` column); the GPU path writes each
      resolution as it finishes.
- [x] Build `Cluster_Stability_Summary.xlsx` only at the end, from the per-resolution
      results, with one sheet per resolution as now. Do not keep the `ExcelWriter`
      open for the whole loop.
- [ ] Optional: resume support. On restart, skip resolutions whose CSV already exists
      and was written with the same `n_iter`, `fraction` and seed. Keep these
      parameters in the CSV or in a sidecar file so that results from a different
      configuration are never reused.
- [x] Decide whether the per-resolution CSVs stay as declared outputs of the rule or
      are only intermediate files. Intermediate for now: not declared outputs.
- [x] Update `docs/cli-and-outputs.md` with the new output files.
