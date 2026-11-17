# IRBF
IR Blocking Foams figures, data, and analysis procedure

## Fit cache

The nested-sampling fits and derived products are cached in `pkls/*.pkl` (~7 GB). They are gitignored, so `git clone` does not include them. To skip refitting, download the latest snapshot from the `pkls-cache` GitHub release:

```
python pkl_cache.py pull              # fetch pickles missing locally (no login needed)
python pkl_cache.py pull --overwrite  # also replace local pickles that differ from the snapshot
python pkl_cache.py push              # replace the snapshot with your local pkls/ (needs write access)
```

The snapshot is overwritten on every push and no history is kept. The caches are fingerprinted (`fit_cache.py`), so a stale pickle is refit rather than reused.
