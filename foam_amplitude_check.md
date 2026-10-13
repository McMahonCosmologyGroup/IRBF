# Foam amplitude check: full Mie with m = 1/n_eff (air cells in the foam medium), wide C prior

If the scattering comes from polymer walls (index n_p, volume fraction phi = (n_eff-1)/(n_p-1)) rather than from air cells of
contrast 1-1/n_eff, the Rayleigh-Gans amplitude ratio is phi (n_p-1)^2 / (1-1/n_eff)^2 ~ (n_p-1)/(n_eff-1). The fitted "cell volume
fraction" f_vol = 2 C d exp(s^2/2) / 3 should then be about that ratio (minus dense-packing suppression).

| Batch | polymer (n_p) | n_eff | phi | f_vol fitted | predicted ratio | fitted / predicted | d [um], s | chi2/dof vs R+M+a_RM |
|---|---|---|---|---|---|---|---|---|
| Styrofoam_batch1 | PS (1.59) | 1.0197 | 0.033 | 14.2 | 31.2 | 0.46 | 18, 0.88 | 0.97 vs 0.92 |
| HD30_batch1 | LDPE (1.514) | 1.0185 | 0.036 | 18.0 | 28.8 | 0.62 | 25, 0.78 | 1.03 vs 1.20 |
| HD30_batch2 | LDPE (1.514) | 1.0142 | 0.028 | 28.1 | 37.2 | 0.75 | 16, 0.91 | 0.80 vs 0.96 |
| HD30_batch3 | LDPE (1.514) | 1.0147 | 0.029 | 19.4 | 35.9 | 0.54 | 14, 0.92 | 0.96 vs 1.09 |
| PP35_batch1 | PP (1.5) | 1.0190 | 0.038 | 4.5 | 27.3 | 0.16 | 16, 0.21 | 1.06 vs 1.06 |
| LD15_batch1 | LDPE (1.514) | 1.0090 | 0.018 | 24.1 | 57.9 | 0.42 | 13, 0.99 | 8.14 vs 15.34 |
| LD15_batch2 | LDPE (1.514) | 1.0063 | 0.012 | 49.4 | 83.1 | 0.59 | 19, 0.84 | 0.62 vs 0.72 |
| LD24_batch1 | LDPE (1.514) | 1.0121 | 0.024 | 18.9 | 43.6 | 0.43 | 15, 0.95 | 0.91 vs 1.19 |
| LD24_batch2 | LDPE (1.514) | 1.0123 | 0.024 | 14.7 | 42.7 | 0.34 | 19, 0.84 | 0.76 vs 0.84 |
| HD80_batch1 | LDPE (1.514) | 1.0378 | 0.073 | 7.4 | 14.7 | 0.51 | 120, 0.03 | 2.83 vs 2.71 |
| HD110_batch1 | LDPE (1.514) | 1.0586 | 0.114 | 5.7 | 9.8 | 0.58 | 294, 0.10 | 0.88 vs 0.97 |
