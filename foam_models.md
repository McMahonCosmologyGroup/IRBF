# Foam alternative models (sandbox; main.ipynb data and likelihood, nlive 1000)

Delta BIC and Delta ln Z (nested sampling) relative to R+M+a_RM fitted the same way. Full Mie: implied cell volume fraction f_vol = 2 C d exp(s^2/2) / 3 for spheres. PP35 data may not be good yet.

## Styrofoam_batch1 (N = 765, n_eff = 1.0197)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 696.9/760 = 0.92 | +0.0 | +0.0 | tan_d=1.87e-05, C=2.64, K_Ray=1.99e+03, K_Mie=2.35e+03, a_RM=0.672 |
| R+M, no coupling | 4 | 780.7/761 = 1.03 | +77.1 | -42.8 | tan_d=2.61e-06, K_R=1.78e+03, C_M=4.83, K_Mie=2.91e+03 |
| R+M+C_RM (physical cut) | 5 | 699.0/760 = 0.92 | +2.1 | -4.1 | tan_d=6.13e-06, K_R=1.61e+03, C_M=6.15, K_Mie=3.47e+03, C_RM=2.49 |
| Full Mie, one population, m free | 5 | 746.5/760 = 0.98 | +49.6 | -24.9 | tan_d=0.000272, C=82.8, d=0.0306, s=0.685, m=0.832; f_vol=0.214 |
| Full Mie, m free, wide C prior | 5 | 739.4/760 = 0.97 | +42.5 | -20.9 | tan_d=0.000265, C=861, d=0.0178, s=0.881, m=0.942; f_vol=1.502 |
| Full Mie, m = 1/n_eff | 4 | 1315.5/761 = 1.73 | +612.0 | -314.4 | tan_d=3.53e-06, C=99.6, d=0.259, s=0.0365; f_vol=1.723 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 739.6/761 = 0.97 | +36.1 | -19.2 | tan_d=0.00026, C=7.86e+03, d=0.0184, s=0.883; f_vol=14.236 |
| Full Mie + Rayleigh | 6 | 726.3/759 = 0.96 | +36.0 | -16.6 | tan_d=0.000208, C=47.2, d=0.0347, s=0.869, m=0.885, K_R=2.09e+03; f_vol=0.159 |
| Full Mie, two populations | 7 | 693.3/758 = 0.91 | +9.7 | -5.6 | tan_d=3.59e-05, C1=47.5, d1=0.0727, C2=0.0246, d2=1.56, s=0.178, m=0.867 |

## HD30_batch1 (N = 340, n_eff = 1.0185)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 402.4/335 = 1.20 | +0.0 | +0.0 | tan_d=2.4e-06, C=0.0169, K_Ray=403, K_Mie=489, a_RM=0.114 |
| R+M, no coupling | 4 | 1063.6/336 = 3.17 | +655.4 | -320.6 | tan_d=6.56e-07, K_R=1.5e+03, C_M=2.82, K_Mie=2.88e+03 |
| R+M+C_RM (physical cut) | 5 | 927.4/335 = 2.77 | +525.0 | -258.0 | tan_d=6.12e-07, K_R=1.47e+03, C_M=4.1, K_Mie=3.53e+03, C_RM=1.17 |
| Full Mie, one population, m free | 5 | 348.6/335 = 1.04 | -53.7 | +36.5 | tan_d=1.79e-06, C=74.2, d=0.0399, s=0.606, m=0.838; f_vol=0.237 |
| Full Mie, m free, wide C prior | 5 | 346.0/335 = 1.03 | -56.3 | +39.4 | tan_d=1.29e-06, C=438, d=0.0185, s=0.846, m=0.914; f_vol=0.773 |
| Full Mie, m = 1/n_eff | 4 | 4431.2/336 = 13.19 | +4023.0 | -2013.8 | tan_d=1.05e-07, C=100, d=0.267, s=0.029; f_vol=1.782 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 345.8/336 = 1.03 | -62.4 | +40.0 | tan_d=1.42e-06, C=8.05e+03, d=0.0246, s=0.785; f_vol=17.965 |
| Full Mie + Rayleigh | 6 | 346.4/334 = 1.04 | -50.2 | +36.8 | tan_d=1.66e-06, C=74.2, d=0.0432, s=0.585, m=0.848, K_R=8.12e+03; f_vol=0.253 |
| Full Mie, two populations | 7 | 348.7/333 = 1.05 | -42.1 | +37.2 | tan_d=1.75e-06, C1=61.4, d1=0.0351, C2=0.00516, d2=0.0705, s=0.583, m=0.852 |

## HD30_batch2 (N = 255, n_eff = 1.0142)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 240.5/250 = 0.96 | +0.0 | +0.0 | tan_d=2.46e-06, C=0.0375, K_Ray=577, K_Mie=535, a_RM=0.104 |
| R+M, no coupling | 4 | 2306.7/251 = 9.19 | +2060.7 | -1016.6 | tan_d=2.26e-07, K_R=1.76e+03, C_M=3.49, K_Mie=3.16e+03 |
| R+M+C_RM (physical cut) | 5 | 1247.9/250 = 4.99 | +1007.4 | -495.0 | tan_d=2.32e-07, K_R=1.65e+03, C_M=4.65, K_Mie=3.85e+03, C_RM=2.18 |
| Full Mie, one population, m free | 5 | 222.9/250 = 0.89 | -17.5 | +26.5 | tan_d=2.95e-05, C=92.6, d=0.0267, s=0.707, m=0.835; f_vol=0.212 |
| Full Mie, m free, wide C prior | 5 | 202.3/250 = 0.81 | -38.2 | +37.0 | tan_d=7.61e-06, C=1.64e+03, d=0.0204, s=0.826, m=0.96; f_vol=3.143 |
| Full Mie, m = 1/n_eff | 4 | 15139.9/251 = 60.32 | +14893.9 | -7444.5 | tan_d=3.56e-08, C=100, d=0.312, s=0.0212; f_vol=2.081 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 200.7/251 = 0.80 | -45.3 | +40.6 | tan_d=8.56e-06, C=1.79e+04, d=0.0156, s=0.908; f_vol=28.102 |
| Full Mie + Rayleigh | 6 | 195.8/249 = 0.79 | -39.1 | +33.9 | tan_d=7.73e-06, C=60.6, d=0.0365, s=0.66, m=0.857, K_R=2.84e+03; f_vol=0.183 |
| Full Mie, two populations | 7 | 209.4/248 = 0.84 | -20.0 | +29.4 | tan_d=1.47e-05, C1=91.7, d1=0.024, C2=0.0113, d2=0.258, s=0.721, m=0.817 |

## HD30_batch3 (N = 425, n_eff = 1.0147)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 458.9/420 = 1.09 | +0.0 | +0.0 | tan_d=1.08e-06, C=0.00517, K_Ray=366, K_Mie=504, a_RM=0.113 |
| R+M, no coupling | 4 | 2894.6/421 = 6.88 | +2429.6 | -1211.0 | tan_d=1.07e-07, K_R=1.9e+03, C_M=1.47, K_Mie=2.98e+03 |
| R+M+C_RM (physical cut) | 5 | 2126.7/420 = 5.06 | +1667.7 | -834.9 | tan_d=9.59e-08, K_R=1.82e+03, C_M=2.32, K_Mie=3.84e+03, C_RM=1.61 |
| Full Mie, one population, m free | 5 | 415.4/420 = 0.99 | -43.5 | +28.6 | tan_d=4.86e-07, C=83, d=0.027, s=0.703, m=0.858; f_vol=0.191 |
| Full Mie, m free, wide C prior | 5 | 405.9/420 = 0.97 | -53.0 | +34.4 | tan_d=3.97e-07, C=1.35e+03, d=0.0119, s=0.959, m=0.948; f_vol=1.692 |
| Full Mie, m = 1/n_eff | 4 | 8143.5/421 = 19.34 | +7678.5 | -3845.7 | tan_d=2.8e-08, C=100, d=0.247, s=0.022; f_vol=1.649 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 405.8/421 = 0.96 | -59.2 | +34.3 | tan_d=4.45e-07, C=1.35e+04, d=0.0141, s=0.918; f_vol=19.375 |
| Full Mie + Rayleigh | 6 | 389.7/419 = 0.93 | -63.1 | +38.1 | tan_d=6.48e-07, C=3.67, d=0.122, s=0.0737, m=0.833, K_R=2.46e+03; f_vol=0.030 |
| Full Mie, two populations | 7 | 415.4/418 = 0.99 | -31.4 | +28.4 | tan_d=6.23e-07, C1=80.7, d1=0.0241, C2=0.0021, d2=0.051, s=0.725, m=0.845 |

## PP35_batch1 (N = 320, n_eff = 1.0190)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 335.0/315 = 1.06 | +0.0 | +0.0 | tan_d=0.000128, C=0.000664, K_Ray=702, K_Mie=3.24e+03, a_RM=0.224 |
| R+M, no coupling | 4 | 263.5/316 = 0.83 | -77.3 | +28.0 | tan_d=0.000102, K_R=4.06e+03, C_M=0.00159, K_Mie=111 |
| R+M+C_RM (physical cut) | 5 | 335.0/315 = 1.06 | -0.0 | -3.5 | tan_d=0.000128, K_R=4.49e+03, C_M=0.00531, K_Mie=3.26e+03, C_RM=0.0505 |
| Full Mie, one population, m free | 5 | 335.3/315 = 1.06 | +0.3 | -1.1 | tan_d=0.000128, C=39, d=0.0206, s=0.202, m=0.855; f_vol=0.055 |
| Full Mie, m free, wide C prior | 5 | 335.1/315 = 1.06 | +0.1 | -0.2 | tan_d=0.000128, C=198, d=0.0154, s=0.203, m=0.892; f_vol=0.208 |
| Full Mie, m = 1/n_eff | 4 | 340.2/316 = 1.08 | -0.5 | -5.7 | tan_d=0.000127, C=81.6, d=0.0523, s=0.102; f_vol=0.286 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 335.1/316 = 1.06 | -5.7 | -0.3 | tan_d=0.000128, C=4.14e+03, d=0.0157, s=0.208; f_vol=4.444 |
| Full Mie + Rayleigh | 6 | 261.9/314 = 0.83 | -67.4 | +29.4 | tan_d=9.9e-05, C=0.002, d=3.93, s=0.144, m=0.884, K_R=4.05e+03; f_vol=0.001 |
| Full Mie, two populations | 7 | 245.8/313 = 0.79 | -77.7 | +34.6 | tan_d=6.12e-05, C1=2, d1=0.068, C2=0.0035, d2=3.91, s=0.292, m=0.892 |

## LD15_batch1 (N = 83, n_eff = 1.0090)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 1196.2/78 = 15.34 | +0.0 | +0.0 | tan_d=3.78e-08, C=0.0145, K_Ray=577, K_Mie=517, a_RM=0.107 |
| R+M, no coupling | 4 | 20073.1/79 = 254.09 | +18872.5 | -9424.8 | tan_d=5.32e-09, K_R=2.35e+03, C_M=1.76, K_Mie=3.6e+03 |
| R+M+C_RM (physical cut) | 5 | 10161.7/78 = 130.28 | +8965.4 | -4479.6 | tan_d=5.71e-09, K_R=2.15e+03, C_M=1.91, K_Mie=3.98e+03, C_RM=2.57 |
| Full Mie, one population, m free | 5 | 833.9/78 = 10.69 | -362.4 | +195.0 | tan_d=1.36e-05, C=99, d=0.0193, s=0.846, m=0.896; f_vol=0.183 |
| Full Mie, m free, wide C prior | 5 | 647.7/78 = 8.30 | -548.5 | +285.5 | tan_d=1.2e-05, C=3.49e+03, d=0.0128, s=0.986, m=0.98; f_vol=4.862 |
| Full Mie, m = 1/n_eff | 4 | 79975.4/79 = 1012.35 | +78774.7 | -39387.3 | tan_d=1.9e-09, C=100, d=0.293, s=0.0204; f_vol=1.956 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 642.9/79 = 8.14 | -557.8 | +295.0 | tan_d=1.19e-05, C=1.72e+04, d=0.0129, s=0.988; f_vol=24.138 |
| Full Mie + Rayleigh | 6 | 313.7/77 = 4.07 | -878.1 | +455.7 | tan_d=4.3e-06, C=70.7, d=0.0305, s=0.816, m=0.935, K_R=3.27e+03; f_vol=0.200 |
| Full Mie, two populations | 7 | 314.5/76 = 4.14 | -872.9 | +448.9 | tan_d=1.7e-06, C1=86.1, d1=0.024, C2=1.4, d2=0.117, s=0.461, m=0.858 |

## LD15_batch2 (N = 170, n_eff = 1.0063)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 119.0/165 = 0.72 | +0.0 | +0.0 | tan_d=1.82e-06, C=0.0116, K_Ray=550, K_Mie=526, a_RM=0.105 |
| R+M, no coupling | 4 | 899.5/166 = 5.42 | +775.3 | -379.9 | tan_d=1.68e-07, K_R=2.34e+03, C_M=1.21, K_Mie=2.97e+03 |
| R+M+C_RM (physical cut) | 5 | 430.1/165 = 2.61 | +311.0 | -151.8 | tan_d=2.19e-07, K_R=2.15e+03, C_M=1.64, K_Mie=3.74e+03, C_RM=2.2 |
| Full Mie, one population, m free | 5 | 110.6/165 = 0.67 | -8.5 | +16.3 | tan_d=1.01e-05, C=86.3, d=0.0195, s=0.801, m=0.87; f_vol=0.155 |
| Full Mie, m free, wide C prior | 5 | 100.6/165 = 0.61 | -18.4 | +21.5 | tan_d=5.88e-06, C=1.91e+03, d=0.0121, s=0.956, m=0.97; f_vol=2.442 |
| Full Mie, m = 1/n_eff | 4 | 6237.5/166 = 37.58 | +6113.4 | -3057.7 | tan_d=4.23e-08, C=99.9, d=0.406, s=0.0258; f_vol=2.708 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 102.1/166 = 0.62 | -22.1 | +21.7 | tan_d=6.2e-06, C=2.72e+04, d=0.0191, s=0.84; f_vol=49.437 |
| Full Mie + Rayleigh | 6 | 86.6/164 = 0.53 | -27.3 | +26.5 | tan_d=2.36e-06, C=22.2, d=0.0376, s=0.691, m=0.882, K_R=3.17e+03; f_vol=0.070 |
| Full Mie, two populations | 7 | 86.3/163 = 0.53 | -22.5 | +25.5 | tan_d=2.2e-06, C1=34.1, d1=0.0444, C2=1.28, d2=0.154, s=0.0793, m=0.863 |

## LD24_batch1 (N = 581, n_eff = 1.0121)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 683.9/576 = 1.19 | +0.0 | +0.0 | tan_d=4.41e-07, C=0.0131, K_Ray=489, K_Mie=492, a_RM=0.11 |
| R+M, no coupling | 4 | 3246.0/577 = 5.63 | +2555.8 | -1272.1 | tan_d=5.67e-08, K_R=1.92e+03, C_M=1.2, K_Mie=2.92e+03 |
| R+M+C_RM (physical cut) | 5 | 2567.6/576 = 4.46 | +1883.8 | -941.0 | tan_d=5.45e-08, K_R=1.86e+03, C_M=1.95, K_Mie=3.8e+03, C_RM=1.48 |
| Full Mie, one population, m free | 5 | 554.7/576 = 0.96 | -129.2 | +74.0 | tan_d=2.57e-06, C=93.3, d=0.0225, s=0.811, m=0.881; f_vol=0.195 |
| Full Mie, m free, wide C prior | 5 | 527.0/576 = 0.91 | -156.9 | +85.7 | tan_d=1.5e-06, C=1.48e+03, d=0.0152, s=0.942, m=0.966; f_vol=2.347 |
| Full Mie, m = 1/n_eff | 4 | 9496.9/577 = 16.46 | +8806.7 | -4407.9 | tan_d=1.65e-08, C=100, d=0.249, s=0.0227; f_vol=1.657 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 525.6/577 = 0.91 | -164.6 | +90.5 | tan_d=1.26e-06, C=1.2e+04, d=0.0149, s=0.953; f_vol=18.880 |
| Full Mie + Rayleigh | 6 | 503.6/575 = 0.88 | -173.9 | +96.8 | tan_d=5e-07, C=25.3, d=0.0769, s=0.497, m=0.927, K_R=2.68e+03; f_vol=0.147 |
| Full Mie, two populations | 7 | 496.3/574 = 0.86 | -174.8 | +89.8 | tan_d=5.39e-07, C1=71.6, d1=0.036, C2=1.84, d2=0.145, s=0.262, m=0.859 |

## LD24_batch2 (N = 255, n_eff = 1.0123)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 210.1/250 = 0.84 | +0.0 | +0.0 | tan_d=1.13e-05, C=0.000225, K_Ray=183, K_Mie=489, a_RM=0.114 |
| R+M, no coupling | 4 | 1042.7/251 = 4.15 | +827.1 | -409.6 | tan_d=1.59e-07, K_R=2.1e+03, C_M=0.867, K_Mie=2.87e+03 |
| R+M+C_RM (physical cut) | 5 | 811.9/250 = 3.25 | +601.8 | -300.6 | tan_d=1.51e-07, K_R=2.02e+03, C_M=1.3, K_Mie=3.6e+03, C_RM=1.34 |
| Full Mie, one population, m free | 5 | 194.1/250 = 0.78 | -16.0 | +16.1 | tan_d=1.52e-06, C=68.9, d=0.023, s=0.762, m=0.857; f_vol=0.141 |
| Full Mie, m free, wide C prior | 5 | 191.2/250 = 0.76 | -18.9 | +17.7 | tan_d=1.49e-06, C=515, d=0.0162, s=0.881, m=0.943; f_vol=0.818 |
| Full Mie, m = 1/n_eff | 4 | 2055.1/251 = 8.19 | +1839.5 | -924.5 | tan_d=6.62e-08, C=99.9, d=0.224, s=0.0264; f_vol=1.493 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 191.1/251 = 0.76 | -24.5 | +18.1 | tan_d=1.35e-06, C=8.04e+03, d=0.0193, s=0.84; f_vol=14.698 |
| Full Mie + Rayleigh | 6 | 185.1/249 = 0.74 | -19.4 | +17.3 | tan_d=2.11e-06, C=8.9, d=0.115, s=0.191, m=0.906, K_R=2.85e+03; f_vol=0.070 |
| Full Mie, two populations | 7 | 186.2/248 = 0.75 | -12.7 | +16.1 | tan_d=1.57e-06, C1=55.8, d1=0.0282, C2=0.0038, d2=0.0823, s=0.687, m=0.872 |

## HD80_batch1 (N = 402, n_eff = 1.0378)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 1075.8/397 = 2.71 | +0.0 | +0.0 | tan_d=2.42e-07, C=0.000519, K_Ray=106, K_Mie=460, a_RM=0.122 |
| R+M, no coupling | 4 | 9511.6/398 = 23.90 | +8429.8 | -4204.2 | tan_d=1.57e-07, K_R=867, C_M=3.82, K_Mie=2.71e+03 |
| R+M+C_RM (physical cut) | 5 | 8641.4/397 = 21.77 | +7565.6 | -3778.9 | tan_d=1.23e-07, K_R=862, C_M=7.88, K_Mie=3.77e+03, C_RM=1.2 |
| Full Mie, one population, m free | 5 | 1167.5/397 = 2.94 | +91.7 | -38.1 | tan_d=1.92e-07, C=95.3, d=0.115, s=0.0306, m=0.876; f_vol=0.729 |
| Full Mie, m free, wide C prior | 5 | 1113.2/397 = 2.80 | +37.4 | -12.5 | tan_d=1.98e-07, C=2.44e+03, d=0.121, s=0.0322, m=0.978; f_vol=19.642 |
| Full Mie, m = 1/n_eff | 4 | 44506.4/398 = 111.83 | +43424.6 | -21716.9 | tan_d=2.39e-08, C=100, d=0.296, s=0.0203; f_vol=1.971 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 1125.9/398 = 2.83 | +44.1 | -13.9 | tan_d=1.96e-07, C=927, d=0.12, s=0.0318; f_vol=7.419 |
| Full Mie + Rayleigh | 6 | 1167.6/396 = 2.95 | +97.8 | -38.7 | tan_d=1.98e-07, C=95.4, d=0.115, s=0.0309, m=0.876, K_R=1.21e+04; f_vol=0.731 |
| Full Mie, two populations | 7 | 1167.4/395 = 2.96 | +103.6 | -38.1 | tan_d=1.93e-07, C1=0.0415, d1=0.0313, C2=95.4, d2=0.115, s=0.0302, m=0.876 |

## HD110_batch1 (N = 69, n_eff = 1.0586)

| Model | k | chi2/dof | Delta BIC | Delta ln Z | parameters (median) |
|---|---|---|---|---|---|
| R+M+a_RM (published form) | 5 | 61.9/64 = 0.97 | +0.0 | +0.0 | tan_d=2.77e-05, C=0.0182, K_Ray=112, K_Mie=168, a_RM=0.112 |
| R+M, no coupling | 4 | 751.4/65 = 11.56 | +685.3 | -330.1 | tan_d=6.73e-06, K_R=403, C_M=23.2, K_Mie=1.78e+03 |
| R+M+C_RM (physical cut) | 5 | 424.7/64 = 6.64 | +362.9 | -173.4 | tan_d=6.41e-06, K_R=373, C_M=68.6, K_Mie=3.34e+03, C_RM=23.8 |
| Full Mie, one population, m free | 5 | 57.6/64 = 0.90 | -4.2 | +17.8 | tan_d=2.3e-05, C=49.6, d=0.28, s=0.0934, m=0.851; f_vol=0.930 |
| Full Mie, m free, wide C prior | 5 | 57.2/64 = 0.89 | -4.7 | +19.3 | tan_d=1.81e-05, C=2.37e+03, d=0.285, s=0.136, m=0.98; f_vol=45.433 |
| Full Mie, m = 1/n_eff | 4 | 1919.7/65 = 29.53 | +1853.6 | -922.1 | tan_d=1.5e-06, C=100, d=0.445, s=0.0261; f_vol=2.970 |
| Full Mie, m = 1/n_eff, wide C prior | 4 | 57.4/65 = 0.88 | -8.7 | +19.1 | tan_d=2.76e-05, C=291, d=0.294, s=0.0957; f_vol=5.732 |
| Full Mie + Rayleigh | 6 | 57.7/63 = 0.92 | +0.0 | +17.9 | tan_d=2.25e-05, C=58.4, d=0.282, s=0.0886, m=0.864, K_R=9.87e+03; f_vol=1.101 |
| Full Mie, two populations | 7 | 57.7/62 = 0.93 | +4.3 | +17.5 | tan_d=2.47e-05, C1=0.109, d1=0.0401, C2=44.1, d2=0.28, s=0.0856, m=0.841 |

