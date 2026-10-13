"""Figure-5-style comparison (arXiv 2607.05003, Fig. 5): data from main.ipynb plotted as in the paper (normalized 5 cm
transmittance, attenuation), with the best fits of the published R+M+a_RM (dashed) and one full Lorenz-Mie population
with m = 1/n_eff and a wide C prior (solid). A second figure decomposes both models for one batch.

    python foam_fig5_compare.py
"""
import contextlib
import io
import json
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import foam_models as fm  # noqa: E402  (chdir to the sandbox root)

MODS = [('rma', 'R+M+a$_{RM}$ (published)', '--'), ('fmr_fix_wc', 'Full Mie, $m=1/n_{eff}$', '-')]


def load_notebook():
    nb = json.load(open('main.ipynb', encoding='utf-8'))
    g = {'__name__': '__main__'}
    plt.show = lambda *a, **k: plt.close('all')
    with contextlib.redirect_stdout(io.StringIO()):
        for i in range(2, 26):
            c = nb['cells'][i]
            if c['cell_type'] != 'code':
                continue
            src = '\n'.join(l for l in ''.join(c['source']).split('\n') if not l.strip().startswith(('%', '!')))
            exec(compile(src, f'main cell {i}', 'exec'), g)
            plt.close('all')
    return g


def best(bname, model):
    res = pickle.load(open(fm.cache(bname, model), 'rb'))
    return res.samples[int(np.argmax(res['logl']))]


def chi2dof(b, model, p, tab):
    f, t, w, err, n = b['f_comb'], b['t_comb'], b['width_comb'], b['err_comb'], float(b['n_eff'])
    mod = np.exp(-w * fm.alpha_model(f, p, n, fm.MODELS[model]['code'], tab))
    return float(np.sum(((t - mod) / err)**2)) / (len(f) - len(p))


def main():
    g = load_notebook()
    batches, comb = g['batches'], g['comb_batches']
    cols = [g['styro_batches'], g['HD30_batches'], [k for k in g['other_zote_keys'] if k != 'HD60 2.03 cm']]
    tab = fm.load_table()
    fgrid = np.linspace(150, 2000, 800)

    fig, axes = plt.subplots(2, 3, figsize=(27, 13), sharex=True, sharey='row')
    lw = 2.5
    rows = []
    for ci, keys in enumerate(cols):
        ax1, ax2 = axes[0, ci], axes[1, ci]
        for bname in keys:
            if bname not in batches or bname not in comb:
                continue
            b = batches[bname]
            for sample in b:
                if '(2)' in sample:
                    continue
                s = b[sample]
                width, mask = s['width'], s['mask']
                ls = ':' if 'batch2' in sample and ci == 2 else '-'
                ax1.plot(s['f'], np.where(mask, np.nan, s['t']**(5 / width)), color=s['color'], lw=lw, ls=ls, zorder=-10)
                ax2.plot(s['f'], np.where(mask, np.nan, s['alpha']), color=s['color'], lw=lw, ls=ls, zorder=-10)
            cb, n = comb[bname], float(comb[bname]['n_eff'])
            fmax = float(np.max(cb['f_comb']))
            fg = fgrid[fgrid <= fmax + 1]
            row = [bname]
            for model, _, mls in MODS:
                p = best(bname, model)
                a = fm.alpha_model(fg, p, n, fm.MODELS[model]['code'], tab)
                kw = dict(color='k' if model == 'rma' else 'crimson', ls=mls, lw=2.5, zorder=10)
                ax1.plot(fg, np.exp(-5 * a), **kw)
                ax2.plot(fg, a, **kw)
                row.append(chi2dof(cb, model, p, tab))
            rows.append(row)

    axes[0, 0].set_xlim(150, 2000)
    axes[0, 0].set_ylim(0, 1.02)
    axes[1, 0].set_ylim(0, 2.5)
    axes[0, 0].set_ylabel('Normalized 5 cm\nTransmittance [-]', fontsize=26)
    axes[1, 0].set_ylabel(r'Attenuation [cm$^{-1}$]', fontsize=26)
    for ax in axes[1, :]:
        ax.set_xlabel('Frequency [GHz]', fontsize=26)
    for ax in axes.flatten():
        ax.grid(True)
        ax.tick_params(labelsize=20)
    for ci, title in enumerate(['Styrofoam', 'Zotefoam HD30', 'Alternate Zotefoams']):
        axes[0, ci].set_title(title, fontsize=28)

    model_h = [mpl.lines.Line2D([0], [1], color='k' if m == 'rma' else 'crimson', ls=ls, lw=3, label=lab)
               for m, lab, ls in MODS]
    chi = {r[0]: r[1:] for r in rows}

    def h(bname, label, ls='-'):
        c = chi[bname]
        return mpl.lines.Line2D([0], [1], color=comb[bname]['color'], lw=3, ls=ls,
                                label=f'{label}  ({c[0]:.2f} / {c[1]:.2f})')

    leg = dict(framealpha=1, edgecolor='black', facecolor='white', fancybox=False)
    axes[1, 0].legend(handles=model_h + [h('Styrofoam_batch1', 'Styrofoam')], loc='upper left', fontsize=17,
                      title=r'$\chi^2$/dof: R+M+a / full Mie', title_fontsize=15, **leg)
    axes[1, 1].legend(handles=model_h + [h(f'HD30_batch{i}', f'HD30 batch {i}') for i in (1, 2, 3)], loc='upper left',
                      fontsize=17, title=r'$\chi^2$/dof: R+M+a / full Mie', title_fontsize=15, **leg)
    other = [('PP35_batch1', 'PP35', '-'), ('LD15_batch1', 'LD15 batch 1', '-'), ('LD15_batch2', 'LD15 batch 2', ':'),
             ('LD24_batch1', 'LD24 batch 1', '-'), ('LD24_batch2', 'LD24 batch 2', ':'), ('HD80_batch1', 'HD80', '-'),
             ('HD110_batch1', 'HD110', '-')]
    axes[1, 2].legend(handles=model_h + [h(*o) for o in other if o[0] in chi], loc='upper right', fontsize=14,
                      title=r'$\chi^2$/dof: R+M+a / full Mie', title_fontsize=13, **leg)
    fig.set_facecolor('white')
    plt.subplots_adjust(wspace=0.07, hspace=0.07)
    fig.savefig('foam_fig5_compare.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    for r in rows:
        print(f'{r[0]:18s} chi2/dof rma {r[1]:7.2f}   fullMie {r[2]:7.2f}')

    with plt.rc_context({'font.size': 13, 'axes.titlesize': 14, 'axes.labelsize': 14, 'xtick.labelsize': 12,
                         'ytick.labelsize': 12, 'legend.fontsize': 9}):
        decompose(comb, tab)


def decompose(comb, tab):
    """Why full Mie works: HD30_batch1 terms, and the efficiency curve shapes."""
    bname = 'HD30_batch1'
    b, n = comb[bname], float(comb[bname]['n_eff'])
    f = np.linspace(150, 2000, 800)
    pr, pm = best(bname, 'rma'), best(bname, 'fmr_fix_wc')
    tan_d, C, KR, KM, a = pr
    loss = 2e7 * np.pi * n * tan_d / fm.C_SI * f
    ray, mie = C * (f / KR)**4, C * fm.q_ada(f, KM)
    cpl = -2 * a * ray * mie / C
    loss_m = 2e7 * np.pi * n * pm[0] / fm.C_SI * f
    sca_m = pm[1] * fm.Q_fm(f, pm[2], pm[3], 1 / n, n, tab)

    fig, ax = plt.subplots(1, 2, figsize=(17, 6.5))
    A = ax[0]
    fa, ta, wa = b['f_comb'], b['t_comb'], b['width_comb']
    ok = ta > 0
    A.plot(fa[ok], -np.log(ta[ok]) / wa[ok], '.', color='0.7', ms=3, label='data, $-\\ln T / w$', zorder=-5)
    A.plot(f, ray, color='tab:blue', ls='--', lw=2, label='R+M+a: Rayleigh $C(f/K_R)^4$')
    A.plot(f, mie, color='tab:green', ls='--', lw=2, label='R+M+a: ADA Mie $C\\,q(f/K_M)$')
    A.plot(f, cpl, color='tab:purple', ls='--', lw=2, label='R+M+a: coupling $-2a\\,\\alpha_R\\alpha_M/C$')
    A.plot(f, loss, color='tab:orange', ls='--', lw=2, label='R+M+a: dielectric $\\alpha_\\delta$')
    A.plot(f, loss + ray + mie + cpl, 'k--', lw=2.5, label='R+M+a total')
    A.plot(f, sca_m, color='crimson', lw=2, label='Full Mie: scattering')
    A.plot(f, loss_m, color='crimson', ls=':', lw=2, label='Full Mie: dielectric $\\alpha_\\delta$')
    A.plot(f, loss_m + sca_m, color='crimson', lw=3, alpha=0.6, label='Full Mie total')
    A.set_xlim(150, 2000)
    A.set_ylim(-1.0, 2.0)
    A.axhline(0, color='k', lw=0.8)
    A.set_xlabel('Frequency [GHz]')
    A.set_ylabel(r'Attenuation [cm$^{-1}$]')
    A.set_title(f'{bname}: model terms (best fits)')
    A.legend(fontsize=9, loc='upper left', ncol=2)
    A.grid(True, alpha=0.3)

    B = ax[1]
    x = np.logspace(-2, 2.5, 600)
    for m, col in [(1 / 1.02, 'crimson'), (1 / 1.06, 'darkred')]:
        q = np.array([fm.q_interp(xx, m, tab) for xx in x])
        B.loglog(x, q, color=col, lw=2.5, label=f'full Mie, $m$ = 1/{1 / m:.2f}')
        rho = 2 * x * abs(m - 1)
        B.loglog(x, 2 - 4 * np.sin(rho) / rho + 4 * (1 - np.cos(rho)) / rho**2, color=col, ls='--', lw=1.5,
                 label=f'ADA ($\\rho = 2x|m-1|$), $m$ = 1/{1 / m:.2f}')
        ray_q = 8 / 3 * x**4 * ((m**2 - 1) / (m**2 + 2))**2
        B.loglog(x, ray_q, color=col, ls=':', lw=1.5, label=f'Rayleigh $\\propto x^4$, $m$ = 1/{1 / m:.2f}')
    B.set_ylim(1e-9, 5)
    B.set_xlabel(r'size parameter $x = \pi d\, n_{eff} f / c$')
    B.set_ylabel(r'scattering efficiency $Q_{sca}$')
    B.set_title('One sphere, low contrast: Rayleigh $x^4$ -> RGD $x^2$ -> saturation')
    B.legend(fontsize=9, loc='lower right')
    B.grid(True, which='both', alpha=0.25)
    fig.tight_layout()
    fig.savefig('foam_mie_terms.png', dpi=200)
    plt.close(fig)
    print(f'{bname}: rma tan_d {tan_d:.2e} C {C:.3g} K_R {KR:.0f} K_M {KM:.0f} a {a:.3f}; '
          f'fullMie tan_d {pm[0]:.2e} C {pm[1]:.3g} d {pm[2] * 1e3:.1f} um s {pm[3]:.2f}')
    print(f'  at 1 THz: rma loss {np.interp(1000, f, loss):.3f}, ray {np.interp(1000, f, ray):.3f}, mie {np.interp(1000, f, mie):.3f}, '
          f'coupling {np.interp(1000, f, cpl):.3f}; fullMie loss {np.interp(1000, f, loss_m):.3f}, sca {np.interp(1000, f, sca_m):.3f}')
    print(f'  at 2 THz: rma ray {np.interp(2000, f, ray):.3f}, coupling {np.interp(2000, f, cpl):.3f}; '
          f'fullMie loss {np.interp(2000, f, loss_m):.3f}, sca {np.interp(2000, f, sca_m):.3f}')


if __name__ == '__main__':
    main()
