import sys, os, json, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from vvlib import *
import c07_easo as c
R, G = c.R, c.GAM
out = []; P = out.append
def soulie(V, th, S, Rg=None):  # code formula, V = bond volume, R2 = R, sqrt(ri rj) = R
    a = -1.1*V**-0.53; b = (-0.148*np.log(V)-0.96)*th*th-0.0082*np.log(V)+0.48; cc = 0.0018*np.log(V)+0.078
    return np.pi*G*R*(np.exp(a*S/R+b)+cc)
def willett(V, th, S):
    Sh = S*np.sqrt(R/(V*R**3)); return 2*np.pi*R*G*np.cos(th)/(1+1.05*Sh+2.5*Sh*Sh)
P("### Pair bridge, equal spheres R = 1 mm, gamma = 0.072 N/m, eta = 0 (Hooke, e = 1)")
P("| V/R^3 | theta | F(0+)/(2 pi R g) | Willett F(0)/(2 pi R g) | max rel diff to code (Soulie) formula | S_c/R measured | Lian (1+theta/2)V^(1/3)/R | F at rupture/(2 pi R g) | F/Willett at S=0.5 S_c | release-baseline max \\|dF\\| [N] |")
P("|---|---|---|---|---|---|---|---|---|---|")
plotdat = {}
for vb in [1e-3, 1e-2]:
    for th in [0.0, 20.0, 40.0]:
        A = load(os.path.join(CASES, "c07_easo", f"pair_release_V{vb}_th{th}", "fs.txt"))
        B = load(os.path.join(CASES, "c07_easo", f"pair_baseline_V{vb}_th{th}", "fs.txt"))
        S = A[:, 2]-A[:, 1]-2*R; F = A[:, 3]
        t = np.radians(th)
        sep = S > 0
        br = sep & (np.abs(F) > 0)
        Sc = S[br].max() if br.any() else np.nan
        i_first = np.where(sep)[0][0]
        Fs = soulie(vb, t, S[br]); err = np.max(np.abs(F[br]-Fs)/Fs)
        lian = (1+t/2)*vb**(1/3)
        k = np.argmin(np.abs(S-0.5*Sc))
        P(f"| {vb:g} | {th:g} | {F[i_first]/(2*np.pi*R*G):.4f} | {np.cos(t):.4f} | {err:.1e} | {Sc/R:.4f} | {lian:.4f} | {F[br][-1]/(2*np.pi*R*G):.4f} | {F[k]/willett(vb, t, S[k]):.3f} | {np.abs(F-B[:, 3]).max():.1e} |")
        plotdat[(vb, th)] = (S[sep], F[sep])
P("\n### Wall bridge (sphere-plane, same gamma, content chosen so that V_bond/R^3 = V if bond = 0.5 V_particle_liquid)")
for vb in [1e-3, 1e-2]:
    A = load(os.path.join(CASES, "c07_easo", f"wall_cap_V{vb}", "fs.txt"))
    S = A[:, 1]-R; F = -A[:, 2]
    br = (S > 0) & (np.abs(F) > 0); Sc = S[br].max()
    i0 = np.where(S > 0)[0][0]
    P(f"V={vb:g}: F(0+) = {F[i0]:.4e} N = {F[i0]/(4*np.pi*R*G):.3f} x sphere-plane theory 4 pi R gamma cos(theta) ; S_c/R = {Sc/R:.4f} (Lian with V: {vb**(1/3):.4f}); implied V_bond/R^3 from S_c = {(Sc/R)**3:.3e}")
P("\n### Viscous force, sphere-plane withdrawing at v = 1e-3 m/s, eta = 1e-3 Pa s (C-19)")
A = load(os.path.join(CASES, "c07_easo", "wall_visc", "fs.txt")); B = load(os.path.join(CASES, "c07_easo", "wall_visc_ref_nocap", "fs.txt"))
S = A[:, 1]-R; Fv = -(A[:, 2]-B[:, 2])  # viscous part = total - capillary-only (same kinematics)
eta, v = 1e-3, 1e-3
for Sq in [2e-6, 5e-6, 1e-5, 2e-5]:
    k = np.argmin(np.abs(S-Sq))
    th = 6*np.pi*eta*R*R*v/S[k]
    P(f"S={S[k]*1e6:.2f} um: F_visc code = {Fv[k]:.4e} N; sphere-plane Reynolds 6 pi eta R^2 v/S = {th:.4e} N; ratio = {Fv[k]/th:.4f}")
A = load(os.path.join(CASES, "c07_easo", "pair_visc", "fs.txt")); B = load(os.path.join(CASES, "c07_easo", "pair_release_V0.01_th0.0", "fs.txt"))
S = A[:, 2]-A[:, 1]-2*R; Fv = A[:, 3]-B[:, 3]
for Sq in [5e-6, 1e-5, 2e-5]:
    k = np.argmin(np.abs(S-Sq)); th = 6*np.pi*eta*(R/2)**2*c.VPULL/S[k]
    P(f"pair S={S[k]*1e6:.2f} um: F_visc code = {Fv[k]:.4e} N; sphere-sphere Reynolds 6 pi eta R*^2 v/S = {th:.4e} N; ratio = {Fv[k]/th:.4f}")
P("\n### Same, minSeparationDistanceRatio = 1e-3 (lubrication branch active)")
A = load(os.path.join(CASES, "c07_easo", "wall_visc_min0.001", "fs.txt")); B = load(os.path.join(CASES, "c07_easo", "wall_visc_ref_nocap_min0.001", "fs.txt"))
S = A[:, 1]-R; Fv = -(A[:, 2]-B[:, 2])
for Sq in [2e-6, 5e-6, 1e-5, 2e-5]:
    k = np.argmin(np.abs(S-Sq)); th = 6*np.pi*eta*R*R*v/S[k]
    P(f"wall S={S[k]*1e6:.2f} um: F_visc code = {Fv[k]:.4e} N; 6 pi eta R^2 v/S = {th:.4e} N; ratio = {Fv[k]/th:.4f}")
A = load(os.path.join(CASES, "c07_easo", "pair_visc_min0.001", "fs.txt")); B = load(os.path.join(CASES, "c07_easo", "pair_nocap_min0.001", "fs.txt"))
S = A[:, 2]-A[:, 1]-2*R; Fv = A[:, 3]-B[:, 3]
for Sq in [5e-6, 1e-5, 2e-5]:
    k = np.argmin(np.abs(S-Sq)); th = 6*np.pi*eta*(R/2)**2*c.VPULL/S[k]
    P(f"pair S={S[k]*1e6:.2f} um: F_visc code = {Fv[k]:.4e} N; 6 pi eta R*^2 v/S = {th:.4e} N; ratio = {Fv[k]/th:.4f}")
txt = "\n".join(out); print(txt); open(os.path.join(LOGS, "c07_easo_table.md"), "w").write(txt)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for i, vb in enumerate([1e-3, 1e-2]):
    for th, col in [(0.0, "C0"), (20.0, "C1"), (40.0, "C2")]:
        S, F = plotdat[(vb, th)]; t = np.radians(th)
        ax[i].plot(S/R, F/(2*np.pi*R*G), "-", color=col, label=f"LIGGGHTS theta={th:g}")
        ss = np.linspace(0, S.max(), 200); ax[i].plot(ss/R, willett(vb, t, ss)/(2*np.pi*R*G), ":", color=col, label=f"Willett 2000 theta={th:g}")
        ax[i].axvline((1+t/2)*vb**(1/3), color=col, lw=0.6, ls="--")
    ax[i].set_title(f"V_bond/R^3 = {vb:g} (dashed: Lian rupture)"); ax[i].set_xlabel("S/R"); ax[i].set_ylabel("F/(2 pi R gamma)"); ax[i].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "c07_easo.png"), dpi=120)
