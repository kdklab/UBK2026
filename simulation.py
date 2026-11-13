"""
UBK2026 simülasyonu: dijital ikiz kovaryatıyla ayarlama (Y1, Y2) ve
sanal kontrol havuzlaması (Y3, Y4), referans t-testine (Y0) karşı.

ADEMP özeti
  DGM: eşzamanlı deney  Y = theta*T + rho'*M + sqrt(1-rho'^2)*e
       tarihsel kontrol  Yh = delta + rho*Mh + sqrt(1-rho^2)*e   (Var = 1)
       rho' = max(rho - eta, 0)
  Estimand: theta (tedavi - kontrol ortalama farkı, kontrol SS birimi)
  Yöntemler: Y0 t-testi (n0/grup); Y1 ANCOVA (ceil(n0*(1-rho^2))/grup);
             Y2 ANCOVA (Y1 + 1/grup); Y3 naif havuzlama (tedavi n0,
             kontrol ceil(n0/2) + Nh tarihsel); Y4 test-sonra-havuzla
             (p >= 0.10 ise havuzla)
  Performans: red oranı (H0: tip I hata, H1: güç), yanlılık, %95 GA kapsamı
"""
import itertools
import math
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 0.05
POOL_ALPHA = 0.10
SEED = 20261005


def n_ttest(theta, power=0.80, alpha=ALPHA):
    """t-testiyle hedef güce ulaşan en küçük grup başı n (kesin, nct ile)."""
    n = 2
    while True:
        df = 2 * n - 2
        tcrit = stats.t.ppf(1 - alpha / 2, df)
        ncp = theta * math.sqrt(n / 2)
        pw = stats.nct.sf(tcrit, df, ncp) + stats.nct.cdf(-tcrit, df, ncp)
        if pw >= power:
            return n, pw
        n += 1


def ttest(a, b):
    """Satır bazında Student t-testi; a ve b: (R, n)."""
    na, nb = a.shape[1], b.shape[1]
    df = na + nb - 2
    sp2 = ((na - 1) * a.var(1, ddof=1) + (nb - 1) * b.var(1, ddof=1)) / df
    se = np.sqrt(sp2 * (1 / na + 1 / nb))
    d = a.mean(1) - b.mean(1)
    p = 2 * stats.t.sf(np.abs(d / se), df)
    q = stats.t.ppf(0.975, df)
    return d, p, d - q * se, d + q * se


def ancova(y1, m1, y0, m0):
    """Satır bazında ANCOVA: y ~ 1 + T + M; T katsayısı için t çıkarımı."""
    R, n = y1.shape
    y = np.concatenate([y1, y0], axis=1)
    m = np.concatenate([m1, m0], axis=1)
    T = np.concatenate([np.ones((R, n)), np.zeros((R, n))], axis=1)
    X = np.stack([np.ones_like(y), T, m], axis=2)
    XtX_inv = np.linalg.inv(np.einsum("rni,rnj->rij", X, X))
    beta = np.einsum("rij,rj->ri", XtX_inv, np.einsum("rni,rn->ri", X, y))
    resid = y - np.einsum("rni,ri->rn", X, beta)
    df = 2 * n - 3
    se = np.sqrt((resid ** 2).sum(1) / df * XtX_inv[:, 1, 1])
    d = beta[:, 1]
    p = 2 * stats.t.sf(np.abs(d / se), df)
    q = stats.t.ppf(0.975, df)
    return d, p, d - q * se, d + q * se


def gen_concurrent(rng, R, n, theta, rho_new):
    """Bir grubun (R, n) verisi ve ikiz tahmini."""
    m = rng.standard_normal((R, n))
    e = rng.standard_normal((R, n))
    y = theta + rho_new * m + math.sqrt(1 - rho_new ** 2) * e
    return y, m


def summarize(method, d, p, lo, hi, theta, animals):
    return dict(
        method=method,
        reject=float(np.mean(p < ALPHA)),
        bias=float(np.mean(d - theta)),
        coverage=float(np.mean((lo <= theta) & (theta <= hi))),
        animals=animals,
    )


def run(reps):
    rng = np.random.default_rng(SEED)
    out = []
    grid = itertools.product(
        [1.0, 1.5],                 # tasarım etki büyüklüğü
        [0.0, 0.3, 0.5, 0.7, 0.9],  # tarihsel korelasyon rho
        [0.0, 0.2, 0.4],            # korelasyon düşüşü eta
        [0.0, 0.25, 0.5, 1.0],      # ortalama kayması delta
        [20, 50],                   # tarihsel kontrol sayısı
        ["H0", "H1"],
    )
    for theta_d, rho, eta, delta, nh, truth in grid:
        theta = 0.0 if truth == "H0" else theta_d
        n0, _ = n_ttest(theta_d)
        n1 = max(2, math.ceil(n0 * (1 - rho ** 2)))
        n2 = n1 + 1
        nc = math.ceil(n0 / 2)
        rho_new = max(rho - eta, 0.0)
        cell = dict(theta_d=theta_d, n0=n0, rho=rho, eta=eta,
                    rho_new=rho_new, delta=delta, nh=nh, truth=truth)

        # Y0: referans t-testi (ikiz kullanılmaz)
        yt, _ = gen_concurrent(rng, reps, n0, theta, rho_new)
        yc, _ = gen_concurrent(rng, reps, n0, 0.0, rho_new)
        out.append({**cell, **summarize("Y0", *ttest(yt, yc), theta, 2 * n0),
                    "n_group": n0})

        # Y1 ve Y2: ikiz kovaryatıyla ANCOVA
        for code, n in (("Y1", n1), ("Y2", n2)):
            yt, mt = gen_concurrent(rng, reps, n, theta, rho_new)
            yc, mc = gen_concurrent(rng, reps, n, 0.0, rho_new)
            out.append({**cell, **summarize(code, *ancova(yt, mt, yc, mc),
                                            theta, 2 * n), "n_group": n})

        # Y3 ve Y4: sanal kontrol (tarihsel kontroller ~ N(delta, 1))
        yt, _ = gen_concurrent(rng, reps, n0, theta, rho_new)
        yc, _ = gen_concurrent(rng, reps, nc, 0.0, rho_new)
        yh = delta + rng.standard_normal((reps, nh))
        pooled = np.concatenate([yc, yh], axis=1)
        d3, p3, lo3, hi3 = ttest(yt, pooled)
        out.append({**cell, **summarize("Y3", d3, p3, lo3, hi3, theta,
                                        n0 + nc), "n_group": n0})
        _, p_pool, _, _ = ttest(yh, yc)
        use = p_pool >= POOL_ALPHA
        d_c, p_c, lo_c, hi_c = ttest(yt, yc)
        out.append({**cell, **summarize(
            "Y4", np.where(use, d3, d_c), np.where(use, p3, p_c),
            np.where(use, lo3, lo_c), np.where(use, hi3, hi_c),
            theta, n0 + nc), "n_group": n0, "pool_rate": float(use.mean())})
    return pd.DataFrame(out)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
    path = sys.argv[2] if len(sys.argv) > 2 else f"results_{reps}.csv"
    t0 = time.time()
    df = run(reps)
    df.to_csv(path, index=False)
    print(f"{len(df)} satır, {reps} tekrar, {time.time() - t0:.1f} sn -> {path}")
