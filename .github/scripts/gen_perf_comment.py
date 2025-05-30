#!/usr/bin/env python3
import json, sys, math, textwrap
from scipy.stats import t

THRESH = 0.03
ALPHA = 0.05


def welch(a, b):
    m1, s1, n1 = a["mean"], a["stddev"], a["runs"]
    m2, s2, n2 = b["mean"], b["stddev"], b["runs"]
    if n1 == 0 or n2 == 0:
        return 1.0
    se = math.hypot(s1 / math.sqrt(n1), s2 / math.sqrt(n2))
    if se == 0:
        return 1.0
    t_stat = (m1 - m2) / se
    df = (s1**2 / n1 + s2**2 / n2) ** 2 / (
        (s1**2 / n1) ** 2 / (n1 - 1) + (s2**2 / n2) ** 2 / (n2 - 1)
    )
    return 2 * (1 - t.cdf(abs(t_stat), df))


def classify(b, p):
    if b is None or p is None:
        return "fail"
    delta = (p["mean"] - b["mean"]) / b["mean"]
    sig = welch(b, p) < ALPHA and abs(delta) >= THRESH
    if not sig:
        return "same"
    return "faster" if delta < 0 else "slower"


def fmt(val):
    return f"{val['mean']:.3f}±{val['stddev']:.3f}" if val else "❌"


def table(rows, extra_p=False):
    if not rows:
        return "_none_\n"
    hd = (
        "| map | main | pr | Δ % |"
        + (" p |" if extra_p else "")
        + "\n|---|---|---|---|"
        + ("---|\n" if extra_p else "\n")
    )
    body = []
    for m, b, p in rows:
        delta = (p["mean"] - b["mean"]) / b["mean"] * 100 if b and p else float("nan")
        ln = f"| `{m}` | {fmt(b)} | {fmt(p)} | {delta:+.1f} |"
        if extra_p:
            ln += f" {welch(b,p):.3f} |"
        body.append(ln)
    return hd + "\n".join(body) + "\n"


data = json.load(open(sys.argv[1]))
groups = {k: [] for k in ("faster", "slower", "same", "fail")}

for e in data:
    cat = classify(e.get("base"), e.get("pr"))
    groups[cat].append((e["map"], e.get("base"), e.get("pr")))

with open(sys.argv[2], "w") as f:
    f.write(
        textwrap.dedent(
            f"""
    ## Performance check (&alpha;={ALPHA:.2f}, ±{THRESH*100:.0f}%)

    ### ❗ Failed maps ({len(groups['fail'])})
    {table(groups['fail'])}

    <details><summary>Speed-ups ({len(groups['faster'])})</summary>

    {table(groups['faster'], extra_p=True)}
    </details>

    <details><summary>Slow-downs ({len(groups['slower'])})</summary>

    {table(groups['slower'], extra_p=True)}
    </details>

    <details><summary>No change ({len(groups['same'])})</summary>

    {table(groups['same'])}
    </details>
    """
        ).strip()
        + "\n"
    )
