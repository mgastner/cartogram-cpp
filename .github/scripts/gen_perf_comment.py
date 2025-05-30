#!/usr/bin/env python3
import json, sys, math
from scipy.stats import t
import html, re

THRESH = 0.05
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


def classify(base, pr):
    if base is None or pr is None:
        return "fail"
    delta = (pr["mean"] - base["mean"]) / base["mean"]
    sig = welch(base, pr) < ALPHA and abs(delta) >= THRESH
    if not sig:
        return "same"
    return "faster" if delta < 0 else "slower"


def fmt(val):
    return f"{val['mean']:.3f}±{val['stddev']:.3f}" if val else "❌"


def wrap(name: str) -> str:
    ZWSP = "\u200b"
    escaped = html.escape(name)
    return escaped.replace("/", "/" + ZWSP).replace("_", "_" + ZWSP)


def fmt_row(m, b, p, show_p):
    if b and p:
        delta = (p["mean"] - b["mean"]) / b["mean"] * 100
        delta_str = f"{delta:+.1f}"
        p_val = f"{welch(b, p):.3f}"
    else:
        delta_str, p_val = "—", ""
    wrapped = wrap(m)
    cell = f'<code title="{html.escape(m)}">{wrapped}</code>'
    row = f"| {cell} | {fmt(b)} | {fmt(p)} | {delta_str} |"
    if show_p:
        row += f" {p_val} |"
    return row


def table(rows, show_p=False):
    if not rows:
        return "_none_\n"

    header = (
        "| map | main | pr | Δ % |"
        + (" p |" if show_p else "")
        + "\n|---|---|---|---|"
        + ("---|\n" if show_p else "\n")
    )

    body = [fmt_row(m, b, p, show_p) for m, b, p in rows]
    return header + "\n".join(body) + "\n"


data = json.load(open(sys.argv[1]))
groups = {k: [] for k in ("faster", "slower", "same", "fail")}

for entry in data:
    cat = classify(entry.get("base"), entry.get("pr"))
    groups[cat].append((entry["map"], entry.get("base"), entry.get("pr")))

report = f"""\
<h2>🚦 Performance Comparison <sup>(α={ALPHA:.2f}, ±{THRESH*100:.0f}%)</sup></h2>

<table>
  <tr>
    <td><strong>🗺️ Total maps</strong></td><td align="right">{len(data)}</td>
    <td><strong>❌ Failed</strong></td><td align="right">{len(groups['fail'])}</td>
    <td><strong>🚀 Speed-ups</strong></td><td align="right">{len(groups['faster'])}</td>
    <td><strong>🐢 Slow-downs</strong></td><td align="right">{len(groups['slower'])}</td>
    <td><strong>⚖️ No change</strong></td><td align="right">{len(groups['same'])}</td>
  </tr>
</table>

---

### ❌ Failures
{table(groups['fail']) if groups['fail'] else '_None 😎_'}  

<details>
  <summary>🚀 Speed-ups ({len(groups['faster'])})</summary>

{table(groups['faster'], show_p=True)}
</details>

<details>
  <summary>🐢 Slow-downs ({len(groups['slower'])})</summary>

{table(groups['slower'], show_p=True)}
</details>

<details>
  <summary>⚖️ No significant change ({len(groups['same'])})</summary>

{table(groups['same'], show_p=True)}
</details>
"""

with open(sys.argv[2], "w") as f:
    f.write(report.strip() + "\n")
