
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PLOTS = Path("plots")

df = pd.read_csv(PLOTS / "dense_upper_notch_step6l.csv")
df = df.sort_values("wavelength_nm").reset_index(drop=True)

old = pd.read_csv("data/step6k_overlap_reference.csv")
band = pd.read_csv("data/band_reference.csv").iloc[0]

eta1550 = float(band["eta1550"])
threshold = float(band["threshold80"])
lower_edge = float(band["lower_edge_nm"])

df["retention_vs_1550"] = df["eta_3xFWHM"] / eta1550

# Validate new 2400-centered source against Step 6K at 2350, 2400, 2450.
overlap = df[df["wavelength_nm"].isin(old["wavelength_nm"])].merge(
    old,
    on="wavelength_nm",
    suffixes=("_step6l", "_step6k")
)

overlap["eta_relative_difference_percent"] = (
    100.0
    * (overlap["eta_3xFWHM"] - overlap["eta3"]).abs()
    / overlap["eta3"].abs()
)

overlap["focus_difference_um"] = (
    overlap["interpolated_focus_y_um"] - overlap["focus_y_um"]
).abs()

overlap.to_csv(PLOTS / "step6k_step6l_overlap_check.csv", index=False)

# Detect all threshold crossings.
crossings = []

for i in range(len(df) - 1):
    a = df.iloc[i]
    b = df.iloc[i + 1]

    ya = float(a["eta_3xFWHM"])
    yb = float(b["eta_3xFWHM"])
    xa = float(a["wavelength_nm"])
    xb = float(b["wavelength_nm"])

    side_a = ya - threshold
    side_b = yb - threshold

    if side_a == 0:
        crossings.append(("exact", xa, xa, xa))
    elif side_a * side_b < 0:
        xcross = xa + (threshold - ya) * (xb - xa) / (yb - ya)

        direction = "downward" if ya > threshold and yb < threshold else "upward"
        crossings.append((direction, xa, xb, xcross))

downward = [c for c in crossings if c[0] == "downward"]
upward = [c for c in crossings if c[0] == "upward"]

first_down_edge = downward[0][3] if downward else None
first_up_edge = upward[0][3] if upward else None

continuous_bw = (
    first_down_edge - lower_edge
    if first_down_edge is not None
    else None
)

notch_width = (
    first_up_edge - first_down_edge
    if first_down_edge is not None and first_up_edge is not None
    else None
)

# Plots.
plt.figure(figsize=(9, 5))
plt.plot(
    df["wavelength_nm"],
    df["eta_3xFWHM"],
    marker="o"
)
plt.axhline(
    threshold,
    linestyle="--",
    linewidth=1.0,
    label="80% of 1550 efficiency"
)
plt.xlabel("Wavelength (nm)")
plt.ylabel("eta_3xFWHM")
plt.title("Step 6L: dense upper-threshold structure")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS / "dense_upper_efficiency_step6l.png", dpi=250)
plt.close()

plt.figure(figsize=(9, 5))
plt.plot(
    df["wavelength_nm"],
    df["retention_vs_1550"],
    marker="o"
)
plt.axhline(0.80, linestyle="--", linewidth=1.0)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Efficiency retention vs 1550")
plt.title("80%-retention crossings near the upper edge")
plt.tight_layout()
plt.savefig(PLOTS / "dense_upper_retention_step6l.png", dpi=250)
plt.close()

plt.figure(figsize=(9, 5))
plt.plot(
    df["wavelength_nm"],
    df["interpolated_focus_y_um"],
    marker="o"
)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Focus y (um)")
plt.title("Dense upper-region focal shift")
plt.tight_layout()
plt.savefig(PLOTS / "dense_upper_focus_step6l.png", dpi=250)
plt.close()

plt.figure(figsize=(9, 5))
plt.plot(
    df["wavelength_nm"],
    df["source_power_relative_to_2400"],
    marker="o"
)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Vacuum source power / 2400-nm value")
plt.title("Local source support")
plt.tight_layout()
plt.savefig(PLOTS / "dense_upper_source_support_step6l.png", dpi=250)
plt.close()

max_eta_overlap = float(overlap["eta_relative_difference_percent"].max())
max_focus_overlap = float(overlap["focus_difference_um"].max())
min_support = float(df["source_power_relative_to_2400"].min())

lines = [
    "Step 6L: dense upper-threshold notch localization",
    "=================================================",
    "",
    "Why this step is different from a simple upper-edge refinement",
    "----------------------------------------------------------------",
    "Step 6K fell below the 80%-retention threshold at 2400 nm but returned",
    "above threshold at 2450 and 2500 nm. Therefore the spectrum is not",
    "monotonic: the 2350-2450 nm region contains a local efficiency dip/notch.",
    "",
    "For the continuous useful band connected to the 1550-nm design point,",
    "the relevant upper edge is the FIRST downward 80% crossing. A later",
    "recovery is reported separately and must not be merged into that continuous",
    "bandwidth.",
    "",
    "Dense results",
    "-------------",
]

for _, r in df.iterrows():
    lines.append(
        f'{r["wavelength_nm"]:.0f} nm: '
        f'focus y={r["interpolated_focus_y_um"]:.3f} um, '
        f'FWHM={r["fwhm_um"]:.3f} um, '
        f'eta_3xFWHM={r["eta_3xFWHM"]:.5f}, '
        f'retention={100*r["retention_vs_1550"]:.2f}%, '
        f'post-lens T={r["postlens_transmission"]:.4f}, '
        f'capture={r["capture_3xFWHM"]:.4f}, '
        f'source support={100*r["source_power_relative_to_2400"]:.1f}%'
    )

lines += [
    "",
    "Overlap consistency with Step 6K",
    "--------------------------------",
    f'max relative eta difference at 2350/2400/2450 nm = {max_eta_overlap:.2f}%',
    f'max focus-position difference at overlap points = {max_focus_overlap:.3f} um',
    f'minimum local source support = {100*min_support:.1f}% of 2400-nm value',
    "",
    "Threshold crossings",
    "-------------------",
    f'80% efficiency threshold = {threshold:.6f}',
    f'lower continuous-band edge from Step 6J = {lower_edge:.2f} nm',
]

if crossings:
    for c in crossings:
        if c[0] == "exact":
            lines.append(f'exact threshold sample at {c[3]:.2f} nm')
        else:
            lines.append(
                f'{c[0]} crossing bracket = {c[1]:.0f}-{c[2]:.0f} nm; '
                f'linear edge estimate = {c[3]:.2f} nm'
            )
else:
    lines.append("no threshold crossing detected in 2350-2450 nm")

lines += [
    "",
    "Continuous-band interpretation",
    "------------------------------",
]

if first_down_edge is not None:
    lines += [
        f'first upper downward crossing = {first_down_edge:.2f} nm',
        f'continuous 80%-retention band connected to 1550 nm = '
        f'{lower_edge:.2f} to {first_down_edge:.2f} nm',
        f'continuous 80%-retention bandwidth = {continuous_bw:.2f} nm',
    ]
else:
    lines += [
        "The first downward upper crossing was not resolved.",
        "Do not quote the continuous upper band edge yet."
    ]

if notch_width is not None:
    lines += [
        f'first recovery crossing = {first_up_edge:.2f} nm',
        f'below-threshold notch width between first down/up crossings = {notch_width:.2f} nm',
    ]

lines += [
    "",
    "Decision",
    "--------",
]

if first_down_edge is not None and max_eta_overlap <= 3.0:
    lines += [
        "The first upper 80%-retention edge is now densely localized with a",
        "source-consistency check. The continuous useful band containing 1550 nm",
        "can be reported using the lower Step 6J edge and the first downward",
        "Step 6L crossing.",
        "",
        "Because the response recovers above threshold at longer wavelength,",
        "describe that recovery as a separate spectral pass region rather than",
        "extending the continuous 1550-centered bandwidth through the notch."
    ]
else:
    lines += [
        "The upper continuous-band edge is not yet sufficiently validated.",
        "Do not finalize the bandwidth."
    ]

lines += [
    "",
    "Scope",
    "-----",
    "This remains the idealized nondispersive n=2, 2D TMz model. The bandwidth",
    "is a numerical property of this model, not a real-material device bandwidth."
]

report = "\n".join(lines)

df.to_csv(PLOTS / "dense_upper_notch_enriched_step6l.csv", index=False)
(PLOTS / "validation_step6l.txt").write_text(report, encoding="utf-8")

print(report)
print("\nOutputs:")
print(" - plots/dense_upper_notch_step6l.csv")
print(" - plots/dense_upper_notch_enriched_step6l.csv")
print(" - plots/step6k_step6l_overlap_check.csv")
print(" - plots/dense_upper_efficiency_step6l.png")
print(" - plots/dense_upper_retention_step6l.png")
print(" - plots/dense_upper_focus_step6l.png")
print(" - plots/dense_upper_source_support_step6l.png")
print(" - plots/validation_step6l.txt")
