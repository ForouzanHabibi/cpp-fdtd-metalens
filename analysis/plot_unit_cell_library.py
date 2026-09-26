
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PLOTS = Path("plots")
df = pd.read_csv(PLOTS / "final_8state_verification_step5l.csv")

for col in ["passes_T", "passes_phase", "passes_energy", "passes_all"]:
    df[col] = (
        df[col]
        .astype(str)
        .str.strip()
        .str.lower()
        .map({
            "true": True,
            "false": False,
            "1": True,
            "0": False
        })
        .fillna(False)
        .astype(bool)
    )

def circ_err(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)

# Temporal stability metrics.
for label in df["label"].unique():
    d = df[df["label"] == label].sort_values("checkpoint_fs")
    if len(d) >= 2:
        phase_shift = circ_err(d.iloc[-1]["phase_deg"], d.iloc[0]["phase_deg"])
        delta_T = d.iloc[-1]["T"] - d.iloc[0]["T"]
        delta_E = d.iloc[-1]["energy_error"] - d.iloc[0]["energy_error"]
        df.loc[d.index, "phase_shift_420_to_600_deg"] = phase_shift
        df.loc[d.index, "delta_T_420_to_600"] = delta_T
        df.loc[d.index, "delta_energy_420_to_600"] = delta_E

df.to_csv(PLOTS / "final_8state_verification_step5l_enriched.csv", index=False)

final = df[df["checkpoint_fs"] == df["checkpoint_fs"].max()].copy()
final = final.sort_values("target_phase_deg").reset_index(drop=True)
final.to_csv(PLOTS / "final_verified_8state_library_step5l.csv", index=False)

# Phase curve.
plt.figure(figsize=(9, 5))
plt.plot(
    final["target_phase_deg"],
    final["target_phase_deg"],
    linestyle="--",
    label="Ideal"
)
plt.plot(
    final["target_phase_deg"],
    final["phase_deg"],
    marker="o",
    label="Verified 600 fs"
)
plt.xlabel("Target phase (deg)")
plt.ylabel("Verified phase (deg)")
plt.title("Final 8-state phase library")
plt.legend()
plt.tight_layout()
plt.savefig(PLOTS / "final_phase_library_step5l.png", dpi=250)
plt.close()

# Transmission.
plt.figure(figsize=(9, 5))
plt.bar(
    [f'{int(x)}°' for x in final["target_phase_deg"]],
    final["T"]
)
plt.axhline(0.80, linestyle="--", linewidth=1.0)
plt.xlabel("Target phase state")
plt.ylabel("Transmission T")
plt.title("Final 8-state transmission")
plt.tight_layout()
plt.savefig(PLOTS / "final_transmission_step5l.png", dpi=250)
plt.close()

# Energy error.
plt.figure(figsize=(9, 5))
plt.bar(
    [f'{int(x)}°' for x in final["target_phase_deg"]],
    final["energy_error"]
)
plt.axhline(0.02, linestyle="--", linewidth=1.0)
plt.xlabel("Target phase state")
plt.ylabel("|R + T - 1|")
plt.title("Final 8-state power-conservation diagnostic")
plt.tight_layout()
plt.savefig(PLOTS / "final_energy_error_step5l.png", dpi=250)
plt.close()

# Temporal phase stability.
plt.figure(figsize=(9, 5))
plt.bar(
    [f'{int(x)}°' for x in final["target_phase_deg"]],
    final["phase_shift_420_to_600_deg"]
)
plt.xlabel("Target phase state")
plt.ylabel("420→600 fs phase shift (deg)")
plt.title("Temporal stability of final library")
plt.tight_layout()
plt.savefig(PLOTS / "final_temporal_stability_step5l.png", dpi=250)
plt.close()

passing = final[final["passes_all"]].copy()

lines = [
    "Step 5L: uniform final verification of the complete 8-state library",
    "==================================================================",
    "",
    "Method",
    "------",
    "All eight phase states are now evaluated under exactly the same final",
    "verification method: fine grid, periodic x boundary, CPML in y, spatially",
    "integrated Poynting flux, far transmission monitor, and a 600-fs run.",
    "",
    "The former 45° geometry (700 x 1600 nm) has been replaced by the Step 5K",
    "verified J2 candidate (450 x 1833.33 nm).",
    "",
    "600-fs state-by-state results",
    "-----------------------------",
]

for _, r in final.iterrows():
    lines.append(
        f'{int(r["target_phase_deg"])}°: '
        f'geometry={r["width_nm"]:.2f} x {r["height_nm"]:.2f} nm, '
        f'phase={r["phase_deg"]:.3f}°, '
        f'phase error={r["phase_error_deg"]:.3f}°, '
        f'T={r["T"]:.6f}, '
        f'R={r["R"]:.6f}, '
        f'energy error={r["energy_error"]:.6e}, '
        f'420->600 phase shift={r["phase_shift_420_to_600_deg"]:.3f}°, '
        f'passes all={bool(r["passes_all"])}'
    )

lines += [
    "",
    "Library summary",
    "---------------",
    f"states passing all project criteria = {len(passing)} / {len(final)}",
    f"minimum T = {final['T'].min():.6f}",
    f"maximum phase error = {final['phase_error_deg'].max():.3f}°",
    f"mean phase error = {final['phase_error_deg'].mean():.3f}°",
    f"maximum energy error = {final['energy_error'].max():.6e}",
    f"maximum 420->600 phase shift = {final['phase_shift_420_to_600_deg'].max():.3f}°",
    "",
    "Decision",
    "--------",
]

if len(passing) == len(final):
    lines += [
        "All eight states pass the current project-level verification criteria.",
        "This library is ready to be used as the discrete phase lookup table for",
        "the next step: metalens phase-profile assignment.",
    ]
else:
    failed = final[~final["passes_all"]]
    fail_labels = ", ".join(
        f'{int(x)}°' for x in failed["target_phase_deg"]
    )
    lines += [
        f"The following states still fail at least one final criterion: {fail_labels}.",
        "Those states should be repaired before using the library for a metalens.",
    ]

summary = "\n".join(lines)
(PLOTS / "validation_step5l.txt").write_text(summary, encoding="utf-8")

print(summary)
print("\nOutputs:")
print(" - plots/final_8state_verification_step5l.csv")
print(" - plots/final_8state_verification_step5l_enriched.csv")
print(" - plots/final_verified_8state_library_step5l.csv")
print(" - plots/final_phase_library_step5l.png")
print(" - plots/final_transmission_step5l.png")
print(" - plots/final_energy_error_step5l.png")
print(" - plots/final_temporal_stability_step5l.png")
print(" - plots/validation_step5l.txt")
