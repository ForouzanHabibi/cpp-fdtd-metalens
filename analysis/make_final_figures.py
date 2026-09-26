from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"

lib = pd.read_csv(DATA / "final_verified_8state_library.csv")
spec = pd.read_csv(DATA / "composite_spectrum.csv")

plt.figure(figsize=(8,5))
plt.plot(lib["state_deg"], lib["state_deg"], linestyle="--", label="Ideal")
plt.plot(lib["state_deg"], lib["phase_deg"], marker="o", label="Verified")
plt.xlabel("Target state (deg)")
plt.ylabel("Verified phase (deg)")
plt.title("Final verified 8-state phase library")
plt.legend()
plt.tight_layout()
plt.savefig(FIG / "final_phase_library.png", dpi=220)
plt.close()

plt.figure(figsize=(10,5))
plt.plot(spec["wavelength_nm"], spec["eta_3xFWHM"], marker="o")
plt.axhline(0.510788, linestyle="--", linewidth=1.0, label="80% threshold")
plt.axvline(1427.50, linestyle=":", linewidth=1.0, label="Lower edge")
plt.axvline(2390.51, linestyle=":", linewidth=1.0, label="First upper edge")
plt.axvline(2434.56, linestyle="-.", linewidth=1.0, label="Recovery")
plt.xlabel("Wavelength (nm)")
plt.ylabel("eta_3xFWHM")
plt.title("Composite spectral robustness")
plt.legend()
plt.tight_layout()
plt.savefig(FIG / "final_spectral_response.png", dpi=220)
plt.close()
