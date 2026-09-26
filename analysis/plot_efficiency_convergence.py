from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

P = Path("plots")
a = pd.read_csv(P / "subpixel_summary_33p333nm_step6g.csv").iloc[0]
b = pd.read_csv(P / "subpixel_summary_25nm_step6g.csv").iloc[0]
binary = pd.read_csv("data/binary_reference_step6e.csv").sort_values("grid_nm", ascending=False).reset_index(drop=True)

def rel(x,y):
    return abs(y-x)/max(abs(y),1e-30)

metrics = [
    ("postlens_transmission", a["postlens_transmission"], b["postlens_transmission"]),
    ("focus_plane_forward_fraction", a["focus_plane_forward_fraction"], b["focus_plane_forward_fraction"]),
    ("eta_1xFWHM", a["eta_1xFWHM"], b["eta_1xFWHM"]),
    ("eta_2xFWHM", a["eta_2xFWHM"], b["eta_2xFWHM"]),
    ("eta_3xFWHM", a["eta_3xFWHM"], b["eta_3xFWHM"]),
    ("capture_3xFWHM", a["capture_3xFWHM"], b["capture_3xFWHM"]),
]
rows=[]
for name,x,y in metrics:
    rows.append({"metric":name,"value_33p333nm":x,"value_25nm":y,"absolute_change":abs(y-x),"relative_change_percent":100*rel(x,y)})
cmp=pd.DataFrame(rows)
cmp.to_csv(P/"subpixel_grid_changes_step6g.csv",index=False)

eta_change=float(cmp[cmp.metric=="eta_3xFWHM"].relative_change_percent.iloc[0])
post_change=float(cmp[cmp.metric=="postlens_transmission"].relative_change_percent.iloc[0])
bin_eta=100*rel(binary.iloc[0].eta_3xFWHM,binary.iloc[1].eta_3xFWHM)
bin_post=100*rel(binary.iloc[0].postlens_transmission,binary.iloc[1].postlens_transmission)
passed=(eta_change<=5.0 and post_change<=5.0)

# plots
plt.figure(figsize=(8,5))
plt.plot([33.333,25.0],[a["eta_3xFWHM"],b["eta_3xFWHM"]],marker="o",label="Subpixel")
plt.plot(binary.grid_nm,binary.eta_3xFWHM,marker="o",label="Binary staircase")
plt.gca().invert_xaxis(); plt.xlabel("Grid size (nm), finer →"); plt.ylabel("eta_3xFWHM"); plt.title("Efficiency convergence: binary vs subpixel"); plt.legend(); plt.tight_layout(); plt.savefig(P/"subpixel_efficiency_convergence_step6g.png",dpi=250); plt.close()

plt.figure(figsize=(9,5))
plt.bar(cmp.metric,cmp.relative_change_percent)
plt.axhline(5.0,linestyle="--",linewidth=1.0,label="5% practical threshold")
plt.ylabel("33.333 → 25 nm relative change (%)"); plt.title("Subpixel power-metric grid sensitivity"); plt.xticks(rotation=25,ha="right"); plt.legend(); plt.tight_layout(); plt.savefig(P/"subpixel_metric_changes_step6g.png",dpi=250); plt.close()

lines=[
"Step 6G: subpixel-averaged efficiency convergence",
"================================================",
"",
"Correction to the Step 6F interpretation",
"-----------------------------------------",
"The actual local eta_3xFWHM changed by only 0.36% when the x-domain was",
"widened from 30 to 40 um. The larger change in total focus-plane power",
"was integrated over different physical transverse widths in the two domains,",
"so it is not an apples-to-apples test of the focusing-efficiency metric.",
"Thus side-boundary clipping is unlikely to explain the ~6.7% binary-grid",
"efficiency difference.",
"",
"Step 6G uses exact rectangle/cell overlap fractions at pillar boundaries and",
"an area-weighted effective epsilon for the TMz Ez update. Both grids use the",
"same physical focus target (25.65 um) and the same 3xFWHM physical window",
"based on FWHM = 1.3752 um.",
"",
"Subpixel results",
"----------------",
f'33.333 nm: post-lens transmission={a["postlens_transmission"]:.6f}, focus-plane fraction={a["focus_plane_forward_fraction"]:.6f}, eta_3xFWHM={a["eta_3xFWHM"]:.6f}, capture={a["capture_3xFWHM"]:.6f}',
f'25 nm: post-lens transmission={b["postlens_transmission"]:.6f}, focus-plane fraction={b["focus_plane_forward_fraction"]:.6f}, eta_3xFWHM={b["eta_3xFWHM"]:.6f}, capture={b["capture_3xFWHM"]:.6f}',
"",
"33.333 -> 25 nm changes",
"-----------------------",
]
for _,r in cmp.iterrows():
    lines.append(f'{r.metric}: absolute change={r.absolute_change:.6f}, relative change={r.relative_change_percent:.2f}%')
lines += [
"",
"Binary-vs-subpixel sensitivity",
"------------------------------",
f'previous binary eta_3xFWHM grid change = {bin_eta:.2f}%',
f'subpixel eta_3xFWHM grid change = {eta_change:.2f}%',
f'previous binary post-lens transmission grid change = {bin_post:.2f}%',
f'subpixel post-lens transmission grid change = {post_change:.2f}%',
"",
"Decision",
"--------",
]
if passed:
    lines += [
        "The primary subpixel efficiency and post-lens transmission both change",
        "by <=5% between 33.333 and 25 nm.",
        "",
        "Use the 25-nm subpixel value as the preferred current 2D efficiency",
        "estimate. The previous binary-staircase value remains a useful reference.",
        "The next device-level study can move to spectral robustness around 1550 nm."
    ]
else:
    lines += [
        "At least one primary subpixel power metric still changes by more than 5%.",
        "Do not start the wavelength sweep yet. The next diagnostic should isolate",
        "source/monitor normalization and temporal-window sensitivity."
    ]
lines += [
"",
"Scope",
"-----",
"This remains a 2D TMz per-unit-out-of-plane-length model. Subpixel averaging",
"improves rasterization of the same idealized rectangles; it does not add",
"substrate dispersion, loss, fabrication roughness, or 3D coupling."
]
report="\n".join(lines)
(P/"validation_step6g.txt").write_text(report,encoding="utf-8")
print(report)
