import pandas as pd
import matplotlib.pyplot as plt

df_in = pd.read_csv(r"C:\Users\6anna\CLionProjects\PlasmaSimulator1D\cmake-build-debug\spectrum_in.csv")
df_out = pd.read_csv(r"C:\Users\6anna\CLionProjects\PlasmaSimulator1D\cmake-build-debug\spectrum_out.csv")


plt.figure(figsize=(8, 5))
plt.plot(df_in["f_over_fL"], df_in["absS2"], lw=2, label="input")
plt.plot(df_out["f_over_fL"], df_out["absS2"], lw=2, label="output")
plt.xlabel(r"$f/f_L$")
plt.ylabel(r"$|S(f)|^2$")
plt.title(r"Flux spectrum")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("flux_spectra_compare.png", dpi=200)
plt.show()
