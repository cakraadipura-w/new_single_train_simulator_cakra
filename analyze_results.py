import pandas as pd

df = pd.read_csv("experiment_results/IS02/E5_results.csv")
metrics = ["hv", "igd", "feasible", "E_best_kWh", "energy_saving_pct", "runtime_s"]
summary = df.groupby("config")[metrics].agg(["mean", "std"])
print("Summary Statistics:")
print(summary)
print("\n")

means = df.groupby("config")[metrics].mean()

def compute_delta(c1, c2):
    if c1 in means.index and c2 in means.index:
        return means.loc[c1] - means.loc[c2]
    return None

delta_c = compute_delta("C_LLM", "C")
delta_d = compute_delta("D_LLM", "D")

print("Deltas (C_LLM - C):")
print(delta_c)
print("\n")
print("Deltas (D_LLM - D):")
print(delta_d)
