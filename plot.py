import pandas as pd
import matplotlib.pyplot as plt

# 1. load the CSV
df = pd.read_csv("positions.csv")

# 2. make a fresh figure
plt.figure(figsize=(8, 6))

# 3. iterate over each body ID and plot its path
for body_id, g in df.groupby("body"):
    plt.plot(g["x"], g["y"], label=f"body {body_id}")

# 4. tidy up
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Planet trajectories")
plt.legend()           # one entry per body
plt.axis("equal")      # keep aspect-ratio 1:1

# 5. save or show
plt.savefig("trajectories.png", dpi=300)
# plt.show()           # uncomment if you prefer an on-screen window
