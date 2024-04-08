import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np

#import values
filename="data/layer_bottom_B_x.txt"

with open(filename, 'r') as file:
    lines = [line.strip() for line in file if not line.startswith('#')]

# Initialize empty lists for each column
X = []
Z = []
B = []

# Iterate through each line and split the values
for line in lines:
    # Split the line into three values using space as the delimiter
    values = line.split()
    
    # Convert each value to the appropriate data type (e.g., float)
    val1, _,_,val2, val3 = map(float, values)
    
    # Append the values to the respective columns
    X.append(val1)
    Z.append(val3)
    B.append(val2)

df = pd.DataFrame({"x":X,"b":B,"z":Z})

# Creazione del grafico bidimensionale
# Creazione del grafico tridimensionale

# Creazione di una griglia per l'interpolazione
xi = np.linspace(min(df["x"]), max(df["x"]), 1000)
zi = np.linspace(min(df["z"]), max(df["z"]), 1000)
xi, zi = np.meshgrid(xi, zi)

# Interpolazione dei valori di B sulla griglia
Bi = griddata((df["x"], df["z"]), B, (xi, zi), method='linear')
print(Bi[999,999])


# Plot della superficie interpolata
fig, ax = plt.subplots(1,2)
ax[0] = fig.add_subplot(111, projection='3d')
ax[0].plot_trisurf(df["x"], df["z"], B, cmap='viridis', edgecolor='k', linewidth=0.1)
ax[0].set_xlabel('X')
ax[0].set_ylabel('Z')
ax[0].set_zlabel('B')

# Aggiunta di un secondo plot con la superficie interpolata
ax[1] = fig.add_subplot(122, projection='3d')
ax[1].plot_surface(xi, zi, Bi, alpha=0.5, cmap='viridis', rstride=100, cstride=100)

plt.show()
