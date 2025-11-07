import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np

# Get the colormap
cmap = plt.cm.twilight_shifted

# Sample 10 colors (for example)
colors = cmap(np.linspace(0, 1, 10))
hex_colors = [to_hex(cmap(i)) for i in np.linspace(0, 1, 10)]


print('Printing color codes:\n')
print(colors)

print('\nPrint hex color codes:\n')
print(hex_colors)
