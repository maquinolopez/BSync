import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def create_sync_timeseries_logo():
    fig, ax = plt.subplots(figsize=(6, 4))

    # Set a white background
    ax.set_facecolor("white")
    fig.patch.set_facecolor('white')

    # Generate a sample time series
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x)
    y2 = np.sin(x - 0.5)  # This is slightly out of phase with y1

    # Plot the time series
    ax.plot(x, y1, color='blue')
    ax.plot(x, y2, color='red', linestyle='--')

    # Annotate to show synchronization
    arrowprops = dict(facecolor='black', edgecolor='black', arrowstyle='->')
    ax.annotate('Synchronizing...', xy=(7, -0.25), xytext=(3, -0.75),
                arrowprops=arrowprops, fontsize=10, ha='center')

    # Add "BSynch" text at the top
    ax.text(5, 1.1, 'BSynch', color='black',
            fontsize=20, ha='center', fontweight='bold')

    # Customize the plot appearance
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Save the image
    plt.tight_layout()
    plt.savefig('logo.png', dpi=300, transparent=False)

    # Crop out whitespace using PIL
    img = Image.open('logo.png')
    img_cropped = img.crop(img.getbbox())
    img_cropped.save('logo.png')


if __name__ == '__main__':
    create_sync_timeseries_logo()
