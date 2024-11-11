import matplotlib.pyplot as plt
import numpy as np
import imageio

# Generate the data
x = np.linspace(0, 2*np.pi, 100)
y = np.abs(np.sin(x))

# Create the frames and save them as images
frames = []
for i in range(1, 21):  # Increase the number of frames (20 frames)
    xtmp = x[:i*5]  # Adjust how much data is plotted to have smoother transitions
    ytmp = y[:i*5]
    
    # Plot the current frame
    plt.plot(xtmp, ytmp)
    plt.xlim(x.min(), x.max())
    plt.ylim(0, 1)
    
    # Save the current frame
    filename = "tmp-frame%02d.png" % i
    plt.savefig(filename)
    plt.close()
    
    # Append the image to the frames list
    frames.append(imageio.imread(filename))

# Create a GIF from the frames with slower speed (duration=1 second between frames)
imageio.mimsave('animation.gif', frames, duration=1)

print("GIF created successfully!")

