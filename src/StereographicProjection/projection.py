from matplotlib import projections
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

theta = np.linspace(0, 2 * np.pi, 100)

nums_lat = 50

z = 0.0
steps_z = 1.0 / nums_lat

for i in range(nums_lat):
    c = np.sqrt(1 - z * z)/(1 + z)
    X = c * np.cos(theta)
    Y = c * np.sin(theta)
    plt.plot(X, Y)
    z += steps_z

phi = np.linspace(0, np.pi, 100)

nums_long = 50

t = 0.0
steps_t = np.pi / nums_long

for i in range(nums_long):
    X = (np.cos(t) * np.cos(phi))/(1 + np.sin(phi))
    Y = (np.sin(t) * np.cos(phi))/(1 + np.sin(phi))
    plt.plot(X, Y)
    t += steps_t

plt.show()