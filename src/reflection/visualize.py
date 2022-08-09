import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt


"""
t = np.linspace(0, np.pi, 100)
s = np.linspace(0, np.pi, 100)
t, s = np.meshgrid(t, s)
x = np.cos(t) * np.cos(s)
y = np.cos(t) * np.sin(s)
z = np.sin(t)
ax = plt.subplot(111, projection='3d')
ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
ax.set_zlabel('z axis') 
ax.quiver(0, 0, 0, 0, 0, 1, color=(0, 0, 1, 0.5)) 
# ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='Blues')
"""

k = np.array([0., 0., 1.])
alpha = 0.1
alpha_2 = 0.01

for i in range(3):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.quiver(0, 0, 0, 0, 0, 1, color=(0, 0, 1, 0.5)) 
    # ax = fig.add_subplot(111, projection='3d')
    wi = np.random.randn(3)
    while np.dot(wi, k) <= 0.0:
        wi = np.random.randn(3)
    # normalize
    wi = wi / np.sqrt(np.sum(wi**2))
    ai, bi, ci = wi
    ax.quiver(0, 0, 0, ai, bi, ci, color=(0,1,0,0.5)) 
    print("wi: ", end="")
    print(wi)
    for j in range(1000):
        x, y = np.random.uniform(0, 1, [2, 1])
        phi = 2.0 * np.pi * y
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        # GGX
        tan_theta_m_2 = alpha_2 * x / (1.0 - x)
        cos_theta = 1.0 / np.sqrt(1.0 + tan_theta_m_2)
        cos_theta_2 = cos_theta * cos_theta
        sin_theta = 1.0 - cos_theta_2
        m = np.array([0., 0., 1.])
        m[0] = cos_phi * sin_theta
        m[1] = sin_phi * sin_theta
        m[2] = cos_theta
        m = m / np.sqrt(np.sum(m**2))
        m = np.transpose(m)
        cos_theta_m_i = np.dot(m, wi)
        h = 2.0 * cos_theta_m_i * m
        wo = h - wi
        ao, bo, co = wo
        if np.dot(wo, k) >= 0.0:
            ax.quiver(0, 0, 0, ao, bo, co, color=(1,0,0,0.5)) 

    figpath = "../../images/reflect_" + str(i+1)
    plt.savefig(figpath)

plt.show()