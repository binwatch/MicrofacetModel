import cv2
import numpy as np

width = 500
height = 500
img = np.zeros((height, width, 3), np.uint8)
img[:, :, 0] = 255
img[:, :, 1] = 255
img[:, :, 2] = 255

xstep = 2 / width
ystep = 2 / height
for i in range(width):
    x = -1.0 + i * xstep
    for j in range(height):
        y = -1.0 + j * ystep
        if (x * x + y * y <= 1.0):
            base = x * x + y * y + 1
            img[i, j, 0] = int((2*x)/(base)*255)
            img[i, j, 1] = int((2*y)/(base)*255)
            img[i, j, 2] = int((1-x*x-y*y)/(base)*255)

cv2.imwrite("./images/color.png", img)