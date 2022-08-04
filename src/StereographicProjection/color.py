import cv2
import numpy as np

width = 500
height = 500
img = np.zeros((height, width, 3), np.uint8)
img[:, :, 0] = 255
img[:, :, 1] = 255
img[:, :, 2] = 255

steps_x = 2 / width
steps_y = 2 / height
for i in range(width):
    x = -1.0 + i * steps_x
    for j in range(height):
        y = -1.0 + j * steps_y
        if ((x * x + y * y) <= (1.0)):
            base = x * x + y * y + 1
            img[i, j, 0] = int(max((1 - x*x - y*y)/(base), 0) * 255)
            img[i, j, 2] = int(((1.0 + (2*x)/(base)) / 2.0) * 255)
            img[i, j, 1] = int(((1.0 + (2*y)/(base)) / 2.0) * 255)

cv2.imwrite("./images/color.png", img)