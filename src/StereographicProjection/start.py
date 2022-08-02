from manimlib import *
from numpy import array

class SteregraphicProjection(Scene):
    def construct(self):
        # axes = Axes((-1, 1), (-1, 1))
        # axes.add_coordinate_labels()
        # self.play(Write(axes, lag_ratio=0.01, run_time=1))

        nums = 10
        step = 1/nums
        z = 0
        for i in range(nums+1):
            circle = Circle(
                radius = 4 * (math.sqrt(1-z*z))/(1+z)
            )
            self.play(ShowCreation(circle))
            z += step
        nums = 10
        step = np.pi/nums
        theta = 0
        for i in range(nums+1):
            line = Line(
                start=array([4 * np.cos(theta), 4 * np.sin(theta)]),
                end=array([-4 * np.cos(theta), -4 * np.sin(theta)])
            )
            self.play(ShowCreation(line))
            theta += step