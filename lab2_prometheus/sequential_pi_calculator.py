#!/usr/bin/env python
import random
from mpi4py import MPI
import socket

radius = 1
random.seed(21596)
NUMBER_OF_POINTS = 10000000


def is_inside_circle(x, y):
    if x * x + y * y < radius * radius:
        return True
    else:
        return False


def pi_aproximator(points_number):
    points_inside_circle = 0
    for i in range(points_number):
        x = random.random()
        y = random.random()
        if is_inside_circle(x, y):
            points_inside_circle += 1
    return 4*(points_inside_circle/points_number)

pi = pi_aproximator(NUMBER_OF_POINTS)
print(pi)

