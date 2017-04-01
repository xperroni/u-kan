#!/usr/bin/env python

from math import cos, sin

from matplotlib import pyplot as pp


def plot_readings(readings, color):
    (x, y) = zip(*readings)
    pp.plot(x, y, color)


def plot(path):
    estimates = []
    ground_truth = []

    for line in open(path):
        values = [float(value) for value in line.split('\t')]
        estimates.append(values[2:4])
        ground_truth.append(values[8:10])

    plot_readings(ground_truth, 'b.')
    plot_readings(estimates, 'r.')

    pp.show()


def main():
    from sys import argv
    plot(argv[1])


if __name__ == '__main__':
    main()
