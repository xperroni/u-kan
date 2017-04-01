#!/usr/bin/env python

from re import search

from matplotlib import pyplot as pp


def plot(path):
    laser = 2
    radar = 3

    e = {
        laser: [],
        radar: []
    }

    for line in open(path):
        match = search(r'e\((\d+)\): ([^\r\n]+)', line)
        if match == None:
            continue

        s = int(match.group(1))
        v = float(match.group(2))
        e[s].append(v)


    pp.plot(e[laser], 'b-')
    pp.plot(e[radar], 'r-')
    pp.plot([5.991] * len(e[laser]), 'b--')
    pp.plot([7.815] * len(e[radar]), 'r--')
    pp.show()


def main():
    from sys import argv
    plot(argv[1])


if __name__ == '__main__':
    main()
