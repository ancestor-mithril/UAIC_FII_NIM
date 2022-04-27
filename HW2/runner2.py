from pathos.pools import ThreadPool as Pool
import subprocess
import random

def read_lines():
    lines = [x.rstrip('\n') for x in open("experiment.txt", "r").readlines()]

    lines1 = lines[:len(lines) //2]
    lines2 = lines[len(lines) //2:]
    return lines1, lines2

def run1(line):
    index, line = line
    print(index)
    f = open(f"experiments/1-{index}.txt", "w")
    args = ["./build/app.exe", "10", "1", ] + line.split() + ["true", "false"]
    subprocess.call(args, stdout=f)
    return 0

def run2(line):
    index, line = line
    print(index)
    f = open(f"experiments/2-{index}.txt", "w")
    args = ["./build/app.exe", "10", "2", ] + line[0].split() + ["true", "false"] + line[1].split() + ["true", "false"]
    subprocess.call(args, stdout=f)
    return 0

def run3(line):
    index, line = line
    print(index)
    f = open(f"experiments/3-{index}.txt", "w")
    args = ["./build/app.exe", "10", "3", ] \
        + line[0].split() + ["true", "false"] + \
            line[1].split() + ["true", "false"] + \
                line[2].split() + ["true", "false"]
    subprocess.call(args, stdout=f)
    return 0

def use1(usedLines):
    with Pool(8) as p:
        results = p.map(run1, enumerate(usedLines))


def use2(usedLines):
    secondLines = usedLines
    random.shuffle(secondLines)
    usedLines = [(x, y) for x, y in zip(usedLines[:len(usedLines)//2], secondLines[len(usedLines)//2:])]
    with Pool(8) as p:
        results = p.map(run2, enumerate(usedLines))

def use3(usedLines):
    secondLines = usedLines
    random.shuffle(secondLines)
    third = usedLines
    random.shuffle(third)
    usedLines = [(x, y, z) for x, y, z in zip(usedLines[:len(usedLines)//2], secondLines[len(usedLines)//2:], third[len(usedLines)//2:])]
    with Pool(8) as p:
        results = p.map(run3, enumerate(usedLines))

use1(read_lines()[0])
