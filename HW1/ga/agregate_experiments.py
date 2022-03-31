import os
import re
import json
# TODO: also do implement helper methods in c++

def agregate_experiments_1():
    file_regex = re.compile("Crossover : (.*)\n" + 
    "Mutation : (.*)\n" +
    "hyperMutations : (.*)\n" +
    "elites : (.*)\n" +
    "selectionPressure : (.*)\n" +
    "crossoverTypes : (.*)\n" +
    "hillclimbings : (.*)\n" +
    "hyperMutationSteps : (.*)\n" +
    "encodingChanges : (.*)\n" +
    "(.*)\n(.*)\n(.*)\n", re.MULTILINE
    )
    experiments_collection = []
    for file in os.listdir("experiments"):
        with open(f"experiments/{file}") as f:
            match = re.match(file_regex,f.read())
        
        if not match:
            print("No match", file)
            continue
        if len(match.groups()) != 12:
            print("Not enough matches", file)
            continue

        current_experiment = dict()
        current_experiment["crossoverProbability"] = match.group(1)
        current_experiment["mutationProbability"] = match.group(2)
        current_experiment["hypermutationRate"] = match.group(3)
        current_experiment["elitesPercentage"] = match.group(4)
        current_experiment["selectionPressure"] = match.group(5)
        current_experiment["crossoverType"] = match.group(6)
        current_experiment["hillclimbingType"] = match.group(7)
        current_experiment["stepsToHypermutation"] = match.group(8)
        current_experiment["encodingChangeRate"] = match.group(9)
        current_experiment["runs"] = 3
        current_experiment["value"] = (float(match.group(10)) + float(match.group(11)) + float(match.group(12)))/3
        experiments_collection.append(current_experiment)

    open("statistics/agregate.json", "w").write(json.dumps(experiments_collection))


# agregate_experiments_1()

def analyze_experiments():
    experiments = json.load(open("statistics/agregate.json", "r"))
    x =  sorted(experiments, key=lambda x: x["value"])
    crossoverProbability = dict()
    mutationProbability = dict()
    hypermutationRate = dict()
    elitesPercentage = dict()
    selectionPressure = dict()
    crossoverType = dict()
    hillclimbingType = dict()
    stepsToHypermutation = dict()
    encodingChangeRate = dict()

    attrs = [
        "crossoverProbability",
        "mutationProbability",
        "hypermutationRate",
        "elitesPercentage",
        "selectionPressure",
        "crossoverType",
        "hillclimbingType",
        "stepsToHypermutation",
        "encodingChangeRate",
    ]
    for i in x[:100]:
        for attr in attrs:
            attr_v = i[attr]
            attr_dict = eval(attr)
            if attr_v in attr_dict:
                attr_dict[attr_v] += 1
            else:
                attr_dict[attr_v] = 1
    print(x[0])
    print(x[1])
    print(x[2])
    print("crossoverProbability")
    print(crossoverProbability)
    print("mutationProbability")
    print(mutationProbability)
    print("hypermutationRate")
    print(hypermutationRate)
    print("elitesPercentage")
    print(elitesPercentage)
    print("selectionPressure")
    print(selectionPressure)
    print("crossoverType")
    print(crossoverType)
    print("hillclimbingType")
    print(hillclimbingType)
    print("stepsToHypermutation")
    print(stepsToHypermutation)
    print("encodingChangeRate")
    print(encodingChangeRate)


def print_statistics():
    import statistics
    f = open("../statistics/1.txt", "w")
    path = "../experiments/20/2"
    for file in os.listdir(path):
        vals = list(map(float, open(f"{path}/{file}", "r").readlines()))
        if len(vals) == 0:
            continue
        vals1 = vals[::2]
        vals2 = list(map(int, vals[1::2]))
        f.write(f"{file}")
        f.write(f" {statistics.mean(vals1):.2f}")
        f.write(f" {statistics.median(vals1):.2f}")
        f.write(f" {statistics.stdev(vals1):.2f}")
        f.write(f" {min(vals1):.2f}")
        f.write(f" {max(vals1):.2f}")
        f.write("\n")
        # f.write(f" {statistics.mean(vals2):.2f}\n")


print_statistics()
