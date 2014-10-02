#!/usr/bin/python3
# Petter Strandmark

import sys

import matplotlib.pyplot as plt

LOOKING_FOR_DATA = 1
READING_DATA = 2
state = LOOKING_FOR_DATA

iterations            = []
objectives            = []
relative_changes_in_x = []
relative_changes_in_y = []
feasibilities         = []

optimal_value = None

for line in sys.stdin.readlines():
	if state == LOOKING_FOR_DATA:
		if "--------" in line:
			state = READING_DATA
			print("Found convergence data.")
		elif "Optimal value:" in line:
			optimal_value = float(line.split(":")[1])
			print("Found optimal value:", optimal_value)
	elif state == READING_DATA:
		try:
			iteration, objective, rel_change_x, rel_change_y, feasibility = \
				[float(n) for n in line.split()]
			iterations.append(iteration)
			objectives.append(objective)
			relative_changes_in_x.append(rel_change_x)
			relative_changes_in_y.append(rel_change_y)
			feasibilities.append(feasibility)
		except:
			state = LOOKING_FOR_DATA

plt.semilogy(iterations,
             feasibilities)
plt.xlabel("Iteration")
plt.title("Feasibility")
plt.show()

if optimal_value:
	plt.semilogy(iterations,
	             [abs(obj - optimal_value) for obj in objectives])
	plt.xlabel("Iteration")
	plt.title("Objective value error")
	plt.show()
