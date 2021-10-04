#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def euler(initial_x, initial_y, step, slope, end):
	current_x = initial_x
	state_history = []
	while( current_x < end+step ): # Can't use <= end b/c floating-point wonkiness
		if( current_x == initial_x ):
			state_history += [(initial_x,initial_y)]
		else:
			last_state = state_history[-1][-1]
			new_estimate = last_state + step*(slope(current_x-step, last_state))
			state_history += [(current_x,new_estimate)]
		current_x += step
	return state_history

def heun(initial_x, initial_y, step, slope, end):
	current_x = initial_x
	state_history = []
	while( current_x < end+step ): # Can't use <= end b/c floating-point wonkiness
		if( current_x == initial_x ):
			state_history += [(initial_x,initial_y)]
		else:
			last_state = state_history[-1][-1]
			euler_estimate = last_state + step*(slope(current_x-step, last_state))
			heun_estimate = last_state + (step/2.0)*(euler_estimate + slope(current_x,euler_estimate))
			state_history += [(current_x,heun_estimate)]
		current_x += step
	return state_history

def plot(beta, step, euler, heun):
	df_e = pd.DataFrame(euler, columns=["t", "infected"])
	df_e["Estimate"] = "Euler"
	df_h = pd.DataFrame(heun, columns=["t", "infected"])
	df_h["Estimate"] = "Heun"
	df = pd.concat([df_e, df_h])
	sns.lineplot(data=df, x="t", y="infected", hue="Estimate")
	plt.title("Beta = %f, Step Size = %f" % (beta,step))
	fname = "sis_%f_%f.pdf" % (beta,step)
	plt.savefig(fname, bbox_inches="tight")
	plt.clf()


if __name__ == "__main__":
	# Example usage
	"""
	def slope(x, y):
		return x + 2*y
	initial = 0

	print(euler(2, 3, 0.1, slope, 2.5))
	print(heun(2, 3, 0.1, slope, 2.5))
	"""

	initial_susceptible = 90
	initial_infected = 10
	initial_y = (initial_susceptible, initial_infected)
	initial_x = 0

	N = 100
	gamma = 0.25
	for beta in [0.03, 0.06, 0.1]:
		for step in [0.01, 0.5, 2.0]:
			def infection_slope(t, current_infected):
				current_susceptible = N-current_infected
				return beta*current_susceptible*current_infected - gamma*current_infected
			euler_estimate = euler(0, initial_infected, step, infection_slope, 50*step)
			heun_estimate = heun(0, initial_infected, step, infection_slope, 50*step)
			plot(beta, step, euler_estimate, heun_estimate)