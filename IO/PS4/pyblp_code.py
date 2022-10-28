import pyblp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pyblp.options.digits = 2
pyblp.options.verbose = False
pyblp.__version__

product_data = pd.read_csv(pyblp.data.NEVO_PRODUCTS_LOCATION)
product_data.head()

X1_formulation = pyblp.Formulation('0 + prices', absorb='C(product_ids)')
X2_formulation = pyblp.Formulation('1 + prices + sugar + mushy')
product_formulations = (X1_formulation, X2_formulation)
product_formulations

agent_data = pd.read_csv(pyblp.data.NEVO_AGENTS_LOCATION)
agent_data.head()

agent_formulation = pyblp.Formulation('0 + income + income_squared + age + child')
agent_formulation

nevo_problem = pyblp.Problem(product_formulations, product_data, agent_formulation, agent_data)
nevo_problem

initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([[ 5.4819,  0,      0.2037,  0     ], [15.8935, -1.2000, 0,       2.6342], [-0.2506,  0,      0.0511,  0     ],[ 1.2650,  0,     -0.8091,  0     ]])
bfgs = pyblp.Optimization('bfgs', {'gtol': 1e-5})
nevo_results = nevo_problem.solve(initial_sigma,initial_pi,optimization=bfgs, method='1s')
nevo_results    

elasticities = nevo_results.compute_elasticities()
diversions = nevo_results.compute_diversion_ratios()

single_market = product_data['market_ids'] == 'C01Q1'
plt.colorbar(plt.matshow(elasticities[single_market]));

plt.colorbar(plt.matshow(diversions[single_market]));

means = nevo_results.extract_diagonal_means(elasticities)

aggregates = nevo_results.compute_aggregate_elasticities(factor=0.1)


plt.hist([means.flatten(), aggregates.flatten()],color=['red', 'blue'],bins=50);
plt.legend(['Mean Own Elasticities', 'Aggregate Elasticities']);

costs = nevo_results.compute_costs()
print("Mean costs: ", np.mean(costs))
print("Median costs: ", np.median(costs))
print("Standard deviation: ", np.std(costs))
plt.hist(costs, bins=50);
plt.legend(["Marginal Costs"]);
plt.show()

prices = nevo_results.compute_prices(costs=costs)
print("Mean prices: ", np.mean(prices))
print("Median prices: ", np.median(prices))
print("Standard deviation: ", np.std(prices))
plt.hist(prices, bins=50);
plt.legend(["Pre-merger prices"]);
plt.show()

margins = nevo_results.compute_markups(costs=costs)
markups = prices - costs
print("Mean markups: ", np.mean(markups))
print("Median markups: ", np.median(markups))
print("Standard deviation: ", np.std(markups))

print("Mean margins: ", np.mean(margins))
print("Median margins: ", np.median(margins))
print("Standard deviation: ", np.std(margins))

plt.hist(markups, bins=50);
plt.legend(["Pre-merger markups"]);
plt.show()
plt.hist(margins, bins=50);
plt.legend(["Pre-merger margins"]);
plt.show()
shares = nevo_results.compute_shares(prices)

hhi = nevo_results.compute_hhi()
profits = nevo_results.compute_profits(costs=costs)
cs = nevo_results.compute_consumer_surpluses()

##post-nabisco merger
product_data['merger_ids'] = product_data['firm_ids'].replace(3, 6)
changed_prices = nevo_results.compute_prices(firm_ids=product_data['merger_ids'],costs=costs)
changed_shares = nevo_results.compute_shares(changed_prices)
changed_margins = nevo_results.compute_markups(changed_prices, costs)
changed_markups = changed_prices - costs
print("Post-Nabisco merger results:")
print("Mean post-merger price change: ", np.mean(changed_prices)- np.mean(prices))
print("Median post-merger price change: ", np.median(changed_prices)-np.median(prices))
print("Standard deviation: ", np.std(changed_prices))
print("Mean post-merger markup change: ", np.mean(changed_markups)-np.mean(markups))
print("Median post-merger markup change: ", np.median(changed_markups) - np.median(markups))
print("Standard deviation: ", np.std(changed_markups))
print("Mean post-merger margin change: ", np.mean(changed_margins) - np.mean(margins))
print("Median post-merger margin change: ", np.median(changed_margins) - np.median(margins))
print("Standard deviation: ", np.std(changed_margins))plt.hist(changed_prices, bins=50)
print("Mean share change",np.mean(changed_shares) - np.mean(shares))
print("Median share change",np.median(changed_shares) - np.median(shares))
plt.legend(["Post-merger prices, Post-Nabisco merger"])
plt.show()
plt.hist(changed_prices - prices, bins=50);
plt.legend(["Price Changes, Post-Nabisco merger"]);
plt.show()
plt.hist(changed_markups, bins=50)
plt.legend(["Post-merger markups, Post-Nabisco merger"])
plt.show()
plt.hist(changed_markups - markups, bins=50);
plt.legend(["Markup Changes, Post-Nabisco merger"]);
plt.show()
plt.hist(changed_margins, bins=50)
plt.legend(["Post-merger margins, Post-Nabisco merger"])
plt.show()
plt.hist(changed_margins - margins, bins=50);
plt.legend(["Margin Changes, Post-Nabisco merger"]);
plt.show()


#gm-quaker merger
product_data['merger_ids'] = product_data['firm_ids'].replace(2, 4)
changed_prices = nevo_results.compute_prices(firm_ids=product_data['merger_ids'],costs=costs)
changed_shares = nevo_results.compute_shares(changed_prices)
changed_margins = nevo_results.compute_markups(changed_prices, costs)
changed_markups = changed_prices - costs
print("GM-Quaker merger results:")
print("Mean post-merger price change: ", np.mean(changed_prices)- np.mean(prices))
print("Median post-merger price change: ", np.median(changed_prices)-np.median(prices))
print("Standard deviation: ", np.std(changed_prices))
print("Mean post-merger markup change: ", np.mean(changed_markups)-np.mean(markups))
print("Median post-merger markup change: ", np.median(changed_markups) - np.median(markups))
print("Standard deviation: ", np.std(changed_markups))
print("Mean post-merger margin change: ", np.mean(changed_margins) - np.mean(margins))
print("Median post-merger margin change: ", np.median(changed_margins) - np.median(margins))
print("Standard deviation: ", np.std(changed_margins))
print("Mean share change",np.mean(changed_shares) - np.mean(shares))
print("Median share change",np.median(changed_shares) - np.median(shares))
plt.hist(changed_prices, bins=50)
plt.legend(["Post-merger prices, GM-Quaker merger"])
plt.show()
plt.hist(changed_prices - prices, bins=50);
plt.legend(["Price Changes, GM-Quaker merger"]);
plt.show()
plt.hist(changed_markups, bins=50)
plt.legend(["Post-merger markups, GM-Quaker merger"])
plt.show()
plt.hist(changed_markups - markups, bins=50);
plt.legend(["Markup Changes, GM-Quaker merger"]);
plt.show()
plt.hist(changed_margins, bins=50)
plt.legend(["Post-merger margins, GM-Quaker merger"])
plt.show()
plt.hist(changed_margins - margins, bins=50);
plt.legend(["Margin Changes, GM-Quaker merger"]);
plt.show()