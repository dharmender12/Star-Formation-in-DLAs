import matplotlib.pyplot as plt
from scipy.stats import bootstrap
import numpy as np
rng = np.random.default_rng()
from scipy.stats import norm
dist = norm(loc=2, scale=4)  # our "unknown" distribution
data = dist.rvs(size=100, random_state=rng)

import matplotlib.pyplot as plt
from scipy.stats import bootstrap



data = (data,)  # samples must be in a sequence
res = bootstrap(data, np.std, confidence_level=0.9,
                random_state=rng)
x = np.linspace(3, 5)
std_sample = np.std(data)
pdf = norm.pdf(x, loc=std_sample, scale=res.standard_error)
fig, ax = plt.subplots()
# ax.hist(res.bootstrap_distribution, bins=25, density=True)
ax.plot(x, pdf)
ax.set_title('Normal Approximation of the Bootstrap Distribution')
ax.set_xlabel('statistic value')
ax.set_ylabel('pdf')
plt.show()