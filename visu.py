import sys

import pandas
import matplotlib.pyplot as plt

# read datafile
filePath = sys.argv[1] if len(sys.argv) > 1 else "data/evrptw_instances/rc208_21.txt"
f = open(filePath, 'r')
cols = list(filter(lambda x : x != '' and x != '\n', f.readline().split(" ")))
values = []
while True:
    line = f.readline()
    if line == "\n":
        break
    values.append(list(filter(lambda x : x != '' and x != '\n', line.split(" "))))
df = pandas.DataFrame(values, columns=cols)
f.close()

print(df)
df.x = df.x.astype(float)
df.y = df.y.astype(float)

# plot the data
# depot
df_depot = df[df.Type == 'd']
plt.scatter(df_depot.x, df_depot.y, c='b', label='depot')
# charging stations
df_stations = df[df.Type == 'f']
plt.scatter(df_stations.x, df_stations.y, c='y', marker='x', label='station')
# customers
df_customers = df[df.Type == 'c']
plt.scatter(df_customers.x, df_customers.y, c='g', label='customer')

plt.legend()
plt.show()
