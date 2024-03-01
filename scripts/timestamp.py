# %%
import pandas as pd

# %%
df = pd.read_csv('./output/timestamp.csv')
df.head()

# %%
#group df by PHASE and get the mean of each group and convert the time in nanoseconds to a seconds

# Calculando a soma por grupo e convertendo para segundos
grouped_sums_seconds = df.groupby('PHASE').sum().apply(lambda x: x / 1e9)

# Ordenando os resultados pela soma dos tempos em segundos
sorted_grouped_sums_seconds = grouped_sums_seconds.sort_values(by='TIME', ascending=True)

# Exibindo o resultado ordenado
print(sorted_grouped_sums_seconds)

