import functools
import time
from Bio import SeqIO
from time import sleep
from mpire import WorkerPool

reference = "/home/joshfactorial/Documents/neat_data/H1N1/H1N1_HA.fa"
reference_index = SeqIO.index(reference, 'fasta')

def time_consuming_function(x: int) -> int:
    sleep(0.1)   # Simulate that this function takes some time to complete
    return x


data = range(100)

start = time.time()
results = [time_consuming_function(x) for x in data]
print(f'list comprehension: {time.time()-start}')

# MPIRE
start = time.time()
with WorkerPool(n_jobs=5) as pool:
    results = pool.map(time_consuming_function, data, progress_bar=True)
print(f'MPIRE: {time.time()-start}')