"""
Run with several values of N, to test.

Time how long each takes.

Verify output: begins with "Hello", ends with "Goodbye" and is actually gzipped.
"""
import gzip
import pathlib

file = "temp.txt.gz"
file = pathlib.Path(file).resolve()
file.touch()

message = "All work and no play makes Homer something something. Go crazy? Don't mind if I do!!"

N = 1000000
with gzip.open(file, 'w') as f:
    f.write("Hello")
for i in range(N):
    with gzip.open(file, 'a') as f:
        f.write(message)
with gzip.open(file, 'a') as f:
    f.write("Goodbye")



