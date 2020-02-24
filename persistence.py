import compress_pickle as pickle
import numpy as np

EXTENSION = "gz" # it is the most fast and good. another choice is lzma, which is 5 times slower, but 1.5 times better

def save(obj, filename):
    pickle.dump(obj, f"{filename}.{EXTENSION}")

def load(filename):
    return pickle.load(f"{filename}.{EXTENSION}")
