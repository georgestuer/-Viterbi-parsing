import numpy as np
import sys


# Read the probability matrix
filename1 = sys.argv[1]
temp = []
with open(filename1, 'r') as f:
    for line in f:
        temp.append(line.strip('\n').split())

for i in range(1, len(temp)):
    temp[i] = list(map(float, temp[i]))

# Read the genome sequence
filename2 = sys.argv[2]
x = ''
with open(filename2, 'r') as f:
    for line in f:
        if line[0] != '>':
            x += line.rstrip().upper()

# probability initialization
init_prob = list(map(float, temp[1]))
e = temp[0][2]  # 'ACGT'
tran_prob = [temp[2][:2], temp[3][:2]]
emis_prob = [temp[2][2:], temp[3][2:]]
transition = {}
emission = {}

# Convert transition prob into dictionary
# e.g. {'+': {'+': 0.999, '-': 0.001}, '-': {'+': 0.01, '-': 0.99}} where "+" indicates A and "-" indicates B
t = ['+', '-']
transition['+'] = {}
transition['-'] = {}
for i in range(len(t)):
    for j in range(len(t)):
        transition[t[i]][t[j]] = transition[t[i]].get(t[j], 0) + tran_prob[i][j]

# Convert emission prob into dictionary
# e.g. {'+': {'A': 0.35, 'C': 0.15, 'G': 0.15, 'T': 0.35}, '-': {'A': 0.15, 'C': 0.35, 'G': 0.35, 'T': 0.15}}
#      where "+" indicates A and "-" indicates B
e = [e[i] for i in range(len(e))] # To convert 'ACGT' into 'A', 'C', 'G', 'T'
emission['+'] = {}
emission['-'] = {}
for i in range(len(t)):
    for j in range(len(e)):
        emission[t[i]][e[j]] = emission[t[i]].get(e[j], 0) + emis_prob[i][j]

NONZERO = 10**(-30)  # To avoid log2(0) in the calculation

# Table initialization
rows, cols = len(t), len(x)
S = np.zeros(shape=(rows, cols), dtype=float)
back = np.zeros(shape=(rows, cols), dtype=int)

# Fill in the first column
for i in range(rows):
    #S[i, 0] = emission[t[i]][x[0]] * init_prob[i]
    S[i, 0] = np.log2(init_prob[i]) + np.log2(emission[t[i]][x[0]])

# Fill the rest of the table
for j in range(1, cols):
    for i in range(rows):
        ep = np.log2(emission[t[i]][x[j]] + NONZERO)
        max, imax = S[0, j - 1] + np.log2(transition[t[0]][t[i]]) + ep, 0
        for i1 in range(1, rows):
            pr = S[i1, j - 1] + np.log2(transition[t[i1]][t[i]]) + ep
            if pr > max:
                max, imax = pr, i1
        S[i, j], back[i, j] = max, imax

# Find final S[k, n] by comparing S[A, n] and S[B, n]
max1, imax1 = S[0, cols-1], 0
for i in range(1, rows):
    if S[i, cols-1] > max1:
        max1, imax1 = S[i, cols-1], i

# Backtrace to find the most likely path
i, path = imax1, [imax1]
for j in range(cols-1, 0, -1):
    i = back[i, j]
    path.append(i)

# Write the position information into .txt file
def write_output(filename, seq):
    cpgs = []
    i = 0
    count_A = 0
    count_B = 0
    while i < len(seq):
        if seq[i] == 0:
            j = i + 1
            while j < len(seq):
                if seq[j] == 1:
                    break
                j += 1
            cpgs.append([i + 1, j, 'state A'])
            count_A += 1
            i = j
        else:
            j = i + 1
            while j < len(seq):
                if seq[j] == 0:
                    break
                j += 1
            cpgs.append([i + 1, j, 'state B'])
            count_B += 1
            i = j
    with open(filename, "w") as f:
        for i in cpgs:
            f.write(str(i[0]) + " " + str(i[1]) + " " + str(i[2]) + "\n")
        f.write("State A " + str(count_A) + "\n")
        f.write("State B " + str(count_B) + "\n")
    return cpgs, count_A, count_B


cpgs, count_A, count_B = write_output('output.txt', path[::-1])
print("The number of segments in state B is ", count_B)
print("See final output in output.txt")





