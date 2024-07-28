import random
# function to read a FASTA file. Returns a string
def read_fasta(filename):
    f = open(filename)
    header = f.readline()
    x = f.read()
    x = x.replace("\n", "")
    return x


# Part A
# Scoring regime
def score(a, b):
    if (a == b):
        return 3
    elif ((a == "-" and b != "-") or (a != "-" and b == "-")):
        return -2
    else:
        return -3


def local_pwalignment(S, T, score=score):
    n = len(S)
    m = len(T)
    d = [[0 for j in range(m + 1)] for i in range(n + 1)]
    max_score = 0
    max_i = 0
    max_j = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d[i][j] = max(0, d[i - 1][j] + score("-", T[j - 1]), d[i][j - 1] + score(S[i - 1], "-"),
                          d[i - 1][j - 1] + score(S[i - 1], T[j - 1]))
            if max_score < d[i][j]:
                max_score = d[i][j]
                max_i = i
                max_j = j
    aligned_s = ""
    aligned_t = ""
    while d[max_i][max_j] != 0:
        if d[max_i][max_j] == d[max_i-1][max_j] + score(S[max_i-1],"-"):
            aligned_s = S[max_i-1] + aligned_s
            aligned_t = "-" + aligned_t
            max_i -= 1
        elif d[max_i][max_j] == d[max_i][max_j-1] + score("-",T[max_j-1]):
            aligned_s = "-" + aligned_s
            aligned_t = T[max_j-1] + aligned_t
            max_j -= 1
        else:
            aligned_s = S[max_i-1] + aligned_s
            aligned_t = T[max_j-1] + aligned_t
            max_j -= 1
            max_i -= 1
    return (max_score,aligned_s,aligned_t)




def permutation_test(S, T):
    original = local_pwalignment(S,T)
    s_lst = list(S)
    t_lst = list(T)
    res = 0
    for i in range(100):
        random.shuffle(s_lst)
        random.shuffle(t_lst)
        S = ''.join(s_lst)
        T = ''.join(t_lst)
        if original<local_pwalignment(S,T):
            res += 1
    return (res/100) < 0.05


if __name__ == '__main__':
    pass  # you can change this section however you would likes

