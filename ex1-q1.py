# Write your code instead of the 'pass' statements
import pandas as pd
import matplotlib.pyplot as plt

def count_mutations_per_sample(path):
    df = pd.read_csv(path,delimiter="\t")
    dict = df.groupby("Sample name")["Mutation ID"].count().astype(int).to_dict()
    return dict


def count_mutation_types(path):
    df = pd.read_csv(path,delimiter="\t")
    dict = df.groupby("Mutation Type")["Mutation Type"].count().astype(int).to_dict()
    return dict


def count_mutations_per_tissue(path):
    df = pd.read_csv(path, delimiter="\t")
    grouped = df.groupby("Primary site")["Mutation Type"].value_counts().unstack(level=0).fillna(0).astype(int).to_dict()
    return grouped

def max_likelihood(d,tissue,mutation):
    n = sum(d[tissue].values())
    k = d[tissue][mutation]
    return k/n


def median_calc(d):
    lst = pd.Series(d.values())
    med = lst.median()
    return med

def plot_histogram (path):
    dict = count_mutations_per_sample(path)
    lst = dict.values
    mean = sum(lst)/len(lst)
    median = median_calc(dict)
    plt.hist(lst, bins=500)
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1, label='Mean: {:.2f}'.format(mean))
    plt.axvline(median, color='green', linestyle='dashed', linewidth=1, label='Median: {:.2f}'.format(median))
    plt.xlim(min(lst)-1,2000)
    plt.legend()
    plt.show()



def extract_above_median_samples(d):
    lst = []
    median = median_calc(d)
    for key in d.keys():
        if d[key] >= median:
            lst.append(key)
    return lst


def count_mutations_per_tissue_2(path, sample_lst):
    df = pd.read_csv(path, delimiter="\t")
    filtered = df[df["Sample name"].isin(sample_lst)]
    grouped = filtered.groupby("Primary site")["Mutation Type"].value_counts().unstack(level=0).fillna(0).astype(int).to_dict()
    return grouped

def max_likely_counter(d):
    count = 0
    for tissue in d.keys():
        if max_likelihood(d,tissue,"intrachromosomal deletion") > 0.2:
            count += 1


if __name__ == '__main__':
    pass  # you can change this main section however you would like



