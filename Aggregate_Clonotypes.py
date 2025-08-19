from pathlib import Path
import pandas as pd
from functools import reduce
from operator import or_

p = Path("./data/BCR")

metadata = pd.read_table(
    "./data/20240502_Tidy_MBL_CLL_metadata_vdj.tsv"
).loc[:, ["patient_id", "pool"]]


def read_clonotypes(file: Path):
    return (
        pd.read_csv(file, index_col=0)["cdr3s_nt"]
        .map(lambda clone: frozenset(clone.split(";")))
        .to_dict()
    )


def build_length_dict(clonotypes: dict[str, frozenset]):
    n_chain2id: dict[int, list] = {}

    for id, chainseq in clonotypes.items():
        n = len(chainseq)
        if n not in n_chain2id:
            n_chain2id[n] = list()
        n_chain2id[n].append(id)
    return n_chain2id


def merge_clonotypes(clonotypes: dict[str, frozenset]):
    n_chain2id = build_length_dict(clonotypes)
    n_chains = sorted(n_chain2id.keys())

    merge_id: dict = {clone: clone for clone in clonotypes}

    for i, n in enumerate(n_chains[:-1]):
        for clone in n_chain2id[n]:
            clone_seq = clonotypes[clone]

            merges = list()
            for m in n_chains[i + 1 :]:
                for clone2 in n_chain2id.get(m, []):
                    n_chain_match = len(clone_seq.intersection(clonotypes[clone2]))
                    if n_chain_match == n:
                        merges.append(clone2)
                if len(merges) == 1:
                    merge_id[clone] = merges[0]
                    break
                elif len(merges) > 1:
                    break

    for current_clone, new_clone in merge_id.items():
        old_clone = current_clone
        while old_clone != new_clone:
            old_clone = new_clone
            new_clone = merge_id[old_clone]
        merge_id[current_clone] = new_clone

    return merge_id


def match_sampleclone_to_mergeclone(all_clonotypes: dict, sample_clonotypes: dict):
    seq2id = {v: k for k, v in all_clonotypes.items()}
    return {s_id: seq2id[seq] for s_id, seq in sample_clonotypes.items()}


for patient, df in metadata.groupby("patient_id"):

    all_samples = {}
    for _, sample in df.iterrows():
        sample_file = p / sample.pool / "clonotypes.csv"
        all_samples[sample.pool] = read_clonotypes(sample_file)

    clonotypes: dict[int, frozenset] = dict(
        enumerate(reduce(or_, map(lambda s: set(s.values()), all_samples.values())))
    )
    aggregated_clonotypes = merge_clonotypes(clonotypes)

    for _, sample in df.iterrows():
        merge_file = p / sample.pool / "aggregated_clonotypes.csv"
        clonotype_sample2merge = match_sampleclone_to_mergeclone(
            clonotypes, all_samples[sample.pool]
        )
        aggregate_clonotype_sample = {
            sample_clone: aggregated_clonotypes[global_clone]
            for sample_clone, global_clone in clonotype_sample2merge.items()
        }
        pd.DataFrame.from_dict(
            aggregate_clonotype_sample, columns=["merged_clone"], orient="index"
        ).reset_index(names=["sample_clone"]).to_csv(merge_file, index=False)
