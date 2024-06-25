import argparse
import pandas as pd
import os
import numpy as np
from scipy import sparse
from umap import UMAP
from plotnine import (ggplot, aes, geom_point,
                      theme_void, scale_size_manual, scale_alpha_manual, ggsave)
import hdbscan
from ortho_classes import synth_clean, Part2
from collections import namedtuple
from itertools import chain



parser = argparse.ArgumentParser()
parser.add_argument("file", help='Clean tRNADB-CE file')
parser.add_argument("amino_acid", help='Amino Acid specified for tRNA generation.')
parser.add_argument("-ac", "--anticodon", help='Optional Anticodon specification (will not change output measurably', default='CTA')
parser.add_argument('-o', '--output_directory', help='Directory to store output files.', default='.')
parser.add_argument('-sf', '--synth_file', help='If supplying a Synth, must supply synth file.')
parser.add_argument('-sn', '--synth_name', help='Name of Synth to find similar Synths for. Must be found in Synth File.', nargs='+')
parser.add_argument('-u', '--usearch', help='Name of usearch .exe file.', default='usearch11.0.667_win32')
parser.add_argument('-d', '--distance', help='BLAST homology distance as cutoff.', default=0.2)
parser.add_argument('-i', '--max_id_dist', help='Maximum number of mutations in defined ID parts between potential tRNAs and the tRNA of the synthetase specified.', type=int)
args = parser.parse_args()

df = pd.read_csv(args.file, usecols=['seq_id', 'gen_id', 'Phylum/Class', 'Species', 'Amino Acid', 'tRNA1-7*',
                                     'tRNA8-9*', 'tRNA10-13*', 'tRNA14-21*', 'tRNA22-25*', 'tRNA26*', 'tRNA27-31*',
                                     'tRNA32-38*', 'tRNA39-43*', 'tRNA44-48*', 'tRNA49-53*', 'tRNA54-60*', 'tRNA61-65*',
                                     'tRNA66-72*', 'tRNA73-76*', 'tRNA14-21* aligned', 'tRNA1-7_66-72*',
                                     'tRNA10-13_22-25*', 'tRNA14-21_54-60*', 'tRNA14-21_54-60* aligned',
                                     'tRNA26_44-48*', 'tRNA27-31_39-43*', 'tRNA49-53_61-65*'],
                     dtype={'Amino Acid': 'category', 'tRNA8-9*': 'category', 'tRNA26*': 'category',
                            'tRNA73-76*': 'category'},
                     engine='c')


def fit_umap(dist_matrix, random_state=None, neighbours=15, min_dist=0.1, spread=1.0):
    reducer = UMAP(
        n_neighbors=neighbours,
        random_state=random_state,
        min_dist=min_dist,
        spread=spread,
        metric='precomputed')
    embedding = reducer.fit_transform(dist_matrix)
    return reducer, embedding


def UMAPPER(file, neighbours=15, min_dist=0.1, spread=1.0):
    distmx = pd.read_csv(file, sep='\t', header=None)

    print('Transforming matrix...')
    diagonal = distmx[0] == distmx[1]
    row = np.concatenate([distmx[0], distmx[1][~diagonal]])
    col = np.concatenate([distmx[1], distmx[0][~diagonal]])
    data = 1 - np.concatenate([distmx[2], distmx[2][~diagonal]])
    distmx_ = sparse.csr_matrix((data, (row, col)), dtype=np.float32)
    new_dist = 1 - distmx_.toarray()

    print('Fitting...')

    r, e = fit_umap(new_dist, neighbours=neighbours, min_dist=min_dist, spread=spread)

    print('Formatting...')
    index = pd.DataFrame()
    index['label'] = [i for i in range(0, len(e))]
    index['umap1'] = e[:, 0]
    index['umap2'] = e[:, 1]
    print(f'{len(e)} entries')
    return index


def UMAPPER_clust(file, neighbours=15, min_dist=0.1, spread=1.0, min_samples=60, min_cluster_size=15, chosen=None):
    distmx = pd.read_csv(file, sep='\t', header=None)

    print('Transforming matrix...')
    diagonal = distmx[0] == distmx[1]
    row = np.concatenate([distmx[0], distmx[1][~diagonal]])
    col = np.concatenate([distmx[1], distmx[0][~diagonal]])
    data = 1 - np.concatenate([distmx[2], distmx[2][~diagonal]])
    distmx_ = sparse.csr_matrix((data, (row, col)), dtype=np.float32)
    new_dist = 1 - distmx_.toarray()
    if chosen:
        nearby = []
        for label in chosen:
            nearby_ = [i for i, d in enumerate(new_dist[label, :]) if d < 0.2]
            nearby += nearby_


    print('Fitting...')

    r, e = fit_umap(new_dist, neighbours=neighbours, min_dist=min_dist, spread=spread)

    print('Clustering...')
    labels = hdbscan.HDBSCAN(
        min_samples=min_samples,
        min_cluster_size=min_cluster_size).fit_predict(e)

    print('Formatting...')
    index = pd.DataFrame()
    index['label'] = [i for i in range(0, len(e))]
    index['umap1'] = e[:, 0]
    index['umap2'] = e[:, 1]
    index['cluster'] = labels
    index['cluster'] = index['cluster'].astype('category')
    index['nearby'] = [True if i in nearby else False for i in range(len(index))]
    #print(f'{len(e)} entries')
    return index


class RangeDict(dict):
    def __getitem__(self, item):
        if not isinstance(item, range):
            for i, key in enumerate(self.keys()):
                if isinstance(key, tuple):
                    if any([item in r for r in key]):
                        return list(self.values())[i]
                if item in key:
                    return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item)


def check(test_str):
    return set(test_str) <= allowed


def dist_id(x, id_parts, dist_dict, trna_ids):
    min_score = min([sum([dist_dict[id_part][trna_id][x[id_part]]
                     for id_part in id_parts])
                     for trna_id in trna_ids])
    return min_score


allowed = set('A' + 'C' + 'T' + 'G' + '_')

df = df[df['Amino Acid'] == args.amino_acid]
df['tRNA32-38*'] = df['tRNA32-38*'].apply(lambda x: x[:2] + args.anticodon + x[-2:] if isinstance(x, str) else x)

# RangeDict structure to return part name based on single canonical base number
base_to_part_2 = RangeDict({(range(1, 8), range(66, 73)): 'tRNA1-7_66-72*', range(8, 10): 'tRNA8-9*',
                            (range(10, 14), range(22, 26)): 'tRNA10-13_22-25*',
                            (range(14, 22), range(54, 61)): 'tRNA14-21_54-60*',
                            (range(26, 27), range(44, 49)): 'tRNA26_44-48*',
                            (range(27, 32), range(39, 44)): 'tRNA27-31_39-43*',
                            range(32, 39): 'tRNA32-38*', (range(49, 54), range(61, 66)): 'tRNA49-53_61-65*',
                            range(73, 77): 'tRNA73-76*'})

# And reverse to return the ranges from the part name
part_to_range_2 = {val: key for key, val in base_to_part_2.items()}
part_to_range_2 = {val: list(range_) if isinstance(range_, range) else list(chain.from_iterable(range_))
                                         for val, range_ in part_to_range_2.items()}

if args.synth_name:
    synth_df = synth_clean(args.synth_file)
    trna_ids = list(synth_df[synth_df.synth.isin(args.synth_name)].trna_id)
    # trna_ids = synth_df[synth_df.synth.isin(args.synth_name)]

# ID elements in the form AlaRS: 2, 3, 5...
id_df = pd.read_excel('identity_elements.xlsx', index_col=0)

# Remove anticodon identity elements - these are changing anyway
id_dict = {index[:3]: {int(val): set() for val in value.strip().split(', ') if int(val) not in {34, 35, 36}}
           for index, value in zip(id_df.index, id_df.iloc[:, 0])}
id_dict = id_dict[args.amino_acid].keys()
id_parts = {base_to_part_2[base] for base in id_dict}
id_parts.update(['tRNA73-76*'])
id_parts = list(id_parts)
id_parts_a = id_parts.copy()
if 'tRNA14-21_54-60*' in id_parts:
    id_parts_a.append('tRNA14-21_54-60* aligned')
id_seq_df = df.loc[:, ['seq_id', 'gen_id', 'Phylum/Class', 'Species']+list(id_parts_a)]
print(f'{len(id_seq_df)} {args.amino_acid} tRNAs')

if args.synth_name:
    chosen_df = id_seq_df[id_seq_df.seq_id.isin(trna_ids)]
    chosen_dict = chosen_df.set_index('seq_id').to_dict('index')

id_seq_df['whole'] = id_seq_df.apply(lambda x: '_'.join([x[part] for part in id_parts]), axis=1)
id_sequ_df = id_seq_df.drop_duplicates('whole')
id_sequ_df['check'] = id_sequ_df.whole.apply(lambda x: check(x))
id_sequ_df = id_sequ_df.drop(id_sequ_df[~id_sequ_df.check].index)
print(f'{len(id_sequ_df)} unique combined identity part sequences')

seqs_dict = {}
part_tuple = namedtuple('Part_tuple', ['seq', 'align'])
for part_type in id_parts:
    if part_type != 'tRNA14-21_54-60*':
        seqs_dict[part_type] = [part_tuple(seq, seq) for seq in
                                id_sequ_df[part_type].unique()
                                if isinstance(seq, str)]
    else:
        # Might be quicker to do this with a zip command?
        # Although that would require subsetting dataframe for unique sequences first, then iterating,
        # so maybe not - point for pipeline optimisation in future perhaps
        seqs_dict[part_type] = [part_tuple(row['tRNA14-21_54-60*'], row['tRNA14-21_54-60* aligned'])
                                for index, row in
                                id_sequ_df.drop_duplicates(subset=['tRNA14-21_54-60*']).iterrows()
                                if isinstance(row['tRNA14-21_54-60*'], str)]

seqs_dict = {part_type: {part_t.seq: Part2(part_t.seq, part_type, args.amino_acid, 'TODO', part_t.align, None)
                         for part_t in part_list
                         if isinstance(part_t.align, str)}
             for part_type, part_list in seqs_dict.items()}

chosen_parts = {}
for part_type in id_parts:
    chosen_parts[part_type] = {}
    id_pos = [base for base in id_dict if base in part_to_range_2[part_type]]
    for trna_id in trna_ids:
        part = [part for seq, part in seqs_dict[part_type].items() if seq == chosen_dict[trna_id][part_type]][0]
        id_seq_dict = {pos: base for pos, base in part.seq_dict.items() if pos in id_pos}
        chosen_parts[part_type].update({trna_id: id_seq_dict})

id_distance_dict = {}
for part_type in id_parts:
    id_distance_dict[part_type] = {}
    id_pos = [base for base in id_dict if base in part_to_range_2[part_type]]
    for trna_id in trna_ids:
        id_distance_dict[part_type][trna_id] = {}
        this_id_dict = chosen_parts[part_type][trna_id]
        for seq, part in seqs_dict[part_type].items():
            count = 0
            for pos in id_pos:
                if this_id_dict[pos] != part.seq_dict[pos]:
                    count += 1
            id_distance_dict[part_type][trna_id][seq] = count

id_sequ_df['dist'] = id_sequ_df.apply(lambda x: dist_id(x, id_parts, id_distance_dict, trna_ids), axis=1)
id_sequ_df = id_sequ_df[id_sequ_df.dist <= args.max_id_dist]

id_sequ_df = id_sequ_df.reset_index()
id_sequ_df = id_sequ_df.drop(columns=['index'])
id_sequ_df = id_sequ_df.reset_index()
id_sequ_df = id_sequ_df.rename(columns={'index': 'label'})
id_seq_df = pd.merge(id_seq_df, id_sequ_df.loc[:, ['whole', 'label']], on='whole')
print(f'{len(id_seq_df)} tRNAs <= {args.max_id_dist} ID elements away.')
print(f'Comprising {len(id_sequ_df)} unique combined identity part sequences.')

with open(f'{args.output_directory}/{args.amino_acid}_id_seqs.fa', 'w+') as f:
    for i, row in id_sequ_df.iterrows():
        f.write(f">{i}\n{row['whole']}\n")

os.system(f'{args.usearch} -calc_distmx {args.output_directory}/{args.amino_acid}_id_seqs.fa \
-tabbedout {args.output_directory}/{args.amino_acid}_id_mat.txt -wordlength 2')

if args.synth_name:
    label_id_dict = {list(id_seq_df[id_seq_df.seq_id == trna_id].label)[0]: trna_id for trna_id in trna_ids}
    dist_clusts = UMAPPER_clust(f'{args.output_directory}/{args.amino_acid}_id_mat.txt',
                                chosen=list(label_id_dict.keys()), neighbours=5)
    id_df_u = pd.merge(id_sequ_df, dist_clusts, on='label')
    df_for_merge = id_df_u.loc[:, ['label', 'umap1', 'umap2', 'cluster', 'nearby']]
    id_df = pd.merge(df_for_merge, id_seq_df, on='label')
    id_df.to_csv(f'{args.output_directory}/{args.amino_acid}_df.csv')
    id_df[id_df.nearby == True].to_csv(f'{args.output_directory}/{args.amino_acid}_nearby_df.csv')
    id_df_u['user_defined'] = id_df_u.label.apply(lambda x: label_id_dict.get(x, '0'))
    sizes = {j: i for i, j in zip([0.01] + [2 for i in range(len(args.synth_name))], ['0'] + trna_ids)}
    shapes = {j: i for i, j in zip(['o'] + ['s' for i in range(len(args.synth_name))], ['0'] + trna_ids)}
    p = (ggplot(id_df_u, aes('umap1', 'umap2', colour='cluster', shape='user_defined', size='user_defined',
                             alpha='nearby')) + theme_void() + geom_point() + \
         scale_size_manual(sizes)) + scale_alpha_manual({True: 1, False: 0.5})#+ scale_shape_manual(shapes) + geom_label(check_overlap=True))
    # p = (ggplot() + geom_point(id_df_u, aes('umap1', 'umap2', colour='cluster'), size=0.1, fill=None) + \
    #      geom_point(id_df_u[id_df_u.user_defined != '0'], aes('umap1', 'umap2', shape='user_defined'), shape='o') + theme_void())

else:
    dist_clusts = UMAPPER_clust(f'{args.output_directory}/{args.amino_acid}_id_mat.txt')
    id_df_u = pd.merge(id_sequ_df, dist_clusts, on='label')
    df_for_merge = id_df_u.loc[:, ['label', 'umap1', 'umap2', 'cluster']]
    id_df = pd.merge(df_for_merge, id_seq_df, on='whole')
    id_df.to_csv(f'{args.output_directory}/{args.amino_acid}_df.csv')
    p = (ggplot(id_df_u, aes('umap1', 'umap2', colour='cluster')) + theme_void() + geom_point())

print(f'{len(id_df[id_df.nearby])} tRNAs passing homology threshold {args.distance}')
print(f'Comprising {len(id_df_u[id_df_u.nearby])} unique combined identity parts.')
ggsave(p, f'{args.output_directory}/{args.amino_acid}_umap.pdf', dpi=300)








