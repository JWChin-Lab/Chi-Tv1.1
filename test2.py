import random
from collections import namedtuple
from test import max_dist_parallel_memo

Part = namedtuple('Part', ['name', 'seq'])

s_seqs1 = [Part(i, ''.join([''.join([random.choice(['A', 'C', 'T', 'G']) for j in range(4)]) + '_' +
           ''.join([random.choice(['A', 'C', 'T', 'G']) for j in range(4)])])) for i in range(30)]


def do(part_list, num_seqs):
    global memo
    memo = {}
    print(max_dist_parallel_memo(part_list, num_seqs, memo))


# def do1():
#     do(s_seqs1, 4)


do(s_seqs1, 4)