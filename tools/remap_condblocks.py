#!/usr/bin/env python3
"""remap_condblocks.py

Tool to remap condition block names from the data files to updated
and more user friendly ones.
"""
import pandas

def remap_cond_blocks(remapping):
    with open('data/cond_blocks.csv') as infile:
        header = infile.readline().strip()
        print(header)
        for line in infile:
            comps = line.strip().split(',')
            if len(comps) != 6:
                raise Exception('wrong number of components')
            comps[4] = remapping[comps[4]]
            print(",".join(comps))

def remap_corem_unique_condblocks(remapping):
    with open('data/20161223.corem.unique.condblocks.csv') as infile:
        header = infile.readline().strip()
        print(header)
        for line in infile:
            comps = line.strip().split(',')
            if len(comps) != 3:
                raise Exception('wrong number of components')
            comps[1] = remapping[comps[1]]
            print(','.join(comps))


if __name__ == '__main__':
    remapping = {}
    with open("data/condition_block_names_new.csv") as infile:
        infile.readline()
        for line in infile:
            comps = line.strip().split(',')
            if len(comps) != 2:
                raise Exception('invalid line !!')
            remapping[comps[0]] = comps[1]

    #remap_cond_blocks(remapping)
    remap_corem_unique_condblocks(remapping)
