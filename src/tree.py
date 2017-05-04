# -*- coding: utf-8 -*-
# Copyright 2014 by Ambuj Kumar, Kimball-Braun lab group, University of Florida.
# All rights reserved. This code is part of the EvoIntNet distribution and governed
# by its license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Bug reports welcome: ambuj@ufl.edu
#
# Performes proteincoevolution analysis by using phylogeny tree data

import glob
import sys

from Bio import Phylo
from Bio.Phylo.PhyloXML import Phylogeny

from scipy.stats import pearsonr


def rem_redundant(pairs): return set((a,b) if a<=b else (b,a) for a,b in pairs)

def unique_pairs(taxa): return rem_redundant([(i, j) for i in taxa for j in taxa if i != j])

def _sort_dict_by_key(d): return sorted(d.items(), key=lambda x: x[0])



def tree_coevol(tree_folder):
    
    """
        @tree_folder - folder name that contains all the phylogeny tree dataset
        
        This function uses pairwise distance obtained from the tree branch lengths
        to perform coevolution test.
        
        """
    
    tree_data = dict()
    tree_objects = glob.glob(tree_folder + "/*.*")
    taxa_dict = dict()
    
    toolbar_width = len(tree_objects)
    print("Importing branch lengths\n")
    for c, tree_file in enumerate(tree_objects):
        tree = Phylo.read(tree_file, "newick")
        phy = Phylogeny.from_tree(tree)
        taxa = [x. name for x in phy.get_terminals('level')]
        taxa_dict[tree_file] = (taxa)
        tree_data[tree_file] = dict()
        taxon_pairs = unique_pairs(taxa)
        for paired_taxa in taxon_pairs:
            tree_data[tree_file][paired_taxa[0] + "-" + paired_taxa[1]] = phy.distance(paired_taxa[0], paired_taxa[1])
        
        p = str((float(c)/toolbar_width)*100)[:5]
        sys.stdout.write("\r%s%%" %p)
        sys.stdout.flush()

    tree_pairs = unique_pairs(tree_objects)
    correl = dict()

    toolbar_width = len(tree_pairs)
    print("Running coevolution test\n")
    for c, Obj in enumerate(tree_pairs):
        if len(tree_data[Obj[0]]) != len(tree_data[Obj[1]]):
            common_taxa = [key for key, val in tree_data[Obj[0]].items() if key in [key2 for key2, val2 in tree_data[Obj[1]].items()]]
        else:
            common_taxa = tree_data[Obj[0]].keys()

        correl[Obj[0].lstrip(tree_folder).lstrip("/") + "-" + Obj[1].lstrip(tree_folder).lstrip("/")] = \
                pearsonr([val for key, val in _sort_dict_by_key(tree_data[Obj[0]]) if key in common_taxa], \
                        [val for key, val in _sort_dict_by_key(tree_data[Obj[1]]) if key in common_taxa])

        p = str((float(c)/toolbar_width)*100)[:5]
        sys.stdout.write("\r%s%%" %p)
        sys.stdout.flush()


    for key, val in correl.items():
        print key, val


