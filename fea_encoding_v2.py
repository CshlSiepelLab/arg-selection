#!/usr/bin/env python3

import os

import dendropy
import numpy as np
import tskit


def encode(nwk_str, no_taxa, der_ls=None):

    dpy_tr = dendropy.Tree.get(data=nwk_str, schema="newick", rooting="default-rooted")

    dpy_tr.calc_node_ages(ultrametricity_precision=1e-01) # assuming branch length in generations
    u_vec = np.zeros(no_taxa) # vector of (internal) node ages, last element always 0

    nd_lab = 0
    for int_nd in dpy_tr.ageorder_node_iter(include_leaves=False, filter_fn=None, descending=True):
        int_nd.label = nd_lab # integer label
        u_vec[nd_lab] = int_nd.age
        nd_lab += 1

    F_mtx = np.zeros((no_taxa-1, no_taxa-1), dtype=int)
    W_mtx = np.zeros((no_taxa-1, no_taxa-1))

    # populate W matrix
    for i in range(no_taxa-1):
        for j in range(i+1):
            W_mtx[i, j] = u_vec[j] - u_vec[i+1]

    for head_node in dpy_tr.preorder_node_iter():
        if head_node.label == 0: # root node
            #print("I am root")
            continue
        if head_node.is_leaf():
            edg_head = no_taxa - 1
        else:
            edg_head = head_node.label
            
        edg_tail = head_node.parent_node.label
        
        #print(edg_head, edg_tail)
        for j in range(edg_tail, edg_head):
            for i in range(j, edg_head):
                F_mtx[i, j] += 1

    if der_ls is None:
        return F_mtx, W_mtx

    R_mtx = np.zeros((no_taxa-1, no_taxa-1), dtype=int) # encoding of derived lineages only
    der_mrca = dpy_tr.mrca(taxon_labels=der_ls)

    for head_node in der_mrca.preorder_iter():
        #print(head_node)
        if head_node == der_mrca: # "root" node
            #print("I am root")
            continue
        if head_node.is_leaf():
            edg_head = no_taxa - 1
        else:
            edg_head = head_node.label
            
        edg_tail = head_node.parent_node.label
        
        #print(edg_head, edg_tail)
        for j in range(edg_tail, edg_head):
            for i in range(j, edg_head):
                R_mtx[i, j] += 1        
    
    return F_mtx, W_mtx, R_mtx
