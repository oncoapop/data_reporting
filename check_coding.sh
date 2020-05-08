#!/bin/bash

# Script to see if all frameshifts lie in known coding genes

# Paths
wd="/home/dyap/ind231"
cd $wd

# Known genes
kgfile="known_genes.csv"

# query file
query="expanded_list_V4.csv"

# mutation
start=`cat $query  | awk -F, '{print $4}'`

