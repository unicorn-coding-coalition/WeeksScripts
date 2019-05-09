#!/usr/bin/env python

# This program will plot distances on a PDB structure from
# a list of pairs read from a file, in the correlation format

#Created by Tom Christy: 9/20/2016
#Edited by Nicole Lama: 2/28/2019

from pymol import cmd,stored
from sys import argv
import numpy as np

pdb = 'D:/Documents2/Rotations/Weeks/RNAFiles/RNaseP/rnasePSupMod.pdb'
corrF = "D:/Weeks/Data/ringRnaseP/rnasep.mi.corr.filt.txt"
		
cmd.reinitialize()
cmd.load(pdb)
#cmd.show('cartoon', 'name p')
cmd.hide('lines')
cmd.color('slate')
cmd.color('black', 'elem O')
pairs = []
pairs = open(corrF, 'r')
pairlines = pairs.readlines()
pairs.close()
pairs= [[int(i.split()[0]), int(i.split()[1])] for i in pairlines]
prefix = "Deletions_"
cmd.select(prefix, "none")
selected = []
stored.list=[]
cmd.iterate("(name O2*)","stored.list.append(resi)")
residues = [int(i) for i in stored.list]
for i in pairs:
	j = i[0]
	k = i[1]
	if j in residues and k in residues:
		already_printed= [] 
		if [j,k] in already_printed:
			print('printed again')
		already_printed.append([j,k])
		for at1 in cmd.index("name O2* in resi "+str(j)):
			for at2 in cmd.index("name O2* in resi "+str(k)):
				dist = cmd.distance(prefix+str(j), "%s`%d"%at1, "%s`%d"%at2)
				if dist < 15:
					cmd.color("firebrick", prefix+str(j))
				elif dist < 20:
					cmd.color("firebrick", prefix+str(j))
				elif dist < 30:
					cmd.color("tv_red", prefix+str(j))
				elif dist < 40:
					cmd.color("deepsalmon", prefix+str(j))
				elif dist > 40:
					cmd.color("salmon", prefix+str(j))
				cmd.set('dash_width',10)
			# create selections containing relevant nts
		if j not in selected:
			selected.append(j)
			cmd.select(prefix, prefix+" + resi "+str(j))
		if k not in selected:
			selected.append(k)
			cmd.select(prefix, prefix+" + resi "+str(k))
		
cmd.show("spheres", "name O2*\' in Link_*")
cmd.set('opaque_background', 'off')
cmd.center()
