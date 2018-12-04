import networkx as nx
import matplotlib.pyplot as plt
import sys
import itertools

#############################
###### file structure ######
############################

#+-- images
#+-- instances
#|   +-- Arabidopsis.txt
#|   +-- smallArabidopsis.txt
#|   +-- Important.txt
#|   +-- smallImportant.txt
#|   +-- Down.txt
#|   +-- smallDown.txt
#|   +-- Up.txt
#|   +-- smallUp.txt
#+-- src
#|   +-- MinSteinerTree.py
#|   +-- ReloadFigure.py
#|   +-- Treeify.py
#|   +-- _edges.txt
#|   +-- _noi.txt
#|   +-- _ndown.txt
#|   +-- _nup.txt
#|   +-- pathways.csv

# This python code requires two arguments: 1) a name for the images produced (they will be stored in the images folder) and 2) an option: -all will print all labels, -imp will print labels for the important nodes only, -down, -up will print labels for the down and up regulated nodes only (respectively), and -not will print labels for all other nodes (other than the important and up/down regulated ones).


# This code requires the produced files below (so the minimum Steiner tree code must have been executed before). It only serves to produce a .csv file (called pathways.csv) which contains all paths from any two leafs in the tree. It also produes three network images for visualization purposes.

# For example python ReloadFigure.py myImage -up will print out three files in the images folder (called myImage1, myImage2, and myImage3) that contain labels for only the up regulated nodes.

# Similarly, python ReloadFigure.py myImage -up -down will print out three files in the images folder (called myImage1, myImage2, and myImage3) that contain labels for both the up and down regulated nodes (but not the rest).



edgeName="_edges.txt"
noiName="_noi.txt"
ndownName="_ndown.txt"
nupName="_nup.txt"


H=nx.read_edgelist(edgeName)

paths=dict(nx.all_pairs_shortest_path(H))
deg=nx.degree_centrality(H)

myPathways=open("pathways.csv", 'w')

leafs=[i for i in H.nodes() if deg[i]*(H.number_of_nodes()-1)==1]

#print(len(leafs))

options=[]
for i in range(2,(len(sys.argv))):
    options.append(sys.argv[i])


for (i,j) in itertools.combinations(leafs, 2):
    for l in paths[i][j]:
        myPathways.write(str(l)+",")
    myPathways.write("\n")
    myPathways.flush()

noi=[]
with open(noiName, 'r') as f:
    for line in f:
        line = line.rstrip()
        noi.append(line)
nblue=[]
for i in H.nodes():
    if str(i) in noi:
        nblue.append(i)
#print(nblue)
ndown=[]
with open(ndownName, 'r') as f:
    for line in f:
        line = line.rstrip()
        ndown.append(line)
nred=[]
for i in H.nodes():
    if str(i) in ndown:
        nred.append(i)
#print(nred)
nup=[]
with open(nupName, 'r') as f:
    for line in f:
        line = line.rstrip()
        nup.append(line)
ngreen=[]
for i in H.nodes():
    if str(i) in nup:
        ngreen.append(i)
#print(ngreen)

myLabels={}
if "-all" in options:
    for i in H.nodes():
        myLabels[i]=str(i)
        myLabels[i]=myLabels[i].replace("3702.","")
else:
    if "-imp" in options:
        for i in nblue:
            myLabels[i]=str(i)
            myLabels[i]=myLabels[i].replace("3702.","")
    if "-down" in options:
        for i in nred:
            myLabels[i]=str(i)
            myLabels[i]=myLabels[i].replace("3702.","")
    if "-up" in options:
        for i in ngreen:
            myLabels[i]=str(i)
            myLabels[i]=myLabels[i].replace("3702.","")
    if "-not" in options:
        for i in H.nodes():
            if i not in nblue and i not in ngreen and i not in nred:
                myLabels[i]=str(i)
                myLabels[i]=myLabels[i].replace("3702.","")

nList=[]
for i in H.nodes():
    nList.append(i)

pos=nx.spring_layout(H, scale=1.0)
posLabels={}
for i in H.nodes():
    posLabels[i]=pos[i]+(0,0.03)
plt.figure(1,figsize=(12,12))
nx.draw(H,pos)
nx.draw_networkx_labels(H, posLabels, myLabels, font_size=12, font_family='sans-serif')
#nx.draw_networkx_edges(H, pos, edgelist=H.edges())
nx.draw_networkx_nodes(H,pos,nodelist=nList,node_color='b')
nx.draw_networkx_nodes(H,pos,nodelist=nred,node_color='r')
nx.draw_networkx_nodes(H,pos,nodelist=ngreen,node_color='g')
nx.draw_networkx_nodes(H,pos,nodelist=nblue,node_color='y')
outName="../images/"+sys.argv[1]+"_1.pdf"
plt.savefig(outName, bbox_inches='tight', dpi=300)
#nx.draw_networkx_labels(H, pos, myLabels, font_size=12, font_family='sans-serif')
plt.show()

pos=nx.shell_layout(H, scale=1.0)
posLabels={}
for i in H.nodes():
    posLabels[i]=pos[i]+(0,0.03)
plt.figure(1,figsize=(12,12))
nx.draw(H,pos)
nx.draw_networkx_labels(H, posLabels, myLabels, font_size=12, font_family='sans-serif')
#nx.draw_networkx_edges(H, pos, edgelist=H.edges())
nx.draw_networkx_nodes(H,pos,nodelist=nList,node_color='b')
nx.draw_networkx_nodes(H,pos,nodelist=nred,node_color='r')
nx.draw_networkx_nodes(H,pos,nodelist=ngreen,node_color='g')
nx.draw_networkx_nodes(H,pos,nodelist=nblue,node_color='y')
outName="../images/"+sys.argv[1]+"_2.pdf"
plt.savefig(outName, bbox_inches='tight', dpi=300)
plt.show()

pos=nx.spectral_layout(H, scale=1.0)
posLabels={}
for i in H.nodes():
    posLabels[i]=pos[i]+(0,0.03)
plt.figure(1,figsize=(12,12))
nx.draw(H,pos)
nx.draw_networkx_labels(H, posLabels, myLabels, font_size=12, font_family='sans-serif')
#nx.draw_networkx_edges(H, pos, edgelist=H.edges())
nx.draw_networkx_nodes(H,pos,nodelist=nList,node_color='b')
nx.draw_networkx_nodes(H,pos,nodelist=nred,node_color='r')
nx.draw_networkx_nodes(H,pos,nodelist=ngreen,node_color='g')
nx.draw_networkx_nodes(H,pos,nodelist=nblue,node_color='y')
outName="../images/"+sys.argv[1]+"_3.pdf"
plt.savefig(outName, bbox_inches='tight', dpi=300)
plt.show()
