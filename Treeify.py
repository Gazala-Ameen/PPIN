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


# The code to generate the tree using networkx was taken from: https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
def hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5,
                  pos = None, parent = None):
    '''If there is a cycle that is reachable from root, then this will see infinite recursion.
        G: the graph
        root: the root node of current branch
        width: horizontal space allocated for this branch - avoids overlap with other branches
        vert_gap: gap between levels of hierarchy
        vert_loc: vertical location of root
        xcenter: horizontal location of root
        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch.'''
    if pos == None:
        pos={}
    pos[root] = (xcenter, vert_loc)
    neighbors = list(G.neighbors(root))
    if parent != None:
        neighbors.remove(parent)
    if len(neighbors)!=0:
        dx = width/len(neighbors)
        nextx = xcenter - width/2 - dx/2
        for neighbor in neighbors:
            nextx += dx
            pos = hierarchy_pos(G, neighbor, width = dx, vert_gap = vert_gap,
                                vert_loc = vert_loc-vert_gap, xcenter=nextx, pos=pos,
                                parent = root)
    return pos

# The code can be run as: python Treeify.py testing AT1G62360.1 -all
# where testing will be the name of the image files produced (in the images folder)
# AT1G62360.1 is the node we are looking for
# and -all signals that all nodes will be labeled

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

gene="3702."+sys.argv[2]

options=[]
for i in range(3,(len(sys.argv))):
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
found=False
for i in H.nodes():
    nList.append(i)
    if str(i)==gene:
        found=True
        geneNode=i




if found:
    pos = hierarchy_pos(H,geneNode)
    plt.figure(1,figsize=(12,12))
    nx.draw(H, pos=pos, with_labels=False)
    nx.draw_networkx_nodes(H,pos,nodelist=nList,node_color='b')
    nx.draw_networkx_nodes(H,pos,nodelist=nred,node_color='r')
    nx.draw_networkx_nodes(H,pos,nodelist=ngreen,node_color='g')
    nx.draw_networkx_nodes(H,pos,nodelist=nblue,node_color='y')
    nx.draw_networkx_labels(H, pos, myLabels, font_size=12, font_family='sans-serif')
    outName="../images/"+sys.argv[1]+"_Tree_"+gene+".pdf"
    plt.savefig(outName, bbox_inches='tight', dpi=300)
    #
    plt.show()
    
    newH=nx.Graph(H)
    newH.remove_node(geneNode)
    allComponents=sorted(nx.connected_component_subgraphs(newH), key = len, reverse=True)
    cnt=1
    for C in allComponents:
        Cnred=[]
        Cngreen=[]
        Cnblue=[]
        for i in C.nodes():
            if str(i) in ndown:
                Cnred.append(i)
            if str(i) in nup:
                Cngreen.append(i)
            if str(i) in noi:
                Cnblue.append(i)
        CmyLabels={}
        for i in C.nodes():
            CmyLabels[i]=str(i)
            CmyLabels[i]=CmyLabels[i].replace("3702.","")

        #print(nred)
        pos=nx.spectral_layout(C, scale=1.0)
        plt.figure(1,figsize=(12,12))
        nx.draw(C,pos)
        nx.draw_networkx_labels(C, pos, CmyLabels, font_size=12, font_family='sans-serif')
        #nx.draw_networkx_edges(H, pos, edgelist=H.edges())
        nx.draw_networkx_nodes(C,pos,nodelist=C.nodes(),node_color='b')
        nx.draw_networkx_nodes(C,pos,nodelist=Cnred,node_color='r')
        nx.draw_networkx_nodes(C,pos,nodelist=Cngreen,node_color='g')
        nx.draw_networkx_nodes(C,pos,nodelist=Cnblue,node_color='y')
        #nx.draw_networkx_labels(C,pos, myLabels, font_size=12, font_family='sans-serif')
        outName="../images/"+sys.argv[1]+"_Subtrees"+str(cnt)+".pdf"
        plt.savefig(outName, bbox_inches='tight', dpi=300)
        #
        #plt.show()
        cnt=cnt+1

else:
    print(sys.argv[2]+" not Found")


