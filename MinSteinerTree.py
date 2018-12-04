import matplotlib.pyplot as plt # needed for printing/visualizing
import sys
import networkx as nx
import heapq as hq

def Prims(G):
    tree=nx.Graph()
    myList=[]
    hq.heapify(myList)
    tryList=list(G.nodes())
    #firstNode=choice(tryList)
    firstNode=tryList[0]
    tree.add_node(firstNode)
    numNodes=G.number_of_nodes()
    for e in G.edges(firstNode, data=True):
        hq.heappush(myList,e)
    while tree.number_of_edges() < numNodes-1:
        minEdge=hq.heappop(myList)
        (u1, u2, minWeight) = minEdge
        if u1 not in tree.nodes():
            for e in G.edges(u1, data=True):
                if e!=minEdge:
                    hq.heappush(myList, e)
        elif u2 not in tree.nodes():
            for e in G.edges(u2, data=True):
                if e!=minEdge:
                    hq.heappush(myList, e)
        else:
            continue
        tree.add_edge(minEdge[0],minEdge[1],weight=minEdge[2])
    return tree

def Steiner(G, noi):
    tree=nx.Graph()
    if len(noi)==0:
        return tree
    elif len(noi)==1:
        tree.add_node(noi[0])
        return tree
    else:
        myHeap = []
        myPaths = {}
        for i in range(len(noi) - 1):
            a=float(i)/float(len(noi)-1)*100
            print("%.1f" % a+"% completed")
            n1 = noi[i]
            for n2  in noi[i+1:]: # take all remaining nodes
                path = nx.bidirectional_dijkstra(G, n1, n2) # check if they are connected
                if path == False: # if they are not, something is wrong!
                    noi.remove(n2) # remove the second node, print an error message, and continue
                    continue # THIS SHOULD NEVER HAPPEN
                nodelist=path[1]
                distance=path[0]
                #print(nodelist)
                pair=[n1,n2]
                pair.sort()
                myPaths["%s%s", pair[0], pair[1]]=nodelist
                hq.heappush(myHeap, (distance, pair))
        while myHeap:
            myItem=hq.heappop(myHeap)
            #print(myItem)
            if myItem[1][0] not in tree or myItem[1][1] not in tree or not nx.has_path(tree, myItem[1][0], myItem[1][1]):
                tree.add_edge(myItem[1][0], myItem[1][1], weight=myItem[0])
        subgraph=nx.Graph()
        for e in tree.edges(data=True):
            pair = [e[0],e[1]]
            pair.sort()
            newList = myPaths["%s%s", pair[0], pair[1]]
            for i in range(len(newList) - 1):
                #w=G[newList[i]][newList[i+1]]
                #print(w)
                #print(newList[i], newList[i+1], G[newList[i]][newList[i+1]]['weight'])
                subgraph.add_edge(newList[i], newList[i+1], weight=G[newList[i]][newList[i+1]]['weight'])
        subgraph = Prims(subgraph)
        return finalTree(subgraph, noi)


def finalTree(graph, noi):
    myList = []
    n1 = noi[0]
    
    myList.append(n1)
    graph = remove(n1, graph, myList, noi)
    return graph

def remove(n, G, myList, noi):
    #print(list(G.neighbors(n)))
    if len(list(G.neighbors(n))) > 1:
        for v in G.neighbors(n):
            if v!=n:
                myList.append(v)
                G = remove(v, G, myList, noi)
    elif len(list(G.neighbors(n))) < 2:
        if n not in noi:
            G.remove_node(n)
    return G


def neigh(cluster, H): # unnecessary (at the moment) function that gets the neighborhood of a node or a set of nodes
    closed=[]
    for i in cluster:
        if i not in closed:
            closed.append(i)
        for j in H[i]:
            if j not in closed:
                closed.append(j)
    return closed


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


#############################
###### setting options #####
############################

#### Filenames #####
fname="../instances/"+"Arabidopsis"+".txt" # The file structure needs to be like the one shown above. You can change this to smallArabidopsis for a small demo.
ImportantFname="../instances/"+"Important"+".txt" # Important filename -- you can use smallImportant for a small demo.
DownFname="../instances/"+"Down"+".txt" # Down regulated filename -- you can use smallDown for a small demo.
UpFname="../instances/"+"Up"+".txt" # Up regulated filename -- you can use smallUp for a small demo.

##### printing
printing=False # Setting to True will print a test.png file in the images folder. Can be used to verify that the program is working.

#################
##### MAIN #####
################
ImportantList=[]
DownList=[]
UpList=[]

readFile=open(ImportantFname, 'r')
for line in readFile:
    line = line.rstrip()
    ImportantList.append("3702."+str(line))

readFile=open(DownFname, 'r')
for line in readFile:
    line = line.rstrip()
    DownList.append("3702."+str(line))

readFile=open(UpFname, 'r')
for line in readFile:
    line = line.rstrip()
    UpList.append("3702."+str(line))

#####################################
##### READING ARABIDOPSIS FILE #####
####################################

G=nx.read_weighted_edgelist(fname)

for (i,j,d) in G.edges(data=True):
    u=int(d['weight'])
    d['weight']=1000-u

print("Done reading network")

allComponents=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)
H=allComponents[0]


nblue=[]
nred=[]
ngreen=[]
noi=[]
for i in H.nodes():
    k=str(i)
    H.node[i]['color']='yellow'
    if k in ImportantList:
        nblue.append(i)
        noi.append(i)
    if k in DownList:
        nred.append(i)
        noi.append(i)
    if k in UpList:
        ngreen.append(i)
        noi.append(i)

print("Finding minimum cost Steiner tree..")
tree=Steiner(H, noi)
print("Writing solution..")

myLabels={}
for i in H.nodes():
    myLabels[i]=""
for i in tree.nodes():
    myLabels[i]=str(i)
    myLabels[i]=myLabels[i].replace("3702.","")

with open("_noi.txt", 'w') as outfile:
    for i in nblue:
        writeLine=str(i)+"\n"
        outfile.write(writeLine)
        outfile.flush()

with open("_ndown.txt", 'w') as outfile:
    for i in nred:
        writeLine=str(i)+"\n"
        outfile.write(writeLine)
        outfile.flush()

with open("_nup.txt", 'w') as outfile:
    for i in ngreen:
        writeLine=str(i)+"\n"
        outfile.write(writeLine)
        outfile.flush()

with open("_edges.txt", 'w') as outfile:
    for (i,j) in tree.edges(data=False):
        writeLine=str(i)+" "+str(j)+"\n"
        outfile.write(writeLine)
        outfile.flush()

if printing:
    pos=nx.spring_layout(H)
    nx.draw(H,pos)
    nx.draw_networkx_labels(H, pos, myLabels, font_size=6, font_family='sans-serif')
    nx.draw_networkx_nodes(H, pos, nodelist=tree.nodes(), node_color='y') # yellow for all other nodes
    nx.draw_networkx_nodes(H,pos,nodelist=nblue,node_color='b') # blue for the nodes in the "important" file
    nx.draw_networkx_nodes(H,pos,nodelist=nred,node_color='r') # red for the nodes in the "down-regulated" file
    nx.draw_networkx_nodes(H,pos,nodelist=ngreen,node_color='g') # green for the nodes in the "up-regulated" file
    nx.draw_networkx_edges(H,pos,edgelist=tree.edges())
    plt.axis('off')
    #plt.show()
    plt.savefig("../images/_test.png", bbox_inches='tight', dpi=300)

print("Terminated.")

