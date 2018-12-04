# PPIN
Protein-Protein interaction network analysis


##########################################################################################
####################################### HOW TO RUN #######################################
##########################################################################################

(-2) Make sure you have a folder structure that looks like the following:

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

Some of the files might look different. For example, you might have a different organism 
text file. Or, you might have chosen to rename your "important", "up", and "down" files 
differently. This is fine. 

You may also be missing files _edges.txt, _noi.txt, _ndown.txt, _nup.txt, pathways.csv. 
These files are automatically generated and populated the first time you run 
MinSteinerTree.py (for the first four files) or ReloadFigure.py/Treeify.py (for pathways).

(-1) Download the organism you are interested in from STRING-DB (https://string-db.org).
Make sure to store the produced text file in the instances folder. In addition to that, 
create files that contain your Down, Up, Important nodes. You can use different names.

(0) Open the python codes and make the following changes. 

For MinSteinerTree.py
=========================================================================================
Update the filename to reflect the locations/names of your files.
Update organism_code to reflect the organism code of your choice.
Select printing=True/False to print figures during execution or not.

For ReloadFigure.py
=========================================================================================
Update organism_code to reflect the organism code of your choice.

For Treeify.py
=========================================================================================
Update organism_code to reflect the organism code of your choice.

(1) Run

python MinSteinerTree.py 

This should produce upon successful execution four files within the src folder, as well
as a test image in the images folder. It will also print a figure, if printing is enabled.

(2) Run

python ReloadFigure.py image_filename [args]

where image_filename is a name of your choice (e.g., MyExperiment) and [args] are optional
arguments including -all or -up -down which decide the nodes that are labeled. 

Example: python ReloadFigure.py myExperiment -up -down

This will print out three files in the images folder (called myExperiment1, myExperiment2, 
and myExperiment3) that contain labels for both the up and down regulated nodes (but not 
the rest).

(3) Run 

python Treeify.py image_filename node_name [args]

where image_filename is a name of your choice (e.g., MyExperiment), node_name is the node
you are interested in using as the root of the tree, and [args] are optional arguments 
including -all or -up -down which decide the nodes that are labeled. 


Example: python Treeify.py myExperiment AT1G62360.1 -all

This will use myExperiment as the base name for image files produced (in the images 
folder), AT1G62360.1 is the node we are looking for, and -all signals that all nodes will 
be labeled.



