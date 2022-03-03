# apply_decision_tree
 Tool to implement a j48 decision tree from WEKA
 
Usage: `apply_decision_tree j48-Decision-Tree-FILE ARFF-File`

Where `j48-Decision-Tree-FILE` is the decisiontree from WEKA and `ARFF-File` is the arff-formatted data from the genomes that you want to test the decision tree on.

The results are outputted to standard-output (the screen) at the moment, so rediect this to a file if necessary

----------------------------------

# To install:

Download the file `apply_decision_tree.c` file and or if you have git installed use the command:
```
git clone https://github.com/ChrisCreevey/apply_decision_tree.git
```
To build apply_decision_tree type:

```
cc apply_decision_tree.c -o apply_decision_tree
```
Copy the executable "apply_decision_tree" to somewhere on your path like "~/bin" to have access to it from anywhere.

--------------------------------



