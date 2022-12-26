# Path association rule mining

This is an implementation of path association rule mining.

The path association rule mining finds regulations of path patterns in a large single graph.

## setup

CMake 3.7 or later

google [densehash](https://github.com/sparsehash/sparsehash)


## dataset
download [Nell and DBpedia](https://github.com/GemsLab/KGist)

download [Pokec](https://snap.stanford.edu/data/soc-pokec.html)

Datasets should include two files:
- edge: a list of triples, two vertex with edge label.
- vertex: a list of attributes associated with vertices

edge
```
v1 l1 u1
v2 p2 u2
...
```

vertex
```
v1 a1 a2 ...
v2 a1 a2 a3 ...
...
```


## run code

1. cmake CMakeLists.txt
2. make
3. ./pathassociation -i {dataset file}

## arguments
Our code has the following parameters. All arguments except for i are optional.

- i X: set X as the input graph
- o X: set X as output file name, default ./result/test
- minsup X: set X as the threshold of relative support, default 0.001
- aminsup X: sex X as the threshold of absolute support
- k X: set X as the maximum length of paths, default 2
- p X: set X as the number of threads, default 1
- pmax: use the maximum number of threads
- nokleene: do not use reachability path patterns, default false
- cr X: set X as the candidate reduction rate, default 1.0 (i.e., exact method without approximation)
- sr X: set X as the sampling rate, default 1.0 (i.e., exact method without approximation)
- es: allow path patterns that attributes of sources are empty, default false
- do: allow rules that are their anticients (resp. consequents) are dominated by their consequents (resp. anticients)
- mink X: set X as the minimum length of paths in rules, default 1
- pmode X: set X as a partition mode (0: balance # of vertices, 1: balance cost, 2: balance cost and # of vertices), default 0
- nooutput: do not output rules, default false
- attributeoutput: output the ids to attribute and label names, default false
- non: use a baseline method, default false
- noextend: do not use optimization methods, default false

example
```
./pathassociation -i ../data/pokec -aminsup 10000 -k 2 -p 32 -sr 0.1 -cr 0.1
```
