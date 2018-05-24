# Python-GSC

This scripts implements a Python2 / 3 version of the Gerstein-Sonnhammer-Chothia algorithm (GSC, Gerstein 1994), a weighing scheme for the leaves of a (phylogenetic) tree. The resulting weights are based on the total length of the branches connecting it to the most internal node of the tree, and on the total number of leaves found in each region of the tree. The final aim is to upweigh those leaves located in scarsely populated regions of the tree or connected to longer branches, and downweigh the others.

The weight of each leaf is computed starting from the most external leaf of the tree and sequentially updated by visiting each level of the tree until the most internal node ("root") is reached. When a leaf *i* is first reached, its weight *w<sub>i</sub>* is initialised to zero. When each internal node *j* is reached, the weight *w<sub>i</sub>* of the *i = 1, ..., s* leaves connected to node *j* are updated according to two possible strategies:
(i) if leaf *i* is a direct child of node *j*,  *w<sub>i</sub>* is updated to *w<sub>i</sub>* + *d<sub>i, j</sub>*, where *d<sub>i, j</sub>* is the length of the branch connecting *j* and *i*;
(ii) if leaf *i* is not a direct child of node *j*, *w<sub>i</sub>* is updated to *w<sub>i</sub>* + *d<sub>k, j</sub>* \* *w<sub>i</sub>* / sum<sub>k</sub>(*w<sub>k</sub>*), where *d<sub>k, j</sub>* is the distance between node *j* and its direct child *k* that is a parent of *i*, and sum<sub>k</sub>(*w<sub>k</sub>*) is the sum of the current weights of all the leaves connected to *k*.

The following example (taken from the paper) illustrates how the algorithm works for a simple ultrametric tree. Here, A, B, C, D are the leaves of the tree; (1), (2), (3) are the internal nodes; and the numbers displayed on the branches of the tree correspond to their length.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/-------------------80------------------&nbsp;D<br>
--&nbsp;(3) |&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/-----------50----------&nbsp;C<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\---30---(2)|<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/--20--&nbsp;A<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\---30---(1)|<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\--20--&nbsp;B<br>

|          | *w<sub>A</sub>*                                | *w<sub>B</sub>*                              | *w<sub>C</sub>*                                 | *w<sub>D</sub>*  |
|-------------------------|----------------------------------|---------------------------------|-----------------------------------|----|
| Initial| 0                                | 0                               | 0                                 | 0  |
| At (1)            | 0 + 20                               | 0 + 20                              | 0                                 | 0  |
| At (2)            | 20 + 30 * 20/40 = 35                  | 20 + 30 * 20/40 = 35                 | 50                                | 0  |
| At (3)            | 35 + 30 * 35 / (35 + 35 + 50) = 43.75  | 35 + 30 * 35 / (35 + 35 + 50) = 43.75 | 50 + 30  * 50 / (35 + 35 + 50)  = 62.5 | 80 |
| Final   | 43.75                            | 43.75                           | 62.5                              | 80 |



Additional details are provided in the original description of the algorithm, which can be found in the Supplementary Materials of *Gerstein M, Sonnhammer EL, Chothia C, Volume Changes in Protein Evolution (1994), J Mol Biol.,* <a href="doi:10.1016/0022-2836(94)90012-4">doi:10.1016/0022-2836(94)90012-4</a>.
