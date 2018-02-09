"""Gerstein-Sonnhammer-Chothia (GSC) weighing algorithm (Gerstein, 1994).

This scripts implements a Python2 / 3 version of the GSC algorithm,
a weighing scheme for the leaves of a (phylogenetic) tree which take into
account both the length of the branches connecting each leaf to the root,
and the total number of leaves found in each region of the tree
(upweighing leaves coming from scarsely populated regions or
connected to longer branches).
The script takes advantage of the existing implementations of tree structures
and level-order traversal from the ete3 Python library.

The complete algorithm is described in the Supplementary Materials of
Gerstein M, Sonnhammer EL, Chothia C, Volume Changes in Protein Evolution (1994),
J Mol Biol., doi:10.1016/0022-2836(94)90012-4

Copyright (c) 2018 Chiara E. Cotroneo
"""
from __future__ import print_function
import ete3
import sys


def testScoring(GSCfunction):
    """Test a function implementing the GSC scoring algorithm (GSCfunction) on the tree used as an example in the paper.
    Note: the weights are assessed as before the normalisation to average=1 step described in the paper."""
    import ete3

    funName = GSCfunction.__name__

    if not callable(GSCfunction):
        sys.exit("The provided input {} is not a function.".format(funName))

    newickTree = "(D:80,(C:50,(A:20,B:20)two:30)three:30);"
    tree = ete3.Tree(newickTree, format=1)

    scores = GSCfunction(tree)

    correctScores = {'A': 43.75, 'C': 62.5, 'B': 43.75, 'D': 80.0}

    correct = 0
    for key in sorted(scores):
        print("Leaf {}: Correct score: {} - Score obtained by function {}: {}".format(key, correctScores[key], funName, scores[key]))
        if correctScores[key] == scores[key]:
            correct += 1
            
    print("{}% scores were correctly predicted ({} out of {})".format(float(correct)/len(scores)*100, correct, len(scores)))


def GSC(t, normalise=True):
    """Take an ete3 tree 't' as input, and returns a dictionary leaf ID : GSC weight with the GSC weights of each leaf of the tree."""
    scores = {}

    if not isinstance(t, ete3.Tree):
        sys.exit("The provided input is not an object of class ete3.Tree().")

    # Traverse the tree level by level, starting from the most external levels
    # up to the root
    for node in reversed(list(t.traverse("levelorder"))):

        children = node.get_children()

        # If the node 'node' is a leaf, its score is initialised to zero
        if len(children) == 0:
            scores[node.name] = 0

        else:

            for child in children:
                # If the *direct* child of node is a leaf, set its score to its
                # distance from node, without weighing
                if len(child.get_children()) == 0:
                    if child.dist != 0:
                        scores[child.name] = child.dist
                    else:
                        # In case there is a polytomy, and the distance between node and
                        # the leaf child is zero, the score is switched from zero to them
                        # smallest system float (epsilon), to prevent the weight of this leaf
                        # from being zero throughout the whole tree
                        scores[child.name] = sys.float_info.min

                # If the child of node is not a leaf, update the scores of the downstream
                # leaves by weighing them
                else:
                    # Find all the leaves downstream of child
                    species = child.get_leaf_names()

                    # Get the cumulative score of all the leaves downstream of child
                    sumWeights = float(sum([scores[name] for name in species]))

                    # For every leaf dowstream of child, update the weight as
                    # their previous score plus the distance between node and child
                    # multiplied by the score of each leaf over the cumulative score
                    # of all the leaves dowstream of child
                    for name in species:
                        # To correct for polytomies and branches of length 0:
                        if child.dist == 0:
                            scores[name] = scores[name] + sys.float_info.min * float(scores[name])/sumWeights

                        else:
                            scores[name] = scores[name] + child.dist * float(scores[name])/sumWeights

    return(scores)


def GSCnormalised(t):
    """Compute the GSC weights of an ete3 tree 't' and normalise them so that the average weight is 1, as described in the paper."""

    if not isinstance(t, ete3.Tree):
        sys.exit("The provided input is not an object of class ete3.Tree().")

    scores = GSC(t)

    totalScore = float(sum([scores[key] for key in scores]))
    nScores = len(scores)

    return({key: scores[key]/totalScore * nScores for key in scores})
