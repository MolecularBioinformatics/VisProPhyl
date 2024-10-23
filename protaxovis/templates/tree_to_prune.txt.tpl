# In this file, you can define which nodes are to be truncated or deleted.
# Empty lines and everything in a line after a hash (#) is ignored.
#
# There are some keywords to switch the modes of action. To switch the
# mode of action, a line must start with a ! directly followed by the
# keyword with nothing else in the line. Possible keywords are:
#
# !pruneafter
# All children of the given node are pruned. That is, the node is kept,
# but all child nodes are detached.
#
# !pruneafterchildren
# The children of the given node are kept but all of their children are
# pruned. In other words, all grandchildren are detached.
#
# !dontprune
# All nodes given here are *not* pruned or detached. However, if a
# higher-level node is detached, this one is deleted as well.
# This can be used to keep a taxon that would otherwise be removed, e.g.
# by !pruneafter.
#
# !detach
# These nodes are directly detached.
#
# !delete
# These nodes are deleted but all children are attached to the parent of
# this node. All !delete instructions will be performed after all other
# instructions were done.
#
# The affected nodes themselves must be given one in each line with
# either the name, the taxid or the combination name^taxid. The name is
# sufficient in most cases, but certain names are ambiguous. When there
# are two nodes with the same name, both will be affected. In these
# cases, write either taxid or name^taxid.
#
# Examples for ambiguous names are (by far not complete):
# - Aphelia (A moth and a plant)
# - Bacteria (!) (The kingdom and an insect genus)
# - Yersinia (Enterobacteria and mantids (insects))
#
# This file is used by `phylotree`.

!pruneafter
Bacteria^2
Archaea^2157

!detach
other_entries^2787854
