# This file is used to configure the tree made with `phylotree`.
# Every line starting with a hash (#) and empty lines are ignored.
#
# FEATURES
# You can chose two or three features to show in the tree, one feature
# per line.
# The first feature written is feature A, the second one feature B and the
# third one feature C (optional)
#
# COUNT
# phylotree can count the occurence of additional features and will
# display the total count of all features listed here (one per line)
# next to the node in the tree.
#
# COLORS
# You can define eight colors here in hexadeximal ("html"), without the
# hash (#; otherwise it would be recognized as comment). The first color
# will be used for pie chart slices that only contain feature A, the
# second color for feature B, and so on (see comments further below).
#
# SCALE
# State the method, the node size should be scale with. Be careful with
# linear as the nodes can become huge!
#
# This file is used by `phylotree`

FEATURES

protein_1	# feature A
prot2		# feature B
someother	# feature C

COUNT

# This will count the occurences of asdf, protB and GapDH in each Node
asdf
protB
GapDH

COLORS

AA0000	# A
00AA00	# B
55AAFF	# C
DDDD00	# A+B
AA00AA	# A+C
3333FF	# B+C
000000	# A+B+C (=all)
c8c8c8	# none

SCALE

log10	# Possible: log10, sqrt, linear
