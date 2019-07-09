=======
 Notes for Overall `parsealignment` Structure
========

Goal
----

Take an alignment of sequences to targets (from `alignparse.target`) and
return information about specified features. 

Inputs
----

Alignment
Features to extract
    Have default return for each feature (say seq or +/-), but also allow
    user to specify what to return


Outputs
----

Features with specified attributes
    Probably as a dataframe

Method
----

Probably have the alignment be a class?
Then parse alignment based on features


Issues to keep in mind
----

Junctions between features could get tricky
Dealing with ambiguous nucleotides
Indels/imperfect alignments?
