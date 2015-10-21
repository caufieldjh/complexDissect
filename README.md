# complexDissect.py
A small utility for comparing sets of protein complexes.

INPUT: 
Two files. Assumed to have names like "complexes*.txt" 
but will take input for names
Each needs to have the following format.
ProteinID	ComplexMembership
ID1	A
ID2	A
ID3	B
ID4	B
ID4	C

One of these files is for the "experimental" set.
The second file is the model set.
Identifiers should be Uniprot IDs but don't have to be as long as both
files use the same type of unique identifier (i.e., not gene names).
IDs will be converted to Uniprot IDs.
Proteins may belong to more than one complex, as shown above, so
interactors may not be unique but lines will be.

OUTPUT:
A file indicating, for each "experimental" complex, conservation
in the model, both within a complex and across the full set.
For example, proteins A and B may both be present in the same complex
in the model set but protein C, though in that complex in the
experimental set, may be in a completely different complex.

Uses ecoli.txt ID conversion file provided by Uniprot/Swiss-Prot.

Will also break down complexes into lists of the format shown above,
so if the data looks like this (e.g., Hu et al 2009, table S7):

CplxID	Complex_members
1	b1234	b1235	b1236

then it can be converted into "long" format.
