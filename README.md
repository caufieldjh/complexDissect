# complexDissect.py
**A small utility for comparing sets of protein complexes.**

**Secondary goal**: predict conservation of complexes across species using
protein orthology.

**Tertiary goals**: 
Allow exploration of specific protein complexes to 
retreive details like best match, size, and functional details.

Search a set of complexes to see which of them have appeared to evolve across a
given set of genomes (as taxids in a pre-defined taxonomy).

Requires several Gb of disk space to store mapping files.

**INPUT**: 

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

"Short" format files can also be parsed. These have the following format:
ComplexID	Proteins
A	ID1	ID2
B	ID3	ID4
C	ID4

Each line contains a complex name/identifier, followed by each of the
components of the complex, separated by tabs.
So if the data looks like this (e.g., Hu et al 2009, table S7):

CplxID	Complex_members
1	b1234	b1235	b1236

then it can be converted into "long" format.

**OUTPUT**:

A file indicating, for each "experimental" complex, conservation
in the model, both within a complex and across the full set.
For example, proteins A and B may both be present in the same complex
in the model set but protein C, though in that complex in the
experimental set, may be in a completely different complex.

Uses ecoli.txt ID conversion file provided by Uniprot/Swiss-Prot.

Output is intended to be used to produce heatmap as seen in Caufield et al. 2015 PLoS Comp Bio.
