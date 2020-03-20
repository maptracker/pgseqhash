### PG SeqHash Oligo Aligner

An oligo-to-genome alignment algorithm capable of exhaustively finding
imperfect matches in a large genomic database even with very short
oligos. Designed to aid prediction of potential off-target activity in
Anti-Sense Oligos (ASOs), where biological activity had been noted
even with two mismatches in an oligo aligning to an intronic region.

The search database is a Postgres DB centered around
polymorphism-permuted 12-mer oligos. The sequence source is the set of
all "pre-mRNA" sequences for several genomes. One sequence is used for
each gene, representing the full genomic sequence spanning from the
"left-most" edge of all transcripts for that gene to the
"right-most". This allows both intronic and exonic sequence to be
considered. In addition, splice boundaries and ORFs for each variant
are also recorded, which allows hits to be assigned to specific
transcripts, as well as allowing classification of hits into exonic,
CDS, splice acceptor, UTR etc.

![DB algorithm][DBalg]

Multiple genomic builds are included in the database. In our
production, these represented a couple human builds (GRCh37 plus
RefSeq Gene), several primates, dog, mouse and rat.

When an oligo (ASO) is submitted for searching, a minimum tiling of
12-mer words is generated for it. The user also specifies the number
of mismatches they wish to allow, the default being 2. Each word is
then fully permuted for all possible mismatch combinations. This
permuted set is then searched in parallel (this is Perl, so forked)
against the database. For speed, some of the permutations will be
discarded if prior, neighboring oligo hits preclude the possibility of
finding a full-oligo match within the requested mismatch limit.

![Search algorithm][SearchAlg]

Results are then provided in three forms:

* Multi-sheet Excel workbook
* Interactive DHTML tables
* Interactive DHTML genomic browser

An example report from the interactive table browser (showing
published oligos from the literature):

![Hit report][Report]

* [BMS Public Disclosure approval](PubD-Disclosure-Approval.md)

[DBalg]: img/SOS-DatabaseAlgorithm.png
[SearchAlg]: img/SOS-SearchAlgorithm.png
[Report]: img/SOS-Overview.png
