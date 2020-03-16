# PgSeqHash - Operating Principles

_PgSeqHash_, also known as _Simple Oligo Search_ or _SOS_, is designed
to detect and report imperfect oligo matches with the following
features and conditions:

* A query (the oligo) is searched against a database (transcipt sequences)
* All matches with less than or equal to `mm` mismatches are
  guaranteed to be found, where the default is `mm = 2`
* All matches (alignments) are __ungapped__
* The minimum oligo size is 12bp
* The default "search space" is presumed to be the combination of both
  pre-mRNA and mRNA transcriptional sequences
* Exonic structure and translation start/end coordinates are used to
  clasify hits (intron, splice acceptor, UTR, exon-exon boundary, etc)

The system is implemented in Perl with a Postgres relational database
that holds both the search databases and caches previously run search
results. This document will not describe the technical details of
these components, but rather the logic they follow to set up the
search database and run the searches.

## Preparing the search database

#### Sequence preparation

To create a search database, the following are needed:

1. __mRNA sequences__. This is the only requirement, the other
   components listed below are optional.
2. __Genomic sequences__, corresponding to the pre-mRNA "footprint" of
   each mRNA. If not provided, then intronic and splice junction hits
   can not be detected or classified.
3. __Coding Sequence boundaries__ - that is, the start and end
   coordinates of the CDS/ORF on protein-coding transcripts. Without
   this information, the program can not differentiate between 5' UTR,
   CDS or 3' UTR regions.

These information can be gleaned from a variety of sources:

1. The [RefSeq Gene][RSG] intiative provides rich information about
   many human transcripts that covers all three sources of
   information. It is implemented as an "alternative" build labeled
   __RSG__ in our system.
2. A [VCF file][VCF] can provide detailed transcript genomic
   coordinates, and can be combined with a genomic fasta sequence file
   to extract underlying mRNA and pre-mRNA sequences. This approach is
   feasible and has been discussed, but not implemented.
3. An mRNA sequence file can be aligned to a genomic sequence file,
   and the resulting alignments used to define pre-mRNA sequences and
   exon coordinates. Aside from RSG, this is how all databases at BMS
   have been built, using a tool called "Genomic Sim4"

A preliminary fasta file is constructed of the provided sequence (this
file is preserved for use in FriendlyBlast, our implementation of
BLAST). The mRNA sequence is redundant to (entirely a subset of) the
genomic sequence, with one exception: Unique sequence at the location
of exon-exon boundaries. The fasta file will then consist of the
following sequences:

1. __pre-mRNA footprints__ - each sequence is effectively a full,
   uninterupted genomic stretch for a single gene, covering all exonic
   and intronic regions
2. __Orphan mRNAs__ - These are "pure" (no intronic sequence) mRNAs for
   transcripts that failed to find a genomic alignment. Rare in human,
   these are unfortunately more common in species with incomplete
   genome assemblies or "noisy" transcriptomes.
3. __Exon-Exon fragments__ - These are short "stub" sequences (30bp on
   each side) capturing mRNA sequence that crosses an exon-exon
   boundary. In some cases, these sequences will include "orphan
   exons" - situations where the mRNA could be aligned to the genome,
   with the exception of one or more exons which failed to align.

For computational efficiency, overlapping genes will be combined into
a single sequence entry. The description line of each entry will
include any original description, plus metadata that describe the
relative [start, end, strand] of the gene(s) covered by the
sequence. The following snippet represents examples of each class of
sequence type taken from `RefSeq RNA vs GRCh37 genome`.

The first entry is a pre-mRNA segment covering CAV3 on the forward
strand and SSUH2 on the reverse.

The second, third and fourth examples are exon-exon fragments: One
shared by CAV3 transcripts NM_001234 and NM_033337 (both following
exon 1) another other found only in CAV3 variant NM_001234 (after exon
2), and a third from OSBPL9 with 51 ("W=53") internal exonic bases
that failed to be aligned to the genome.

Finally, an orphan mRNA, XM_006710173.1, is included following its
failure to align to GRCh37.

```
>3.GRCh37:8661086-8788448.FR [CAV3 SSUH2] /LOC51066="[SSUH2] 1-125641 R 100.0% 8 RNAs. ssu-2 homolog (C. elegans)" /LOC859="[CAV3] 114401-127363 F 99.5% 2 RNAs. Caveolin 3" /len='127363'
TTGCACTTTATGCCTAAGTGCAGTATTTATTGTGCATTAGGGCCAGATATTTGTGTGTATCCTTGAACATTTCAGTAATACTCAGATCGTCAAGAGGCTC
... etc etc ...

>3.GRCh37:8775676^8787212.F [CAV3] P=30 ON=3.GRCh37:8661086-8788448.FR NM_001234=1 NM_033337=1
CCCAAGAACATTAACGAGGACATAGTCAAggTGGATTTTGAAGACGTGATCGCAGAGCCT
>3.GRCh37:8787559^8787656.F [CAV3] P=30 ON=3.GRCh37:8661086-8788448.FR NM_001234=2
GTGGTGCTGCGGAAGGAGGTCTAAAGCCAggGACTGCTCCATACCCCATGATGGAGCACA
>1.GRCh37:52082893^52179675.F [OSBPL9] P=30 W=53 ON=1.GRCh37:52082546-52344609.FR XM_006710324=1 XM_006710325=1
TACAATGCAGGACTGCTCTCCTACTACACgtccaaggacaaaatgatgagaggctctcgcagaggatgtgttagactcagacCCGTGATGCTGATGAGCG
AGAGAAGTGGA

>XM_006710173.1 PREDICTED: Homo sapiens glycine-rich cell wall structural protein 1.0-like (LOC101930481), transcript variant X2 /gene='LOC101930481' /orphan='reported absent from Homo sapiens Genome Build GRCh37' /taxa='Homo sapiens'
gagaatagatattttctcatttttgcttgtctgattttccaagtttaaaaaataattcacatgtacagcttctacaaagggtatcaaaaatcatccatga
... etc etc ...
```

Not shown in the descriptions are the finer details of exon structure
coordinates and CDS boundaries. That information is preserved while
the search database is constructed.

#### "Word" extraction and RDBS loading

The fasta file is now broken into identical-length overlapping
"words". The default length (which must be consistent in a given
database) is 12bp, chosen to be modestly unique in a genome but also
short enough to accomodate small oligos. Oligos smaller than this
value can NOT be searched against the database. All words on the
forward strand of each sequence are considered eg:

```
acggcgtgtggaaggtgagctacacca  <- Sequence in fasta file
acggcgtgtgga   <- First word
 cggcgtgtggaa   <- Second word
  ggcgtgtggaag   etc
```

Special considerations:

* The reverse strand is not considered - both strands will be searched
  by generating reverse complement queries from the oligo.
* Sequence case is ignored. That is, "a" = "A"
* The presence of any _masking_ characters (at the moment only "x") will
  cause the word to be ignored (not entered into the database)
* _Ambiguity_ characters are expanded. For example, the nucleotide "R"
  indicates "A or G" while Y indicates "C or T" if a word
  `acgRcgtgYgga` was encountered, four entries would be made in the
  database for that one unique genomic coordiate, corresponding to the
  words `acgAcgtgCgga`, `acgAcgtgTgga`, `acgGcgtgCgga` and
  `acgGcgtgTgga`
  * This expansion has an upper limit (default 50 words) to prevent
    "explosions" that might occur if a word has many ambiguous
    characters (pure N-stretches are the most severe, with 4^12 =
    16,777,216 words being possible, but some immunologic loci have
    significant expansion issues as well due to exceedingly high
    variant density).

Aside from the special exclusions mentioned above (masked or exploding
ambiguous words) every word is recorded in the database with the
following information:

* `GenomeBuild` (eg "GRCh37" or "MacFas5"), the public token
  associated with the genome assembly
* `SubjectID` (eg "1.GRCh37:154929503-154946959.FR" or "NM_005960.1"),
  an identifier for the sequence the word was extracted from
* `SubjectPos`, the position of the word, recorded at the first base
  and counting from 1 on the SubjectID
* _Note: Internally, the example values above are what would be
  reported to a human user, but are represented internally in the
  database by normalized integer IDs)_

#### Recording gene structure and metadata

In addition to storing sequence information in the form of words, gene
structure and annotations are also recorded in the database. Basic
metadata captured are:

* Gene ID ("LOC859")
* Species ("Homo sapiens") (redundant to Genome Build)
* Symbols ("DUSP27), both preferred (generally the "official" symbol)
  and aliases.
* Description ("EF-hand calcium binding domain 2")
* Human orthologue information (LOC23418 | CRB1 | 97.51%) in the form
  of GeneID, primary symbol, and sequence similarity (generally the
  "percent total identity" for the best alignment between the protein
  isoforms available to the two loci).

Additionally, the "structure" of the gene is also included. This can
encompass:

* Splice variants are collected for each gene
* Exon structure on the genome, or exon/exon boundaries for orphan
  RNAs (though often not available for orphans) for each splice
  variant
* CDS (ORF) start/stop coordinates if the transcript is protein coding
* Alignment quality in the form of "percent total identity" for
  RNA-genome alignments

* "percent total identity" = `2 * (number identical positions) /
  (lengthSeq1 + lengthSeq2)`

## Searching a Query Sequence

Any nucleotide query can be searched against the database, so long as
it is at least as long as the wordsize used to populate the database
(this value is referenced as `$wordsize` below). The initial search is
to identify word matches in the database against the query, which
takes the following steps:

#### Determination of "raw" word hits

1. The query is broken into _non-overlapping_ words of length
   `$wordsize`. The 'last' word may be overlapping, in order to extend
   'neatly' to the end of the query. For example, a 25bp oligo will
   have three words: [1,12], [13,24], [14,25].
2. Each word is then "expanded" in the same way that database
   sequences were - ambiguity characters are replaced by all
   combinatoric permutations to generate unique, non-ambiguous
   words. Unlike database generation, __all__ query permutations are
   kept. In practice, queries generally do not have ambiguities, so
   this "expansion" should still result in only a single word.
3. An additional permutation is performed to account for
   mismatches. When searching, the user provides a maximum tolerable
   number of mismatches, which is by default `2`. All possible
   permutations that include that number of mismatches are then
   included on the ambiguity-permuted word set.
4. The exapanded words are then searched for hits in the database. If
   a hit is found, the `SubjectID`, `SubjectPos` and `MM` (number of
   mismatches) are recorded.
5. This process is performed "left-to-right". The first word will
   always be searched in full. All subsequent words will have their
   hits compared to surviving hits from the prior word(s). For a hit
   on a 'new' word to be kept, it must:
   1. Be from an already noted `SubjectID`. For example, if the second
      word finds a hit to 'NM_123456.3' but that SubjectID was not
      seen in the first word, then the hit is discarded.
   2. Be exactly `$wordsize` bases away from a previously seen
      `SubjectPos`. For example, Word1 and Word2 both hit
      'NM_123456.3', but the positions were '12' and '302'. Those
      positions are not "perfectly adjacent", so the hit is discarded.
   3. Have the total `MM` count not exceed the limit. For example,
      Word1 hits 'NM_123456.3' at position 12 with MM=1 and Word2 hits
      'NM_123456.3' at position 24 with MM=2. Even though the hits are
      on the same sequence, and are adjacent, the total number of
      mismatches needed to build the hit is now 3, so the hit is
      discarded.
   4. The process described above is "looking backwards" - the new
      hits are discarded if they are not supported by "previous
      hits". However, the process also "looks forward" - For __every
      word__, there must be a chain of adjacent hits that meet the
      three above criteria for the query to be counted as "hitting"
      the subject.
6. Any subjects meeting the above criteria are recorded with
   `SubjectID`, `SubjectPos`, `HitStrand` and `MM`.
7. The above process is then repeated, but using the reverse
   complement of the query.

#### Mapping raw hits to genes

Presuming that at least one raw hit survived the process above, the
hit is then interpreted in the context of gene structure. Because the
database excludes purely intergenic sequence, each hit will overlap
(in some way) one _or more_ genes. For each gene, calculate:

1. `GeneStrand` is calculated relative to the query. In most cases the
   genes reside on the forward strand of the `Subject`; However, if a
   subject has two or more overlapping genes, the gene may have been
   from the reverse strand. `GeneStrand` is then the product of
   `HitStrand` and `SubjectStrand`. That is, a hit that was from the
   -1 (reverse complement) strand of the query, and is overlapping a
   gene that is on the -1 strand of the subject, is then a +1 (-1 *
   -1) strand for the gene. This calculation will generally eliminate
   some hits (and some genes) because the default search parameters
   include a filter that only allows -1 strand hits between the query
   and genes.
2. Two categories of position are calculated for the hit:
   * Genomic cooridnates, presuming the hit was not an orphan RNA
   * RNA coordinates, which will be of the form `145..159` for a hit
     that falls within an exon, or `546^547` (falling between RNA
     bases 546 and 547) for hits in introns, or more complex
     specifications for esoteric hits (like splice junctions)
3. The `HitType` is calculated based on where in the gene structure
   the query lands. This can be values like 'Intron' or 'UTR3' - see
   below for a description of reported Impact types.

#### Impact types

The following query:gene types are reported by the system:

* Intron :: The oligo is predicted to only fall in introns
* Exon :: The oligo covers an exon, at least in one variant. If
  CDS/ORF coordinates were available, will be broken out into a
  subtype:
  * UTR5 (5' UTR)
  * CDS
  * UTR3 (3' UTR)
* ExonExon :: The oligo spans an exon-exon boundary
* FullSpan :: The oligo COMPLETELY covers an exon - presumably this is
  a tiny exon or a huge oligo
* SpliceJunction :: The oligo covers an exon-intron boundary. It will
  not be reported directly, but will rather have one of the following
  subtypes:
  * SpliceDonor
  * SpliceAcceptor
* RNA :: The oligo hits an mRNA (not a pre-mRNA). It is exonic, but we
  can not tell if it crosses an Exon-Exon boundary or not
* CrypticExon :: The oligo is within exonic sequence that does not
  appear to be represented on the reference genome build. That is,
  there is RNA sequence that appears to be absent in the gDNA.
* FirstExon :: The oligo is MOSTLY in the first exon, but a little bit
  hangs off the 5' end
* LastExon :: The oligo is MOSTLY in the last exon, but a little bit
  hangs off the 3' end

It is important to note that a query may have _multiple_ types
assigned to it, as it might overlap multiple transcripts or even
multiple genes.

#### Caching hits

The process of finding and categorizing hits is somewhat
laborious. For this reason, all query hit results are cached in the
database for rapid recovery later (it is presumed that most queries
will be submitted to the database multiple times). While this is
primarirly a performance detail, it should be kept in mind if the
database is modified (for example, by adding new genome builds to
it). In such circumstances it would be desirable to clear out cached
results (a manual SQL process).

#### Reporting hits

Several structures and formats are available for reporting
hits. Structurally, hits may be reported one of three ways:

* Summary - For each Query+GenomeBuild pair the _count_ of perfect, 1
  mismatch and 2 mismatch hits are reported.
* Overview - Similar to summary, but a condensed list of each gene is
  reported, broken out as a list of gene symbols for each of six
  categories: {Exon vs Intron} x {0MM, 1MM, 2MM}
* Details - Every hit is reported in detail; Each record reports:
  * Query
  * SubjectID
  * GenomeBuild + Species
  * GeneID
  * GeneSymbol
  * HumanID
  * HumanSymbol
  * QuerySequence
  * SubjectSequence
  * SubjectPosition
  * RnaPosition
  * HitStrand
  * GeneStrand
  * ImpactType
  * MM
  * AM (number of matches that fell on an ambiguous nucleotide)

These results can be reported in one of several formats:

* Dynamic (sortable, filterable) HTML tables. Default format for
  searches performed from the human-driven web form
* Excel worksheet
* TSV table. Default format for command line searches
* Simple text report
* Rich JSON structure. Designed for web service integration with other
  tools.

[RSG]: https://www.ncbi.nlm.nih.gov/refseq/rsg/
[VCF]: https://en.wikipedia.org/wiki/Variant_Call_Format
