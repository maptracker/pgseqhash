#!/stf/biobin/perl -w

my $isBeta;
BEGIN {
    # Utilize patched BioPerl libraries
    use lib '/stf/biocgi/tilfordc/patch_lib';
    
    # Allows usage of beta modules to be tested:
    my $progDir = join(' ', $0, `pwd`);
    if ($progDir =~ /(working|perllib)/) {
        require lib;
        import lib '/stf/biocgi/tilfordc/perllib';
        $isBeta = 1;
    } else {
        $isBeta = 0;
    }
    $| = 1;
    print '';
}

use strict;
use BMS::BmsArgumentParser;
use BMS::MapTracker::AccessDenorm;
use BMS::Utilities::SequenceUtilities;
use BMS::MapTracker::GenAccService;
use BMS::SequenceLibraryFinder;
use BMS::FastaMeta;

use Bio::DB::Fasta;
use Bio::SeqIO;
use JSON;
use DB_File;

my $args = BMS::BmsArgumentParser->new
    ( -nocgi        => $ENV{HTTP_HOST} ? 0 : 1,
      -nonredundant => 0,
      -block        => 100,
      -paramalias => {
          limit     => [qw(lim)],
          verbose   => [qw(vb)],
          capture   => [qw(extract select)],
          progress  => [qw(prog)],
          fork      => [qw(forknum numfork)],
          output    => [qw(out outfile export)],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
      });

$args->shell_coloring();


my $age         = $args->val('age') || '8 Sep 2016 4pm';
my $ad          = BMS::MapTracker::AccessDenorm->new
    ( -age => $age, -cloudage => $age);
my $mt          = $ad->tracker();
my $slf         = BMS::SequenceLibraryFinder->new();
my $su          = BMS::Utilities::SequenceUtilities->new();
my $limit       = $args->val('limit');
my $clobber     = $args->val(qw(clobber));
my $progCnt     = 1000;
my $veryBigNum  = 99999999999999999999999;

my $forkNum     = $args->val('fork') || 30;
my $nsr         = $ad->namespace_token( $args->val('rnans') || 'RSR' );
my $nsg         = $ad->namespace_token( $args->val('genens') || 'LL' );
my $block       = $args->val(qw(block)) || 100;
my $xxSize      = 30;
my $xxMod       = $xxSize - 1;
my $noSeq       = $args->val(qw(noseq));
my $rnaInput    = $args->val(qw(rna transcript));
my $dnaFa       = $args->val(qw(dna genome));
my $outF        = $args->val('output');
my $nonRedun    = $args->val(qw(nonredundant));
my $allowOrphan = $args->val('allowOrphan');
my ($fans)      = $slf->stnd_ns($nsr eq 'RSR' ? 'RefSeq' : 'Ensembl');
my %specMeta    = map { $_ => 1 } qw(id v sym desc src);

if (my $capRe = $args->val('capture')) {
    &capture($capRe);
}
if (my $id = $args->val('test')) {
    &test($id);
    exit;
}

if (my $bReq = $args->val('build')) {
    my @dnas = split(/[\n\r]/, `ls -1 /genomes01/*.$bReq.fa`);
    if ($#dnas == 0) {
        $dnaFa = $dnas[0];
    } elsif ($#dnas == -1) {
        $args->death("Failed to find genome for build '$bReq'");
    } else {
        $args->death("Multiple genome for build '$bReq'", @dnas);
    }
}

$args->death("You must provide the genome file with -dna",
             "Or provide a build with -build") unless ($dnaFa);

my $shrtDna = $dnaFa; $shrtDna =~ s/.+\///;
my $sdbNm   = "$shrtDna.13.13.btdb";
my $sdb     = $mt->get_searchdb( -name => $sdbNm );

$args->death("Failed to find genomic search DB in MapTracker",
             $sdbNm, "You may need to first perform NOSHA searches")
    unless ($sdb);


my ($absTag, %anchor, $taxa, $build, %collected,
    $bioMeta);

if ($shrtDna =~ /^([^\.]+)\.genome\.([^\.]+)\.fa$/) {
    ($taxa, $build) = ($slf->stnd_taxa($1), $2);
    #$taxa =~ s/_/ /g;
    #substr($taxa, 0, 1) = uc(substr($taxa, 0, 1));
    $absTag = sprintf("%s Genome Build %s", $taxa, $build);
}

unless ($rnaInput) {
    my @rnas = $slf->get_libraries( -taxa => $taxa,
                                    -type => 'rna',
                                    -ns   => $fans, );
    if ($#rnas == 0) {
        $rnaInput = $rnas[0];
    } elsif ($#rnas == -1) {
        $args->death("Failed to find transcripts for NS=$fans, TAXA=$taxa");
    } else {
        $args->death("Multiple transcripts for NS=$fans, TAXA=$taxa", @rnas);
    }
}
$args->death("You must provide the RNA file with -rna") unless ($rnaInput);

if (!$outF && $taxa && $build && $fans) {
    $outF = sprintf("%s %s pre-mRNA %s", lc($taxa), $fans, $build);
    $outF =~ s/\s+/_/g;
}

unless ($outF) {
    $args->msg("[!]","Please provide an output filename with -output");
    exit;
}

my $shrtRna = $rnaInput; $shrtRna =~ s/.+\///;

my $fdb    = Bio::DB::Fasta->new($dnaFa);
my $rdb    = Bio::DB::Fasta->new($rnaInput);
my $num    = 0;

if (0) {
    my $orid = "REFSEQN:XM_011531503.1";
    my $sl   = $rdb->length($orid);
    my $seq  = uc($rdb->seq($orid) || "");
    warn sprintf("%s : len = %d, length(seq) = %d", $orid, $sl, length($seq||""));
    die;
}

my $msgF = $shrtDna . "-preMrnaMessages.txt";
my $posF = join('-vs-', $shrtRna, $shrtDna);
$posF   .= "-NonRedun" if ($nonRedun);
my $rnaF = join('-OrphanRNA-', $shrtRna, $shrtDna);
my $orpF = "$shrtRna-NoGene";
$rnaF    =~ s/\.fa$//;
$outF    =~ s/\.fa$//;
if ($limit) {
    $posF   .= "-LIMIT";
    $rnaF   .= "-LIMIT";
    $orpF   .= "-LIMIT";
    $outF   .= "-LIMIT";
}
my $gasBase = $posF . "-GAS-";
$posF   .= ".tsv";
$rnaF   .= ".fa";
$orpF   .= ".fa";

if (my $pReq = $args->val(qw(posfile))) {
    $args->death("-posfile '$pReq' does not exist") unless (-e $pReq);
    $rnaF  = "";
    $orpF  = "";
    $posF  = $pReq;
    $outF .= "-CustomPos";
} elsif ($clobber || $limit)  {
    unlink($posF) if (-e $posF);
    unlink($rnaF) if (-e $rnaF);
}
$outF   .= ".fa";

open(MSGS, ">$msgF") || $args->death
    ("Failed to open message file", $msgF, $!);

my $fh    = *STDOUT;

my $smFmt    = '%15s : %s';

&msg_list('IN', {
    Build => $build,
    Taxa  => $taxa,
    NS    => $fans,
    'Build Set' => $absTag
}, [sprintf($smFmt, "DNA", $dnaFa),sprintf($smFmt, "RNA", $rnaInput)] );

my $outLines = {
    Output          => $outF,
    Messages        => $msgF,
    Positions       => $posF,
    'Orphaned mRNA' => $rnaF,
    'RNAs w/o Loci' => $orpF,
};

&process($rnaInput);

close MSGS;
&msg_list('OUT', $outLines );



sub msg_list {
    my ($tag, $hash, $list) = @_;
    $list ||= [];
    foreach my $k (sort keys %{$hash}) {
        if (my $v = $hash->{$k}) { push @{$list}, sprintf($smFmt, $k, $v); }
    }
    $args->msg("[$tag]", @{$list});
}

sub process {
    if ($args->val(qw(trial))) {
        $args->msg("[-]", "Trial run, just checking parameters");
        return;
    }
    if (&_is_up_to_date($posF, $rnaInput) && !$limit) {
        $args->msg("[*]", "Using existing position file", $posF);
    } else {
        &make_pos_file($rnaInput);
    }
    &read_pos_file();
}

sub make_pos_file {
    my $file = shift;
    open(POS, ">$posF") || $args->death
        ("Failed to write position file", $posF, $!);
    open(MRNA, ">$rnaF") || $args->death
        ("Failed to write mRNA file", $rnaF, $!);
    open(NOG, ">$orpF") || $args->death
        ("Failed to write gene-less RNA file", $orpF, $!);
    my $instream = Bio::SeqIO->new( -file   => $file,
                                    -format => 'fasta');
    $args->msg("[>]", "Creating genomic position file", $posF,
               "Source: $file", "mRNAs: $rnaF", "RNAs w/o genes: $orpF");
    &create_biological_metadata( $rnaInput );

    &read_biological_metadata( $rnaInput);

    while ( my $bs = $instream->next_seq( ) ) {
        $num++;
        &extract_rna($bs);
        last if ($limit && $num >= $limit);
    }
    close POS;
    close MRNA;
    close NOG;

    my $sF = $posF."-sort";
    # my $cmd = "sort -S 10G -k 1,1 -k 2,2nr -k 3,3 \"$posF\" > \"$sF\"";
    my $cmd = "sort -S 10G -k 1,1b \"$posF\" > \"$sF\"";
    $args->msg("[+]","Sorting position data...", $cmd);
    system($cmd);
    system("mv \"$sF\" \"$posF\"");
    $args->msg("[+]","Done");
}

sub _is_up_to_date {
    my $file = shift;
    return 0 unless (-s $file);
    my $mod = -M $file;
    foreach my $other (@_) {
        next unless (-s $other);
        my $omod = -M $other;
        return 0 if ($omod < $mod);
    }
    return 1;
}

sub create_biological_metadata {
    my $file  = shift;
    my $mfile = "$file.dbfile";
    if (&_is_up_to_date($mfile, $file)) {
        return $mfile;
    }
    unlink($mfile) if (-e $mfile);
    open(RFA, "<$file") || $args->death("Failed to read fasta", $file, $!);
    my (@ids, %data);
    my $rNum = 0;
    $args->msg("[<]", "Parsing fasta file", $file);
    while (<RFA>) {
        if (/^>(.+)/) {
            my $id = $1;
            $id =~ s/[\n\r]+$//;
            my $desc = "";
            if ($id =~ /^(\S+)\s+(.+?)$/) {
                ($id, $desc) = ($1, $2);
            }
            my $raw = $id;
            $id =~ s/^.+\://;
            if ($id =~ /^(.+)\.(\d+)$/) {
                $id = $1;
                $data{$id}{v} = $2;
            }
            push @ids, $id;
            $data{$id}{id}   = $id;
            $data{$id}{desc} = &_clean_desc($desc);
            $data{$id}{src}  = $raw unless ($raw eq $id);
            $rNum++;
        }
    }
    close RFA;
    my $tsv   = "$mfile.tsv";
    $args->msg("[<]","Recovering metadata for $rNum RNA IDs", $tsv);
    unless (&_is_up_to_date($tsv, $file)) {
        $args->msg("[-]","Refreshing RNA metadata...");
        unlink($tsv);
    }
    my $rows  = &cached_data
        ( -ids      => \@ids,
          -ns1      => $nsr, -ns2      => $nsg,
          -keepnull => 1,    -scramble => 1,
          -cols     => 'termin,termout,sym,taxa,desc',
          -output   => $tsv );

    my %g2r;
    my $noGene = 0;
    foreach my $row (@{$rows}) {
        my ($rid, $gid, $sym, $tax, $gdesc) = @{$row};
        $sym = &_clean_sym( $sym );
        $data{$rid}{sym}   = $sym;
        $data{$rid}{taxa}  = $tax;
        if ($gid) {
            $data{$rid}{gene} = $gid;
            $data{$gid}{sym}  = $sym;
            $data{$gid}{taxa} = $tax;
            $data{$gid}{desc} = &_clean_desc($gdesc);
            $g2r{$gid}{$rid}  = 1;
        } else {
            $noGene++;
        }
    }
    my $gNum = 0;
    while (my ($gid, $rH) = each %g2r) {
        $data{$gid}{rna} = [ sort keys %{$rH} ];
        $gNum++;
    }

    $args->msg("[*]","Building DB_File metadata index for RNAs", $mfile);
    my $btinfo = new DB_File::BTREEINFO;
    $btinfo->{'cachesize'} = $args->val('cachesize') || 250000000;
    $btinfo->{'psize'}     = $args->val('pagesize')  || 8192;

    my %dbf;
    tie(%dbf, 'DB_File', $mfile, O_CREAT|O_RDWR, 0644, $DB_BTREE ) ||
        $args->death("Failed to tie DBFILE hash", $mfile, $!);
    # Some file-level metadta:
    my $globalData  = {
        FASTA_FILE     => $file,
        GDNA_SOURCE    => $dnaFa,
        RNA_SOURCE     => $rnaInput,
        PREPARE_TIME   => &_nice_date(),
        RNA_COUNT      => $rNum,
        GENE_COUNT     => $gNum,
        RNA_NO_GENE    => $noGene,
        TAXA           => $taxa,
        BUILD          => $build,
        BUILD_SET      => $absTag,
        RNA_NAMESPACE  => $fans,
        GENE_NAMESPACE => $nsg,
    };
    $dbf{GLOBAL_DATA} = encode_json( $globalData );

    foreach my $id (sort keys %data) {
        $data{id} = $id;
        my $json  = encode_json( $data{$id} );
        $dbf{$id} = $json;
    }
    untie %dbf;
    $args->msg("[>]", "Generated metadata file for $rNum RNAs and gNum genes",
               $mfile, $tsv);
    return $mfile;
}

sub read_biological_metadata {
    return $bioMeta if ($bioMeta);
    my $file  = shift;
    my $mfile = "$file.dbfile";
    my %dbf;
    tie(%dbf, 'DB_File', $mfile, O_RDONLY, 0644, $DB_BTREE ) ||
        $args->death("Failed to tie DBFILE hash", $mfile, $!);
    $args->msg("[+]","RNA metadata file loaded", $mfile);
    return $bioMeta = \%dbf;
}

sub biometa {
    my $id = shift;
    my $rv = {};
    if (my $bmTxt = &bioJson( $id )) {
        $rv = decode_json( $bmTxt );
    }
    $rv->{id} ||= $id;
    return $rv;
}

sub bioJson {
    my $id = shift;
    if ($id) {
        $id =~ s/\.\d+$//;
        return $bioMeta->{$id} || "";
    }
    return "";
}

sub msg {
    my @m = @_;
    # $args->msg(@m);
    print MSGS join("\t", @m)."\n";
}

sub extract_rna {
    my ($bs)   = @_;
    my $rid    = &_clean_id( $bs->display_id() );
    my $unv    = &_unv_id( $rid );
    my $rmeta  = &biometa( $rid );
    # die $args->branch($rmeta);
    my $geneID = $rmeta->{gene};
    unless ($geneID) {
        # We were not able to associate this RNA with a gene
        # Ignore it
        $rmeta->{ignored} = "No locus found for this RNA";
        print NOG &meta_fasta( $bs, $rmeta);
        return;
    }

    my @maps = sort { $b->score() <=> $a->score() } 
    $mt->get_mappings( -sdb  => $sdb, -name => $rid );

    my $bestSc = 0;
    my @kept;
    if ($#maps == -1) {
        # No maps found for this RNA
        my @knownAbsent;
        if ($absTag) {
            @knownAbsent = $mt->get_edge_dump( -name1 => $rid,
                                               -name2 => $absTag,
                                               -type  => 'is absent from');
        }
        $rmeta->{orphan} = $#knownAbsent == -1 ? 
            'no locations found' : "reported absent from $absTag";
        $rmeta->{build} = $build;
        # Record the raw RNA to the file
        print MRNA &meta_fasta( $bs, $rmeta);
        return;
    } else {
        $bestSc = $maps[0]->score();
        foreach my $map (@maps) {
            # Keep all maps with the top score
            last if ($map->score() < $bestSc);
            push @kept, $map;
        }
    }

    my $mnum = $#kept + 1;
    if ($mnum > 1) {
        &msg("[$num]", "$rid : $mnum maps at $bestSc");
    } elsif (!($num % $progCnt)) {
        $args->msg("[$num]", "$rid");
    }
    my $sym = $rmeta->{sym} || $geneID;
    
    foreach my $map (@kept) {
        my $chrSeq = $map->other_seq( $rid );
        my $str    = $map->strand() < 0 ? 'R' : 'F';
        my @locs   = sort { $a->[0] <=> $b->[0] } 
        $map->locations_for_seq( $chrSeq );
        my $cid    = $chrSeq->name();
        my ($s,$e) = ($locs[0][0], $locs[-1][1]);

        my $show   = $cid;
        if ($show =~ /^([^\.]+)\.([^\.]+)\.([^\.]+)\.([^\.]+)$/) {
            my ($tx, $typ, $id, $bld) = ($1, $2, $3, $4);
            if ($typ eq 'chromosome') {
                $id = "0$id" if ($id =~ /^\d$/);
                $show = "$id.$bld";
            } else {
                $show = $id;
            }
        }
        # my ($es, $ee, $rs, $re) = @{$loc};
        my $locTxt = join('|', map { join(',', @{$_}) } @locs);
        my $key    = sprintf("%s:%012d.%012d", $show, $s, $e);
        my @row    = ($key, $str, $unv, $sym, $mnum, $bestSc, $cid, $locTxt);
        print POS join("\t", @row)."\n";
    }
}

sub meta_fasta {
    my ($bs, $meta) = @_;
    my $id = $meta->{id};
    if (my $v = $meta->{v}) {
        $id .= ".$v";
    }
    my $desc = $meta->{desc} || "";
    if (my $sym = $meta->{sym}) {
        $desc = "[$sym]".($desc ? " $desc" : "");
    }
    foreach my $sk (sort keys %{$meta}) {
        unless ($specMeta{$sk}) {
            my $v  = $meta->{$sk};
            if (defined $v && $v ne '') {
                $desc .= " " if ($desc);
                $desc .= sprintf("/%s='%s'", $sk, $v);
            }
        }
    }
    my $seq = $bs->seq();
    my $sl  = CORE::length($seq);

    my $rv = ">$id";
    $rv   .= " ".$desc if ($desc);
    $rv   .= "\n";
    for (my $p = 0; $p < $sl; $p += $block) {
        $rv .= substr($seq, $p, $block)."\n";
    }
    return $rv;
}

sub _nice_date {
    my $dt   = `date +'%Y-%b-%d %H:%M:%S'`;
    $dt =~ s/[\n\r\s]+$//;
    return $dt;
}

sub _clean_id {
    my $id = shift || "";
    # Remove SeqStore container
    # REFSEQN:NM_001145054.1
    $id =~ s/^.+\://;
    return $id;
}

sub _unv_id {
    my $id = shift || "";
    # Remove version number
    $id =~ s/\.\d+$//;
    return $id;
}

sub _clean_sym {
    my $sym = shift;
    if ($sym) {
        $sym =~ s/[\*~]$//;
    }
    return $sym;
}

sub _clean_desc {
    my ($d, $retFlag) = @_;
    return "" unless ($d);
    $d =~ s/( \[LOC\d+\]| \[[NP|XP]_\d+\])+\s*$//;
    $d =~ s/\/taxid=\d+/ /;
    $d =~ s/\,?\s*mRNA/ /;
    $d =~ s/[\s\t\n\r\"]+/ /g;
    $d =~ s/\s+/ /g;
    $d =~ s/[\.\s\,]+$//;

    my $foundTag;
    while ($d =~ /(\/([^=]+)=(\S+))/) {
        my ($r, $k, $v) = ($1, $2, $3);
        $foundTag ||= {};
        $foundTag->{$k}{$v} = 1;
        $d =~ s/\Q$r\E//g;
    }

    if ($foundTag) {
        # Deal with found tags
        while (my ($k, $vs) = each %{$foundTag}) {
            if ($retFlag) {
                # Set tags to a hash
                map { $retFlag->{$k}{$_} = 1 } keys %{$vs};
            } else {
                # Add the tags back in to the description
                $d .= " " if ($d);
                $d .= sprintf("/%s='%s'", $k, join(',', sort keys %{$vs}));
            }
        }
    }
    return $d;
}

sub metadata {
    my ($req, $ns) = @_;
    $req = uc($req || "");
    unless (defined $anchor{$req}) {
        my ($sym, $desc);
        $sym  = &_clean_sym( $ad->best_possible_symbol
                             ( $req, $ns, "best short" ) || "");
        $desc = &_clean_desc( $ad->description
                              ( -id => $req, -ns => $ns ) || "");

        my ($id, $seq) = $ad->standardize_id($req, $ns);
        $anchor{$req} = {
            id   => $id,
            ns   => $ns,
            sym  => $sym,
            desc => $desc,
        };
    }
    return $anchor{$req};
}

sub all_metadata {
    my $cols = shift;
    my %byNS;
    foreach my $col (@{$cols}) {
        foreach my $ans (sort keys %{$col->{anc}}) {
            foreach my $anc (sort keys %{$col->{anc}{$ans}}) {
                $byNS{$ans}{$anc}++;
            }
        }
    }
    my @idBits;
    foreach my $ns (sort keys %byNS) {
        my @ids   = keys %{$byNS{$ns}};
        push @idBits, scalar(@ids)." $ns";
        my $file  = $outLines->{"Annot $ns"} = $gasBase . $ns . ".tsv";
        unlink($file);
        my $rows  = &cached_data
            ( -ids => \@ids, -ns1 => $ns, -mode => 'simple',
              -cols => 'termout,sym,desc', -output => $file );
        foreach my $row (@{$rows}) {
            my ($id, $sym, $desc) = @{$row};
            $anchor{uc($id)} = {
                id   => $id,
                ns   => $ns,
                sym  => $sym,
                desc => $desc,
            };
        }
        map { $anchor{uc($_)} ||= {
                id   => $_,
                ns   => $ns,
                sym  => "",
                desc => "",
            } } @ids;
    }
    $args->msg("[+]", "Anchors used: ".join(' ,', @idBits));
}

sub cached_data {
    my $isBeta = 1;
    my $gas = BMS::MapTracker::GenAccService->new
        ( -fork       => $forkNum,
          -ignorecase => 1,
          @_,
          -format  => 'tsv',
          -age      => $age,
          -cloudage => $age,
          -verbose  => 1,
          -scramble => 1, @_ );
    $gas->use_beta( $isBeta );
    my $rows = $gas->cached_array( $clobber );
    return $rows;
}


sub read_pos_file {

    &read_biological_metadata( $rnaInput);
    open(POS, "<$posF") || $args->death
        ("Failed to open position file", $posF, $!);
    $args->msg("[<]", "Reading genomic position file", $posF);

    if ($outF) {
        open(OUT,">$outF") || $args->death
            ("Failed to open output", $outF, $!);
        $fh = *OUT;
    }

    my $mfile  = $outF . ".dbfile";
    $args->msg("[*]","Building DB_File metadata index for output", $mfile);
    my $btinfo = new DB_File::BTREEINFO;
    $btinfo->{'cachesize'} = $args->val('cachesize') || 250000000;
    $btinfo->{'psize'}     = $args->val('pagesize')  || 8192;

    my %dbf;
    tie(%dbf, 'DB_File', $mfile, O_CREAT|O_RDWR, 0644, $DB_BTREE ) ||
        $args->death("Failed to tie DBFILE hash", $mfile, $!);
    # Some file-level metadta:
    $dbf{SOURCE_FILE}  = $outF;
    $dbf{GENOME_FILE}  = $dnaFa;
    $dbf{PREPARE_TIME} = &_nice_date();

    my $data  = { show => "", e => 0, s => 0 };
    my $count = 0;
    my $lnum  = 0;
    my $done  = {};
    while (<POS>) {
        s/[\n\r]+$//;
        my ($key, $str, $rid, $sym, $mnum, $bestSc, $cid, $locTxt) 
            = split(/\t/);
        my ($show, $se) = split(':', $key);
        my ($s, $e) = map { $_ + 0 } split(/\./, $se);

        if ($data->{e} < $s || $e < $data->{s} ||
            $data->{show} ne $show) {
            # This segment does not overlap with the prior one
            # Or it is a different anchor (chromosome)

            my $oldId = &process_data( $data, \%dbf, $done );
            $count++ if ($oldId);
            $data = {
                show  => $show,
                cid   => $cid,
                s     => $s,
                e     => $e,
                prior => $oldId,
                str   => {},
                rna   => {},
            };
        } else {
            # Overlaps with prior position. Extend the end position
            $data->{e} = $e if ($data->{e} < $e);
            $data->{s} = $s if ($data->{s} > $s);
        }
        $data->{str}{$str}++;
        my $unv   = &_unv_id( $rid );
        my $rtarg = $data->{rna}{$unv} ||= {
            id  => $unv,
            num => $mnum,
            loc => [],
        };
        my @locs = map { [ split(',', $_) ] } split(/\|/, $locTxt);
        push @{$rtarg->{loc}}, [ $bestSc, $str, \@locs ];
    }
    $count++ if ( &process_data( $data, \%dbf, $done ) );

    &add_orphans( \%dbf, $done );

    my @rids = keys %{$done->{rna}};
    my @orph = keys %{$done->{orph}};
    my @gids = sort keys %{$done->{gene}};
    foreach my $gid (@gids) {
        my @ids   = @{$done->{gene}{$gid}};
        my $gmeta = &biometa( $gid );
        $gmeta->{gdna} = \@ids;
        if (my $gunk = $done->{gunk}{$gid}) {
            $gmeta->{unaligned} = $gunk;
        }
        my $json   = encode_json( $gmeta );
        $dbf{$gid} = $json;
    }

    my $gNum = $#gids + 1;
    my $rNum = $#rids + 1;
    my $oNum = $#orph + 1;
    $dbf{GENOME_COUNT}  = $count;
    $dbf{RNA_COUNT}     = $rNum;
    $dbf{ORPHAN_COUNT}  = $oNum;
    $dbf{GENE_COUNT}    = $gNum;

    untie %dbf;

    my @ua = keys %anchor;
    $args->msg("[+]", 
               "        gDNA: $count",
               "        Gene: $gNum",
               "Anchored RNA: $rNum",
               "Orphaned RNA: $oNum" );
    
    
    if ($outF) {
        close OUT;
    } else {
        $outF = "STDOUT";
    }

}

sub add_orphans {
    my ($dbf, $done) = @_;
    return unless (-s $rnaF);
    open(MRNA, "<$rnaF") || $args->death
        ("Failed to read mRNA file", $rnaF, $!);
    while (<MRNA>) {
        print $fh $_;
        if (/^>(\S+)/) {
            my $rmeta = &biometa($1);
            my $rid   = $rmeta->{id};
            my $unv   = &_unv_id( $rid );
            unless ($done->{orph}{$unv}++) {
                $rmeta->{orphan} = 1;
                my $json       = encode_json( $rmeta );
                $dbf->{$unv}   = $json;
            }
            my $gid = $rmeta->{gene} || 'UNK';
            if ($done->{gene}{$gid}) {
            } else {
                $done->{gene}{$gid} ||= [];
            }
            push @{$done->{gunk}{$gid}}, $unv;
        }
    }
    close MRNA;
}

sub process_data {
    my ($data, $dbf, $done) = @_;
    my $show = $data->{show};
    return "" unless ($show);
    my $cid  = $data->{cid};
    $show    =~ s/^0//;
    my ($s, $e) = ($data->{s}, $data->{e});
    my $strTxt  = join('', sort keys %{$data->{str}});
    my $id      = sprintf("%s:%d-%d.%s", $show, $s, $e, $strTxt );
    
    my (%genes, %exonExon);
    # This is going to be the serailized data structure in the DB_File
    # that describes the gene structure:
    my $jData = {
        id  => $id,
        s   => $s,
        e   => $e,
        chr => $cid,
    };
    my $chrText = join(' ', $show, $s, $e);
    while (my ($rid, $rdata) = each %{$data->{rna}}) {
        # We will store versioning information in the RNA metadata
        # Record RNAs as just the unversioned ID
        # This will be a problem if multiple versions of an RNA are in 
        # the fasta file. That would be silly. So do not do that.
        $rid =~ s/\.\d+$//;

        my $rmeta  = &biometa( $rid );
        my $geneID = $rmeta->{gene} || 'UNK';
        my $gt = $genes{$geneID} ||= {
            id  => $geneID,
            s   => $veryBigNum,
            e   => 0,
            sc  => 0,
            rna => {},
            str => {},
        };
        $gt->{rna}{$rid} = 1;
        my $orid = $rmeta->{src};
        unless ($orid) {
            warn $args->branch($rmeta);
            $args->death("Failed to recover original RNA from fasta file",
                         "Accession: $rid", );
                         
        }
        my $sym  = $rmeta->{sym};
        foreach my $ld (@{$rdata->{loc}}) {
            my ($sc, $str, $locs) = @{$ld};
            my ($rs, $re) = ($locs->[0][0], $locs->[-1][1]);
            $gt->{sc} = $sc if ($gt->{sc} < $sc);
            $gt->{s}  = $rs if ($gt->{s} > $rs);
            $gt->{e}  = $re if ($gt->{e} < $re);
            $gt->{str}{$str} = 1;
            my $linds = $#{$locs};
            for my $i (0..$linds) {
                my $loc = $locs->[$i];
                my ($gs,$ge,$rs,$re) = @{$loc};
                if ($i) {
                    # Capture exon-exon boundaries
                    my ($pgs,$pge,$prs,$pre) = @{$locs->[$i-1]};
                    my $eeTarg = $exonExon{"$pge $gs $str"} ||= {
                        s => $pge,
                        e => $gs,
                        str => $str,
                        src => [],
                    };
                    my ($us, $ue, $num) = ($str eq 'F') ?
                        ($pre, $rs, $i) : ($re, $prs, $linds - $i + 1);
                    push @{$eeTarg->{src}}, [$rid, $orid, $sym,$us, $ue, $num];
                }
            }
            my @gc = ($str eq 'R') ? reverse @{$locs} : @{$locs};
            # Record the full coordinate data for this RNA in the segment:
            push @{$gt->{hsps}}, {
                str => $str,
                sc  => $sc,
                rna => $rid,
                hsp => \@gc,
            };
        }
    }
    my (@descs, %uGene);
    foreach my $gt ( sort { $a->{s} <=> $b->{s} || 
                                $a->{id} cmp $b->{id} } values %genes) {
        my $gid    = $gt->{id};
        my $gmeta  = &biometa( $gid );
        # Normalize the strand to a string
        my $str    = $gt->{str} = join('', sort keys %{$gt->{str}});
        # Normalize RNA IDs to an array:
        my @rids   = sort keys %{$gt->{rna}};
        $gt->{rna} = \@rids;

        # Make sure that the underlying RNA/gene metadata are copied over:
        push @{$done->{gene}{$gid}}, $id;
        foreach my $rid (@rids) {
            unless ($done->{rna}{$rid}++) {
                my $bmTxt = &bioJson( $rid );
                $dbf->{$rid} = $bmTxt;
            }
        }
        my ($gs, $ge, $sc) = map { $gt->{$_} } qw(s e sc);
        my $sym = $gmeta->{sym};

        my $desc = "";
        $desc .= "[$sym] " if ($sym);
        $uGene{$gid} ||= $sym || $gid;

        # Note relative coordinates of the gene:
        $desc .= sprintf("%d-%d %s %.1f%% ", 
                         $gs - $s + 1, $ge - $s + 1, $str, $sc);

        my $rnum = $#rids + 1;
        if ($rnum == 0) {
            $desc .= "$rnum RNAs? ";
        } else {
            $desc .= sprintf("%d RNA%s. ", $rnum, $rnum == 1 ? '' : 's');
        }
        if (my $gd = $gmeta->{desc}) {
            $desc .= $gd;
        }
        $desc =~ s/\"+//g;
        $desc =~ s/\s+$//;
        push @descs, "/$gid=\"$desc\"";
        push @{$jData->{genes}}, $gt;
    }

    my @exEx = sort { $a->{s} <=> $b->{s} } values %exonExon;
    foreach my $xx (@exEx) {
        # Record the Exon-Exon boundary sequences
        my ($s, $e, $str) = ($xx->{s}, $xx->{e}, $xx->{str});
        my $xxBase = sprintf("%s:%d^%d", $show, $s, $e);
        my @seqs;
        foreach my $src (@{$xx->{src}}) {
            my ($rid, $orid, $sym, $rl2, $rr1, $num) = @{$src};
            unless ($orid) {
                die $args->branch($src);
            }
            my $w = $rr1 - $rl2 + 1;
            unless ($w == 2) {
                # $args->msg("[?]","$rid Exon-Exon = $rl2-$rr1");
            }
            my $sl  = $rdb->length($orid);
            unless ($sl) {
                $args->msg_once("[?]","Failed to recover RNA size for $orid");
                next;
            }
            my $rl1 = $rl2 - $xxMod;
            $rl1    = 1 if ($rl1 < 1);
            my $rr2 = $rr1 + $xxMod;
            $rr2    = $sl if ($rr2 > $sl);
            my $seq = uc($rdb->seq($orid, $rl1, $rr2) || "");
            unless ($seq) {
                $args->msg_once("[?]","Failed to recover sequence $rl1-$rr2 for $orid");
                next;
            }
            my $sd;
            for my $i (0..$#seqs) {
                if ($seqs[$i]{seq} eq $seq) {
                    # There is already an entry that matches this sequence
                    # Hopefully this is the most common case when two
                    # or more RNAs overlap a junction
                    $sd = $seqs[$i];
                    last;
                }
            }
            unless ($sd) {
                $sd = {
                    seq => $seq,
                    pos => $rl2 - $rl1 + 1,
                    w   => $w,
                };
                push @seqs, $sd;
            }
            $sd->{sym}{$sym} = 1 if ($sym);
            $sd->{ex}{$rid} = $num;
        }
        for my $j (0..$#seqs) {
            my $sd   = $seqs[$j];
            my $xxid = $xxBase . ($#seqs == 0 ? '' : '.'.($j+1)).".F";
            my $pos  = $sd->{pos};
            my $ex   = $sd->{ex};
            my $w    = $sd->{w};
            my @syms = sort keys %{$sd->{sym}};
            my $xdat = {
                id  => $xxid,
                str => $str,
                s   => $s,
                e   => $e,
                par => $id,
                pos => $pos,
                ex  => $ex,
                w   => $w,
            };
            my @descs;
            push @descs, "[".join(' ', @syms)."]" unless ($#syms == -1);
            push @descs, "P=$pos";
            push @descs, "W=$w" unless ($w == 2);
            push @descs, "ON=$id";
            foreach my $rid (sort keys %{$ex}) {
                push @descs, "$rid=$ex->{$rid}";
            }
            my $seq = $sd->{seq};
            my $sl  = CORE::length($seq);
            substr($seq, $pos - 1, $w) = lc(substr($seq, $pos - 1, $w));
            printf( $fh ">%s %s\n", $xxid, join(' ', @descs));
            unless ($noSeq) {
                for (my $p = 0; $p < $sl; $p += $block) {
                    print $fh substr($seq, $p, $block)."\n";
                }
            }
            my $json      = encode_json( $xdat );
            $dbf->{$xxid} = $json;
        }
    }
    
    if ($#exEx > -1) {
        $jData->{exex} = $#exEx + 1;
    }
    # Write the metadata index for the main gDNA:
    my $json    = encode_json( $jData );
    $dbf->{$id} = $json;

    my @ug     = sort values %uGene;
    my $seq    = $fdb->seq($cid, $s, $e) || "";
    my $expect = $e - $s + 1;
    my $sl     = CORE::length($seq);
    if ($sl != $expect) {
        unshift @descs, "/error='Incorrect number of bases recovered'";
        $args->msg("[!]","$id recovers $sl bases, not $expect");
    }
    unshift @descs, "[".join(' ', @ug)."]";
    push    @descs, "/len='$sl'";

    print $fh ">$id ".join(' ', @descs)."\n";
    unless ($noSeq) {
        for (my $p = 0; $p < $sl; $p += $block) {
            print $fh substr($seq, $p, $block)."\n";
        }
    }
    return $id;
}

sub test {
    my $id     = shift;
    my $file   = $args->val(qw(db fasta));
    my $faMeta = BMS::FastaMeta->new(-fasta => $file );
    $args->msg("[+]","Metadata interface", $faMeta->fasta, $faMeta->dbfile);
    if (my $obj = $faMeta->fetch($id)) {
        print $obj->to_text();
        print $args->branch($obj->data());
    } else {
        $args->msg("[?]","Unknown request '$id'");
    }
    die "TESTING";
}

sub capture {
    my ($re) = @_;
    my @dbs = $args->each_split_val(qw(file fasta db));
    if ($#dbs == -1) {
        $args->msg("[?]", "Can not extract sequences",
                   "Provide one or more databases using -db");
        exit;
    }
    my $outFile = $args->val('output') || "$re.fa";
    my $mfile   = $outFile . ".dbfile";
    unlink($mfile) if (-e $mfile);
    open(OUT, ">$outFile") || $args->death
        ("Failed to write fasta", $outFile, $!);
    $args->msg("[+]", "Capturing sequences with '$re'");
    my $debug = $args->val('debug');
    my $ofh = *STDOUT;
    my %dbCap;
    foreach my $db (@dbs) {
        $args->msg("[<]", $db );
        my $meta = BMS::FastaMeta->new(-fasta => "$db.dbfile" );
        unless (open(IN, "<$db")) {
            $args->warn("Failed to read fasta", $db, $!);
            next;
        }
        my $cap = 0;
        while (<IN>) {
            if (/^>(\S+)/) {
                my $id = $1;
                if (/$re/i) {
                    $cap = 1;
                    $args->msg("[+]", $id);
                    my $obj = $meta->fetch($id);
                    # die $args->branch($obj->data()) if ($id eq 'NM_001290095.1');
                    my $sj  = $obj->supporting_json();
                    while (my ($key, $json) = each %{$sj}) {
                        if (my $prior = $dbCap{$key}) {
                            if ($prior ne $json) {
                                $args->msg('[?]',"Inconsistent metadata for $key", $prior, $json);
                            }
                        } else {
                            $dbCap{$key} = $json;
                        }
                    }
                } else {
                    $cap = 0;
                }
            }
            print OUT $_ if ($cap);
        }
        close IN;
    }
    close OUT;
    
    my (%dbf, @msg);
    tie(%dbf, 'DB_File', $mfile, O_CREAT|O_RDWR, 0644, $DB_BTREE ) ||
        $args->death("Failed to tie DBFILE hash", $mfile, $!);
    while (my ($key, $json) = each %dbCap) {
        next unless ($json);
        $dbf{$key} = $json;
        push @msg, sprintf("%15s: %s", $key, $json);
    }
    $args->msg("[DEBUG]", @msg) if ($debug);
    untie %dbf;
    $args->msg("Done", $outFile, $mfile);
    exit;
}

