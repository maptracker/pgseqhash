#!/stf/biobin/perl -w

BEGIN {
    # Allows usage of beta modules to be tested:
    my $prog = $0; my $dir = `pwd`;
    if ($prog =~ /working/ || $dir =~ /working/) {
	warn "\n\n *** This is Beta Software ***\n\n";
	require lib;
	import lib '/stf/biocgi/tilfordc/perllib';
    }
    require lib;
    # import lib '/stf/biocgi/tilfordc/released/Bio/SeqIO';
    require BMS::ErrorInterceptor;
    my $foo = BMS::ErrorInterceptor->new();
    $foo->ignore_error('Replacement list is longer than search list');
}

=head1 Usage

 alias lsh /stf/biocgi/tilfordc/working/bmscat/loadSeqHash.pl

=head2 Loading gene metadata

Set the genomic boundaries and exonic structure of genes:

 lsh -genemeta /path/to/file.fa

Read a TSV file to capture aliases associated with genes:

 lsh -genealias /path/to/metadata.tsv

=cut
    
use strict;

use BMS::ArgumentParser;
use BMS::PgSeqHash;
use BMS::FastaMeta;
use BMS::ForkCritter;
use BMS::Utilities::SequenceUtilities;
use Bio::SeqIO;
use JSON;

my $args = BMS::ArgumentParser->new
    ( -testmode  => 1,
      -verbose   => 1,
      -limit     => 0,
      -nocgi     => 1,
      -wordsize  => 12,
      -progress  => 120,
      -fork      => 30,
      -prefix    => 2,
      -maxexpect => 0,
      -maxword   => 0,
      -maxambig  => 50,
      -benchmark => 1,
      -tempdir   => '/scratch/PgSeqHash',
      -paramalias => {
          wordsize  => [qw(ws wordlen wordlength)],
          limit     => [qw(lim)],
          verbose   => [qw(vb)],
          progress  => [qw(prog)],
          fork      => [qw(forknum numfork)],
          benchmark => [qw(bench dobench showbench)],
          fasta     => [qw(input db file)],
          tempdir   => [qw(tmpdir tmp temp)],
          wordsize  => [qw(ws)],
          maxambig  => [qw(maxexpand)],
          trial     => [qw(istrial dotrial)],
          genemeta  => [qw(dogene gene meta)],
          split     => [qw(splitter separator)],
          rsg       => [qw(refseqgene)],
          clobber   => [qw(overwrite)],
          rebuild   => [qw(builddb)],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
      });


my $limit         = $args->val('limit')    || 0;
my $progress      = $args->val('progress') || 300;
my $vb            = $args->val('verbose');
my $fnum          = $args->val('fork')      || 15;
my $tmpdir        = $args->val('tempdir')   || "/tmp/PgSeqHash/";
my $ws            = $args->val('wordsize')  || 12;
my $inputFile     = $args->val('fasta');
my $maxWord       = $args->val('maxword')   || 0;
my $maxExpect     = $args->val('maxexpect') || 0;
my $maxExpand     = $args->val('maxambig');
my $isTrial       = $args->val('trial')     || 0;
my $clobber       = $args->val('clobber')   || 0;
my $block         = 100;
my $bigNum        = 999999999999;
my $smlNum        = $bigNum * -1;
my $addDetail     = $args->val('adddetail');
my $tmpLoc        = "/scratch";
my $su            = BMS::Utilities::SequenceUtilities->new();

$tmpdir =~ s/\/+$//;
$args->assure_dir($tmpdir);
$args->shell_coloring();

my $strMap = { F => 1, R => -1, 1 => 1, -1 => -1 };
my $lastRpt = 0;
my ($fc, $psh, $hitSTH, $wrdSTH, $whSTH, $wsrcSTH, 
    $globals, $srcID, $maxExpNum, $outputFastaFh);

&build_all();
&load_refseq_gene();
&load_gene_aliases( $args->val('genealias') );
&process_file();
# &show_bench();

sub show_bench {
    my $bench = $args->val('benchmark');
    return "" unless ($bench);
    my $txt = $su->show_bench( -shell => 1 );
    warn $txt;
}

sub global_info {
    my $file = shift;
    return undef unless ($file);
    my $short   = $file;
    $short      =~ s/^.*\///;
    my $workDir = sprintf("%s/%s", $tmpdir, $short);
    $args->assure_dir($workDir);
    my $meta = {
        fasta => $file,
        short => $short,
        work  => $workDir,
    };
    return $meta;
}

sub load_refseq_gene {
    my $dir = $args->val('rsg');
    return unless ($dir);
    my $pattern = 'refseqgene.\d+.genomic.gbff.gz';
    my @rsg = $args->read_dir( -dir => $dir, -keep => $pattern );
    
    $args->death("No RefSeqGene entries found in directory:", $dir,
                 "Expecting files of format:", $pattern) if ($#rsg == -1);
    $args->msg("[<]", scalar(@rsg)." RefSeqGene files found:", @rsg,`date`);

    $args->death("To load RefSeqGene you must provide the path to the simplified fasta file with -fasta") unless ($inputFile);
    if (-e $inputFile && -s $inputFile && !$clobber) {
        $args->msg("Existing output file will not be overwritten",
                   "Pass -clobber to replace it with newly parsed data",
                   $inputFile);
    } else {
        $args->msg("[>]", "Sequences extracted will be written to:",
                   $inputFile);
        open($outputFastaFh, ">$inputFile") || $args->death
            ("Failed to write RefSeqGene to fasta", $inputFile, $!);
        my $weirdF = "$inputFile-Issues.txt";
        open(WEIRD, ">$weirdF") || $args->death
            ("Failed to write Issues file", $weirdF, $!);
        select((select(WEIRD), $| = 1)[$[]); # Bob's call to autoflush
    }

    $args->ignore_error("import is deprecated");

    $psh = &dbh();
    foreach my $rf (@rsg) {
        my $rfh;
        open($rfh, "gunzip -c '$rf'|") || $args->death
            ("Failed to gunzip RefSeqGene file", $rfh, $!);
        $args->msg("[<]", $rf);
        &_process_rsg_file($rfh);
        close $rfh;
    }
    if ($outputFastaFh) {
        close $outputFastaFh;
        $args->msg("[>]", "New fasta file created", $inputFile);
    }
    $args->msg("Finished.",`date`);
    exit;
}

sub _process_rsg_file {
    my ($rfh) = @_;
    my $sio = Bio::SeqIO->new( -format => 'genbank', -fh => $rfh, );
    my %amCache;
    my $bsnum = 0;
    # We will not gather information from some genbank tags:
    my $ignoreTag = { 
        map { $_ => 1 } 
        qw(source STS misc_difference enhancer LTR misc_feature precursor_RNA) };
    while (my $bs = $sio->next_seq()) {
        my ($acc, $vers) = ($bs->accession_number, $bs->version);
        unless ($vers) {
            $args->msg("[?]","Skipping unversioned accession $acc!");
            next;
        }
        my $vacc = "$acc.$vers";
        $args->msg("[+]", $vacc) unless ($bsnum % 100);
        $bsnum++;
        my %data = ( var => {} );

        $args->death("Multiple transcripts for some genes",
                     "Need to get exon data from mRNA features",
                     "Do ncRNAs have different tags?",
                     "Should record overlapping gene footprints as tagval",
                     "Need to make fasta, with variant ambiguities",
                     "Need to index sequence") if (0);
        
        foreach my $feat ($bs->get_SeqFeatures()) {
            my $pt = $feat->primary_tag();
            if ($pt eq 'variation') {
                my ($s, $e) = ($feat->start(), $feat->end());;
                # Only consider SNPs
                next unless ($s == $e);
                # Only consider SNPs with a dbSNP ID:
                my %idH;
                foreach my $xr (&_tag_vals($feat, 'db_xref')) {
                    if ($xr =~ /^dbSNP:(\d+)$/) {
                        $idH{$1} = 1;
                    }
                }
                my ($id) = sort { $a <=> $b } keys %idH;
                next unless ($id);
                my $targ = $data{var}{$s} ||= [ {}, {} ];
                $targ->[1]{$id} = 1;
                foreach my $allele (map {uc($_)} &_tag_vals($feat, 'replace')) {
                    if ($allele =~ /^[A-Z]$/) {
                        $targ->[0]{$allele} = 1;
                    }
                }
                
            } elsif ($pt eq 'mRNA' || $pt eq 'ncRNA' || $pt eq 'misc_RNA') {
                # This is where we are going to gather the main alignments
                my $sym = &_gene_symbol($feat);
                next unless ($sym);
                push @{$data{rnaid}{$sym}}, $feat;
                if ($pt eq 'mRNA') {
                    # Set up to capture the next CDS
                    push @{$data{mrna}{$sym}}, $feat;
                }
            } elsif ($pt eq 'exon') {
                # The exon features describe the "Primary" gene stored
                # by this entry
                my @genes = &_tag_vals($feat, 'gene');
                if ($#genes == -1) {
                    &note_weird("[!]\t$vacc\t\tExon without gene");
                    next;
                } elsif ($#genes != 0) {
                    &note_weird("[!]\t$vacc\t".join(',',@genes)."\tExon with multiple genes");
                    next;
                }
                my $sym  = uc($genes[0]);
                # We were initially presuming the exon entries were the main
                # data source
                # Now we just parse them to determine what the primary gene
                # in this NG entry is.
                $data{exons}{$sym}++;
            } elsif ($pt eq 'gene') {
                my $sym = &_gene_symbol($feat);
                next unless ($sym);
                my %gidH;
                foreach my $xr (&_tag_vals($feat, 'db_xref')) {
                    if ($xr =~ /^GeneID:(\d+)$/) {
                        $data{locid}{$sym}{$1}++;
                    }
                }
            } elsif ($pt eq 'CDS') {
                my $sym = &_gene_symbol($feat);
                next unless ($sym);
                # Sigh. CDS features do not include an internal link back
                # to the relevant transcript ID. We will ASSUME that this CDS
                # should go to the earliest un-claimed mRNA
                if (my $rfeat = shift @{$data{mrna}{$sym} || []}) {
                    $rfeat->{_CAT_CDS} =[ $feat->start(), $feat->end() ];
                } else {
                    &note_weird("[!]\t$vacc\t$sym\tCDS entry could not find mRNA");
                }
            } elsif ($ignoreTag->{$pt}) {
                # Ignore these
                $ignoreTag->{$pt}++;
            } else {
                $args->msg_once("[?]","Unanticipated tag '$pt'");
            }
            # warn "   $pt\n";
        }
        foreach my $sym (keys %{$data{locid} || {}}) {
            my @u = sort {$a <=> $b} keys %{$data{locid}{$sym}};
            if ($#u == 0) {
                $data{locid}{$sym} = $u[0];
            } else {
                &note_weird("[!]\t$vacc\t$sym\tMultiple GeneIDs\t".
                            join(',', @u));
                $data{locid}{$sym} = "";
            }
        }

        my @syms = keys %{$data{exons} || {}};
        if ($#syms == -1) {
            &note_weird("[!]\t$vacc\t\tNo primary gene entries" );
            next;
        }
        my $allSym = join(',', @syms);
        if ($#syms != 0) {
            &note_weird("[?]\t$vacc\t$allSym\tMultiple gene entries" );
        }
        my $sbjID = $psh->sbj_id( $vacc, $srcID );

        # Nearly all entries have a single gene represented
        # A few have multiple, though
        # eg: NG_012989.1     PEG3,ZIM2       Multiple primary gene entries
        foreach my $sym (@syms) {
            my @rfeats = @{$data{rnaid}{$sym} || []};
            delete $data{rnaid}{$sym};
            my $lid = $data{locid}{$sym};
            unless ($lid) {
                &note_weird("[!]\t$vacc\t$sym\tNo GeneID" );
                next;
            }
            # These scores will not really always be 100%
            # We are not provided with alignment quality information, but
            # will take these locations as 'canonical' and report them as 100:
            my $score = 100;
            my @rnaJsons;
            my ($firstL, $lastR);
            foreach my $rfeat ( @rfeats ) {
                my ($rid, $iderr) = &_unique_tagval( $rfeat, 'transcript_id' );
                if ($iderr) {
                    &note_weird("[!]\t$vacc\t$sym\t$iderr" );
                    next;
                }
                if ($rfeat->strand() < 0) {
                    &note_weird("[!!]\t$vacc\t$sym\tPrimary RNA on opposite strand!\t$rid");
                    next;
                }
                my @locs  = sort { $a->[0] <=> $b->[0] } &_location( $rfeat );
                $firstL = $locs[0][0] if
                    (!defined $firstL || $firstL > $locs[0][0]);
                $lastR = $locs[-1][1] if
                    (!defined $lastR || $lastR < $locs[-1][1]);
                my $prior = 0;
                my @hspJson;
                foreach my $loc (@locs) {
                    my ($l, $r) = @{$loc};
                    my $s = $prior + 1;
                    my $e = $prior = $s + $r - $l - 2;
                    push @hspJson, "[$l,$r,$s,$e]";
                }
                
                my $rjson = sprintf('["%s",%s,[%s]', $rid, $score,
                                    join(',', @hspJson));
                if ($rfeat->primary_tag() eq 'mRNA') {
                    if (my $cds = $rfeat->{_CAT_CDS}) {
                        $rjson .= sprintf(",null,null,[%d,%d]", @{$cds});
                    } else {
                        &note_weird("[!]\t$vacc\t$sym\tNo CDS associated with mRNA\t$rid");
                    }
                }
                $rjson .= ']';
                push @rnaJsons, $rjson;
            }
            unless (defined $firstL) {
                &note_weird("[!]\t$vacc\t$sym\tNo coordinates defined" );
                next;
            }
            
            my $gid   = $psh->gene_id( "LOC".$lid );
            # push @geneLoad, [ $gid, $l, $r, 1, $sc, $json];
            my $gjson = sprintf('{"rna":[%s]}', join(',', @rnaJsons));
            $psh->set_gene_loc_details
                ( $sbjID, $gid, $firstL, $lastR, 1, $score, $gjson );
        }

        # Note the non-primary loci. The primary ones should have been removed
        # from {rnaid} in the loop above.
        my $overlap;
        foreach my $osym (sort keys %{$data{rnaid}}) {
            my $olid = $data{locid}{$osym};
            next unless ($olid);
            $overlap ||= {};
            my $targ = $overlap->{"LOC$olid"} ||= {
                sym  => $osym,
            };
            foreach my $ofeat ( @{$data{rnaid}{$osym}} ) {
                my @locs  = &_location( $ofeat );
                my ($rid, $iderr) = &_unique_tagval( $ofeat, 'transcript_id' );
                if ($rid) {
                    my $rtarg = $targ->{rnas} ||= {};
                    $rtarg->{$rid} = \@locs;
                    $targ->{str} ||= $ofeat->strand() + 0;
                };
            }
        }

        my %dbSnp;
        my $sd = $outputFastaFh ? $bs->seq() : "";
        my $varNum = 0;
        while (my ($pos, $varDat) = each %{$data{var}}) {
            my $ind = $pos - 1;
            my ($charH, $idH) = @{$varDat};
            my @alls = sort keys %{$charH};
            # Only consider variants with two or more alleles:
            next if ($#alls < 1);
            my $chk  = uc(substr($sd, $ind, 1));
            unless ($charH->{$chk}) {
                &note_weird("[?]\t$vacc\t$allSym\tReference allele not represented in variant\t$pos\t$chk vs ".join('/',@alls)."\t".join(' ',map { "rs$_"} sort { $a <=> $b } keys %{$idH}));
                push @alls, $chk;
            }
            if ($sd) {
                # Update the sequence
                my $akey  = join('',@alls);
                my $ambig = $amCache{$akey} ||= $su->ambiguous_code( \@alls );
                substr($sd, $ind, 1) = $ambig;
            }
            $dbSnp{$pos + 0} = [ sort { $a <=> $b } map { $_ + 0 } keys %{$idH} ];
            $varNum++;
        }
        
        if ($outputFastaFh) {
            my $slen = CORE::length($sd);
            print $outputFastaFh ">$vacc [$allSym]";
            if (my $d = $bs->desc()) {
                print $outputFastaFh " $d";
            }
            if ($varNum) {
                printf($outputFastaFh " /dbSNPcount=%d", $varNum);
            }
                    
            print $outputFastaFh "\n";
            for (my $i = 0; $i < $slen; $i += $block) {
                print $outputFastaFh substr($sd, $i, $block);
                print $outputFastaFh "\n";
            }
        }
        $psh->clear_tag( $sbjID, "OtherGenes" );
        $psh->set_tag( $sbjID, "OtherGenes", encode_json($overlap) )
            if ($overlap);
    }
    undef $sio;
    my @ignored;
    foreach my $it (sort { $ignoreTag->{$b} <=> 
                               $ignoreTag->{$a} } keys %{$ignoreTag}) {
        if (my $num = $ignoreTag->{$it} - 1) {
            push @ignored, "$it : $num";
        }
    }
    $args->msg("[-]","Ignored tags:", @ignored) unless ($#ignored == -1);
}

sub note_weird {
    my ($msg) = @_;
    return unless ($outputFastaFh);
    print WEIRD "$msg\n";
}

sub _unique_feature {
    my ($array, $type) = @_;
    return (undef, "no $type entries") if (!$array || $#{$array} == -1);
    return (undef, "multiple $type entries") if ($#{$array} > 0);
    return ($array->[0]);
}

sub _unique_tagval {
    my ($feat, $tag) = @_;
    my %valH = map { $_ => 1 } &_tag_vals($feat, $tag);
    my @u    = keys %valH;
    return ("", "No $tag entries") if ($#u == -1 || !$u[0]);
    return ("", "Multiple $tag values: ".join(' + ', @u)) if ($#u > 0);
    return $u[0];
}

sub _gene_symbol {
    my $feat = shift;
    my @genes = &_tag_vals($feat, 'gene');
    return $#genes == 0 ? uc($genes[0]) : "";
}

sub _tag_vals {
    my ($feat, $tag) = @_;
    return $feat->has_tag($tag) ? $feat->get_tag_values( $tag ) : ();
}

sub _location {
    my $feat = shift;
    my @rv;
    foreach my $loc ($feat->location->each_Location()) {
        push @rv, [ $loc->start - 1, $loc->end + 1 ];
    }
    return @rv;
}

sub initialize {
    $fc->ignore_error('Replacement list is longer than search list');
    $psh     = &dbh();
    # $hitSTH  = $psh->_set_hit_sth();
    $wrdSTH  = $psh->_set_word_sth();
    # $whSTH   = $psh->_set_wordhit_sth();
    $wsrcSTH = $psh->_set_srchit_sth();
}

sub finalize {
    
}

sub dbh {
    my $psh    = BMS::PgSeqHash->new();
    $srcID   ||= $psh->src_id( $inputFile, $ws );
    return $psh;
}

sub build_all {
    return unless ($args->val('rebuild'));
    $args->msg("[.]", "Rebuilding database schema");
    my $psh = BMS::PgSeqHash->new();
    my $dbh = $psh->dbh();
    $dbh->make_all();
    $psh->disconnect();
    $args->msg("[.]", "Finished");
}


sub get_db {
    return BMS::PgSeqHash->new();
}

sub process_file {
    if (my $gm = $args->val('genemeta')) {
        $gm = $args->val('fasta') if ($gm eq '1' && $args->val('fasta'));
        &load_gene_meta($gm);
        return;
    }
    return undef unless ($inputFile);
    my $short   = $inputFile;
    $short      =~ s/^.*\///;
    my $workDir = sprintf("%s/%s", $tmpdir, $short);
    $workDir   .= "-LIMIT$limit" if ($limit);
    $args->assure_dir("$workDir/words");
    $args->msg("[+]", "Working directory",$workDir);
    $globals = &_make_word_file( $inputFile, $workDir);
    if ($args->val('profile')) {
        &profile_hits();
    } else {
        &load_words();
    }
    # die $args->branch($meta);
}

sub profile_hits {
    my $expect    = $globals->{expect};
    my $dir       = $globals->{dir};
    my $countFile = "$dir/Counts.tsv";
    $args->msg("[<]","Building hit profile", $countFile);
    unless (-s $countFile) {
        $args->msg("[!]","Can not profile - count file is absent");
        exit;
    }
    open(CF, "<$countFile") || $args->death
        ("Failed to read count file", $countFile, $!);
    my @counts;
    my $tot = 0;
    while (<CF>) {
        s/[\n\r]+$//;
        my ($cnt, $num) = split(/\t/);
        $counts[$cnt] += $num;
        $tot          += $num * $cnt;
    }
    close CF;
    my $x = 0;
    my $runSum = 0;
    my @msg     = ("Total of $tot oligo locations",
                   "  Expect #Repeat  %Excl (#Excluded");
    my $priExcl = -1;
    for my $cnt (1..$#counts) {
        my $num  = $counts[$cnt] || 0;
        $runSum += $cnt * $num;
        my $fold = $cnt / $expect;
        if ($fold >= $x || $cnt == $#counts) {
            my $excl = $tot - $runSum;
            unless ($excl == $priExcl) {
                my $percExcl = 100 * $excl / $tot;
                push @msg, sprintf("  %5dx %7d %6.2f (%d)", 
                                   $fold, $cnt, $percExcl, $excl);
            }
            $priExcl = $excl;
            $x += $x >= 100 ? 100 : 10;
        }
        # last if ($fold >= 1000);
    }
    $args->msg("[+]", @msg);
    exit;
}

sub load_words {
    my $dir  = $globals->{dir};
    $fc = BMS::ForkCritter->new
        ( -init_meth   => \&initialize,
          -finish_meth => \&final_build,
          -input_type  => 'array',
          -method      => \&_parse_word_file,
          -limit       => $limit,
          -progress    => $progress,
          -verbose     => $vb );
    $fc->reset();
    $fc->input( $globals->{files} );

    my $num    = $#{$globals->{files}} + 1;
    my @msg    = ("Loading $num word files into DB");
    my $expect = $globals->{expect};
    push @msg, sprintf("Expecting %s hits per word in %s",
                       $expect, $su->size_with_unit( $globals->{size} ));
    # $maxExpNum = $maxExpect ? int(0.999 + $expect * $maxExpect) : 0;
    my $exclMsg;
    $maxExpNum = 0;
    if ($maxExpect) {
        $maxExpNum = int(0.999 + $expect * $maxExpect);
        $exclMsg = "> $maxExpNum (${maxExpect}x)";
    }
    if ($maxWord) {
        if ($maxWord > $maxExpNum) {
            $maxExpNum = $maxWord;
            $exclMsg = "> $maxExpNum";
        }
    }
    if ($exclMsg) {
        if ($addDetail) {
            push @msg, "**ONLY** loading if hits $exclMsg";
        } else {
            push @msg, "Excluding oligos $exclMsg hits";
        }
    } elsif ($addDetail) {
        $args->msg("Request to -adddetail, but neither -maxexpect or -maxword are set","Exitting.");
        exit;
    }
    my ($countFile);
    if ($globals->{limit} || $isTrial) {
        if (my $limit = $globals->{limit}) {
            push @msg, "These data were limited to $limit entries";
        } else {
            push @msg, "Trial mode, only counting";
        }
        $countFile = "$dir/Counts.tsv";
        $fc->output_file( 'COUNT', "$countFile.temp" );
    }
    $args->msg("[<]", @msg);

    my $exclFile = "$dir/Excluded.tsv";
    $fc->output_file( 'EXCL', $exclFile );
    if (my $failed = $fc->execute( $fnum )) {
        $args->death("$failed child processes failed!");
    }
    @msg = ("Processing finished");
    if ($countFile) {
        my $temp = "$countFile.temp";
        if (-s $temp) {
            open(TEMP, "<$temp") || $args->death
                ("Failed to read count file", $temp, $!);
            my %counts;
            while (<TEMP>) {
                s/[\n\r]+$//;
                my ($cnt, $num) = split(/\t/);
                $counts{$cnt} += $num;
            }
            close TEMP;
            open(COUNT, ">$countFile") || $args->death
                ("Failed to write count file", $countFile, $!);
            print COUNT &_write_counts( \%counts );
            close COUNT;
            push @msg, "Counts: $countFile";
        } else {
            $args->msg("[?]","Failed to find temp count file", $temp);
        }
    }
    if (-s $exclFile) {
        # Sort with most common words at top
        my $temp = "$exclFile.temp";
        my $cmd = sprintf("sort -nr -k2 -T \"%s\" \"%s\" > \"%s\"",
                          $tmpLoc, $exclFile, $temp);
        system($cmd);
        if (-s $temp) {
            system("mv \"$temp\" \"$exclFile\"");
        } else {
            $args->msg("[?]", "Failed to sort excluded file", $exclFile);
        }
        push @msg, "Excluded: $exclFile";
    }
    $args->msg("[>]", @msg);
}

sub _write_counts {
    my ($cH) = @_;
    my $txt = "";
    foreach my $count (sort { $a <=> $b } keys %{$cH}) {
        my $num = $cH->{$count};
        $txt .= "$count\t$num\n";
    }
    return $txt;
}

sub final_build {
    if (my $cH = $globals->{wdCnt}) {
        $fc->write_output("COUNT",  &_write_counts( $cH ) );
    }
}

sub _parse_word_file {
    my $file = shift;
    my $dir  = $globals->{dir};
    my $path = "$dir/words/$file";
    my $srt  = "$path.sort";
    unless (-s $srt) {
        $args->death("Failed to find word file", $path) unless (-s $path);
        my $sz  = "1G";
        my $cmd = sprintf("sort -k1 -T \"%s\" -S %s \"%s\" > \"%s\"",
                          $tmpLoc, $sz, $path, $srt);
        system($cmd);
        if (-s $srt) {
            unlink($path);
        } else {
            $args->death("Failed to generate sorted file", $srt, $cmd);
        }
    }
    open(SRT, "<$srt") || $args->death("Failed to read sorted hits", $srt, $!);
    my $data;
    my $method = ($globals->{limit} || $isTrial) ? 
        \&_summarize_word_tidbit : \&_load_word_tidbit;
    my $done = 0;
    my $loadLimit = $args->val('loadlimit');
    while (<SRT>) {
        s/[\n\r]+$//;
        my ($wrd, $sid, $pos, $isNrm) = split(/\t/);
        if (!$data || $data->{word} ne $wrd) {
            $data = &{$method}( $data, $wrd );
            $done++;
            last if ($loadLimit && $done > $loadLimit);
        }
        if ($pos < 0) {
            $args->msg("[!!]","Illegal postion for $wrd in $sid : $pos");
        } else {
            push @{$data->{sbj}{$sid}}, $pos + 1;
            $data->{count} += $isNrm;
        }
    }
    &{$method}( $data );
    close SRT;
}

sub _summarize_word_tidbit {
    my ($data, $wrd) = @_;
    if ($data) {
        $globals->{wdCnt}{ $data->{count} }++;
        &_exclude_word($data);
    }
    return &_new_data( $wrd );    
}

sub _load_word_tidbit {
    my ($data, $wrd) = @_;
    if ($data) {
        my @hitbits;
        my $skipDetails = &_exclude_word($data);
        my $count       = $data->{count};
        my $word        = $data->{word};
        if ($addDetail) {
            # Backfilling previously excluded hit detail
            # Ignore data that would have been excluded
            # (assume it has been loaded earlier)
            return &_new_data( $wrd ) unless ($skipDetails);
            # Unset the flag so the repeats can get added in:
            $skipDetails = 0;
            # $args->msg("[-]","Backfill $word : $count");
        }
        unless ($skipDetails) {
            # Make note of every location
            my @sids  = sort { $a <=> $b } keys %{$data->{sbj}};
            for my $s (0..$#sids) {
                my $sid = $sids[$s];
                # Separate multiple sbj_id values with zeros:
                push @hitbits, 0 if ($s);
                # Add the sbj_id
                push @hitbits, $sid;
                # Add all the hit locations for that subject:
                push @hitbits, @{$data->{sbj}{$sid}};
            }
        }
        my $wid = $wrdSTH->get_single_value( $word );
        $wsrcSTH->execute($wid, $srcID, $count, \@hitbits );
        # $args->msg("[-]","wrd_id = $wid '$word' : $count ".($#hitbits+1));
        if (my $err = $wsrcSTH->err()) {
            print $err;
        }
    }
    return &_new_data( $wrd );
}

sub _new_data {
    my ($wrd) = @_;
    return {
        word  => $wrd,
        sbj   => {},
        count => 0,
    };
}

sub _exclude_word {
    return 0 unless ($maxExpNum);
    my $data  = shift;
    my $count = $data->{count};
    if ($count > $maxExpNum) {
        my $word = $data->{word};
        $fc->write_output("EXCL", "$word\t$count\n");
        return $count;
    } else {
        return 0;
    }
}

sub _make_word_file {
    my ($file, $wdDir) = @_;
    my $meta = sprintf("%s/meta", $wdDir);
    $meta .= ".json";
    if (-s $meta) {
        my $txt = "";
        open(META, "<$meta") || $args->death
            ("Failed to read metadata", $meta, $!);
        while (<META>) {
            $txt .= $_;
        }
        close META;
        my $rv = decode_json($txt);
        return $rv;
    }

    my $psh = &dbh();
    my $amFile = sprintf("%s/Ambiguity.tsv", $wdDir);
    open(AMF, ">$amFile") || $args->death
        ("Failed to make ambiguity file", $amFile, $!);
    my $prefixLen = 5;
    my $output    = {
        dir    => $wdDir,
        files  => [],
        fh     => [],
        size   => 0,
        masked => 0,
        exon   => 0,
        normal => 0,
        srcid  => $srcID,
        limit  => $limit,
        prefix => $prefixLen,
    };
    my %fhs;
    $args->msg("[<]", "Scanning input database", $file);
    my $instream = Bio::SeqIO->new( -file   => "<$file",
                                    -format => 'fasta');
    if (!defined $maxExpand) {
        $maxExpand = 2 ** $ws;
        $args->msg("[-]","Allowing at most $maxExpand variants per word",
                   "(Two-allele SNP at every base)");
    } elsif ($maxExpand) {
        $args->msg("[-]","User-defined limit of $maxExpand variants per word");
    }
    while ( my $bs  = $instream->next_seq( ) ) {
        my $acc     = $bs->display_id();
        my $isNrm   = $acc =~ /\^/ ? 0 : 1;
        my $seq     = uc($bs->seq());
        my $len     = CORE::length($seq);
        my $sbjID   = $psh->sbj_id( $acc, $srcID, $len );
        my $masked  = 0;
        my $mWords  = 0;
        while ($seq =~ /([NX]{3,})/) {
            my $nRun = $1;
            my $mLen = CORE::length($nRun);
            my $mask = 'x' x $mLen;
            $seq =~ s/$nRun/$mask/;
            $masked += $mLen;
        }
        if ($isNrm) {
            $output->{normal}++;
            $output->{size}   += $len - $masked;
            $output->{masked} += $masked;
        } else {
            $output->{exon}++;
        }
        my $wLen = $len - $ws;
        my $ambig;
        my %words;
        for (my $i = 0; $i <= $wLen; $i++) {
            my $baseWord = substr($seq, $i, $ws);
            # Ignore word with masked characters
            if ($baseWord =~ /x/) {
                $mWords++;
                next;
            }
            my @expand = $su->expand_ambiguous_sequence($baseWord);
            my $exNum  = $#expand + 1;
            if ($maxExpand && $exNum > $maxExpand) {
                print AMF join("\t", $baseWord, $sbjID, $i, $exNum, "MASKED")."\n";
                $mWords++;
                next;
            } elsif ($exNum > 2) {
                print AMF join("\t", $baseWord, $sbjID, $i, $exNum)."\n";
            }
            foreach my $word (@expand) {
                my $pfx = substr($word, 0, $prefixLen);
                my $fh  = $fhs{$pfx} ||= &_prefix_file($pfx, $output);
                print $fh join("\t", $word, $sbjID, $i, $isNrm)."\n";
            }
        }
        if (time - $lastRpt > $progress) {
            # Note progress
            my $num = $output->{normal} || 1;
            $args->msg("[$num]", $acc);
            $lastRpt = time;
        }
        last if ($limit && $output->{normal} >= $limit);
    }
    close AMF;
    foreach my $fh (values %fhs) {
        close $fh;
    }
    my $sz = $output->{size};
    my $expect = $output->{expect} = int(0.5 + 1000 * $sz / (4 ** $ws)) / 1000;
    $args->msg("[+]", "$sz bases total vs $ws bp words = $expect expected hits");
    open(META, ">$meta") || $args->death
        ("Failed to write metadata", $meta, $!);
    print META encode_json( $output );
    close META;
    &_set_database_meta( $output ); 
    return $output;
}

sub load_gene_meta {
    my $passed = shift;
    $inputFile = $passed if (-e $passed);
    return unless ($inputFile);
    my $mfile = "$inputFile.dbfile";
    unless (-s $mfile) {
        $args->msg("[?]", "Failed to find database file for metadata input", $mfile);
        return;
    }
    my $meta = BMS::FastaMeta->new( -fasta => $inputFile);
    my $psh  = &dbh();
    open(FASTA, "<$inputFile") || $args->death
        ("Failed to read fasta file (metadata extraction)", $inputFile, $!);
    $args->msg("[<]", "Parsing fasta metadata file", $mfile);
    my $fanum = 0;
    while (<FASTA>) {
        next unless (/^>(\S+)/);
        my $sname = $1;
        $fanum += &load_subject_id( $sname, $psh, $meta );
        last if ($limit && $fanum >= $limit);
        $args->msg("[+]", $sname) unless ($isTrial || $fanum % 1000);
    }
    close FASTA;
    foreach my $sid ($args->each_split_val('id')) {
         &load_subject_id( $sid, $psh, $meta );
    }
    $args->msg("[+]","Metadata parse complete");
}

sub load_subject_id {
    my ($sid, $psh, $meta) = @_;
    return 0 unless ($sid);
    my $subj  = $meta->fetch($sid);
    my $sbjID = $psh->sbj_id( $sid, $srcID );
    unless ($subj) {
        $args->msg("[?]","Failed to find metadata for $sid");
        return 0;
    }
    # if ($sid =~ /\^/) { next;} else { die $args->branch( $subj->data());}
    my $sdata = $subj->data();



    return 0 unless ($sdata->{ex});



    my @geneLoad;
    $args->msg("[+]", $sid) if ($isTrial);

    if (my $genes = $sdata->{genes}) {
        # We will report left coordinates relative to subject start,
        # 1-indexed. So the value left to the of the subject start is zero.
        # $leftMod will be subtracted from all left coordinates
        my $leftMod  = $sdata->{s} || 1;
        # $rightMod will be subtracted from all right coordinates
        # The position to the right of the last base in the subject will
        # be N+1 for a subject of length N
        my $rightMod = $leftMod - 2;
        foreach my $gdat (@{$genes}) {
            my $gacc = $gdat->{id};
            next unless ($gacc);
            my $gid = $psh->gene_id( $gacc );
            unless ($gid) {
                $args->msg_once("Failed to recover gene_id for '$gacc'");
                next;
            }

            my %byStr;
            foreach my $hspDat ( @{$gdat->{hsps} || []} ) {
                # Normally we expect the gene to be on only one strand
                # It is possible that INVERTED duplications could have
                # the same gene aligning to both strands, however
                my $str  = $meta->numeric_strand( $hspDat->{str} );
                my $sc   = $hspDat->{sc};
                my $racc = $hspDat->{rna};
                my $robj = $meta->fetch($racc);
                my $rdat = $robj->data();
                if (my $v = $rdat->{v}) {
                    # Version the ID
                    $racc .= ".$v" unless ($racc =~ /\.$v$/);
                }
                my $targ = $byStr{$str} ||= {
                    rng => [ $bigNum, $smlNum ],
                    rna => [],
                    sc  => 0,
                };
                $targ->{sc} = $sc if ($sc && $sc > $targ->{sc});
                my (@hspJson, $notes);
                # We will record HSPs in genome order - this will allow
                # easier analysis of hits
                my @hspDats = sort {$a->[0] <=> $b->[0]} @{$hspDat->{hsp}};
                my $eNum    = $#hspDats + 1;
                for my $hi (0..$#hspDats) {
                    my $hd  = $hspDats[$hi];
                    # Convert Genomic start,end to relative left,right:
                    my ($l, $r) = ( $hd->[0] - $leftMod, 
                                    $hd->[1] - $rightMod );
                    my @hsp = ($l, $r);
                    $targ->{rng}[0] = $l if ($targ->{rng}[0] > $l);
                    $targ->{rng}[1] = $r if ($targ->{rng}[1] < $r);
                    if (my $rs = $hd->[2]) {
                        # RNA Coordinates as well
                        if (my $re = $hd->[3]) {
                            my $diff = $re - $rs;
                            if ($diff < 10) {
                                my $en = $str < 0 ? $eNum -$hi : $hi +1;
                                $notes ||= {};
                                if ($diff < 0) {
                                    $notes->{'FlippedCoordinates'}{"Exon$en"}++;
                                    ($rs, $re) = ($re, $rs);
                                } else {
                                    $diff++;
                                    $notes->{'SuspiciousShortExon'}{"$en:${diff}bp"}++;
                                }
                            }
                            push @hsp, ($rs, $re);
                        }
                    }
                    push @hspJson, sprintf('[%s]', join(',', @hsp));
                }
                # For TANDEM duplications, we could have the same RNA
                # aligned more than once to the same strand. We will
                # record information as an array rather than a hash
                # keyed to RNA ID
                my $rjson = sprintf(' ["%s",%s,[%s]', $racc, $sc || 'null',
                                    join(',', @hspJson));
                if ($notes) {
                    my @bits;
                    foreach my $k (sort keys %{$notes}) {
                        push @bits, sprintf("/%s=\'%s\'", $k, join(',', sort keys %{$notes->{$k}}));
                    }
                    $rjson .= sprintf(',null,"%s"', join(' ', @bits));
                }
                $rjson .= ']';
                push @{$targ->{rna}}, $rjson;
            }
            # Now set up the gene information load
            foreach my $str (sort { $b <=> $a } keys %byStr) {
                my $targ = $byStr{$str};
                my ($l, $r) = @{$targ->{rng}};
                if ($l == $bigNum) {
                    # Huh. Should not happen
                    $args->msg("[?]","Failed to parse coordinate range for $gacc [$str] in $sid");
                    next;
                }
                # For some reason the FastaMeta files have duplicate 
                # RNA entries in some places
                my %u = map { $_ => 1 } @{$targ->{rna}};
                my $json = sprintf('{"rna":['."\n".'%s ]}',join
                                   (",\n", sort keys %u));
                push @geneLoad, [ $gid, $l, $r, $str, $targ->{sc}, $json ];
            }
        }
    } elsif (my $exons = $sdata->{ex}) {
        # This is an exon-exon boundary
        
        my $l  = $sdata->{pos};
        my $w  = $sdata->{w} || 2;
        my $r  = $l + $w - 1;
        # This is the strand relative to the genome
        #  (RNA strand of subject for ExonExon entries is always 1):
        my $xs = $strMap->{ $sdata->{str} || "" } || 0;
        # FastaMeta is not recording absolute L/R RNA coordinates,
        # so our little 'HSP' will only have subject l,r coords:
        my $base  = '{"type":"ExEx",';
        if (my $par = $sdata->{par}) {
            # Note the parent object that the ExonExon entry resides on:
            if (my $parID = $psh->sbj_id( $par, $srcID )) {
                $base .= sprintf('"par":%d,', $parID);
            } else {
                $args->msg_once("[?]","Failed to find ExonExon parent '$par'");
            }
        } else {
            $args->msg("[?]","ExonExon entry '$sid' lacks a parent");
        }
        $base    .= sprintf('"p":%d,"w":%d,"s":%d,"rna":['."\n",
                            $l, $w, $xs);
        my %rnas = %{$sdata->{ex} || {}};
        my @accs = sort keys %rnas;
        my %gids;
        for my $ri (0..$#accs) {
            my $racc = $accs[$ri];
            # Hm. Neglected to capture score in FastaMeta for Exon-Exon
            if (my $robj = $meta->fetch($racc)) {
                if (my $gobj = $robj->gene()) {
                    # Ok, we were able to link that RNA ID to a Gene
                    my @xxhsp = ($l, $r);
                    # I did not capture RNA coordinates or scores
                    # for the ExEx entries
                    # Try to rescue from match-up to gDNA
                    my $sc;
                    foreach my $gdna ($robj->each_gdna) {
                        my %cap = ($sdata->{s} => [], $sdata->{e} => []);
                        # warn $args->branch($gdna) if ($sid eq '1.GRCh37:15905^15906.F');
                        my $isRC = ($gdna->{STR} || 0) < 0 ? 1 : 0;
                        foreach my $hsp (@{$gdna->{HSP} || []}) {
                            # warn $args->branch($hsp) if ($sid eq '1.GRCh37:15905^15906.F');
                            for my $hi (0..1) {
                                if (my $arr = $cap{$hsp->[$hi]}) {
                                    my $grab = $isRC ? ($hi ? 0 : 1) : $hi;
                                    if (my $rcoord = $hsp->[$grab + 2]) {
                                        push @{$arr}, $rcoord;
                                    }
                                }
                            }
                            # warn $args->branch($hsp);
                        }
                        my ($f1, $f2) = values %cap;
                        if ($#{$f1} == 0 && $#{$f2} == 0) {
                            push @xxhsp, sort { $a <=> $b } ($f1->[0], $f2->[0]);
                            $sc = $gdna->{SC};
                        }
                    }
                    $sc ||= 'null';
                    my $hjson = sprintf("[[%s]]", join(',', @xxhsp));

                    my $num = $rnas{$racc} || 0;
                    if (my $v = $robj->val('v')) {
                        # Version the ID
                        $racc .= ".$v" unless ($racc =~ /\.$v$/);
                    }
                    push @{$gids{$gobj->id()}}, sprintf
                        (' ["%s",%s,%s,%d]', $racc, $sc, $hjson, $num );
                }
            }
        }
        foreach my $gacc (sort keys %gids) {
            my $sc   = undef;
            if (my $gid  = $psh->gene_id( $gacc )) {
                my $json = $base;
                $json .= join(",\n", @{$gids{$gacc}});
                $json .= ' ]}';
                push @geneLoad, [ $gid, $l, $r, 1, $sc, $json];
            } else {
                $args->msg_once("[?]","Failed to recover gene_id for '$gacc'");
            }
        }
    } elsif ($sdata->{orphan}) {
        my $len = $subj->length();
        unless ($len) {
            $args->msg_once("Failed to recover RNA length for '$sid'");
            next;
        }
        if (my $gacc = $sdata->{gene}) {
            my $gid = $psh->gene_id( $gacc );
            unless ($gid) {
                $args->msg_once("Failed to recover gene_id for '$gacc'");
                next;
            }
            my $json = sprintf('{"type":"RNA","rna":['."\n".
                               ' ["%s",100,[[0,%d,1,%d]] ]]}', $sid, $len+1, $len);
            push @geneLoad, [ $gid, 0, $len+1, 1, 100, $json];
        } else {
            $args->msg_once("Failed to recover gene for orphan '$sid'");
            next;
        }
    } else {
        $args->msg_once("[!]", "No metadata loaded for '$sid'");
        next;
    }

    if ($isTrial) {
        my $msg = "";
        foreach my $gl (@geneLoad) {
            $msg .= sprintf("### gene_id = %d [%d..%d] {%d} %s\n%s\n", map { defined $_ ? $_ : "" } @{$gl});
        }
        print $msg;
    } else {
        foreach my $gl (@geneLoad) {
            $psh->set_gene_loc_details( $sbjID, @{$gl} );
        }
    }
    return 1;
}

sub _set_database_meta {
    my ($globals) = @_;
    return unless ($globals);
    my $psh = &dbh();
    my %get = (
        size   => "SearchSize",
        normal => "GeneEntries",
        exon   => "ExonJunctions",
        masked => "MaskedCharacters",
        expect => "ExpectedHits",
        );
    my $bld;
    if ($bld = $args->val('build')) {
        $args->msg("[+]", "Manual definition of build as '$bld'");
    } elsif (my $first = `head -n1 "$inputFile"`) {
        # Try to extract the genome build
        # Expecting entries like:
        # 1.GRCh37:12227^12613.1.F
        # 1.GRCh37:11731-29370.FR
        if ($first =~ /^(\S+)/) {
            my $acc = $1;
            if ($acc =~ /^([^:]+):/) {
                # Remove trailing coordinates and strand
                $acc = $1;
            }
            if ($acc =~ /\.([^\.]+)$/) {
                # Cope with scaffold and contig designations
                # canis_familiaris.scaffold.JH374193.1.CanFam3:1507-2200.R
                # Will mangle and builds with internal '.' characters
                $bld = $1;
                
            }
        }
    }
    $psh->set_tag( $srcID, "GenomeBuild", $bld) if ($bld);
    while (my ($hk, $val) = each %{$globals}) {
        if (my $key = $get{$hk}) {
            $psh->set_tag( $srcID, $key, $val);
        }
    }
}

sub aliases_for_db {
    my $dbAliBase = "PgSeqHash_Aliases";
    if ($limit) {
        $dbAliBase .= "-LIMIT";
    } 
    my $dbAliFile = "$dbAliBase.tsv";
    if (!$limit && -s $dbAliFile) {
        $args->msg("[+]","Using existing database alias file", $dbAliFile);
        return $dbAliFile;
    }
    
    my $psh    = &dbh;
    my $dbh    = $psh->dbh();
    my $getAll = $dbh->prepare( -sql => "SELECT acc, taxa FROM gene",
                                -limit => $limit);
    # We will trust the taxa as given in DB
    my $txFile = sprintf("%s-Accessions.tsv", $dbAliBase);
    my (@ids, %isHuman, @noTax);
    if (!$limit && -s $dbAliFile) {
        open(TAXACC, "<$txFile") || $args->death
            ("Failed to read accessions", $!, $txFile);
        my $head = <TAXACC>;
        my $num = 0;
        while (<TAXACC>) {
            s/[\n\r]+$//;
            my ($id, $tax) = split(/\t/);
            push @ids, $id;
            if ($tax) {
                $isHuman{uc($id)} = 1 if ($tax eq 'Homo sapiens');
            } else {
                push @noTax, $id;
            }
            $num++;
        }
        close TAXACC;
        $args->msg("[<]", "$num IDs in existing DB accession file",
                   $txFile);
    } else {
        my $rows = $getAll->selectall_arrayref(  );
        open(TAXACC, ">$txFile") || $args->death
            ("Failed to write accessions", $!, $txFile);
        print TAXACC "Accession\tTaxa\n";
        foreach my $row (@{$rows}) {
            my ($id, $tax) = @{$row};
            push @ids, $id;
            if ($tax) {
                $isHuman{$id} = 1 if ($tax eq 'Homo sapiens');
            } else {
                push @noTax, $id;
            }
            print TAXACC join("\t", $id, $tax)."\n";
        }
        close TAXACC;
        $args->msg("[>]", "Processing ".($#{$rows}+1)." genes from database",
                   $txFile);
    }

    # Get all symbols
    my $symFile  = sprintf("%s-Symbols.tsv", $dbAliBase);
    unlink($symFile) if ($limit && -e $symFile);
    my $symRows  = &cached_data
        ( -idfile     => $txFile,
          -istable    => 1,
          -hasheader  => 1,
          -ns2        => 'SYM',
          -scramble   => 1,
          -nullscore  => -1,
          -sort       => 'termin',
          -cols       => 'termin,termout,score',
          -warn       => 0,
          -output     => $symFile );
    my %syms;
    foreach my $row (@{$symRows}) {
        my ($id, $sym, $sc) = @{$row};
        if ($id && $sym) {
            my $targ = $syms{uc($id)} ||= {};
            $targ->{$sym} = $sc if 
                (!defined $targ->{$sym} || $targ->{$sym} < $sc);
        }
    }
    $args->msg("[+]","Parsed ".($#{$symRows}+1)." gene symbols", $symFile);

    # Get descriptions
    my $descFile  = sprintf("%s-Description.tsv", $dbAliBase);
    unlink($descFile) if ($limit && -e $descFile);
    my $descRows  = &cached_data
        ( -idfile     => $txFile,
          -istable    => 1,
          -hasheader  => 1,
          -mode       => 'desc',
          -scramble   => 1,
          -nullscore  => -1,
          -nonull     => 1,
          -sort       => 'termin',
          -cols       => 'termin,desc',
          -output     => $descFile );
    my %descs;
    foreach my $row (@{$descRows}) {
        $descs{uc($row->[0])} = $row->[1];
    }
    $args->msg("[+]","Parsed ".($#{$descRows}+1)." descriptions", $descFile);

    # Get orthologues
    my $orthFile  = sprintf("%s-Orthologues.tsv", $dbAliBase);
    unlink($orthFile) if ($limit && -e $orthFile);
    my $orthRows  = &cached_data
        ( -idfile     => $txFile,
          -istable    => 1,
          -hasheader  => 1,
          -sort       => 'termin',
          -ns2        => 'orth',
          -scramble   => 1,
          -nullscore  => -1,
          -keepbest   => 1,
          -int        => 'homo sapiens',
          -intns      => 'TAX',
          -nonull     => 1,
          -aspercent  => 1,
          -minscore   => 5,
          -cols       => 'termin,termout,symout,score',
          -output     => $orthFile );
    my %orths;
    foreach my $row (@{$orthRows}) {
        my ($id, $hid, $sym, $sc) = @{$row};
        $orths{uc($id)} = [$hid, $sym, $sc] if ($id && $hid);
    }
    $args->msg("[+]","Parsed ".($#{$orthRows}+1)." orthologues", $orthFile);
    
    # Write the final file
    open(ALIFILE, ">$dbAliFile") || $args->death
        ("Failed to write alias file", $dbAliFile, $!);
    print ALIFILE join("\t", qw(Accession Symbol HumanAccession HumanSymbol HumanMatch Aliases Description))."\n";
    foreach my $id (@ids) {
        my $ucid = uc($id);
        my $targ = $syms{$ucid} ||= {};
        my @syms = sort { $targ->{$b} <=> $targ->{$a} ||
                              $a cmp $b } keys %{$targ};
        my $sym = shift @syms || "";
        my $ali = join(',', @syms) || "";
        my ($hid, $hsym, $hsc) = $isHuman{$ucid} ?
            ($id, $sym, 100) : @{$orths{$ucid} || ['','','']};
        my @row = ($id, $sym,  $hid, $hsym, $hsc, $ali, $descs{$ucid} || "");
        print ALIFILE join("\t", @row)."\n";
    }
    $args->msg("[>]","Alias file written", $dbAliFile);
    close ALIFILE;
    if ($#noTax != -1) {
        my $taxIds   = sprintf("%s-NoTax.tsv", $dbAliBase);
        open(NOTAX, ">$taxIds") || $args->death
            ("Failed to write taxa-less ID file", $taxIds, $!);
        map { print NOTAX "$_\n" } @noTax;
        close NOTAX;
        my $taxFile  = sprintf("%s-Species.tsv", $dbAliBase);
        unlink($taxFile) if ($limit && -e $taxFile);
        my $taxRows  = &cached_data
            ( -idfile     => $taxIds,
              -sort       => 'termin',
              -ns2        => 'tax',
              -scramble   => 1,
              -keepbest   => 1,
              -nonull     => 1,
              -minscore   => 1,
              -cols       => 'termin(Accession),termout(Species)',
              -output     => $taxFile );
        $args->msg("[+]",($#noTax+1)." IDs lacked species, ".($#{$taxRows}+1).
                   " were found", "You can load this file separately",$taxFile);
    }

    die;
    return $dbAliFile;
}

sub cached_data {
    # Needs GenAcc available:
    require BMS::MapTracker::GenAccService;
    my $isBeta = 1;
    my $age = $args->val('age') || '8 Sep 2016 4pm';
    my $gas = BMS::MapTracker::GenAccService->new
        ( -fork       => $args->val('fork') || 16,
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


sub load_gene_aliases {
    my $file = shift;
    return unless ($file);
    if ($file eq 'db') {
        # We want to get aliases for all accessions in the database
        $file = &aliases_for_db;
    }
    open(FILE, "<$file") || $args->death("Failed to read aliases", $file, $!);
    $args->msg("[<]","Loading gene metadata from file", $file);
    my $psh    = &dbh;
    my $dbh    = $psh->dbh();
    my $headRow = <FILE>;
    $headRow =~ s/[\n\r]+$//;
    my @head = split(/\t/, $headRow);
    my ($symCol, $descCol);
    my @loaders;
    for my $h (1..$#head) {
        my $cname = lc($head[$h]);
        my $sth;
        if ($cname =~ /human/) {
            if ($cname =~ /(sym)/) {
                $sth = $dbh->named_sth("Set Human Symbol");
            } elsif ($cname =~ /(gene|loc|acc)/) {
                $sth = $dbh->named_sth("Set Human Gene");
            } elsif ($cname =~ /(score|match)/) {
                $sth = $dbh->named_sth("Set Human Score");
            }
        } elsif ($cname =~ /(sym)/) {
            $sth = $dbh->named_sth("Set Gene Symbol");
        } elsif ($cname =~ /(tax|species)/) {
            $sth = $dbh->named_sth("Set Gene Taxa");
        } elsif ($cname =~ /(desc)/) {
            $sth = $dbh->named_sth("Set Gene Description");
        } elsif ($cname =~ /(alias)/) {
            my $aliSth = $dbh->named_sth("Set Gene Aliases");
            my $splitter = $args->val('split') || ',';
            $args->msg("[+]", "Aliases in Column #$h : $head[$h]",
                       "Will be split on '$splitter'");
            push @loaders, sub {
                my $row = shift;
                my $av  = $row->[$h];
                return unless ($av);
                my @alis = split(/\s*$splitter\s*/, $av);
                map { s/^\s+//; s/\s+$//; } @alis;
                $aliSth->execute( \@alis, $row->[0] );
            }
        }
        if ($sth) {
            $args->msg("[+]", "Capturing Column #$h : $head[$h]");
            push @loaders, sub {
                my $row = shift;
                my $val = $row->[$h];
                $val = undef if (defined $val && $val eq '');
                $sth->execute( $val, $row->[0] );
            };
        }
    }
    $args->death("[?]","Failed to recognize any columns",
                 "Please make sure the file has a header row",
                 "Recognized columns:","Symbol, Taxa, Description, Aliases")
        if ($#loaders == -1);
    my $num = 0;
    while (<FILE>) {
        s/[\n\r]+$//;
        my @row = split(/\t/);
        if (my $acc = $row[0]) {
            if (my $gid = $psh->gene_id( $acc )) {
                $row[0] = $gid;
                map { &{$_}( \@row ) } @loaders;
                $num++;
                last if ($limit && $num >= $limit);
                $args->msg("[$num]", $acc) unless ($num % 1000);
            }
        }
    }
    close FILE;
    $args->msg("Loaded aliases for $num genes","Exitting");
    exit;
}

sub _prefix_file {
    my ($pfx, $output) = @_;
    my $dir  = $output->{dir};
    my $file = "$pfx.tsv";
    push @{$output->{files}}, $file;
    my $path = "$dir/words/$file";
    my $fh;
    open($fh, ">$path") || $args->death("Failed to make word file", $path, $!);
    return $fh;
}


=head2 PG Space benchmarks

 relname  | tablespace | kilopages | kilotuples |  mega_xid  |  disk   
----------+------------+-----------+------------+------------+---------
 wordsrc  |   GRCh37   |   466.843 |      16716 |  26.167199 | 16 GB
 wordsrc  |  +MacFas5  |   914.332 |      32826 |   42.58231 | 31 GB
 
From the prior doomed-to-fail 2D-array implementation, and with 50x masking

 relname  | tablespace | kilopages | kilotuples |  mega_xid  |  disk   
----------+------------+-----------+------------+------------+---------
 ++ CanFam3:
 word     |            |    25.519 |      16046 | 102.476613 | 1837 MB
 wordsrc  |            |    87.604 |      10803 |  27.516984 | 3009 MB
 ++ BgiIR:
 wordsrc  |            |   172.179 |      21712 |  38.805488 | 5930 MB
 word     |            |    25.519 |      16046 | 113.765117 | 1858 MB
 ++ BgiCR:
 wordsrc  |            |   245.478 |      32122 |  49.253143 | 8595 MB
 word     |            |    25.519 |      16046 | 124.212772 | 1860 MB
 ++ GRCh37
 wordsrc  |            |   320.421 |      41749 |  59.004248 | 11 GB
 word     |            |    25.519 |      16046 | 133.963877 | 1866 MB
 ++ GRCm38 Rnor_5 MacFas5 BgiCE
 wordsrc  |            |   591.768 |      78276 |  96.229563 | 20 GB
 word     |            |     26.73 |      16769 | 100.000705 | 1873 MB

=head2 Location savings when excluding oligos over 50x expected hit count

             Expect #Repeat  %Excl (#Excluded
     GRCh37     50x    4033   7.19 (97308711)
     GRCm38     50x    3730   4.58 (57303353)
     CanFam3    50x    3108   5.29 (55095764)
     Rnor_5     50x    3177   3.48 (37607621)
     MacFas5    50x    3803   7.10 (90472770)
     BgiCE      50x    3447   5.71 (65983181)
     BgiCR      50x    2698   6.70 (60529103)
     BgiIR      50x    2819   7.32 (69136013)

=cut
