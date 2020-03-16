#!/stf/biobin/perl -w
my $isBeta;
BEGIN {
    # Allows usage of beta modules to be tested:
    my $prog = $0; my $dir = `pwd`;
    require lib;
    if ($prog =~ /working/ || $dir =~ /working/) {
	# warn "\n\n *** This is Beta Software ***\n\n";
	import lib '/stf/biocgi/tilfordc/perllib';
        $isBeta = 1;
    }
    import lib '/stf/biocgi/tilfordc/patch_lib';
    # import lib '/stf/biocgi/tilfordc/released/Bio/SeqIO';
}

=head1

=head2 Dump oligos to fasta file

  psql -P 'format=Unaligned' -P 't' -c "select fasta from v_fasta order by daterun desc" seqhash > foo.fa

=cut

use strict;
use Bio::SeqIO;
use Bio::PrimarySeq;
use JSON;

use BMS::ArgumentParser;
use BMS::CleanSeqString;
use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::ColorUtilities;
use BMS::ForkCritter;
use BMS::Utilities::Benchmark;
use BMS::Utilities::Serialize;
use BMS::ExcelHelper;
use BMS::FastaMeta;
use BMS::PgSeqHash;
use BMS::PgTaskManager;

my $args = BMS::ArgumentParser->new
    ( -nocgi     => $ENV{HTTP_HOST} ? 0 : 1,
      -testmode  => 1,
      -verbose   => 1,
      -strand    => 'R',
      -limit     => 0,
      -cache     => 500,
      -prefix    => 6,
      -basedir   => 'Ensembl',
      -progress  => 120,
      -mismatch  => 2,
      -noorphan  => 1,
      -fork      => 5,
      -exclreq    => &_temp_exclude(),
      -paramfile  => [ 'BMS::PgSeqHash' ],
      -paramalias => {
          strand     => [qw(str)],
          limit      => [qw(lim)],
          verbose    => [qw(vb)],
          project    => [qw(proj)],
          progress   => [qw(prog)],
          mismatch   => [qw(mis)],
          fork       => [qw(forknum numfork)],
          projbase   => [qw(projectbase basename)],
          benchmark  => [qw(bench dobench showbench)],
          start      => [qw(s)],
          end        => [qw(e)],
          minentropy => [qw(minent)],
          srchid     => [qw(srch_id)],
          htmltmp    => [qw(htmltemp tmphtml temphtml)],
          bgnd       => [qw(background)],
          genereq    => [qw(gene rnaacc rnaid rnaname)],
          exclreq    => [qw(exclude)],
          addcol     => [qw(addcols addcolumn addcolumns)],
          output     => [qw(out outfile)],
          pathtourl  => [qw(path2url)],
          mismatch   => [qw(mm maxmiss)],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
          xxxx => [qw()],
      });

if ($args->val('showinc')) {
    # Show the modules being used
    foreach my $mod (sort keys %INC) {
        print "$mod\t$INC{$mod}\n";
    }
    &do_exit();
}

my $outFH = *STDOUT;
my $formId = 'mainform';

my $user  = $ENV{'REMOTE_USER'} || $ENV{'LDAP_USER'} || $ENV{'USER'} || 
    $ENV{'LOGNAME'} || "USER";

my ($fc, %doneFmt, $preMrnaFiles, $doneSomething, %nck,
    %prefixCache, $nowMeta, %metaCache, %doneDB, %fixSeq,
    $forkPsh, $forkPgtm, $srcCache, $resultHandler, $htmlTmp, $cdbCnt,
    $finalURL, $jsonStructure, $extDivDone, $dhtmlajax,
    $filterNotes, $defaultCols, $foundMeta, $pgtm
    );

my $usedCols = {};

$args->xss_protect('all');
my @protectedArgs = 
    qw(cxurl toolurl errormail htmltmp favicon errorfile adminmail
       stylesheet javascript confstyles dbuser pathtourl
       imagedir urlmap pgport pghost tiddlywiki);
$args->default_only(@protectedArgs);
my $webCall = $args->val('webcall');
if ($webCall || $args->val('dhtmlajax')) {
    $dhtmlajax = {};
    $args->set_mime( -mail => $args->val('errmail'),
                     -mime => 'application/json' );
    $SIG{__DIE__}  = sub {
        my $msg = shift;
        push @{$dhtmlajax->{error}},"FATAL: $msg";
        &do_exit();
    };
    $SIG{__WARN__}  = sub {
        my $msg = shift;
        push @{$dhtmlajax->{error}}, $msg;
    };
}

my $norun       = $args->val('norun');
my $nocgi       = $webCall ? 1 : $dhtmlajax ? 0 : ($args->val('nocgi') || 0);
my $svcReq      = $args->val('service') || 0;
my $bgndReq     = $args->val('bgnd')    || 0;
my $clean       = $args->val('clean');
my $dumpSql     = $args->val('dumpsql') || 0;
my $tempDir     = $args->val('htmltmp') || "/tmp/PgSeqHash";
$tempDir        =~ s/\/$//;
$args->assure_dir($tempDir);
my $doHTML      = $nocgi ? 0 : $norun ? 1 : ($clean || $dhtmlajax) ? 0 : 1;
my $bgndCnt     = 0;
my $webSite     = "http://bioinformatics.bms.com";
my $genericName = 'UserSequence';
my $headerPrfx  = '###';
my $filtPrfx    = '@FILTER@ ';
my $tempSfx     = '.temp.tsv';
my $path2url    = $args->val('pathtourl');
my $psh         = BMS::PgSeqHash->new( );
my $usePsh      = $psh;
my $jobId       = $args->val('jobid');

if ($path2url) {
    my @reqs = ref($path2url) ? @{$path2url} : ($path2url);
    my @pairs;
    foreach my $req (@reqs) {
        if ($req =~ /^\s*(.+?)\s+=\s+(.+?)\s*$/) {
            push @pairs, [$1, $2];
        }
    }
    if ($#pairs != -1) {
        $args->url_callback( sub {
            my $path = shift;
            foreach my $pair (@pairs) {
                my ($sdir, $url) = @{$pair};
                if ($path =~ /^\Q$sdir\E/) {
                    my $rv = $path;
                    $rv =~ s/^\Q$sdir\E/$url/;
                    return $rv;
                }
            }
            return undef; });
    }
}

if ($jobId) {
    # This is a specialized task interacting with a PgTaskManager job entry
    $pgtm = BMS::PgTaskManager->new
        ( -jobid    => $jobId,
          -instance => $usePsh->instance(),
          -dbuser   => $usePsh->dbuser(),
          -dbpass   => $usePsh->dbpass(),);
    # Share the URL callback with PGTM - it will need it for doing
    # JSON calls against local files:
    $pgtm->url_callback( $args->url_callback() );
    &get_job_status() if ($args->val('getstatus'));
    my $meta = $pgtm->job_metadata();
    my $wdir = $meta ? $meta->{dir} : undef;
    unless ($wdir) {
        my $failFile = "$tempDir/genericFailures.txt";
        if (open FAIL, ">>$failFile") {
            print FAIL "Failed to identify job_id = $jobId\n";
            close FAIL;
        }
        &do_exit();
    }
    # Redirect STDOUT and STDERR to files
    my $outF = "$wdir/stdout.txt";
    my $errF = "$wdir/stderr.txt";
    open(STDOUT, ">>$outF");
    open(STDERR, ">>$errF");
    warn "## This file captures 'Standard Error'\n";
    warn "## Most entries will be innocuous messages rather than errors\n";
    warn "\n";
    warn "$$ was given job ID $jobId";
    system("ln -sf \"$errF\" /stf/biocgi/tilfordc/working/stderr.txt");
}

my @seqParams = qw(oligo oligos seq input file sequence id query);


# my $tl = '/tmp/SOS_Timeline.txt'; system("date >> $tl; echo '  $nocgi $bgndReq $tickPath' >> $tl;");

# die $args->branch(\%INC);

my $debugDump = "/tmp/ParamDump.param";
unlink($debugDump) if (-e $debugDump);

if ($dhtmlajax || $webCall) {
    
} elsif ($pgtm) {
    # Do not color output, it will be going to a text file
} elsif ($nocgi) {
    $args->shell_coloring();
} else {
    $args->set_mime( -mail => $args->val('errmail'),
                     -mime => $clean ? 'text' : 'html' );
}

# if (open(DMP, ">$debugDump")) { print DMP $args->to_text(); close DMP; chmod(0666, $debugDump); }

&START_HTML();

# &test_nck();

my $iter        = 0;
my $lastRpt     = 0;
my $fmtIter     = 0;
my $ln10        = log(10);
my $ln4         = log(4);
my $lodBand     = 2;
my $rcTok       = "<span class='rctok'>&#x21BB;</span>";
my $css         = BMS::CleanSeqString->new();
my $renameSth   = $usePsh->prepare("SELECT obj_id FROM tagval WHERE tag = ?");
my $metaPrfx    = $usePsh->user_col_prefix();

# $usePsh->dbh->make_all() if ($args->val('rebuild'));
# my $foo = get_pgtm(); $foo->dbh->make_all(); die;


my @clsPriority = qw(Exon ExonExon RNA FirstExon LastExon ExonInDel Intron);
my %clsPri      = map { $clsPriority[$_] => $_ + 1 } (0..$#clsPriority);

my %strNames = ( F  => "Forward strand",
                 R  => "Reverse strand",
                 FR => "both strands");

my %buildColors = ( GRCh37  => 'lime',
                    MacFas5 => 'green',
                    BgiIR   => 'green',
                    BgiCR   => 'green',
                    BgiCE   => 'green',
                    Rnor_5  => 'brown',
                    GRCm38  => 'brown',
                    CanFam3 => 'orange',
                    RSG     => 'lime');

my $stReq     = lc($args->val('strand')  || 0);
my $maxMiss   = $args->val('mismatch')   || 0;
my $limit     = $args->val('limit')      || 0;
my $vb        = $args->val('verbose')    || 0;
my $prog      = $args->val('progress')   || 300;
my $forkNum   = $args->val('fork');
my $redo      = $args->val('redoexcel');
my $keepScan  = $args->val('keepscan');
my $noRna     = $args->val(qw(norna));
my $noOrphan  = $args->val(qw(noorphan));
my $clobber   = $args->val('clobber');
my $projBase  = $args->val('projbase');
my $mode      = lc($args->val('mode') || "");
my $frmReq    = lc($args->val('format') || "");
my $fileReq   = $args->val('output') || "";
my $geneStr   = $stReq =~ /^(-1|r)/ ? -1  : $stReq =~ /^(\+?1|f)/ ? 1   : 0;
my $chkStrand = $stReq =~ /^(-1|r)/ ? 'R' : $stReq =~ /^(\+?1|f)/ ? 'F' : 'FR';

my $fileEnd   = '******';
my $bnch      = BMS::Utilities::Benchmark->new();
my $cu        = BMS::Utilities::ColorUtilities->new();
my $su        = BMS::Utilities::SequenceUtilities->new();
my $ser       = BMS::Utilities::Serialize->new();

my $head = ['OligoID', 'Symbol', 'HumanSymbol', 'Build', 'Mis', 'HitID', 'GeneID', 'RnaID', 'Start', 'End', 'Str', 'Exonic', 'GenomeFootprint', 'OligoSeq', 'HitSeq', 'OligoDesc', 'GeneDesc' ];
my $tsvRows = 0;
my $unBuild = "UnkBld";
my $unTaxa  = "Unknown species";
my $unGene  = "UnkGene";
my $fakeSym = 'zzzzzzz';
my $outFh   = *STDOUT;
my $proj    = $args->val('project') || ""; $proj =~ s/\/+$//;
my $projSht = $proj; $projSht =~ s/^.+\///;

my $outF    = $args->val(qw(output out)) || "";
$outF       =~ s/\.[\.]{2,6}$//;
my @seqs    = ();

my @dbs; #     = &get_dbs(); OLD CODE
my $search;#  = &get_query(); # OLD CODE

my $srcBldSth  = $usePsh->named_sth('Get source ID for build');
my $subNameSth = $usePsh->named_sth('Get sbj_id for name');
my ($rnaReq, $gidReq, $srcReq, $subReq, $keepReq) = &_process_gene_request
    ( 'genereq', "Require" );

my ($rnaExcl, $gidExcl, $srcExcl, $subExcl, $exclReq) = &_process_gene_request
    ( 'exclreq', "Exclude" );

my @filters = (-rna        => $rnaReq,
               -geneids    => $gidReq,
               -srcids     => $srcReq,
               -subjects   => $subReq,
               -norna      => $rnaExcl,
               -nogeneids  => $gidExcl,
               -nosrcids   => $srcExcl,
               -nosubjects => $subExcl,
               -strand     => $geneStr,
               -dumpsql    => $dumpSql,
               -mm         => $maxMiss);

#OLD:
my $qsLU = 
    $chkStrand eq 'F' ? { F => 'F', R => 'R', FR => 'FR' } :
    $chkStrand eq 'R' ? { F => 'R', R => 'F', FR => 'FR' } : 
{ F => 'FR', R => 'FR', FR => 'FR' };

my %okStr;
if ($chkStrand =~ /F/) {
    $okStr{1}  = { FR => 1, F => 1 };
    $okStr{-1} = { FR => 1, R => 1 };
}
if ($chkStrand =~ /R/) {
    $okStr{1}  = { FR => 1, R => 1 };
    $okStr{-1} = { FR => 1, F => 1 };
}


# die $args->branch(\@dbs);
my $misPick = sub {
    my $v = shift;
    if (!$v) {
        return 'zeromiss';
    } elsif ($v == 1) {
        return 'onemiss';
    } elsif ($v == 2) {
        return 'twomiss';
    } elsif ($v > 2) {
        return 'threemiss';
    }
};

my $strCB = sub {
    my $v = shift; 
    return ! $v ? 'strUnk' : $v < 0 ? 'rev' : 'fwd';
};

my $formatters = {
    Symbol   => sub {
        my ($v, $row) = @_;
        return $row->{isPseudo} ? 'pseudosym' : 'sym';
    },
    Build   => sub {
        my ($v, $row) = @_;
        return $v;
    },

    'Ent1'   => \&_entropy_fmt,
    'Ent2'   => \&_entropy_fmt,
    'Ent3'   => \&_entropy_fmt,

    Taxa     => sub { return 'species'; },
    Species  => sub { return 'species'; },
    OligoSeq => sub { return 'mono'; },
    OligoID  => sub { return 'oligoid'; },
    HitSeq   => sub { return 'mono'; },
    Str      => $strCB,
};


my $htmlFormatters = {
    Sequence => sub {
        return sprintf("<td class='align'>%s</td>", shift);
    },
    MM => \&_mismatch_html_cell,
    Type => \&_self_classed_html_cell,
};


my $queries   = &get_oligos();

my ($fileSfx, $format);
my $isAuto = $args->val('isauto');
if (!$frmReq || $frmReq =~ /auto/) {
    # This is a flag to help configure the web gui after AJAX calls
    $frmReq = "";
    $isAuto = 1;
}
if ($frmReq =~ /(excel|xls)/ || $fileReq =~ /\.(xlsx?)$/i) {
    $format  = 'Excel';
    $fileSfx = 'xlsx';
} elsif ($frmReq =~ /(text|txt)/ || $fileReq =~ /\.(te?xt)$/i) {
    $format = 'Text';
    $fileSfx = 'txt';
} elsif ($frmReq =~ /(json)/ || $fileReq =~ /\.(json)$/i) {
    $format = 'JSON';
    $fileSfx = 'json';
} elsif ($frmReq =~ /(htm)/ || $fileReq =~ /\.(html?)$/i) {
    $format = 'HTML';
    $fileSfx = 'html';
} elsif ($frmReq =~ /(tsv)/ || $fileReq =~ /\.(tsv)$/i || $nocgi) {
    $format = 'TSV';
    $fileSfx = 'tsv';
} else {
    $format = 'JSON';
    $fileSfx = 'json';
}
$frmReq = $isAuto ? '' : $frmReq ? $format : "";

if ($mode =~ /det/) {
    $mode = 'Details';
    $defaultCols = $format eq 'Excel' ?
        [qw(QueryID Build Symbol HumanSymbol MM AM Strand Type GenomeFootprint RnaFootprint NumRNA ExonNum RnaStart RnaEnd QuerySeq SubjectSeq GeneAccession HumanAccession HumanScore GeneDescription)] : 
        [qw(OligoID Symbol HumanSymbol Build MM AM Strand HitType GenomeFootprint NumRNA QuerySeq SubjectSeq GeneAccession HumanAccession HumanScore GeneDescription GeneID)];
} elsif ($mode =~ /over/) {
    $mode = 'Overview';
    $defaultCols = [qw(OligoID GenomeBuild Ex0 Ex1 Ex2 Int0 Int1 Int2 
GeneIDsEx0 GeneIDsEx1 GeneIDsEx2 GeneIDsInt0 GeneIDsInt1 GeneIDsInt2 Species DatabasePath)];
} elsif ($mode =~ /gen/) {
    $mode   = 'Genomic View';
    $defaultCols = [];
} else {
    $mode = 'Summary';
    $defaultCols = [qw(OligoID OligoSeq Build MM0 MM1 MM2 Species DatabasePath)];
}

my $addCols = $args->val('addcol');
if ($addCols) {
    # These are extra columns to add to the default
    my @list = ref($addCols) ? @{$addCols} :
        split(/\s*[\n\r\t\,]+\s*/, $addCols);
    $addCols = \@list;
}
if (my $userCols = $args->val('cols')) {
    # The user is specifying a particular column order
    my @list = ref($userCols) ? @{$userCols} :
        split(/\s*[\n\r\t\,]+\s*/, $userCols);
    $defaultCols = \@list;
}

my $oligoMeth = 
    $mode eq 'Summary'  ? \&oligo_summary : \&run_oligo;

&service($svcReq) if ($svcReq);

&test();
&make_oligos();

my @startmsg = ("Checking ".$strNames{$chkStrand},
                "Allowed Mismatches: $maxMiss");
push @startmsg, "Limiting to $limit hits per DB" if ($limit);
&msg("[+]",@startmsg) unless ($doHTML || $dhtmlajax);


my $excelMeta = {


    
    'Source'     => [15, "For specificity scans, the transcript ID the oligos were derived from"],
    'Ent1'       => [8, "Single nucleotide entropy, a measure of how diverse the oligo is at 1bp resolution. Low entropy means low complexity in the oligo."],
    'Ent2'       => [8, "Dinucleotide entropy, a measure of how diverse the oligo is at 2bp resolution"],
    'Ent3'       => [8, "Trinucleotide entropy, a measure of how diverse the oligo is at 3bp resolution"],

    '<Mis'      => [5, "The minimum number of mismatches observed"],
    'Mis>'      => [5, "The maximum number of mismatches observed"],

    
    'Hits'      => [8, "The total number of hits observed within a database"],
    'Score'     => [8, "A crude measure of how 'hitable' a database was - it is a function of both the number of hits, and how few mismatches those hits had. This is just used to order the databases in reports so 'more interesting' ones are shown first"],
    
    'Path'      => [40, "The full path of the database on the file system"],
    'TotSymbols' => [5, "The total number of symbols observed for an oligo within a particular species. Should be equal to or smaller than TotLoci. Loci without symbols are generally the 'speculative' genes", 'cen'],
    'TotLoci' => [5, "The total number of loci observed for an oligo within a particular species. These represent distinct genes. Note that one pre-mRNA region may contain two or more genes", 'cen'],
    'AllSymbol' => [10, "On the symbol summary page, this column simply shows the specific capitalization of observed symbols (since case standards vary across species)", 'cen'],
};

if ($args->val(qw(list))) {
    my $files = &_premrna_files();
    my @all = sort keys %{$files};
    my @msg = map { sprintf("%10s : %s", $_, $files->{$_}) } @all;
    &msg("Detected pre-mRNA files:", @msg,"To use all of them:",
               "-db ".join(',', @all));
} elsif ( $pgtm ) {
    &run_job( );
} else {
    &run();
}

my $fsmsg = "";
my @fsb = sort keys %fixSeq;
foreach my $fsb (sort keys %fixSeq) {
    $fsmsg .= "Build $fsb\n";
    foreach my $res ( sort { $fixSeq{$fsb}{$b} <=>
                                 $fixSeq{$fsb}{$a} } keys %{$fixSeq{$fsb}}) {
        $fsmsg .= sprintf("   %03d %s\n", $fixSeq{$fsb}{$res}, $res);
    }
}

warn "Corrupted Hit sequences fixed:\n$fsmsg" if ($fsmsg);
&HTML_INTERFACE();

if (my $txt = &show_bench()) {
    warn $txt;
}

&do_exit();

sub prepare_for_json {
    $args->set_mime( -mail => $args->val('errmail'),
                     -mime => 'json' );
    $jsonStructure = {};
}

sub get_job_status {
    return "" unless ($pgtm);
    &prepare_for_json();
    my $statJson = encode_json( $pgtm->job_status_json() );
    print $statJson."\n";
    my $tmp = "/stf/biocgi/tilfordc/working/stat.json"; if (open(STAT,">$tmp")) { print STAT $statJson; close STAT; }
    exit;
}

sub crude_debug {
    return;
    my $txt = shift;
    return unless ($txt);
    my $dbFile = "/tmp/SOS_Timeline.txt";
    unless ($cdbCnt++) {
        unless (-e $dbFile) {
            system("touch $dbFile; chmod 0666 $dbFile");
        }
        my $dt = `date`; $dt =~ s/[\n\r]+$//;
        system(sprintf("echo -e '\n\n##< %d >## %s' >> %s", $$, $dt, $dbFile));
    }
    system(sprintf("echo -e '[%2d] %s' >> %s", $cdbCnt, $txt, $dbFile));
}



sub make_oligos {
    my $req = $args->val('makeoligo');
    # for my $i (1..10) { my $tst = "ACTG" x $i; warn sprintf("%6s %s\n", &_seq_entropy( $tst), $tst) }; die;

    return unless ($req);
    my ($sr, $er) = ($args->val('start'), $args->val('end'));
    my $len = $args->val('length') || 15;
    my $ofh = *STDOUT;
    my ($file, $tsv, $entf);
    my @files;
    if ($outF) {
        $file  = $outF;
        $file .= ".fa" unless ($file =~ /\.fa$/);
        $args->assure_dir($file, 'isfile');
        $ofh = undef;
        open($ofh, ">$file") || $args->death
            ("Failed to write fasta file", $file, $!);
        my $tfile = $file . ".tsv";
        open($tsv, ">$tfile") || $args->death
            ("Failed to write oligo TSV file", $tfile, $!);
        print $tsv join("\t", "Sequence", "1-mer Entropy", "2-mer Entropy", "3-mer Entropy", "Position", "Source", "Name")."\n";
        push @files, ($file, $tfile);
    }

    my $minEnt = $args->val('minentropy');
    &msg("[+]","Rejecting oligos with entropy < $minEnt");
    my $css = BMS::CleanSeqString->new();
    foreach my $bs ( $css->clean_bioseq($req) ) {
        my $s = $sr || 1;
        my $e = $er || $bs->length() - $len + 1;
        my $id = $bs->display_id();
        my $seq = $bs->seq();
        for my $p ($s..$e) {
            my $sseq  = $su->revcom( substr($seq, $p - 1, $len) );
            my $name  = sprintf("%s:%d..%d", $id, $p, $p + $len - 1);
            my @ents  = map { &_seq_entropy( $sseq, $_) } (1, 2, 3);
            my $block = sprintf(">%s /pos=%d /len=%d /entropy=%s,%s,%s\n%s\n",
                                $name, $p, $len, @ents, $sseq);
            my ($worst) = sort {$a <=> $b } @ents;
            if ($minEnt && $worst < $minEnt) {
                if ($file && !$entf) {
                    my $ef = "$file.rejected";
                    push @files, $ef;
                    open($entf, ">$ef") || $args->death
                        ("Failed to write oligo reject file", $ef, $!);
                }
                print $entf $block if ($entf);
                $name = "FILTERED [$worst < $minEnt] $name";
            } else {
                print $ofh $block;
            }
            print $tsv join("\t", $sseq, @ents, $p, $id, $name)."\n";
        }
    }
    unless ($#files == -1) {
        close $ofh;
        close $tsv;
        close $entf if ($entf);
        &msg("[>]","Oligos generated", @files);
    }
    &do_exit();
}

sub _seq_entropy {
    my $seq  = uc(shift || "");
    my $sz   = shift || 1;
    $seq     =~ s/[^ACTG]+//g;
    my $slen = CORE::length($seq);
    my $denom = $slen - $sz + 1;
    my %u;
    for (my $i = 0; $i <= ($slen - $sz); $i++) {
        my $bit = substr($seq, $i, $sz);
        $u{$bit}++;
    }
    my $entropy = 0;
    foreach my $count (values %u) {
        my $frac  = $count / $denom;
        $entropy -= $frac * log($frac);
    }
    # What is the maximum entropy for this size?
    my $maxVars = 4 ** $sz;
    $maxVars = $denom if ($maxVars > $denom);
    my $maxEnt = $maxVars ? 0 - log(1/$maxVars) : 0;
    $entropy /= $maxEnt if ($maxEnt);
    return int(0.5 + 1000 * $entropy) / 1000,
}

sub testoligo {
    my $seqs = shift;


    die "Need to normalize with newer code";
    
    $doneSomething = 1;
    my @tsv;
    &init();
    $resultHandler = ($nocgi || $clean) ? die "no" : die "not this either";
    foreach my $bs (@{$seqs}) {
        &msg("[+]", sprintf("%s : %s", $bs->display_id(), $bs->seq()));
        my ($srchid, $isNew) = &search_id_for_bioseq( $bs );
        my $flat = &{$oligoMeth}( [$srchid, $isNew] );
    }
    print join("", map { join("\t", @{$_})."\n" } @tsv)."\n";
    if (my $sb = &show_bench()) {
        warn $sb;
    }
    &do_exit();
}

sub testlocation {
    my $req = shift;
    &init();
    foreach my $loc (split(/[\n\r\s]+/, $req || "")) {
        next unless ($loc);
        my ($sid, $s, $e, $str) = split(',', $loc);
        $str ||= 1;
        my $sname;
        if ($sid =~ /^\d+$/) {
            $sname = $usePsh->subject_name($sid);
        } else {
            $sname = $sid;
            $sid = $usePsh->sbj_id( $sname );
        }
        my $genes = $usePsh->subject_position_to_genes( $sid, $s, $e, $str );
        warn $args->branch($genes);
    }
    &do_exit();
}

sub print_dhtmlajax {
    return 0 unless ($dhtmlajax);
    if (exists $dhtmlajax->{files} && $dhtmlajax->{files}) {
        my $urls = $dhtmlajax->{urls} = {};
        while (my ($key, $path) = each %{$dhtmlajax->{files}}) {
            if (my $url = $args->path2url( $path )) {
                $urls->{$key} = $url;
            }
        }
    }
    print encode_json( $dhtmlajax );
    return 1;
}

sub do_exit {
    &note_filters() unless (defined $dhtmlajax && !$dhtmlajax);
    if (&print_dhtmlajax()) {

    } elsif ($doHTML) {
        my $fh = *STDOUT;
        print $fh &contact_html();
        print "<script>document.location='$finalURL'</script>" if ($finalURL);
        print "</body></html>\n";
    }
    exit;
}

sub test {
    if (my $seqs = &bioseqs_for_param('testoligo')) {
        &testoligo($seqs);
        &do_exit();
    } elsif (my $loc = $args->val('location')) {
        &testlocation($loc);
        &do_exit();
    }
    
    my $id = $args->val('test');
    return unless ($id);
    if ($#dbs == -1) {
        if (my $db = &_db_for_ID( $id )) {
            push @dbs, $db;
        } else {
            &msg("[!]","To run a test, please provide at least one -db");
            &do_exit();
        }
    }
    foreach my $db (@dbs) {
        my $file = $db->{file};
        my $faMeta = BMS::FastaMeta->new(- fasta => $file );
        unless ($faMeta) {
            &msg("[!!]","Not a metadata-supported file", $file);
            next;
        }
        if (my $obj = $faMeta->fetch($id)) {
            &msg("[+]","Metadata interface",
                       $faMeta->fasta, $faMeta->dbfile);
            my $data = $obj->data();
            print $args->branch($data);
            print join("\n", "", $obj->json(), "", "");
            print $obj->to_text();
            my ($s, $e) = ($args->val('start'), $args->val('end'));
            $s ||= 1;
            $e ||= $s + 99;
            &msg("[-]", "Sequence recovery $s..$e",
                       $obj->subseq($s, $e) || "NONE FOUND");
        } else {
            &msg("[?]","Unknown request '$id'");
        }
    }
    &do_exit();
}

sub _premrna_files {
    my $req = shift;
    unless ($preMrnaFiles) {
        $preMrnaFiles = {};
        foreach my $path (split(/[\n\r]/, `ls -1 /gcgblast/*_pre-mRNA_*.fa`)) {
            if ($path =~ /pre-mRNA_(\S+).fa$/) {
                $preMrnaFiles->{uc($1)} = $path;
            }
        }
    }
    return $preMrnaFiles->{uc($req)} if ($req);
    return wantarray ? sort keys %{$preMrnaFiles} : $preMrnaFiles;
}

sub run_job {
    return unless ($pgtm);
    my $jid = $pgtm->job_id();
    &msg("[#]","Child manager process $$ is taking on job_id = $jid" );
    my $meta = $pgtm->job_metadata();
    my $wdir = $meta->{work};
    my $outFile;
    if ($format eq 'JSON') {
        # These files will be made on-the-fly
        $resultHandler = \&task_to_json;
    } else{ 
        $resultHandler = \&store_data_in_file;
        $outFile       = "$wdir/Output.tsv";
    }
    
    $fc = BMS::ForkCritter->new
        ( -inputtype      => 'pgtm',
          -method         => \&pgtm_method,
          -progress       => 0,
          -init_meth      => \&init_pgtm_fork,
          -finish_meth    => \&finish_pgtm,
          -verbose        => $vb );
    $meta->{time}{start} = $fc->nice_date();
    $fc->output_file( 'OUT', ">$outFile" ) if ($outFile);
    my $failed = $fc->execute( $forkNum );
    # Get a fresh copy of the metadata, in case it was changed
    $meta = $pgtm->job_metadata();
    if ($failed) {
        push @{$meta->{error}}, "$failed jobs failed to properly run";
    }
    if (my $remain = $pgtm->total_remaining_tasks()) {
        push @{$meta->{error}}, "$remain tasks were never attempted";
    }
    if (my $orph = $pgtm->orphaned_tasks()) {
        push @{$meta->{error}}, "$orph tasks were reserved but never finished";
    }
   
    # NEED TO MANAGE TEMP FILE AND USER FILE HERE
    if ($format eq 'JSON') {
        # Nothing needs to be done here.
        # The results are JSON structures stored on the server
        # The PGTM tasks will be updated with the URL to each
        # That URL will be picked up by the AJAX monitoring Javascript
        # The results will be pulled into the framework via settimeout()
    } else {
        # The children need to be assembled into a final file
        my $out = $fileReq || $meta->{out};
        &post_process( $out, $outFile );
        delete $meta->{wait};
        $meta->{path} = $out;
    }
    $meta->{time}{end} = $fc->nice_date();
    $pgtm->change_job_metadata( $meta );
    exit;
}

sub get_pgtm {
    my $obj = $usePsh || $psh;
    return BMS::PgTaskManager->new
        ( -instance => $obj->instance(),
          -dbuser   => $obj->dbuser(),
          -dbpass   => $obj->dbpass(),);
}

sub init_pgtm_fork {
    # Give each forked child its own PgTaskManager object:
    &init();
    $forkPgtm = BMS::PgTaskManager->new
        ( -jobid    => $jobId,
          -instance => $usePsh->instance(),
          -dbuser   => $usePsh->dbuser(),
          -dbpass   => $usePsh->dbpass(),);
    $usePsh->pgtm( $forkPgtm );
    $fc->input_args( [-pgtm => $forkPgtm] );
    my $remain = $forkPgtm->total_remaining_tasks();
    if ($remain) {
        &msg("[#]", "Child $$ is instantiated, $remain task".
             ($remain == 1 ? '' : 's')." available");
    } else {
        &msg("[-]", "Child $$ not needed, all tasks taken by other children");
    }
}

sub pgtm_method {
    # This is the child-fork method to perform a task
    my $taskid = shift;
    my $compCode = -2;
    if ($taskid) {
        &msg("[+]", "Child $$ reserves ".$forkPgtm->task_token($taskid));
        my $meta = $forkPgtm->task_metadata( $taskid );
        if ($meta && $meta->{seq}) {
            my ($srchid, $isNew) = &search_id_for_bioseq( $meta->{seq});
            &{$oligoMeth}([$srchid, $isNew]);
            $compCode = -1;
        } else {
            $compCode = -3;
        }
    } else {
        &msg("[?]", "Child $$ was supposed to have a task, but it is not there!");
    }
    $forkPgtm->complete_task( $taskid, $compCode );
}

sub task_to_json {
    my ($data) = @_;
    my $obj  = $forkPgtm || $pgtm;
    my $jm   = $obj->job_metadata();
    my $wdir = $jm->{work};
    # warn $args->branch($pgtm->{TASK_IDS});
    my $tid  = $obj->task_id();
    my $file = "$wdir/Task_$tid.json";
    my $meta = $obj->task_metadata();
    &make_json($data, $file);
    $meta->{json} = $file;
    $obj->change_task_metadata( $meta, $tid );
}

sub finish_pgtm {
    my $obj  = $forkPgtm || $pgtm;
    return unless ($obj);
    if (my $remain = $obj->total_remaining_tasks()) {
        &msg("[!!]","Child $$ is exiting but reports $remain remaining tasks!");
    }
}

sub setup_background {
    my @uniqSeq = sort keys %{$queries};
    my $uNum = $#uniqSeq + 1;
    unless ($uNum) {
        if ($dhtmlajax) {
            $dhtmlajax->{alert} = "You need to specify at least one sequence";
            &do_exit();
        }
        return;
    }
    $pgtm = BMS::PgTaskManager->new
        ( -instance => $usePsh->instance(),
          -dbuser   => $usePsh->dbuser(),
          -dbpass   => $usePsh->dbpass(),);
    # $pgtm->make_tables(); die "Made schema";
    my $wdir      = &_work_dir();
    
    my $seqFile   = "$wdir/queries.fa";
    my $paramFile = "$wdir/settings.param";
    my $outDir    = "$wdir/output";
    my $files = {
        STDERR     => "$wdir/stderr.txt",
        STDOUT     => "$wdir/stdout.txt",
        Parameters => $paramFile,
        WorkDir    => "$wdir/",
        OligoSeqs  => $seqFile,
    };
    $dhtmlajax->{files} = $files if ($dhtmlajax);
    $args->set_param('format', $format);
    #if (($dhtmlajax || $doHTML) && $format ne 'Excel') {
    #    # We will want output written in JSON for recovery via AJAX
    #    $args->set_param('format', 'JSON');
    #    $format = 'JSON';
    #}
    my $outFile   = ($format eq 'JSON') ? undef : &pick_outfile();
    if (my $err = &write_sequences( $seqFile )) {
        $args->err("Failed to set up sequence file for search", $seqFile, $err);
        return;
    }
    # warn "$outDir = ".$args->path2url( $outDir );
    &msg("[>]", "STDERR: $wdir/stderr.txt", "STDOUT: $wdir/stdout.txt");

    $args->assure_dir( $outDir );
    my $meta = {
        dir   => $wdir,
        seq   => $seqFile,
        work  => $outDir,
        param => $paramFile,
    };
    $meta->{wait} = "Assemble $format output" unless ($format eq 'JSON');
    $meta->{out}  = $outFile if ($outFile);
    my $jid  = $pgtm->new_job
        ( -token => "$mode x $uNum",
          -meta  => $meta, );
    # We have normalized all sequence params to one file:
    $args->set_param('query', $seqFile);
    # Make sure we do not have duplicate entries from alias parameters:
    my %ignQry = map { $_ => 1 } @seqParams;
    delete $ignQry{'query'};
    
    # -jobid is what will trigger the background execution in the child
    $args->set_param('jobid', $jid );
    $args->set_param('isAuto', 1) if ($isAuto);

    open(PARAM, ">$paramFile") || $args->death
        ("Failed to serialize search parameters", $paramFile, $!);
    my @paramArgs =
        ( -comments => ["Parameter file for simpleOligoSearch.pl",
                        "Accessing job_id = $jid in pgtm_job"],
          -blockquote => [qw(EXCLREQ INPUT)],
          -ignore => [@protectedArgs, keys %ignQry, 
                      qw(background junk paramfile norun dhtmlajax)], );
    print PARAM $args->to_text( @paramArgs );
    close PARAM;
    $dhtmlajax->{params} = $args->to_hash(@paramArgs) if ($dhtmlajax);

    foreach my $seq (@uniqSeq) {
        my $sidIsNew = &prepare_seq_for_searching( $seq );
        my ($srchid, $isNew) = @{$sidIsNew};
        my $label;
        my @qry = @{$queries->{$seq} || []};
        foreach my $bs (@qry) {
            my $name = $bs->display_id();
            unless ($name =~ /^User/) {
                $label = $name;
                last;
            }
        }
        if ($label) {
            $label .= " (";
            # Note if two or more requests have the same sequence
            $label .= "+$#qry, " if ($#qry > 0);
            $label .= "$seq)";
        } else {
            $label = $seq;
        }
        $pgtm->queue_task
            ( -weight => $isNew ? 5 : 50,
              -token => sprintf("%s %s", $isNew ? "Run" : "Recover", $label),
              -meta => {
                  seq => $seq,
              });
    }
    $pgtm->change_job_task_count( $uNum );
    my $kid;
    $ENV{FORCE_NOCGI} = "Force NOCGI";
    my @cmd = ($0, '-nocgi', '-valuefile', $paramFile);
    if ($kid = fork) {
        # Parent does not need to do anything
        #&msg("Parent $$ passes job_id = $jid to $kid", $wdir)
        #    if ($nocgi);
        # The child will take care of the ajax
        $dhtmlajax = 0;
    } elsif (defined $kid) {
        #my $gkid = 0;
        #delete $ENV{PERL5LIB};
        if ($nocgi) {
            # From command line, wait for manager to finish
            my $rv = system @cmd;
            die "Job manager fails with code $rv" if ($rv);
        } else {
            # When run in a web context, we need to close STD handles,
            # launch manager, then exit.
            &END_JSON( $jid, $paramFile );
            # &END_HTML();
            close STDERR;
            close STDOUT;
            close STDIN;
            exec @cmd;
        }
        exit;
    } else {
        $args->death("Failed to fork off background child process!");
    }
    if ($nocgi) {
        # Wait for children to finish, then show output
        &msg("[-]","Parent [$$] is waiting for Job manager [$kid] ...");
        waitpid($kid, 0);
        my $err = 0;
        my $exit_value = $? >> 8;
        &msg("[!!]","CAUTION - Child process exits with $exit_value")
            if ($exit_value);
        &msg("[>]","Job manager [$kid] finished",
                   $outFile ? $outFile : $wdir);
    }
    &do_exit();
}

sub msg {
    my @msg = @_;
    if ($dhtmlajax || $webCall) {
        shift @msg if ($msg[0] =~ /^\[/);
        push @{$dhtmlajax->{msg}}, \@msg;
    } else {
        $args->msg(@msg);
    }
}

sub END_JSON {
    my ($jid, $paramFile) = @_;
    if ($dhtmlajax) {
        $dhtmlajax->{jobid} = $jid;
        $dhtmlajax->{alert} = "Job $jid: Started";
        &print_dhtmlajax();
        
    } else {
        print "<div id='exterr'></div>\n";
        if ($jid) {
            print "<script>monitor_pgtm( $jid )</script>\n";
        } else {
            &msg("[??]","CODE ERROR",
                 "No job_id was set, can not process request!");
        }
        print "<a class='setpop' href='simpleOligoSearch.pl?jobid=0&norun=1&paramfile=$paramFile'>Change search settings</a> - <a class='setpop' href='simpleOligoSearch.pl'>New Search</a><br />\n" if ($paramFile);
    }
}

sub write_sequences {
    my $file = shift;
    my @bss = sort { uc($a->display_id()) cmp
                         uc($b->display_id()) } 
    map { @{$_} } values %{$queries};
    return "No sequences" if ($#bss == -1);
    my $sio = Bio::SeqIO->new( -file => ">$file", -format => 'fasta' );
    foreach my $bs (@bss) {
        $sio->write_seq( $bs );
    }
    return 0;
}

sub run {
    return if ($args->val('norun'));
    my @uniqSeq = sort keys %{$queries};
    if ($mode eq 'Genomic View') {
        return &run_gene();
    } elsif ($#uniqSeq == -1) {
        # No oligos defined
        if ($rnaReq || $gidReq || $subReq) {
            # But there is one or more genes requested...
            return &run_gene();
        } elsif ($dhtmlajax) {
            $dhtmlajax->{alert} = "You need to specify at least one sequence";
            &do_exit();
        }
        &msg("[?]", "Can not run oligo search",
                   "Provide one or more oligo sequences using -oligo")
            unless ($doHTML);
        return;
    }
    my @srchids;
    foreach my $seq (@uniqSeq) {
        if (my $sidIsNew = &prepare_seq_for_searching( $seq )) {
            push @srchids, $sidIsNew;
        }
    }
    my $snums = $#srchids + 1;
    unless ($snums) {
        &msg("[!!]","No valid sequences to search");
        return;
    }
 
    unless ($webCall) {
        return &setup_background() if ($dhtmlajax || $bgndReq || $doHTML);
    }

    my $outFile;
    if ($format eq 'JSON') {
        if ($#srchids > 0) {
            &msg("[!]", "Can not generate JSON for multiple search_ids");
            return;
        } else {
            $resultHandler = \&make_json;
        }
    } elsif ($doHTML) {
        $resultHandler = undef;
        return &html5_run_oligo( \@srchids );
    } else {
        $outFile = &pick_outfile();
    }
    
    &msg("[<]","Processing $snums distinct sequences");
    $fc = BMS::ForkCritter->new
        ( -inputtype      => 'array',
          -method         => $oligoMeth,
          -progress       => $prog,
          -init_meth      => \&init,
          -finish_meth    => \&final,
          -verbose        => $webCall ? 0 : $vb );
    $fc->reset();
    $fc->input(\@srchids);
    if ($outFile) {
        $fc->output_file( 'OUT', ">$outFile$tempSfx" );
    }
    if (my $failed = $fc->execute( $forkNum )) {
        $args->death("$failed jobs failed to properly convert");
    }
    my $file = &post_process( $outFile );
    if (!$file) {
        &msg("[!]","Failed to recover file", "Should have been at: $outFile")
            if ($outFile);
    } elsif ($format eq 'Excel') {
        my @msg = ("Your excel file is ready", $file);
        if ($doHTML) {
            if (my $url = $finalURL = $args->path2url( $file )) {
                print "<script>document.location = '$url'</script>\n";
                push @msg, "<a href='$url'>Link</a>";
            }
        }
        &msg("[OUT]", @msg );
    } else {
        &show_stdout( $file );
    }
}

sub prepare_seq_for_searching {
    my $seq = shift;
    return undef unless ($seq);
    my $bss   = $queries->{uc($seq)};
    unless ($bss) {
        &msg("[??]","Request for unexpected sequence '$seq'");
        return undef;
    }
    my $bs    = $bss->[0];
    my ($srchid, $isNew) = &search_id_for_bioseq( $bs );
    if ($srchid) {
        my ($seq, $name, $descr) = $usePsh->search_details( $srchid );
        foreach my $obs (@{$bss}) {
            my ($oname, $oseq) = ($obs->display_id(), $obs->seq());
            unless ($oname eq $genericName) {
                $usePsh->set_tag( $srchid, "Other ID", $oname)
                    if (uc($oname) ne uc($name));
                $usePsh->set_tag( $srchid, "Other Seq", $oseq)
                    if ($oseq ne $seq);
                $usePsh->set_tag( $srchid, "Case $oseq", $oname);
            }
        }
        return [$srchid, $isNew];
    } else {
        &msg("[!!]", "Failed to set search for oligo", $seq);
    }
    return undef;
}


sub pick_outfile {
    my $outFile;
    my $force = shift;
    if ($format eq 'JSON' && !( $force || $args->val('forcefile'))) {
        # JSON is always to STDOUT unless otherwise requested
        return "";
    }
    if ($format eq 'Excel') {
        $resultHandler = \&store_data_in_file;
        $oligoMeth     = \&run_oligo;
        $outFile       = &_output_file();
    } else {
        # Command line
        $outFile       = &_output_file();
        $resultHandler = \&store_data_in_file;
    }
    return $outFile;
}

sub run_gene {
    unless ($rnaReq || $gidReq || $subReq) {
        &msg('[?]',"At least one gene, RNA or genomic subject must be provided to perform a gene search");
        return;
    }
    my $hits = $usePsh->recover_search_hits( @filters );
    &filter_hits($hits);
    my $file = &pick_outfile();
    return &generate_output( $hits, $file );
}

sub _output_file {
    my $file = $fileReq;
    if (!$file || lc($file) eq 'stdout') {
        my $tdir = &_work_dir();
        my $base = sprintf("%s/SeqHashOutput", $tdir);
        my $iter = 0;
        while (1) {
            $file = $base.($iter ? ".$iter" : "") . ".$fileSfx";
            last unless (-e $file);
            $iter++;
        }
    } elsif (-d $file) {
        if ($fileSfx) {
            $file = "$file/Output.$fileSfx";
        } else {
            &msg("[?]", "Output specified as directory, but a file suffix was not determined:");
            $file = "$file/Output";
        }
    }
    return $file;
}

sub post_process {
    my ($file, $tmp) = @_;
    if (!$file) {
        if ($tmp) {
            $file = &_output_file();
        }
        return unless ($file);
    }
    $tmp ||= "$file$tempSfx";
    unless (-e $tmp) {
        &msg("[?]","Failed to find output", $tmp);
        return;
    }
    my $data = {
        cols => [],
        flat => 1,
    };
    my $rows = $data->{rows} = [];
    my $colMap;
    my $jsonCol = {};
    if (open(TSV, "<$tmp")) {
        while (<TSV>) {
            s/[\n\r]+$//;
            if (/^\Q$headerPrfx\E(.+)/) {
                # This is a header row
                # Due to reassembly of forked files, there may be multiple
                # header rows, and potentially with different column order
                my @src = split(/\t/, $1);
                $colMap = {};
                for my $si (0..$#src) {
                    my $oi = $usePsh->data_index( $data, $src[$si], 'add');
                    unless ($oi == -1) {
                        $colMap->{$si} = $oi;
                    }
                }
                $jsonCol = $usePsh->json_encoded_columns( $data );
            } elsif (/^\Q$filtPrfx\E(.+)/) {
                my $filts = decode_json( $1 );
                push @{$data->{filters}}, @{$filts};
            } else {
                my @src  = split(/\t/, $_);
                my @row;
                push @{$rows}, \@row;
                while (my ($si, $oi) = each %{$colMap}) {
                    $row[$oi] = $src[$si];
                }
                while (my ($ind, $nullCB) = each %{$jsonCol}) {
                    $row[$ind] = decode_json( $row[$ind] || &{$nullCB} );
                }
            }
        }
        close TSV;
        # data_index() will have set 'need' flags, we need to clear them:
        delete $data->{need};
        return &generate_output( $data, $file );
    } else {
        &msg("[!!]","Failed to read output", $file, $!);
        return;
    }
}

sub generate_output {
    my ($data, $file) = @_;

    if ($format eq 'JSON') {
        # &_prepare_expanded_columns( $data );
        $file ||= &pick_outfile( 'force' );
        my $jf  = &make_json( $data, $file);
        my $url = &feed_ext_json( $jf );
        &msg("[DEBUG]", "JSON file generated", $jf, $url);
        return $file;
     }
    
    if ($mode eq 'Overview') {
        # &prebranch($data);
        $data = $usePsh->compact_hits( $data );
    }
    &set_query_column_info( $data );

    if ($format eq 'Excel') {
        &make_excel( $data, $file );
    } elsif ($format eq 'Text') {
        &make_txt( $data, $file);
    } elsif ($format eq 'HTML') {
        &make_html( $data, $file);
    } else {
        &make_tsv( $data, $file);
    }
    return $file;
}

sub feed_ext_json {
    my $file = shift;
    return 0 unless ($file);
    my $url = $args->path2url( $file );
    unless ($url) {
        &msg("Could not convert file path into URL", $file);
        return 0;
    }
    return $url if ($nocgi);
    my $extjs = &ext_setup_script( $formId );

    print <<EOF;
$extjs
<script>
    pgtm_basic_ajax_launch('$url', null, sos_fetch_resource);
</script>
EOF

    return $url;
}

sub show_stdout {
    my $file = shift;
    return unless ($file);
    return unless (lc($fileReq || "") eq 'stdout' && -s $file);
    if (open(FILE, "<$file")) {
        while (<FILE>) {
            print $_;
        }
        close FILE;
    } else {
        &msg("[!]","Failed to print output to STDOUT", $!, $file);
    }
}

sub _header_row {
    my ($colNames, $inds) = @_;
    $inds ||= [ 0..$#{$colNames} ];
    my @rv = map { defined $_ ? $_ : "" } map { $colNames->[$_] } @{$inds};
    map { s/^\Q$metaPrfx\E// } @rv;
    return @rv;
}

sub make_tsv {
    my ($data, $file) = @_;
    unless (open(TSV, ">$file")) {
        &task_err("Failed to write TSV output", $file, $!);
        return;
    }
    my ($colNames, $inds) = &_prepare_expanded_columns( $data );
    # die $args->branch($data);
    print TSV join("\t", &_header_row($colNames, $inds))."\n";
    #warn $args->branch($defaultCols);
    # print $args->branch( -ref => $data, -skipkey => [qw(-1 1)]);
    my $reformat = $usePsh->column_format_callbacks($data, 'text', $colNames);
    foreach my $row (@{$data->{rows}}) {
        my @newRow;
        for my $i (0..$#{$inds}) {
            my $src   = $inds->[$i];
            my $val   = $row->[$src];
            if (my $cb = $reformat->{$src}) {
                $val = &{$cb}( $usePsh, $val );
            }
            $val = "" unless (defined $val);
            push @newRow, $val;
        }
        print TSV join("\t", @newRow)."\n";
    }
    close TSV;
    &report_file($file, "TSV file created");
}

sub make_html {
    my ($data, $file) = @_;
    my $fh;
    if ($file) {
        unless (open($fh, ">$file")) {
            &task_err("Failed to write HTML output", $file, $!);
            return;
        }
    } else {
        $fh = *STDOUT;
    }
    my ($colNames, $inds) = &_prepare_expanded_columns( $data );
    my ($cnCB) = &_get_class_name_formatters( $inds, $colNames );
    print $fh &_html_header() unless ($clean);
    print $fh "<table class='tab'><tbody>\n";
    print $fh "  <tr>".join('', map { "<th>$_</th>" } 
                            &_header_row($colNames, $inds))."</tr>\n";
    my $reformat = $usePsh->column_format_callbacks($data, 'text', $colNames);
    my $rows     = $data->{rows};
    # warn $args->branch({ raw => $data->{cols}, inds => $inds, names => $colNames});
    for my $r (0..$#{$rows}) {
        my $row = $data->{rows}[$r];
        # warn $args->branch($row);
        my @newRow;
        for my $i (0..$#{$inds}) {
            my $src   = $inds->[$i];
            my $val   = $row->[$src];
            if (my $cb = $reformat->{$src}) {
                $val = &{$cb}( $usePsh, $val );
            }
            my $cls = &{$cnCB->[$i]}($usePsh, $val, $data, $src, $r);
            $val = "" unless (defined $val);
            push @newRow, sprintf("<td%s>%s</td>", $cls ? " class='$cls'" : "", $val);
        }
        print $fh "  <tr>".join('', @newRow )."</tr>\n";
    }
    print $fh "</tbody>\n";
    # printf($fh " <caption><span class='seq'>%s</span> <span class='date'>Run on %s</span></caption>\n", $seq, $date);
    print $fh "</table>\n";
    print $fh "</body></html>\n" unless ($clean);
        
    if ($file) {
        close $fh;
        &report_file($file, "HTML table created");
    }
}

sub task_err {
    my $obj = $forkPgtm || $pgtm;
    if ($obj && $obj->task_id()) {
        my $errTxt = join('; ', map { defined $_ ? $_ : '-UNDEF-' } @_);
        my $meta = $obj->task_metadata();
        push @{$meta->{error}}, $errTxt;
        $obj->change_task_metadata( $meta );
    } else {
        &msg("[!!]", @_ );
    }
}

sub make_json {
    my ($data, $file) = @_;
    my $fh;
    if ($file) {
        unless (open($fh, ">$file")) {
            &task_err("Failed to write json output", $file, $!);
            return undef;
        }
    } else {
        $fh = *STDOUT;
    }
    &set_query_column_info();
    if ($mode eq 'Overview') {
        # print $args->branch($data); die
        &note_genes( $data );
        $data = $usePsh->compact_hits( $data );
    } elsif ($mode eq 'Details') {
        &note_genes( $data );
    }
    # warn $#{$data->{rows}} . " rows total";
    my ($colNames, $inds) = &_prepare_expanded_columns( $data );
    my $reformat = $usePsh->column_format_callbacks($data, 'json', $colNames);
    my $json     = $usePsh->data_to_sencha( -data     => $data,
                                            -inds     => $inds,
                                            -reformat => $reformat );
    $json->{svctype} = $svcReq;
    print $fh encode_json( $json );
    if ($file) {
        close $fh;
        &report_file($file, "JSON data structure generated");
        return $file;
    }
    return undef;
}

sub report_file {
    my $file = shift;
    my @out = ("No file path provided");
    if ($file) {
        @out = ($file);
        if (my $url = $args->path2url( $file )) { 
            push @out, $nocgi ? $url : "<a targ='json' href='$url'>$url</a>"; 
        }
        if (my $msg = shift) { unshift @out, $msg; }
    }
    &msg("[OUT]", @out );
}

sub note_genes {
    my $data = shift;
    return unless ($data);
    my $ind = $usePsh->data_index($data, 'GeneDetails');
    if ($ind != -1) {
        my $rows = $data->{rows}  || [];
        my $targ = $data->{genes} ||= {};
        foreach my $row (@{$rows}) {
            if (my $gd = $row->[$ind]) {
                foreach my $gid (keys %{$gd}) {
                    $targ->{$gid} ||= [ $usePsh->gene_info( $gid ) ];
                }
            }
        }
        return;
    }
    $ind = $usePsh->data_index($data, 'GeneID');
    if ($ind != -1) {
        my $rows = $data->{rows}  || [];
        my $targ = $data->{genes} ||= {};
        foreach my $row (@{$rows}) {
            if (my $gid = $row->[$ind]) {
                $targ->{$gid} ||= [ $usePsh->gene_info( $gid ) ];
            }
        }
        return;
    }
}

sub _prepare_expanded_columns {
    my ($data, $colReq, $addReq) = @_;
    $addReq   = $addCols unless (defined $addReq);
    $colReq ||= $defaultCols;
    my ($colNames, $inds, $expand) = 
        $usePsh->custom_column_order( $data, $colReq, $addReq );
    # Expand the columns, if needed
    foreach my $expCB (@{$expand}) {
        my ($cb, $cols) = ($expCB->{cb}, $expCB->{cols});
        &{$cb}($usePsh, $data, $cols);
    }
    # warn $args->branch($defaultCols);
    #for my $i (0..$#{$colNames}) {
    #    $colNames->[$i] =~ s/^\Q$metaPrfx\E// if ($colNames->[$i]);
    #}
    return ($colNames, $inds);
}

sub make_txt {
    my ($data, $file) = @_;
    unless (open(TXT, ">$file")) {
        &task_err("Failed to write TXT output", $file, $!);
        return;
    }
    #warn $args->branch($defaultCols);
    my ($colNames, $inds) = &_prepare_expanded_columns( $data );
    #warn $args->branch($data->{cols});
    my $reformat = $usePsh->column_format_callbacks($data, 'text', $colNames);
    my @head = &_header_row($colNames, $inds);
    my @wids = map { CORE::length($_) } @head;
    my @show;
    foreach my $row (@{$data->{rows}}) {
        my @newRow;
        push @show, \@newRow;
        for my $i (0..$#{$inds}) {
            my $src   = $inds->[$i];
            my $val   = $row->[$src];
            if (my $cb = $reformat->{$src}) {
                $val = &{$cb}( $usePsh, $val );
            }
            $val = "" unless (defined $val);
            push @newRow, $val;
            my $len   = CORE::length($val);
            $wids[$i] = $len if ($wids[$i] < $len);
            # warn "[$c] $len -> $wids[$c]\n";
        }
    }
    # &prebranch(\@show); die;
    my $fmt = '| '.join(' | ', map { '%-'.$_.'s' } @wids). ' |'."\n";
    my $bar = sprintf($fmt, map {''} @wids);
    $bar =~ s/ /\-/g;
    $bar =~ s/\|/\+/g;
    
    print TXT $bar;
    printf(TXT $fmt, @head);
    print TXT $bar;
    foreach my $row (@show) {
        printf(TXT $fmt, @{$row});
    }
    print TXT $bar;
    if (my $filts = $data->{filters}) {
        print TXT "Filtered by: ".join(', ', @{$filts})."\n";
    }
    close TXT;
    &report_file($file, "Text file created");
    # warn $fmt;
    
    return;
}

sub make_excel {
    my ($data, $file) = @_;
    my $eh = BMS::ExcelHelper->new( $file );
    &excel_formats($eh);
    my $sname = "Hits";
    my ($colNames, $inds) = &_prepare_expanded_columns( $data );
    my $reformat = $usePsh->column_format_callbacks($data, 'excel', $colNames);
    my ($ehFmtCB, $wids ) = &_get_class_name_formatters( $inds, $colNames );
    my @head = &_header_row($colNames, $inds);
    $eh->sheet( -name    => $sname,
                -freeze  => [1],
                -width   => $wids,
                -columns => \@head, );

    map { $usedCols->{$_}++ } @head;
    &_make_excel_hits($eh, $sname, undef, $data,
                      $inds, $reformat, $ehFmtCB);
    &_make_excel_dbsum( $eh, $data );
    &_make_excel_ng($eh, $data);

    &help_sheet($eh, $data);
    $eh->close();
    my $eFile = $eh->file_path();
    
}

sub _get_class_name_formatters {
    my ($inds, $colNames) = @_;
    my (@wids, @formatters);
    my $colDet = $usePsh->details('columns');
    # die $args->branch($colDet);
    for my $i (0..$#{$inds}) {
        my $src   = $inds->[$i];
        my $cn    = $colNames->[$src];
        my $stnd  = $usePsh->standard_data_name( $cn );
        my $wid   = $usePsh->get_metadata('columns', $cn, 'width') || 10;
        push @wids, $wid;
        if (my $efcb = $usePsh->get_metadata('columns', $cn, 'classNameCB')) {
            # &msg("$cn has callback");
            push @formatters, $efcb;
        } else {
            #warn "$cn = ".$args->branch($colDet->{$stnd});
            push @formatters, sub { return undef; };
        }
    }
    return (\@formatters, \@wids);
}

sub _make_excel_hits {
    my ($eh, $sname, $srcRows, $data, $inds, $reformat, $ehFmtCB) = @_;
    $srcRows ||= $data->{rows};
    # warn $args->branch({inds => $inds, data => $srcRows});
    for my $r (0..$#{$srcRows}) {
        my (@vals, @fmts);
        for my $i (0..$#{$inds}) {
            my $c      = $inds->[$i];
            my $val    = $srcRows->[$r][$c];
            if (my $cb = $reformat->{$c}) {
                $val = &{$cb}( $usePsh, $val );
            }
            push @vals, $val;
            push @fmts, &{$ehFmtCB->[$i]}
            ($usePsh, $val, $data, $c, $r, $srcRows);
        }
        $eh->add_row($sname, \@vals, \@fmts );
    }
}

sub _make_excel_ng {
    my ($eh, $data) = @_;
    # Skip this step unless the user has requested specific genes
    return unless ($keepReq);
    my $bldInd = $usePsh->data_index( $data, 'GenomeBuild', 'add');
    my $fpInd  = $usePsh->data_index( $data, 'GenomicFootprint', 'add');
    return if ($fpInd == -1 || $bldInd == -1);
    # die $args->branch($data->{need});
    my $ngRecover;
    foreach my $row (@{$data->{rows}}) {
        # Find RefSeqGene rows (Build = RSG)
        my $bld = $row->[$bldInd] || "";
        next unless ($bld eq 'RSG');
        my $foot = $row->[$fpInd] || "";
        if ($foot =~ /^(NG[^:]+):(.+)$/) {
            $ngRecover ||= {};
            $ngRecover->{$foot} = $1;
        }
    }
    return unless ($ngRecover);
    
    my @cols = ('BMS Number', "ASO Sequence", "GenomeStart", "GenomeEnd",
                "CDS/UTR", "Exon/Intron", 
                "Symbol", "Strand", "MM", "AM", "GenomeSeq", "HitType");
    push @cols, @{$foundMeta} if ($foundMeta);
    &msg("[#]","Adding RSG Sheet");
    my ($colNames, $inds) = &_prepare_expanded_columns( $data, \@cols, 0 );
    my $symInd = $usePsh->data_index( $data, 'Symbol');
    my (%ngData, %syms);
    foreach my $row (@{$data->{rows}}) {
        my $gid = $row->[$fpInd] || "";
        if (my $ng = $ngRecover->{$gid}) {
            push @{$ngData{$ng}}, $row;
            if (my $sym = $row->[$symInd]) {
                $syms{$ng}{$sym}++;
            }
        }
    }
    my @withSym;
    foreach my $ng (keys %ngData) {
        my $sH    = $syms{$ng} || {};
        my @syms  = sort {$sH->{$b} <=> $sH->{$a} } keys %{$sH};
        my $sname = (join(' ', @syms) || $ng). " RSG";
        push @withSym, [$sname, $ng];
    }
    
    # warn $args->branch({ colnames => $colNames, inds => $inds, aRow => $data->{rows}[0]});
    foreach my $ngd (sort {$a->[0] cmp $b->[0]} @withSym) {
        my ($sname, $ng) = @{$ngd};
        my $rows = $ngData{$ng};
        # die $args->branch({ rows => $rows->[0], cols => $data->{cols}, inds => $inds});
        my $reformat = $usePsh->column_format_callbacks($data, 'excel', $colNames);
        my ($ehFmtCB, $wids) = &_get_class_name_formatters( $inds, $colNames );
        my @head = &_header_row($colNames, $inds);
   
        $eh->sheet( -name    => $sname,
                    -freeze  => [1],
                    -width   => $wids,
                    -columns => \@head, );

        map { $usedCols->{$_}++ } @head;
        &_make_excel_hits($eh, $sname, $rows, $data,
                           $inds, $reformat, $ehFmtCB);
    }
}

sub html5_run_oligo {
    my $srchids = shift;
    &init();
    if (!$srchids || $#{$srchids} == -1) {
        print "<i>No searches requested</i><br />\n";
        return;
    }
    my $svc = $mode eq 'Details' ? 'search_detail' : 
        $mode eq 'Overview' ? 'search_overview' : 'search_summary';
    my @ajaxCalls;
    print &ext_setup_script( $formId );
    print "<div class='box'><b>Pending Tasks</b>\n";
    foreach my $srchDat (@{$srchids}) {
        my ($srchid, $isNew) = @{$srchDat};
        my $id     = "srch$srchid";
        my $params = &ajax_call( service => $svc,
                                 srchid  => $srchid );
        push @ajaxCalls, [$params, $id];
        my ($seq, $name, $descr, $date) = $usePsh->search_details( $srchid );
        print "<div class='box' id='$id'>\n";
        print "<b>$seq</b> : ".($isNew ? "Running new search" : "Recovering data from $date");
        print "</div>";
    }
    print "</div>\n";
    print "<div id='exterr'></div>\n";
    print "<script>\n".join("\n", map { 
        "  launch_service($_->[0], '$_->[1]');"} @ajaxCalls ).
            "\n</script>\n";
    # &prebranch(\%ENV);
}

sub service {
    my $svc    = lc(shift || "");
    my $srchid = $args->val('srchid');
    $args->death("Unused method service() called");
    # used to do OLDservice...
    &do_exit();
}

sub ajax_call {
    my %params = ( mismatch => $maxMiss,
                   strand   => $geneStr,
                   genereq  => $keepReq,
                   exclreq  => $exclReq,
                   bgnd     => ++$bgndCnt,
                   @_ );
    return encode_json( \%params );
}

sub summary_html_for_srchid {
    my $srchid = shift;
    if (my $err = &start_or_wait_for_search( $srchid )) {
        return "<p class='err'>$err</p>\n";
    }
    my ($seq, $name, $descr, $date, $mm) = $usePsh->search_details( $srchid );
    return  "<p class='err'>srch_id = $srchid not found</i>\n" unless ($seq);
    my $dat = $usePsh->search_summary(  $srchid, $maxMiss, $geneStr );
    my $html = "<table class='tab'>\n";
    $html .= " <tbody>\n";
    my @head = ('Build', (map { "MM$_" } (0..$maxMiss)), "Database");
    $html .= "<tr>".join('', map { "<th>$_</th>" } @head)."</tr>\n";
    my $tot = 0;
    foreach my $srcid (sort { $a <=> $b } keys %{$dat}) {
        my $info = $usePsh->source_info( $srcid );
        $html .= sprintf("<tr><th>%s</th>", $info->{GenomeBuild} || "");
        for my $mm (0..$maxMiss) {
            my $num = $dat->{$srcid}{$mm} || 0;
            $tot += $num;
            $html .= sprintf("<td class='cen'>%s</td>", $num || "");
        }
        $html .= sprintf("<td class='db'>%s</td></tr>\n", $info->{Path} || "?");
    }
    $html .= "</tbody>\n";
    $html .= sprintf(" <caption><span class='seq'>%s</span> <b>%d hits</d> <span class='date'>Run on %s</span></caption>\n", $seq, $tot, $date);
    $html .= "</table>\n";
    # $html .= $args->branch($dat);
    return $html;
}

sub summary_grid_for_srchid {
    my $srchid = shift;
    my $hash;
    if (my $err = &start_or_wait_for_search( $srchid )) {
        $hash = &ext_error( $err );
        return encode_json( $hash );
    }
    my ($seq, $name, $descr, $date, $mm) = $usePsh->search_details( $srchid );
    unless ($seq) {
        $hash = &ext_error( "$srchid not found" );
        return encode_json( $hash );
    }
    my $dat = $usePsh->search_summary( -searchid => $srchid,
                                        @filters );
    my @head = ('oligoid', 'oligoseq', 'build');
    push @head, map  { "mm$_" } (0..$maxMiss);
    push @head, 'database';
    my @rows;
    foreach my $srcid (sort { $a <=> $b } keys %{$dat->{$srchid}}) {
        my $info = $usePsh->source_info( $srcid );
        my $row = {
            oligoid  => $name,
            oligoseq => $seq,
            database => $info->{Path},
            build    => $info->{GenomeBuild},
        };
        for my $mm (0..$maxMiss) {
            $row->{"mm$mm"} = ($dat->{$srchid}{$srcid}{$mm} || 0) + 0;
        }
        push @rows, $row,
    }
    return encode_json( {
        cols => \@head,
        rows => \@rows,
        svctype => 'summary_grid',
                        });
}

sub ext_error {
    my $err = shift;
    my $hash = {
        type => 'error',
        text => $err,
    };
    return $hash;
}

sub overview_grid_for_srchid {
    my $srchid = shift;
    if (my $err = &start_or_wait_for_search( $srchid )) {
        my $hash = &ext_error( $err );
        print encode_json( $hash );
        return;
    }
    my $flat = &recover_old_search( $srchid );
    return &overview_grid_for_rows( $flat, $srchid );
}

sub overview_grid_for_rows {
    my ( $flat, $srchid ) = @_;
    my $compact = $usePsh->compact_hits( $flat, $srchid );
    my @rows;
    my @head = qw(oligoid oligoseq build ex0 int0 ex1 int1 ex2 int2 ex0sym int0sym ex1sym int1sym ex2sym int2sym);
    foreach my $srchid  (sort { $a <=> $b } keys %{$compact}) {
        my ($seq, $name, $descr, $date, $mm) = 
            $usePsh->search_details( $srchid );
        my $locSTH  = $usePsh->named_sth('Get Gene Information');
        my $srcH = $compact->{$srchid};
        my @allBuilds = sort { $a <=> $b } keys %{$srcH};
        my @allEts    = qw(Ex Int);
        my @allMMs    = (0..3);
        my $tot  = ($#allBuilds + 1) * ($#allEts + 1) * ($#allMMs + 1);
        my $done = 0;
        foreach my $srcid (@allBuilds) {
            my %geneCache;
            my $sh = $usePsh->source_info( $srcid );
            my $bld = $sh->{GenomeBuild};
            my $row = {
                oligoid  => $name,
                oligoseq => $seq,
                build    => $bld,
            };
            $bld ||= "UnkBuild";
            push @rows, $row;
            foreach my $et (@allEts) {
                my $exTrg = $srcH->{$srcid}{$et};
                for my $mm (0..3) {
                    $done++;
                    next unless ($exTrg);
                    my $gidH = $exTrg->[$mm];
                    next unless ($gidH);
                    $pgtm->change_task_status
                        ({task => "Create grid: $bld",
                          subtask => [["$et MM$mm", ($done-1)/$tot]]}) if ($pgtm);
                    my @gids = keys %{$gidH};
                    $row->{lc($et).$mm} = $#gids + 1;
                    my @locs;
                    foreach my $gid (@gids) {
                        my $ginfo = $geneCache{$gid} ||= 
                            $locSTH->selectall_arrayref( $gid ) || [];
                        my ($lid, $sym, $desc) = @{$ginfo->[0] || []};
                        next unless ($lid);
                        $sym =~ s/[\*\~]+$// if ($sym);
                        push @locs, [$lid, $sym, $desc];
                    }
                    @locs = sort { ($a->[1] || $fakeSym) cmp ($b->[1] || $fakeSym)
                                       || $a->[0] cmp $b->[0] } @locs;
                    $row->{lc($et).$mm.'sym'} = \@locs;
                }
            }
        }
    }
    return encode_json( {
        cols => \@head,
        rows => \@rows,
        svctype => 'overview_grid',
                        });
}

sub start_or_wait_for_search {
    my $srchid = shift;
    return "No search ID provided" unless ($srchid);
    return "srch_id = '$srchid' is malformed" unless ($srchid =~ /^\d+$/);
    my $ti = time;
    while (1) {
        my $stat = $usePsh->search_status( $srchid );
        if (defined $stat) {
            if ($stat eq '-1') {
                # Search is done
                last;
            } elsif ($stat !~ /^\d+$/) {
                # Should never happen, but we will be backticking later
                # so checking for safety
                return "Serious error - search status for srch_id = '$srchid' returns non-integer result";
            } elsif ($stat == 0) {
                # Need to run
                &run_new_search( $srchid );
                last;
            } elsif ($stat =~ /^\d+$/) {
                # $stat represents the PID for a running job
                my @found = split(/[\n\r]+/, 
                                  `ps -ef | grep $stat | grep -v grep`);
            }
        } else {
            return "srch_id = $srchid does not appear valid";
        }
        my $waiting = int(0.5 + 10 * (time - $ti) / 60) / 10;
        if ($waiting > 15) {
            return "srch_id = $srchid waited $waiting minutes with no result";
        }
        $pgtm->change_task_status
            ({task => "Waiting for Search $stat to finish",
             subtask => [["Wait", 0]]}) if ($pgtm);
        sleep(5);
    }
    return 0;
}

sub search_id_for_bioseq {
    my $bs = shift;
    my $isNew = 0;
    my $srchid;
    $srchid   = $usePsh->get_srch_ids( $bs, $maxMiss ) unless ($clobber);
    if ($srchid) {
        my $stat = $usePsh->search_status( $srchid );
        $isNew = 1 if (defined $stat && !$stat);
    } else {
        $isNew = 1;
        $srchid = $usePsh->set_srch_id( $bs, $maxMiss);
    }
    return ($srchid, $isNew);
}

sub oligo_summary {
    my $srchDat = shift;
    my ($srchid, $isNew) = @{$srchDat};
    &run_new_search( $srchid ) if ($isNew);
    my $sum = $usePsh->search_summary( -srchid => $srchid, @filters );
    &{$resultHandler}( $sum ) if ($resultHandler);
    return $sum;
    # warn $args->branch($sum);
}

sub run_oligo {
    my $srchDat = shift;
    my ($srchid, $isNew) = @{$srchDat};
    my $flat = $isNew ? 
        &run_new_search( $srchid ) : &recover_old_search( $srchid );
    return &{$resultHandler}( $flat ) if ($resultHandler);
    return 0;
}

sub run_new_search {
    my $srchid = shift;
    my ($seq) = $usePsh->search_details( $srchid );
    unless ($seq) {
        &task_err("Failed to recover sequence for srch_id = $srchid");
        return [];
    }
    $usePsh->search_status( $srchid, $$);
    &msg("[+]","Running new search for $seq") unless ($doHTML);
    $usePsh->set_tag( $srchid, "Mismatch", $maxMiss || 0);
    my $hits  = $usePsh->hits_for_sequence( -seq => $seq, -mm => $maxMiss );
    my $flat  = $usePsh->flatten_hits( $hits, $srchid );
    if ($#{$flat->{rows}} == -1) {
        &msg("[??]","Your search found no hits");
        $usePsh->clear_search( $srchid );
    } else {
        # warn $args->branch($flat);
        my $wrote = $usePsh->store_search_hits( $flat );
        if ($wrote > 0) {
            $usePsh->search_status( $srchid, -1);
        } else {
            &msg("[??]","Failed to store hits!");
            $usePsh->clear_search( $srchid );
        }
    }
    &filter_hits( $flat );
    return $flat;
}

sub recover_old_search {
    my $srchid = shift;
    my ($seq, $name, $desc, $dt) = $usePsh->search_details( $srchid );
    unless ($seq) {
        &task_err("Failed to recover sequence for srch_id = $srchid");
        return [];
    }
    &msg("[+]","Recovering data for $seq generated on $dt") unless ($doHTML);

    my $flat = &recover_hits( $srchid );
    &filter_hits( $flat );
    # &prebranch($flat); die;
    return $flat;
}

sub recover_hits {
    my $srchid = shift;
    return undef unless ($srchid);
    my $found = $usePsh->recover_search_hits
        ( -searchid   => $srchid,
          @filters );
    return $found;
}

sub filter_hits {
    my $hits = shift;
    return $usePsh->filter_flat_hits
        ( -hits       => $hits,
          @filters );
}

sub expand_oligo_info {
    my $self    = shift; # This will be a PgSeqHash object
    my $data    = $self->_assure_flattened_hits( shift );
    my $request = shift;
    my $srcCol  = 'SearchID';
    my @baseCol = qw(OligoAccession OligoSeq OligoDescription);
    # Add in any metadata columns, if present:
    map { push @baseCol, $_ if (/^\Q$metaPrfx\E/) } @{$request || []};
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $request, $srcCol, \@baseCol,
          sub {
              my $col = shift;
              if ($col eq 'OligoAccession') {
                  return sub { return shift->display_id() };
              } elsif ($col eq 'OligoSeq') {
                  return sub { return shift->seq() };
              } elsif ($col eq 'OligoDescription') {
                  return sub { return shift->desc() };
              } elsif ($col =~ /^\Q$metaPrfx\E(.+)/) {
                  my $tag = $1;
                  # warn  "TAG FOUND FOR '$tag'";
                  return sub {
                      my $bs = shift;
                      my @arr = @{$bs->{CAT_TAG}{$tag} || []};
                      return join(',', sort @arr) || "";
                  };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $input    = $data->{rows};
    my $output   = $data->{rows} = [];

    foreach my $row (@{$input}) {
        my $sr = $row->[$srcInd];
        my @srcs = (0);
        if (my $srchid = $row->[$srcInd]) {
            @srcs = @{$queries->{$srchid} || []};
            if ($#srcs == -1) {
                if (my $bs = $usePsh->bioseq_for_search_id( $srchid )) {
                    @srcs = ($bs);
                } else {
                    @srcs = (0);
                }
            }
        }
        foreach my $bs (@srcs) {
            my $copy = $#srcs == 0 ? $row : [ @{$row} ];
            push @{$output}, $copy;
            if ($bs) {
                while (my ($ind, $cb) = each %{$needed}) {
                    $copy->[$ind] = &{$cb}( $bs );
                }
            }
        }
    }
    $self->bench_end();
    return;
}

sub write_data_as_json {
    my $data = shift;
    my $rows    = $data->{rows};
    return if ($#{$rows} == -1);
    
}

sub store_data_in_file {
    my $data = shift;
    my $rows    = $data->{rows};
    return if ($#{$rows} == -1);
    my $jsonCol = $usePsh->json_encoded_columns( $data );
    # &prebranch($data); die;
    
    my @head  = @{$data->{cols}};
    my $headtxt = join("\t", @head);
    $fc->write_output('OUT', "$headerPrfx$headtxt\n");
    # &prebranch($data);
    foreach my $row (@{$data->{rows}}) {
        my @copy = @{$row};
        while (my ($ind, $nullCB) = each %{$jsonCol}) {
            $copy[$ind] = encode_json( $copy[$ind] || &{$nullCB} );
        }
        map { $_ = "" unless (defined $_) } @copy;
        map { s/[\t\r\n]+/ /g } @copy;
        $fc->write_output('OUT', join("\t", @copy)."\n");
    }
    if (my $filts = $data->{filters}) {
        $fc->write_output('OUT', $filtPrfx.encode_json($filts)."\n");
    }
    if (my $obj  = $forkPgtm || $pgtm) {
        my $meta = $obj->task_metadata();
        $meta->{done} = 1;
        $obj->change_task_metadata( $meta );
    }
}

sub set_query_column_info {
    my ($data) = @_;
    return if ($usePsh->{SOS_QCI_SET}++);
    # Find what meta tags, if any, were in the submitted queries:
    my %meta;
    while (my ($seqKey, $bsArr) = each %{$queries || {}}) {
        # Also associate the metadata with the searchID
        my ($srchid) = &search_id_for_bioseq( $bsArr->[0] );
        $queries->{$srchid} = $bsArr;
        foreach my $bs (@{$bsArr}) {
            if (my $ct = $bs->{CAT_TAG}) {
                map { $meta{$_}++ } keys %{$ct};
            }
        }
    }
    my @mc = sort keys %meta;
    $foundMeta = $#mc == -1 ? undef : [ map { "$metaPrfx$_" } @mc ];
    my $det = $usePsh->details('columns');
    foreach my $mcol (@mc) {
        $addCols ||= [];
        my $dCol = $metaPrfx . $mcol;
        push @{$addCols}, $dCol;
        $usePsh->reset_detail($dCol, 'expandCB', \&expand_oligo_info, 'columns');
        $usePsh->reset_detail($dCol, 'CBcol', $mcol, 'columns');
        $usePsh->reset_detail($dCol, 'desc', "User-provided column", 'columns');
    }
    # die $args->branch($queries);

    # Redefine the expansion method for oligo name, sequence and description:
    foreach my $col (qw(OligoAccession OligoSeq OligoDescription)) {
        $usePsh->reset_detail($col, 'expandCB', \&expand_oligo_info, 'columns');
        #$usePsh->reset_detail($col, 'foobar', \&expand_oligo_info, 'columns');
        #die $args->branch($det->{$col});
    }
}

sub _subject_base_class_hash {
    my %rv;
    # Upper case bases are normal, unclassed
    map { $rv{$_} = '' } qw(A C G T);
    # Upper case ambiguities are special:
    map { $rv{$_} = 'mmamb' } split('', 'KMRYSWBVHDXN-');
    # Lower case nucleotide letters are mismatched:
    map { $rv{$_} = 'mmseq' } split('', 'acgtkmryswbvhdxn');
    # Gap
    $rv{'-'} = 'mmgap';
    return \%rv;
}

sub _basic_html_cell {
    my $val = shift;
    return sprintf("<td>%s</td>", defined $val ? $val : "");
}

sub _mismatch_html_cell {
    my $val = shift;
    my $col = !defined $val ? 'mmunk' :
        $val == 0 ? 'mm0' :
        $val == 1 ? 'mm1' :
        $val == 2 ? 'mm2' : 'mm3';
    return sprintf("<td class='%s'>%s</td>", $col, defined $val ? $val : "");
}

sub _self_classed_html_cell {
    my $val = shift;
    $val = "" unless (defined $val);
    return sprintf("<td class='%s'>%s</td>", lc($val), $val);
}

sub show_bench {
    my $bench = $args->val('benchmark');
    return "" unless ($bench);
    return if (!$doneSomething);
    my $txt = $bnch->show_bench( $doHTML ? (html => 1) : (-shell => 1 ));
    return $txt;
}


sub init {
    $forkPsh  = $usePsh = BMS::PgSeqHash->new( );
    $srcCache = {};
}

sub subject_while_forked {
    my ($sbjId, $srcId) = @_;
    my $sc = $srcCache->{$srcId}; 
    unless ($sc) {
        $sc = $srcCache->{$srcId} = $usePsh->source_info( $srcId );
        my @keys = sort keys %{$sc};
        my $path = $sc->{Path};
        $sc->{FastaMeta} = BMS::FastaMeta->new(-fasta => $path );
    }
    my $fm = $sc->{FastaMeta};
    return $fm->fetch( $sbjId );
}

sub final {
    if (my $txt = &show_bench()) {
        unless ($txt =~ /\b0\.\d+ seconds/) {
            &msg("[DEBUG]", "Forking Benchmarks");
            print STDERR $txt;
        }
    }
}

sub project_base {
    unless ($projBase) {
        $projBase = `date +'%Y-%m-%d.%H:%M:%S'`;
        $projBase =~ s/[\n\r]+$//;
    }
    return $projBase;
}

sub project_file {
    return "" unless ($proj);
    my ($name, $prfx, $sfx) = @_;
    return "" unless ($name);
    my $path = $proj . '/';
    $path .= $prfx if ($prfx);
    $path .= $name;
    $path .= $sfx if ($sfx);
    $args->assure_dir($path, 'isFile');
    return $path;
}

sub oligo_file {
    my ($oid) = @_;
    my $hDir = &project_file( "html/oligos" );
    unless ($hDir) {
        $hDir = $outF || "";
        $hDir =~ s/[^\/]+$//;
        $hDir ||= "./";
        $hDir .= "html/oligos";
    }
    $args->assure_dir($hDir);
    return "$hDir/$oid.html";
}

sub oligo_url {
    my $oid = shift;
    my $file = &oligo_file($oid);
    return $args->path2url( $file );
}

sub gene_file {
    my ($ll) = @_;
    my $hDir = &project_file( "html/genes" );
    unless ($hDir) {
        $hDir = $outF || "";
        $hDir =~ s/[^\/]+$//;
        $hDir ||= "./";
        $hDir .= "html/genes";
    }
    $args->assure_dir($hDir);
    return "$hDir/$ll.html";
}

sub gene_url {
    my $ll = shift;
    my $file = &gene_file($ll);
    return $args->path2url( $file );
}

sub set_out_fh {
    my $file = shift;
    unless ($file) {
        &msg_once("[>]","Output to STDOUT");
        return *STDOUT;
    }
    my $oFH;
    open($oFH, ">$file") || $args->death
        ("Failed to open output", $file, $!);
    select((select($outFh), $| = 1)[$[]); # Bob's call to autoflush
    &msg_once("[>]", "Output to:", $file);
    return $oFH;
}

sub record {
    $bnch->bench_start();
    my ($canInd, $mm, $db, $ssInd, $id, $subSeq) = @_;
    my ($qDat, $str) = @{$search->{cand}[$canInd]};
    my $qlen   = $qDat->{len};
    my $subEnd = $ssInd + $qlen;
    my $subObj = $nowMeta->fetch($id);


    # &msg_once("[<]", $subObj->to_text() );
    if (time - $lastRpt > 60 * 15) {
        # Note progress every 15 minutes
        &msg_once("[+]", $subObj->id());
        $lastRpt = time;
    }
    
    if ($qlen < $search->{maxlen}) {
        # Need to trim of right padding off subject oligo
        $subSeq = substr($subSeq, 0, $qlen);
    }
    if ($str eq 'R') {
        # Need to reverse complement the subject oligo
        $subSeq = $su->revcom( $subSeq );
    }
    my ($s, $e) = ($ssInd + 1, $subEnd);
    my $hit = BMS::FastaMeta::Hit->new
        ( -queryseq  => $qDat->{seq},
          -queryid   => $qDat->{id},
          -querydesc => $qDat->{desc},
          -pad       => $mm,
          -start     => $s,
          -end       => $e,
          -strand    => $str,
          -subject   => $subObj, );

    # die "$ssInd = ".$args->branch($hit);

    my @rows = $hit->row_data( $chkStrand );
    &tsv_out( \@rows );
    $bnch->bench_end();
}

sub locs_to_rows {
    $bnch->bench_start();
    
    $bnch->bench_end();
}

sub tsv_out {
    my $rows = shift;
    my $rn   = $#{$rows} + 1;
    return unless ($rn);
    unless ($tsvRows) {
        print $outFh join("\t", @{$head}) . "\n";
    }
    $tsvRows += $rn;
    foreach my $rh (@{$rows}) {
        my @row = map { defined $_ ? $_ : "" } map { $rh->{$_} } @{$head};
        print $outFh join("\t", @row)."\n";
    }
}

sub get_oligos {
    $bnch->bench_start();
    my %rv;
    my $num = 0;
    foreach my $param (@seqParams) {
        foreach my $bs (&bioseqs_for_param($param)) {
            $num++;
            my $key = uc( $bs->seq() );
            push @{$rv{$key}}, $bs;
            my $desc = " ".($bs->desc() || "");
            while ($desc =~ /( \/([a-z0-9_]+)=\'([^\']+)\')/i ||
                   $desc =~ /( \/([a-z0-9_]+)=\"([^\"]+)\")/i ||
                   $desc =~ /( \/([a-z0-9_]+)=([^\"\'\s]+))/i) {
                my ($rep, $k, $v) = ($1, $2, $3);
                $desc =~ s/\Q$rep\E/ /;
                $bs->{CAT_TAG} ||= {};
                push @{$bs->{CAT_TAG}{$k}}, $v;
            }
        }
    }
    # if ($num) { &prebranch(\%rv); &do_exit(); }
    $bnch->bench_end();
    return \%rv;
}

sub bioseqs_for_param {
    my $param = shift;
    my @rv;
    my $byName = $usePsh->named_sth('Get search details for name');
    if (my $req = $args->rawval($param)) {
        my @reqs = ref($req) ? @{$req} : ($req);
        foreach my $r (@reqs) {
            next unless ($r);
            my $num = 0;
            foreach my $bs ( $css->clean_bioseq($r, 'split',
                                                \&_fetch_by_name) ) {
                my $seq = uc($bs->seq() || "");
                next unless ($seq);
                &_rename_bioseq( $bs ) if ($bs->display_id() eq $genericName);
                push @rv, $bs;
                $num++;
            }
            # &prebranch({ param => $param, req => $r, rv => \@rv});
        }
    }
    return wantarray ? @rv : $#rv == -1 ? undef : \@rv;
}

sub _fetch_by_name {
    my $id    = shift;
    my $clean = $id;
    $clean    =~ s/^\s+//;
    $clean    =~ s/\s+$//;
    my $found = $usePsh->search_details_by_name( $id );
    my (@rv, %uniq);
    foreach my $details (@{$found}) {
        my ($seq, $name, $desc) = @{$details};
        $uniq{uc($seq)} ||= [ $seq, $desc ];
    }
    foreach my $dat (values %uniq) {
        push @rv, Bio::PrimarySeq->new
            ( -seq        => $dat->[0],
              -alphabet   => 'DNA',
              -desc       => $dat->[1],
              -display_id => $clean, );
    }
    if ($#rv == -1) {
        &msg("[?]", "Failed to find a sequence named '$id'");
    }
    return @rv;
}

sub _rename_bioseq {
    my $bs     = shift;
    my $seq   = $bs->seq();
    my $srchid;
    my @found = sort $renameSth->get_array_for_field("Case $seq");
    if ($#found != -1) {
        $srchid = $found[0];
    } else {
        $srchid = $usePsh->get_srch_ids( $seq );
    }
    if ($srchid) {
        my ($seq, $name, $descr) = $usePsh->search_details( $srchid );
        $bs->display_id( $name );
        $bs->desc( $descr );
    }
}

sub _db_for_ID {
    my $id = shift;
    my $rv;
    if ($id && $id =~ /^[^\.]+\.([^:]+)\:/) {
        ($rv) = &_add_db( $1 );
    }
    return $rv;
}

sub _add_db {
    my $file = shift;
    return () unless ($file);
    unless (-s $file) {
        my $chk = &_premrna_files($file);
        $file = $chk if ($chk);
    }
    return () if ($doneDB{$file}++);
    $bnch->bench_start();
    {
        # Make sure the index is calculated
        my $temp = BMS::FastaMeta->new( -fasta => $file );
    }
    my $short = $file;
    $short =~ s/.+\///;
    # $short = "[!] $short" unless ($file =~ /^\/gcgblast\//);
    my $db = {
        file  => $file,
        short => $short,
        hitcount => 0,
    };
    &_set_db_fh( $db );

    my ($tax, $ns, $bld);
    if ($short =~ /(\S+)_([^_]+)_pre-mRNA_(\S+)\.fa/) {
        ($tax, $ns, $bld) = ($1, $2, $3);
    } elsif (my $dat = $db->{now}) {
        if ($dat->[0] =~ /^.+\.(\S+):\d+\-\d+[\+\-]1 /) {
            $bld = $1;
        }
    }
    if ($tax) {
        $tax =~ s/_/ /g;
        substr($tax, 0, 1) = uc(substr($tax, 0, 1));
    }
    $db->{taxa}  = $tax || $unTaxa;
    $db->{build} = $bld || $unBuild;
    $db->{ns}    = $ns  || "UNK";
    
    # Set up output paths
    if ($proj) {
        # Project folder to collect results over time
        $db->{out} = &project_file
            ( &project_base(), "tsv/",  sprintf(".%05d", $$));
    } elsif ($outF) {
        # Single specified output file
        $db->{out} = $outF;
    }
    if ($db->{out}) {
        $db->{out} .= $bld ? "-$bld" : "-".++$iter;
        $db->{out} .= "-LIMIT$limit" if ($limit);
        $db->{tsv}  = $db->{out} . ".tsv";
        $db->{meta} = &_make_db_meta( $db );
    }
    return ($db);
}

sub _make_db_meta {
    my $db = shift;
    return undef unless ($db);
    my $file = $db->{meta};
    unless ($file) {
        if (my $tsv = $db->{tsv}) {
            $file = "$tsv.meta";
        }
    }
    return undef unless ($file);
    return $file if (-s $file && !$clobber);
    open(META, ">$file") || $args->death
        ("Failed to create database meta file", $file, $!);
    print META "# simpleOligoSearch database meta file\n";
    foreach my $k (sort keys %{$db}) {
        my $v = $db->{$k};
        next if (!defined $v || ref($v));
        printf(META "%s : %s\n", $k, $v);
    }
    close META;
    return $file;
}

sub _read_db_meta {
    my $tsv = shift || "";
    return $metaCache{$tsv} if (defined $metaCache{$tsv});
    return $metaCache{$tsv} = 0 unless ($tsv);
    my $file = "$tsv.meta";
    return $metaCache{$tsv} = 0 unless (-s $file);
    $bnch->bench_start();
    my $db = $metaCache{$tsv} = {};
    open(META, "<$file") || $args->death
        ("Failed to read database meta file", $file, $!);
    while (<META>) {
        s/[\n\r]+$//;
        if (/([^:]+?)\s+:\s+(.+?)\s*$/) {
            my ($k, $v) = ($1, $2);
            $db->{$k} = $v;
        }
    }
    close META;
    if (my $file = $db->{file}) {
        if (-s "$file.dbfile") {
            my $fm = BMS::FastaMeta->new( -fasta => $file );
            $db->{FM} = $fm;
        }
    }
    $bnch->bench_end();
    return $db;
}

sub _close_fh {
    my $db = shift;
    if (my $fh = $db->{fh}) {
        close $fh;
    }
    delete $db->{fh};
}

sub _set_db_fh {
    my $db   = shift;
    if (my $fh = $db->{fh}) {
        return $fh;
    }
    my $file = $db->{file};
    delete $db->{now};
    my $fh;
    open($fh, "<$file") || $args->death
        ("Failed to open fasta database", $file, $!);

    $db->{fh} = $fh;
    while (<$fh>) {
        if (/^>(.+)/) {
            $db->{now} = [ $1, 0, ""];
            last;
        }
    }
    $args->death("Failed to identify fasta header in file", $file)
        unless ($db->{now});
    return $fh;
}


sub html_styles {
    $bnch->bench_start();
    my $allBuild = shift || [];
    my $com = '/* http://www.codeproject.com/Tips/166266/Making-Text-Upside-down-using-CSS */';
    my $rv = <<EOF;
 .ol  { color:green; }
 .fwd { color: green;  font-style: italic; }
 .mat { color:\#ccc;}
 .err { color:pink; }
 .mis { background-color:yellow; }
 .ext { background-color:lime; }
 .bld { background-color:\#fac; }
 .rev { color: blue; font-style: italic; }
 .nul { color: \#888; font-style: italic; font-weight: bold;}
 .tax { color:brown; font-style: italic;  }
 .seq { font-family: monospace; color: #060; }
.rctok { color: red; font-size: 0.8em; padding-left: 3px; vertical-align: middle; }

 .bigsym   { font-weight: bold; font-size: 1.3em; color: #360; text-align: center; font-family: sans-serif; }
 .linkid   { font-weight: bold; font-size: 1.3em; color: #603; }
 .note     { color:#066; font-style: italic;  }
 .normid   { font-weight: bold; color: #603; }
 .bigid    { font-weight: bold; font-size: 2em; color: #603; }
 .file     { color: #93f; font-size: 0.7em; font-family: monospace; }
 .rnaonly  { color: brown; font-style: italic; font-size: 0.8em; }
 .exonexon { color: purple; font-style: italic; font-size: 0.8em; }
 .multihit { color: red; font-style: italic; font-size: 0.8em; }
 .variant  { color: #309; font-weight: bold; }

.tab, .tab table
{ border-collapse: collapse; }
.tab th
{ background-color: #ffc; }
.tab td
{ background-color: #fff; }
.tab th, .tab td
{ border: #fc9 solid 1px; padding: 2px;
  empty-cells: show; vertical-align: top; }

$com
.revcom {
    color : blue;
   filter:  progid:DXImageTransform.Microsoft.BasicImage(rotation=2);
   ms-filter: "progid:DXImageTransform.Microsoft.BasicImage(rotation=2)";
   -moz-transform: rotate(-180deg);  /* FF3.5+ */
   -o-transform: rotate(-180deg);  /* Opera 10.5 */
   -webkit-transform: rotate(-180deg);  /* Safari 3.1+, Chrome */
   position: absolute; 
}

EOF
    foreach my $bld (@{$allBuild}) {
        my $hex = $cu->pastel_text_color($bld);
        $rv .= sprintf(".%s { background-color: %s }\n", $bld, $hex);
    }
    $bnch->bench_end();
    return $rv;
}

sub db_stats {
    $bnch->bench_start();
    my ($allRows, $dbLU) = @_;
    foreach my $row (@{$allRows}) {
        next if ($row->{SKIP});
        my $gnm  = $row->{GenomeFootprint} || "";
        my $rbld = $row->{Build} || "";
        $rbld = "" if ($rbld eq $unBuild);
        if ($gnm =~ /^.+\.([^\.]{4,}):/) {
            my $bld = $1;
            if ($rbld =~ /^.+\.\Q$bld\E/) {
                # if ($bld =~ /^(scaffold|contig|\d+\.)/) {
                # Misparsed builds for contigs
                &msg_once("[+]", "Correcting build '$rbld' to '$bld'");
                $row->{Build} = $bld;
            }
            my $db = $dbLU->{$bld} ||= {
                build => $bld,
                file  => "",
            };
            my $mm = $row->{Mis} || 0;
            if (!defined $db->{Locs}{$gnm} || $db->{Locs}{$gnm} > $mm) {
                $db->{Locs}{$gnm} = $mm;
            }
        }
    }

    my @dbObjs = values %{$dbLU};
    foreach my $db (@dbObjs) {
        my $lh       = $db->{Locs} || {};
        my @locs     = keys %{$lh};
        my $sc       = 0;
        my $rpt      = $db->{report} ||= {
            Build    => $db->{build} || $unBuild,
            Database => $db->{short} || "",
            Path     => $db->{file}  || "",
            Taxa     => $db->{taxa}  || $unTaxa,
            NS       => $db->{ns}    || "",
        };
        $rpt->{Hits}  = $#locs + 1;
        map { $sc += 1 / (1 + $lh->{$_}) } @locs;
        $rpt->{Score} = int(0.5 + 10 * $sc) / 10;
    }
    my @allBuild = map { $_->{Build} }  sort { 
        $b->{Score} <=> $a->{Score} ||
            $b->{Hits} <=> $a->{Hits} ||
            $a->{Taxa} cmp $b->{Taxa} } map { $_->{report} } @dbObjs;
    for my $i (0..$#allBuild) {
        my $bld = $allBuild[$i];
        if (my $db = $dbLU->{ $bld }) {
            $db->{report}{Rank} = $i + 1;
        }
    }
    # die $args->branch(\@allBuild);
    $bnch->bench_end();
    return \@allBuild;
}

sub _make_excel_dbsum {
    my ($eh, $data) = @_;
    my $srcRows = $data->{rows};
    my $srcInd = $usePsh->data_index( $data, 'SourceID');
    return if ($srcInd == -1);
    my @cols = ("Build", "Species");
    # Are there summarized mismatch count columns?
    my ($mms, $mmInd);
    for my $i (0..5) {
        my $col = "MM$i";
        my $ind = $usePsh->data_index( $data, $col);
        last if ($ind == -1);
        $mms ||= [];
        push @{$mms}, $ind;
        push @cols, $col;
    }
    # Is there a specific mismatch column?
    unless ($mms) {
        my $ind = $usePsh->data_index( $data, 'MisMatch');
        unless ($ind == -1) {
            push @cols, map { "MM$_" } (0..2);
            $mmInd = $ind;
        }
    }
    my $summary = {};
    for my $r (0..$#{$srcRows}) {
        my $row = $srcRows->[$r];
        my $srcid = $row->[$srcInd];
        next unless ($srcid);
        my $targ = $summary->{$srcid} ||= {
            SourceID => $srcid,
            Order => 0,
        };
        if (defined $mmInd) {
            my $mm = $row->[$mmInd];
            if (defined $mm && $mm <= 2) {
                $targ->{"MM$mm"}++;
            }
        } elsif ($mms) {
            for my $j (0..$#{$mms}) {
                if (my $num = $row->[$mms->[$j]]) {
                    $targ->{"MM$j"} += $num;
                }
            }
        }
    }
    my @weights = (1, 0.3, 0.1);
    foreach my $targ (values %{$summary}) {
        my $srcid = $targ->{SourceID};
        my $info = $usePsh->source_info( $srcid );
        map { $targ->{$_} ||= $info->{$_} } ('Species', 'WordSize');
        $targ->{Build} = $info->{GenomeBuild};
        $targ->{DatabasePath} = $info->{Path};
        for my $j (0..$#weights) {
            if (my $num = $targ->{"MM$j"}) {
                $targ->{Order} += $weights[$j] * $num;
            }
        }
    }
    # my @icols = qw(Path Entries Characters WordSize);
    push @cols, 'Order';
    push @cols, qw(DatabasePath WordSize);
    my ($ehFmtCB, $wids ) = &_get_class_name_formatters( [0..$#cols], \@cols );
    my $reformat = $usePsh->column_format_callbacks
        ({ cols => \@cols }, 'excel', \@cols);
    
    my @head  = @cols;
    my $sname = "DBs";
    $eh->sheet( -name    => $sname,
                -freeze  => [1],
                -width   => $wids,
                -columns => \@head, );
    map { $usedCols->{$_}++ } @head;
    foreach my $targ ( sort { $b->{Order} <=> $a->{Order} ||
                                  $a->{Build} cmp $b->{Build} }
                       values %{$summary}) {
        my @row  = map { $targ->{$_} } @cols;
        my (@vals, @fmts);
        for my $c (0..$#row) {
            my $val    = $row[$c];
            if (my $cb = $reformat->{$c}) {
                $val = &{$cb}( $usePsh, $val );
            }
            push @vals, $val;
            push @fmts, &{$ehFmtCB->[$c]}
            ($usePsh, $val);
        }
        $eh->add_row($sname, \@vals, \@fmts );
    }
}

sub nchoosek {
    my ($n, $k) = @_;
    $n ||= 0;
    $k ||= 0;
    unless ($nck{$n}{$k}) {
        if ($n && $k && $n > $k) {
            my $v = 1;
            map { $v *= $_ } (($n - $k + 1) .. $n);
            map { $v /= $_ } (1..$k);
            $nck{$n}{$k} = $v;
        } else {
            $nck{$n}{$k} = 1;
        }
    }
    return $nck{$n}{$k};
}

sub _hit_obj_for_row {
    my ($targ) = @_;
    if (my $fm = $targ->{DBMETA}{FM}) {
        if (my $sid = $targ->{HitID}) {
            if (my $sobj = $fm->fetch( $sid )) {
                return ($sobj);
            } else {
                if ($sid =~ /^(.+?)\, /) {
                    my $first = $1;
                    if (my $sobj = $fm->fetch( $first )) {
                        &msg_once("[-]","Using $first for multiple hits '$sid'");
                        return ($sobj);
                    }
                }
                return (undef, "Failed to recover HitID $sid");
            }
        } else {
            return (undef, "No HitID available");
        }
    } else {
        return (undef, "No FastaMeta object");
    }
}

sub _set_sorter {
    my ($targ) = @_;
    my $cbit = $targ->{GenomeFootprint} ||"";
    if ($cbit =~ /^([^\.]{1,2})\.\S+\:(\d+)[^\d]+(\d+)/) {
        my ($chr, $s, $e) = ($1, $2, $3);
        # $targ->{Start} = $s;
        # $targ->{End}   = $e;
        $chr = '0'.$chr if ($chr =~ /^\d$/);
        $cbit = sprintf("%2s %010d %010d", $chr, $s, $e);
    } else {
        # warn $targ->{GenomeFootprint}."\n";
    }
    my $sym = &_symbolSorter($targ->{Symbol}, $targ);
    $targ->{Sort} = sprintf("%-20s %d %-15s %s", $sym, 
                            $targ->{Mis}, uc($targ->{Build}), $cbit);
}

sub _symbolSorter {
    my ($sym, $targ) = @_;
    my $rv  = $sym;
    if (!$rv) {
        $rv = $fakeSym;
    } elsif ($rv =~ /(Rik|orf\d+)$/ ||
             $rv =~ /^(LOC|Gm|LINC|RGD)\d+$/) {
        # Sketchy genes
        $targ->{PredGene} = 1 if ($targ);
        $rv = substr($fakeSym,1) . $rv;
    } elsif ($rv =~ /^([^\d]+)(\d+)$/) {
        $rv = sprintf("%s%03d", $1, $2);
    } elsif ($rv =~ /^([^\d]+)(\d+)([^\d]+)$/) {
        $rv = sprintf("%s%03d%s", $1, $2, $3);
    }
    return uc($rv);
}

sub make_canvasXpress {
    $bnch->bench_start();
    my ($byLoc, $allBuild, $scanHits) = @_;
    my $locIndFile = &gene_file('index');
    my $inFh;
    if (open($inFh,">$locIndFile")) {
        print $inFh "<html><head><title>Gene index";
        print $inFh " for $projSht" if ($proj);
        print $inFh "</title><style>\n";
        print $inFh &html_styles( $allBuild );
        print $inFh "</style>\n";
        print $inFh "<link rel='shortcut icon' href='/biohtml/images/geneChart_16x16.png'>\n";
        print $inFh "</head><body>\n";
        print $inFh "<p class='note'>Each gene that was hit by at least one oligo is listed below. You can also view the <a href='../oligos/index.html'>Oligo Summary</a> directly.</p>\n";
        print $inFh "<table class='tab'><tbody>\n";
        print $inFh " <tr><th>".join('</th><th>', ('Locus Accession', 'Symbol', 'Build', 'Oligo Hits', 'Description'))."</th></tr>\n";
    } else {
        &msg("[!]","Failed to create oligo index file", $locIndFile, $!);
    }
    my @locDat = sort {
        $b->{sorter} <=> $a->{sorter}
        || ($a->{sym} || $fakeSym ) cmp ($b->{sym} || $fakeSym) 
        } values %{$byLoc};
    foreach my $ld (@locDat) {
        &_make_cx_file( $ld, $inFh, $allBuild, $byLoc, $scanHits );
    }
    if ($inFh) {
        print $inFh "</tbody></table>\n";
        print $inFh "</body></html>\n";
        close $inFh;
        &msg("[>]","Gene index : $locIndFile");
        if ($locIndFile =~ /(.+)\/html\/genes\/([^\/]+)$/) {
            # Make symlink
            my $sl = "$1/GeneIndex.html";
            # system("ln -sf \"$locIndFile\" \"$sl\"");
        }
    }
    $bnch->bench_end();
}

sub _make_cx_file {
    $bnch->bench_start();
    my ($dat, $inFh, $allBuild, $byLoc, $scanHits)= @_;
    my $ll   = $dat->{ll};
    my $sym  = $dat->{sym};
    my $desc = $dat->{desc} || "";

    # Organize the information by build
    my %builds;
    foreach my $row (@{$dat->{rows}}, @{$dat->{scanrows}}) {
        next if ($row->{SKIP});
        my $bld  = $row->{Build};
        my $hID  = $row->{HitID};
        my ($obj, $err) = &_hit_obj_for_row( $row );
        if ($obj) {
            if (my $par = $obj->val('par')) {
                if (my $pobj = $obj->db->fetch($par)) {
                    $obj = $pobj;
                } else {
                    &msg_once("Failed to recover parent '$par' for $hID");
                }
            }
            $hID = $obj->id();
        }
        my $targ = $builds{$bld}{$hID} ||= {
            obj  => $obj,
            rows => [],
            scan => {},
        };
        if ($row->{IS_SCAN}) {
            my $oid = $row->{OligoID};
            my $src = $row->{scanKey};
            my $off = ($obj && $obj->can('offset')) ? $obj->offset() : 0;
            my $fp  = $row->{GenomeFootprint};
            $targ->{scan}{$src}{$oid}{$fp}++;
        } else {
            push @{$targ->{rows}}, $row;
        }
    }
    my @bbs = sort keys %builds;
    my $bnum = $#bbs + 1;

    my @countBits;
    foreach my $mm (sort { $a <=> $b } keys %{$dat->{counts}}) {
        my $mc = !$mm ? 'ext' : $mm <= 2 ? "mm$mm" : "mm3";
        my $num = $dat->{counts}{$mm};
        my $tok = sprintf(" <span class='%s'>%d &times; %s</span>", 
                          $mc, $num, $mm ? "mm$mm" : "Perfect");
        push @countBits, $tok;
    }
    my $countTxt = join(', ', @countBits) || "?No counts found?";

    my @buildBits;
    foreach my $bld (@bbs) {
        my $btok = sprintf("<span class='%s'>[%s]</span>", $bld, $bld);
        push @buildBits, $btok;
    }
    my $buildTxt = join(' ', @buildBits) || "";
    
    my $file = &gene_file($ll);
    unless (open(CX, ">$file")) {
        &msg("[!]", "Failed to write CX file", $file, $!);
        printf($inFh "<tr><th>%s</th><td class='bigsym'>%s</td><td>%s</td><td>%s</td><td class='err'>%s</td></tr>\n", $ll, $sym || "", $buildTxt, $countTxt, "$desc<br />Failed to write details!") if ($inFh);
        $bnch->bench_end();
        return "";
    }
    if ($inFh) {
        my $lnk = sprintf("<a href='%s.html'>%s</a>", $ll, $ll);
        printf($inFh "<tr><th class='linkid'>%s</th><td class='bigsym'>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n", $lnk, $sym || "", $buildTxt, $countTxt, $desc)
    }
    print CX &_cx_header( $sym ? "$sym : $ll" : $ll, $allBuild );
    if (my $ov = $dat->{overlap}) {
        my @olls = sort keys %{$ov};
        print CX "<b>Perfectly overlaps:</b>";
        foreach my $ol (@olls) {
            my $od = $byLoc->{$ol};
            my $show = $od->{sym} || $ol;
            printf(CX " <a href='%s.html'>%s</a>", $ol, $show);
        }
        print CX "<br />\n";
    }
    print CX "<b>Normalized hit counts:</b> $countTxt<br />\n";
    my %mmMod;
    foreach my $bld (@bbs) {
        my $btok = sprintf("<span class='%s'>[%s]</span></h3>\n", $bld, $bld);
        print CX "<h3>Genome Build $btok</h3>\n";
        foreach my $hID (sort keys %{$builds{$bld}}) {
            my @rows = @{$builds{$bld}{$hID}{rows}};
            my $obj  = $builds{$bld}{$hID}{obj};
            my @errs;

            # Get hit description
            my %sdH  = map { ($_->{GeneDesc} || "") => 1 } @rows;
            delete $sdH{""};
            # my $desc = join(' // ', sort keys %sdH);

            my $olTrack = {
                name => "Oligo Matches",
                type      => "box",
                height    => 12,
                connect   => "true",
            };
            
            my $bestStr     = 1;
            my $coordMapSub = sub { return @_; };
            my ($gdna, $exTrack, $otherGenes, $rnaTrack);
            if (!$obj) {
                push @errs, "Failed to find metadata for $hID";
                $desc ||= "";
            } elsif ($obj->type() eq 'RNA') {
                $desc ||= $obj->desc() || "";
                if (my $v = $obj->variant()) {
                    $desc .= sprintf(" <span class='variant'>Variant %s</span>", $v);
                }
                $bestStr = 1;
                # die $obj->{JSON};
                if (my $len = $obj->length()) {
                    my $rname = &rna_sym_name( $ obj );
                    $rnaTrack = {
                        name      => "$rname : An RNA not aligned to genome",
                        type      => "box",
                        height    => 10,
                        connect   => "true",
                        data      => [{
                            hideName => 1,
                            dir      => 'right',
                            fill     => '#990033',
                            data     => [ [1, $len] ],
                        }],
                    };
                }
            } else {
                $desc    ||= $obj->desc() || "";
                my $gLen   = $obj->length();
                my @genes  = $obj->each_gene();
                my @rnaCls = ([],[]);

                # Find out which strand to use
                my (%strs, @offCoords);
                foreach my $gene (@genes) {
                    my @r = $gene->each_rna();
                    if ($gene->id() eq $ll) {
                        my $se = $gene->se();
                        push @offCoords, ($se->[0], $se->[1]);
                        # The gene of interest
                        push @{$rnaCls[0]}, @r;
                        foreach my $rna (@r) {
                            # die $args->branch($rna);
                            $strs{$rna->strand()}++;
                        }
                    } else {
                        # Another gene on the same segment
                        push @{$rnaCls[1]}, @r;
                    }
                }
                ($bestStr) = sort { $strs{$b} <=> $strs{$a} ||
                                        $b <=> $a } keys %strs;
                $bestStr ||= 1;
                if ($bestStr < 0) {
                    # Find the right-most genomic coordinate for our gene
                    my ($offset) = sort { $b <=> $a } @offCoords;
                    $offset++;
                    $coordMapSub = sub {
                        my @rv = map { $offset - $_ } @_;
                        return reverse @rv;
                    };
                } else {
                    # Find the left-most genomic coordinate for our gene
                    my ($offset) = sort { $a <=> $b } @offCoords;
                    $offset--;
                    $coordMapSub = sub {
                        my @rv = map { $_ - $offset } @_;
                        return @rv;
                    };
                }
                for my $ri (0..$#rnaCls) {
                    my @rnas = @{$rnaCls[$ri]};
                    next if ($#rnas == -1);
                    my $track;
                    my $trackData = {
                        type      => "box",
                        height    => 10,
                        connect   => "true",
                    };
                    if ($ri) {
                        $track = $otherGenes ||= {
                            name      => "Overlapping Genes",
                            type      => "box",
                            height    => 10,
                            connect   => "true",
                          };
                    } else {
                        $track =  $exTrack ||= {
                            name      => "Exonic Regions",
                            type      => "box",
                            height    => 10,
                            connect   => "true",
                        };
                    }
                    foreach my $rna (@rnas) {
                        my $acc = $rna->id();
                        my $str = $rna->strand() * $bestStr;
                        my $gene = {
                            name    => &rna_sym_name( $rna ),
                            dir     => $str < 0 ? 'left' : 'right',
                            outline => '#cc9900',
                            fill    => $ri ? '#999999' : '#aaffaa',
                        };
                        my @hsps = @{$rna->{HSP} || []};
                        my @hdata;
                        foreach my $hsp (@hsps) {
                            my ($ss, $se, $rs, $re) = @{$hsp};
                            my @h = &{$coordMapSub}( $ss, $se );
                            push @hdata, \@h;
                        }
                        @hdata = sort { $a->[0] <=> $b->[0] } @hdata;
                        if ($str < 0) {
                            # die $args->branch(\@hdata);
                            # @hdata = reverse @hdata;
                            # die $args->branch(\@hdata);
                        }
                        $gene->{data} = \@hdata;
                        push @{$track->{data}}, $gene;
                    }
                }
            }
                
            print CX "<h4>$hID";
            print CX " $rcTok" if ($bestStr < 0);
            print CX " - $desc" if ($desc);
            print CX "</h4>\n";

            # Check for scan data
            my (%scanTracks, %allScanPos);
            my $scanDat = $builds{$bld}{$hID}{scan} || {};
            while (my ($src, $oHash) = each %{$scanDat}) {
                my $shs = $scanHits->{Sources}{$src};
                unless ($shs) {
                    &msg("Failed to find scanning data for '$src'");
                    next;
                }
                while (my ($oid, $seHash) = each %{$oHash}) {
                    my $sho = $shs->{pos}{$oid};
                    unless ($sho) {
                        &msg("Failed to find scanning data for '$oid' in '$src'");
                        next;
                    }
                    my $len = $sho->{Len};
                    my $sfx = sprintf("%03d %s", $len, $src);
                    
                    my @poses;
                    foreach my $fp (keys %{$seHash}) {
                        my ($ft, $str, $gdata, $data) = &_map_footprint
                            ( $fp, $bestStr, $coordMapSub);
                        my $pos  = $data->[0][0];
                        my $gpos = $gdata->[0][0];
                        # &msg_once("[DEBUG]", "$fp => $pos vs ".join(",", map {join("..",@{$_}) } @{$data}));
                        push @poses, $pos;
                        if (my $prior = $allScanPos{$pos}) {
                        } else {
                            $allScanPos{$pos} = {
                                pos  => $pos,
                                gpos => $gpos,
                                sz   => {},
                            };
                        }
                        my $asp = $allScanPos{$pos};
                        my $aspsz = $asp->{sz}{$len} ||= {
                            Len => $len,
                        };
                        $aspsz->{fp} = $fp;
                        $aspsz->{sho} = $sho;
                    }
                    foreach my $row (@{$sho->{rows}}) {
                        my $bld = $row->{Build};
                        my $mm  = $row->{Mis};
                        my $fp  = $row->{GenomeFootprint};
                        my $scDat = $scanTracks{$bld}{$sfx} ||= {
                            Build  => $bld,
                            Source => $src,
                            Len    => $len,
                            pos    => {},
                        };
                        foreach my $pos (@poses) {
                            my $pDat = $scDat->{pos}{$pos} ||= {
                                pos => $pos,
                                mm  => {},
                            };
                            $pDat->{mm}{$mm}{$fp}++;
                        }
                    }
                }
            }
            my @asp = sort {$a <=> $b } keys %allScanPos;

            # Get exon coordianates and oligos
            my %olRngs;
            foreach my $row (@rows) {
                my $oid    = $row->{OligoID};
                unless ($row->{GenomeFootprint}) {
                    push @errs, "No footprint defined for $oid on $hID";
                    next;
                }
                # Remove leading object name
                my ($ft, $str, $gdata, $data) = &_map_footprint
                    ( $row->{GenomeFootprint}, $bestStr, $coordMapSub);
                my $dir    = $str < 0 ? 'left' : 'right';
                my $hkey   = join("\t", $ft, $str);
                my $rngDat = $olRngs{$hkey} ||= {
                    se  => [ $row->{Start}, $row->{End} ],
                    dir => $dir,
                    ft  => $ft,
                    data => $data,
                    gdat => $gdata,
                };
                my $oDat = $rngDat->{oligos}{$oid} ||= {
                    id => $oid,
                };
                $oDat->{mm}{ $row->{Mis} || 0 }++;
                if (my $oseq = $row->{OligoSeq}) {
                    $oDat->{seq}{uc($oseq)} ||= $oseq;
                }
                if (my $oHack = $row->{OligoDesc}) {
                    # Extract out CanvasXpress tags from description
                    while ($oHack =~ /(CX:([a-z]+)=\'([^\']+)\')/i) {
                        my ($r, $k, $v) = ($1, $2, $3);
                        $oDat->{tags}{$k} = $v;
                        $oHack =~ s/\Q$r\E/ /g;
                    }
                }
            }

            # Add the oligos:
            my $genomicPad  = 10;
            my @allData = map {[$_->{gpos} - $genomicPad,
                                $_->{gpos} + $genomicPad]} values %allScanPos;
            foreach my $rngDat (sort { $a->{se}[0] <=> $b->{se}[0] ||
                                     $a->{se}[1] <=> $b->{se}[1] }
                              values %olRngs) {
                my @data = @{$rngDat->{data}};
                foreach my $se (@{$rngDat->{gdat}}) {
                    my ($s, $e) = @{$se};
                    push @allData, [$s - $genomicPad, $e + $genomicPad];
                }
                my $se = $rngDat->{se};
                my $lr = $rngDat->{dir};
                foreach my $oDat (values %{$rngDat->{oligos}}) {
                    my @mms = sort { $a <=> $b } keys %{$oDat->{mm}};
                    my $oid = $oDat->{id};
                    push @errs, "Multiple mismatch values for $oid ".
                        join('-', @{$se}) if ($#mms > 0);
                    my $mm  = $mms[0] || 0;
                    my $col = !$mm ? '#00ff00': $mm == 1 ? '#ff9900':'#ff0000';
                    my @seqs = sort values %{$oDat->{seq}}; 
                    push @errs, "Multiple sequence values for $oid ".
                        join('-', @{$se}) if ($#seqs > 0);
                    my $cobj = {
                        data    => [ sort { $a <=> $b } @data ],
                        outline => $col,
                        fill    => $col,
                        dir     => $lr,
                        name    => $oid, #join(',', sort keys %{$info->{ids}}),
                        sequence => $seqs[0],
                        outline  => '#ffffff',
                    };
                    # Add description-based CX overrides:
                    while (my ($k, $v) = each %{$oDat->{tags} || {}}) {
                        $cobj->{$k} = $v;
                    }
                    push @{$olTrack->{data}}, $cobj;
                }
            }

            if ($obj) {
                my $ot = $obj->type();
                my @gSegs;
                # Add genomic segments
                @allData = sort { $a->[0] <=> $b->[0] ||
                                      $b->[1] <=> $a->[1]} @allData;
                if (my $seed = shift @allData) { push @gSegs, $seed; }
                foreach my $se (@allData) {
                    if ($gSegs[-1][1] + 1 >= $se->[0]) {
                        # Overlaps with last segment, extend end
                        $gSegs[-1][1] = $se->[1] if ($gSegs[-1][1] < $se->[1]);
                    } else { 
                        push @gSegs, $se;
                    }
                }
                foreach my $gse (@gSegs) {
                    my ($gs, $ge) = @{$gse};
                    if (my $seq = $obj->subseq_genomic( $gs, $ge)) {
                        # warn "$hID [$gs..$ge] = $seq\n";
                        # $seq = $obj->revcom($seq) if ($bestStr < 0);
                        $gdna ||= {
                            name      => $ot eq 'RNA' ? 
                                "RNA Sequence" : "Genomic Context",
                            type      => "box",
                            height    => 12,
                            connect   => "true",        
                        };
                        my ($s, $e) = &{$coordMapSub}( $gs, $ge );
                        my $seg = {
                            data     => [[$s, $e]],
                            hideName => 1,
                            dir      => $bestStr < 0 ? 'left' : 'right',
                            sequence => $seq,
                        };
                        push @{$gdna->{data}}, $seg;
                        
                    } else {
                        push @errs, "Failed to recover genomic sequence for $hID:$gs..$ge";
                    }
                }
            }
            

            # Add the genomic segments
            my $cxData = {};
            my $tracks = $cxData->{tracks} = [];
            foreach my $trk ($rnaTrack, $exTrack, $gdna, $olTrack) {
                push @{$tracks}, $trk if ($trk && $trk->{data});
            }

            # Scan plots
            my @scanCols = ( 'rgba(0,100,0,0.6)',
                             'rgba(255,153,0,0.6)',
                             'rgba(255,0,0,0.6)' );

            my $scanCommon = {
                type      => 'line',
                autowidth => 1,
                height    => 20,
                hideFeatureNames => 1,
                outline   => \@scanCols,
                # staggered => 0,
                featureStaggered => 1,
            };
            my @posClusters;
            foreach my $pos (@asp) {
                my $asp   = $allScanPos{$pos};
                next unless ($asp);
                if ($#posClusters == -1 || 
                    $posClusters[-1][-1] < $pos -1
                    # || $#{$posClusters[-1]} > 81 # CX issue with long clusters
                    ) {
                    # Start a new cluster
                    push @posClusters, [$pos];
                } else {
                    # Extend a cluster
                    push @{$posClusters[-1]}, $pos;
                }
            }
            
            if ($#posClusters != -1) {
                # die $args->branch(\@posClusters);
                my ($minLen) = sort { $a <=> $b } map { $_->{Len} } map { values %{$_} } values %scanTracks;
                
                # Add entropy plot
                my @datas;
                foreach my $pClus (@posClusters) {
                    my $data = {
                        offset => $pClus->[0],
                        width  => $#{$pClus} + 1,
                        data   => [[],[],[]],
                    };
                    push @datas, $data;
                    foreach my $pos (@{$pClus}) {
                        my $asp   = $allScanPos{$pos};
                        my $aspsz = $asp->{sz}{$minLen} || {};
                        my $sho   = $aspsz->{sho};
                        my $seq   = $sho->{OligoSeq};
                        for my $nt (1..3) {
                            my $ent = &_seq_entropy( $seq, $nt);
                            push @{$data->{data}[$nt-1]}, $ent;
                            $sho->{"Ent$nt"} = $ent;
                        }
                    }
                }
                my $track = {
                    %{$scanCommon},
                    name      => "Oligo Entropy",
                    names     => ['1bp', '2bp', '3bp' ],
                    bumpSpace => 0,
                    data      => \@datas,
                    setMaxY   => 1,
                    #featureStaggered => 1,
                };
                push @{$tracks}, $track;
            }

            my %doneScanBld;
            foreach my $bld (@{$allBuild}, keys %scanTracks) {
                next if ($doneScanBld{$bld}++);
                my $trkdat = $scanTracks{$bld};
                next unless ($trkdat);
                foreach my $sf (sort keys %{$trkdat}) {
                    my $scDat = $trkdat->{$sf};
                    my $len   = $scDat->{Len} || 1;
                    my $src   = $scDat->{Source};
                    my $sctn  = sprintf("Specificity Scan : %s vs. %s (%dbp)",
                                        $bld, $src, $len );
                    my @datas;
                    foreach my $pClus (@posClusters) {
                        my $data = {
                            offset => $pClus->[0],
                            width  => $#{$pClus} + 1,
                            data   => [[],[],[]],
                        };
                        push @datas, $data;
                        foreach my $pos (@{$pClus}) {
                            my $asp   = $allScanPos{$pos};
                            my $pDat  = $scDat->{pos}{$pos} || {};
                            my $aspsz = $asp->{sz}{$len} || {};
                            my $sho   = $aspsz->{sho};
                            my $exRow =  {
                                Source          => $src,
                                Symbol          => $sym,
                                Build           => $bld,
                                GenomeFootprint => $aspsz->{fp},
                            };
                            map { $exRow->{$_} = $sho->{$_} } 
                            qw(Start End OligoSeq OligoID Ent1 Ent2 Ent3);
                            push @{$scanHits->{Excel}}, $exRow;
                            for my $mm (0..2) {
                                # Powers of three because for mismatches each site
                                # has only three variants (one is a match)
                                my $mod = $mmMod{$len}{$mm} ||= !$mm ? 1 : 
                                    (3 ** $mm) * &nchoosek($len, $mm);

                                # Number of unique locations
                                my @uniq = keys %{$pDat->{mm}{$mm} || {}};
                                my $num  = $#uniq + 1;
                                my $val  = int(0.5 + 1000 * $num / $mod) / 1000;
                                push @{$data->{data}[$mm]}, $val;
                                $exRow->{$mm."MM"} = $num if ($num);
                            }
                        }
                    }
                    my $track = {
                        %{$scanCommon},
                        name    => $sctn,
                        names   => ['0MM', '1MM', '2MM' ],
                        data    => \@datas,
                        setMaxY => 3.5,
                    };
                    push @{$tracks}, $track;
                }
            }
            
            # Add related genes
            push @{$tracks}, $otherGenes if ($otherGenes);


            my $cid = "$bld-$hID";
            print CX &_cx_script($cxData, $cid);
            print CX "<span class='err'>ERRORS: ".join(' / ', @errs).
                "</span><br />\n" unless ($#errs == -1);
        }
    }
    my $mismatchMod = "";
    foreach my $bp (sort { $a <=> $b } keys %mmMod) {
        foreach my $mm (sort { $a <=> $b } keys %{$mmMod{$bp}}) {
            next unless ($mm);
            $mismatchMod .= sprintf
                ("%d bp, %d mismatch%s: <b>1/%d</b><br />\n",
                 $bp, $mm, $mm == 1 ? "" : "es", $mmMod{$bp}{$mm});
        }
    }
    print CX "<h4>Specificity Scan Scaling Factors</h4>\n<i>To accomodate higher hit counts for higher mismatches, the specificity plots have been scaled by the following values</i><br />$mismatchMod" if ($mismatchMod);
    print CX "&larr;<a href='index.html'>Return to Summary</a><br />\n";
    print CX "</body></html>\n";
    $bnch->bench_end();
    return {
        path => $file,
    };
}

sub _map_footprint {
    my ($ft, $bestStr, $coordMapSub) = @_;
    $ft =~ s/^.+\://;
    my $str    = 1;
    if ($ft =~ /(.+)\[([\-\+]*1)\]$/) {
        ($ft, $str) = ($1, $2 + 0);
    }
    $str *= $bestStr;
    my (@gdata, @data);
    foreach my $fb (split(',', $ft)) {
        my ($s, $e);
        if ($fb =~ /^\d+$/) {
            $s = $e = $fb;
        } elsif ($fb =~ /^(\d+)\.\.(\d+)$/) {
            ($s, $e) = ($1, $2);
        } else {
            next;
        }
        push @gdata, [$s, $e];
        if ($coordMapSub) {
            ($s, $e) = &{$coordMapSub}( $s, $e );
            push @data, [$s, $e];
        }
    }
    return ($ft, $str, \@gdata, \@data);
}

sub rna_sym_name {
    my $obj = shift;
    return "" unless ($obj);
    my @dbits;
    if (my $sym = $obj->sym()) {
        push @dbits, $sym;
    }
    if (my $var = $obj->variant()) {
        push @dbits, "Var.$var";
    }
    push @dbits, $obj->id();
    return join(" ", @dbits);
}

sub _cx_script {
    my ($cxData, $cid) = @_;
    my $cxSer = $ser->obj_to_json($cxData, 0, { basicArray => 1 });
    my $w      = 1000;
    my $h      = 600;
    my $rv = "<canvas id='$cid' width='$w' height='$h'></canvas>\n";
    $rv .= <<EOF;
<script>
new CanvasXpress('$cid', $cxSer,
 { "autoExtend": "true",
   "debug": null,
   "featureNameFontColor": "rgb(29,34,43)",
   "featureStaggered": "true",
   "filterSkipNullKeys": 1,
   "graphType": "Genome",
   "infoTimeOut": 300000,
   "margin": 2,
   "sequenceFill": "#cccccc",
   "trackNameFontColor": "rgb(256,100,0)",
   "trackNameFontSize": 10,
   "wireColor": "rgba(29,34,43,0.1)",
   "subtracksMaxDefault": 1000,
   "xAxisTickColor": "rgb(29,34,43)" });
</script>
<div style='clear:both'></div>
EOF

    return $rv;
}

sub _cx_file_path {
    my ($dir, $id) = @_;
    return "$dir/$id.html";
}

sub _cx_header {
    my ($lab, $allBuild) = @_;
    my $cxHead = &cx_header();
    my $sty    = &html_styles( $allBuild );
    my $rv = <<EOF;
<html>
 <head>
  <link rel="shortcut icon" href="/biohtml/images/geneDetail_16x16.png">
  <title>$lab Oligo Report</title>
  <style>
$sty
  </style>
$cxHead
 </head><body>
<h2>$lab</h2>
EOF

    if ($isBeta) {
    }
    return $rv;
}

our $cxBoilerHeader;
sub cx_header {
    unless ($cxBoilerHeader) {
        my $site  = "http://xpress.pri.bms.com";
        my $debug = $args->val(qw(debug debugcx cxdebug));
        $debug = 1 unless (defined $debug);
        my $url   = "$site/JAVASCRIPT/canvas/js/canvasXpress.min.js";
        $cxBoilerHeader = 
            '  <meta http-equiv="X-UA-Compatible" content="chrome=1">'."\n";
        if ($debug) {
            $url =~ s/\.min\./\.debug\./;
            $cxBoilerHeader .= "  <script type='text/javascript' src='http://xpress.pri.bms.com/JAVASCRIPT/canvas/js/canvasXpress.public.min.js'></script>\n";
        }
        my $ieKludge = $url;
        $ieKludge =~ s/\/[^\/]+$//;
        $cxBoilerHeader .= <<EOF;
  <!--[if IE]>
    <script type='text/javascript' src='$ieKludge/flashcanvas.js'></script>
  <![endif]-->
  <script type='text/javascript' src='$url'></script>
EOF

    }
    return $cxBoilerHeader;
}

sub query_summary_sheet {
    $bnch->bench_start();
    my ($eh, $allRows, $dbByBuild) = @_;
    # die $args->branch($dbByBuild);
    my @head    = ("BMS Number", "Species", "TotSymbols", "TotLoci");
    my $hOff    = $#head + 1;
    my $cntWid  = 5;
    my $geneWid = 20;
    my @wid     = map { $cntWid } @head;
    $wid[0] = 15;
    $wid[1] = 20;

    my %byQuery;
    my $maxMis = 0;
    foreach my $row (@{$allRows}) {
        next if ($row->{SKIP});
        my $oid = $row->{OligoID};
        my $bld = $row->{Build} || $unBuild;
        my $db  = $dbByBuild->{$bld} || {};
        my $tax = $db->{taxa} || $unTaxa;
        my $dat = $byQuery{$oid}{$tax} ||= {
            'BMS Number' => $oid,
            'Species'    => $tax,
            OligoSeq     => $row->{OligoSeq},
            Len          => CORE::length($row->{OligoSeq}),
            Loci         => {},
            Syms         => {},
        };
        my $mm   = $row->{Mis};
        $maxMis  = $mm if ($maxMis < $mm);
        if (my $ll = $row->{GeneID}) {
            my $ldat = $dat->{Loci}{$ll} ||= {
                id => $ll,
                mm => 999,
                ex => {},
            };
            my $ex     = $row->{Exonic} || "";
            my $exFlag = ($ex eq 'Intron') ? "Int" : "Ex";
            
            $ldat->{mm} = $mm if ($ldat->{mm} > $mm);
            $ldat->{ex}{$exFlag} = $mm if (!defined $ldat->{ex}{$exFlag} ||
                                           $ldat->{ex}{$exFlag} > $mm);
            my $sym = $row->{Symbol};
            if (!$ldat->{sym} && $sym ne $ll) {
                $ldat->{sym} = $sym;
                $dat->{Syms}{uc($sym)} ||= $sym;
            }
        }
    }

    my @useHead = @head;
    my %hitTokMap;
    for my $mm (0..$maxMis) {
        foreach my $isEx (1, 0) {
            my $tok     = ($isEx ? 'Ex' : 'Int').$mm;
            my $part    = $isEx ? "Exon" : "Intron";
            my $geneCol = $hitTokMap{$tok} = sprintf
                ("%s %d Genes", $part, $mm);
            $formatters->{$tok} = sub {
                my $v = shift;
                return $v ? 'cenbold' : 'cengray';
            };
            $formatters->{$geneCol} = sub {
                my $v = shift;
                return $v ? 'symblock' : '';
            };
            push @useHead, ($tok, $geneCol);
            push @wid, ($cntWid, $geneWid);
            my $niceMM = 'mismatch' . ($mm == 1 ? '' : 'es');
            $excelMeta->{$tok} = [ $cntWid, "Count of distinct genes where the oligo aligned in an $part with $mm $niceMM"];
            $excelMeta->{$geneCol} = [ $cntWid, "A list of genes where the oligo aligned in an $part with $mm $niceMM. Symbols are used where available; if no symbol is provided, the locus ID is used instead."];

        }
    }

    my @oids = sort keys %byQuery;
    if ($#oids == -1) {
        $bnch->bench_end();
        return;
    }

    foreach my $ocol (qw(Len OligoSeq)) {
        push @useHead, $ocol;
        push @wid, $excelMeta->{$ocol}[0];
    }
    my %lu = map { $useHead[$_] => $_ } (0..$#useHead);
    my $sname = "Oligo Summary";
    $eh->sheet( -name    => $sname,
                -freeze  => [1,2],
                -cols    => \@useHead,
                -width   => \@wid);

    # Exon #mm | Exon # Genes | Intron #mm | Intron # Genes
    foreach my $oid (@oids) {
        foreach my $taxa (sort keys %{$byQuery{$oid}}) {
            my $dat  = $byQuery{$oid}{$taxa};
            my %sSort;
            foreach my $sym (values %{$dat->{Syms}}) {
                my $srt = &_symbolSorter( $sym );
                $sSort{$srt} ||= $sym;
            }
            my @uSym = map { $sSort{$_} } sort keys %sSort;
            $dat->{TotSymbols} = ($#uSym + 1) || "";
            my @uLoc = values %{$dat->{Loci}};
            $dat->{TotLoci} = ($#uLoc + 1) || "";
            map { delete $dat->{$_} } qw(Syms Loci);
            my %hitCols;
            foreach my $ldat (@uLoc) {
                while (my ($tok, $mm) = each %{$ldat->{ex}}) {
                    push @{$hitCols{$tok.$mm}}, $ldat->{sym} || $ldat->{id};
                }
                #my $isEx = $ldat->{ex};
                #my $mm   = $ldat->{mm};
                #my $tok  = ($isEx ? 'Ex' : 'Int').$mm;
                #push @{$hitCols{$tok}}, $ldat->{sym} || $ldat->{id};
            }
            while (my ($tok, $arr) = each %hitCols) {
                $dat->{$tok} = $#{$arr} + 1;
                my %sSort;
                foreach my $sym (@{$arr}) {
                    my $srt = &_symbolSorter( $sym );
                    $sSort{$srt} ||= $sym;
                }
                my @uSym = map { $sSort{$_} } sort keys %sSort;
                $dat->{ $hitTokMap{$tok} } = join(' ', @uSym);
            }
            my (@row, @fmt);
            while (my ($cn, $ind) = each %lu) {
                my $val = $dat->{$cn};
                if ($hitTokMap{$cn}) {
                    $val ||= 0;
                } elsif ($cn eq 'BMS Number') {
                    if (my $url = &oligo_url($val)) {
                        $val = ['url', $url, $val];
                    }
                }
                $row[$ind] = $val;
                if (my $cb = $formatters->{$cn}) {
                   $fmt[$ind] = &{$cb}($val);
               }
            }
            $eh->add_row_explicit($sname, \@row, \@fmt);
        }
        $eh->blank_row($sname);
    }
    $bnch->bench_end();
}

sub help_sheet {
    $bnch->bench_start();
    my ($eh, $data) = @_;
    my $sname = "Help and Glossary";
    $eh->sheet( -name    => $sname,
                -width   => [30, 60, 60, 60] );

    $eh->add_row($sname, ["If a row is too long to read, sometimes double-clicking the cell will wrap it for you."]);
    $eh->blank_row($sname);

    my @cols = sort keys %{$usedCols};
    unless ($#cols == -1) {
        &_help_header_row($eh, $sname, "Column Glossary");
        $eh->add_row($sname, ["Column","Description"], ['bold', 'bold']);
        foreach my $cn (@cols) {
            next unless ($cn);
            my $uCn = $cn;
            $uCn =~ s/^\Q$metaPrfx\E//;
            my $desc = $usePsh->get_metadata('columns', $cn, 'desc');
            $eh->add_row($sname, [$uCn, $desc], ['bold']);
        }
        $eh->blank_row($sname);
    }

    my $htxt = &_help_boilerplate();
    foreach my $line (split(/[\n\r]/, $htxt)) {
        next if ($line =~ /^\#/);
        if ($line =~ /^\s*$/) {
            $eh->blank_row($sname);
        } elsif ($line =~ /^!(.+)/) {
            &_help_header_row($eh, $sname, $1);
        } elsif ($line =~ /^(.+)\s*::\s*(.+)/) {
            my ($c1, $c2) = ($1, $2);
            my $fmt = 'bold';
            if ($c2 =~ /^(.+)\s+\[([^\]]+)\]\s*$/) {
                ($c2, $fmt) = ($1, $2);
            }
            $eh->add_row($sname, [$c1, $c2],  [$fmt]);
        } else {
            $eh->add_row($sname, [$line]);
        }
    }

    $bnch->bench_end();
}

sub _help_boilerplate {
    my @hc = map { 4 ** (($_ + 1) * $lodBand) } (0..2);
    my $htxt = <<EOF;
!Worksheet Overview
Hits :: This is the most detailed worksheet. Each row is a hit of one oligo against one database entry
# Oligo Summary :: Each row is one oligo in one species. Tallies up the genes hit by exonic/intronic, and how many mismatches were present
# Symbols :: Each row is one gene symbol. Columns represent the different genome builds searched, and report the specific oligos that hit the gene. "mm2" means "two mismatches". Note that genes may not always share the same symbol across species.
# Found Seqs :: Each column is a query oligo. The rows show specific sequences found across all databases, with [##] showing how many times that sequence was seen.
# Query Oligos :: Basic information (ID, sequence, description) for the oligos that were searched.
DBs :: Basic information (file name, build, species, total number of observed hits) for the databases that were searched.

!About pre-mRNAs
The pre-mRNAs are 'cookie cutter' extracted from the genome from the boundaries of RNA-to-genome alignments (Sim4 algorithm)
All the overlapping RNA splice variants from the same loci are clustered together into one pre-mRNA.
The furthest "left" coordinate and furthest right coordinate from all the variants are used to "cut out" the pre-mRNA from the genome.
If the RNA aligned to the -1 strand of the reference genome, then the sequence is reversed complemented. This way all the pre-mRNAs will be "in the right direction" (5' to 3').
If two loci overlap and have the EXACT SAME left and right coordinates then they will be combined into one pre-mRNA. These are almost always highly-conserved family members.
Only the best genome alignment is used for each locus. However, in some cases a locus has two or more genomic locations with identical top scores.
In these cases, each of the top-scored genomic locations will be extracted as a separate pre-mRNA. These cases are also almost always highly conserved gene families.

!About Hit IDs
7.GRCm38:81347027..81454758.R :: This is a pre-mRNA on Chromosome 7 in genome build GRCm38, taken from coordinates 81347027-81454758. There is at least one gene on the reverse (R) strand. Segments that have genes on both strand will be suffixed 'FR'
4.GRCm38:47472458^47473939.F :: A '^' character indicates an discontinuity in the genome. In this case, the  entry is representing the spliced Exon-Exon boundary, having removed the intron from 47472459 to 47473938. These are little stub entries needed to allow matches to oligos spanning an mRNA splice junction.
NM_001038609.2 :: This is an mRNA (non-genomic). NM transcripts are well-supported protein coding
XM_006725265.1 :: Another protein-coding mRNA, but XM transcripts are more speculative (more likely to be false)
NR_024583.1 :: This is a non-coding, but well-supported RNA
XR_248730.2 :: Also a non-coding RNA, but XR transcripts are speculative

!About Genome Builds
A 'build' is a specific reference assembly of a genome, usually represented by a short 'token'.
The build token abbreviates the species, and indicates the version number
As a species is refined and finished, the difference between subsequent versions is generally small.
Early genome builds can have significant missing data and assembly errors.
Human builds are GRCh##, with 37 being the current (2014) build. It is exceptionally good. Older human builds were NCBI##
RSG is a special build representing RefSeqGene, a public gDNA/pre-mRNA data set put together by the NCBI.
Mouse builds are GRCm##, 38 most current. GRCm38 is a decent assembly, but still has some significant errors. Older mouse builds were NCBIM##
Rat is Rnor_##, Rnor_5 being current. It is not a terribly good assembly
Dog is CanFam##, CanFam3 being current. Also not a great assembly.
Cyno has been assembled by both BGI (BgiCE = Crab Eating) and Wash U (MacFas5). The Wash U build is superior but still has major holes.
Rhesus has two BGI assemblies, BgiIR (Indian Rhesus) and BgiCR (Chinese Rhesus). These assemblies have lots of holes.

!About Symbols and Accessions
Symbols (eg CAV3) are friendly terms to use when discussing genes, are easy to remember, and often suggest basic biological function.
Accessions (eg LOC859, 3.GRCh37:142720366-142779567+1) are cold and ugly, and convey no human information by themselves.
However, symbols are notoriously ill-defined, very hard for software to properly recognize, and may even change over time.
Whenever possible, when recording information in your notebooks, please be sure to include an accession along with the symbol.
This will make future analysis MUCH easier, and reduce the chance of error.

!Various table markup
Intron :: The oligo is predicted to only fall in introns [intron]
Exon :: The oligo covers an exon, at least in one variant [exon]
ExonExon :: The oligo spans an exon-exon boundary [exonexon]
FullSpan :: The oligo COMPLETELY covers an exon - presumably this is a tiny exon or a huge oligo [fullspan]
SpliceJunction :: The oligo covers an exon-intron boundary [splicejunction]
RNA :: The oligo hits an mRNA (not a pre-mRNA). It is exonic, but we can not tell if it crosses an Exon-Exon boundary or not [rna]
CrypticExon :: The oligo is within exonic sequence that does not appear to be represented on the reference genome build. That is, there is RNA sequence that appears to be absent in the gDNA. [crypticexon]
FirstExon :: The oligo is MOSTLY in the first exon, but a little bit hangs off the 5' end [firstexon]
LastExon :: The oligo is MOSTLY in the last exon, but a little bit hangs off the 3' end [lastexon]

0 :: No mismatches between the oligo and the target [MM0]
1 :: One mismatch between the oligo and the target [MM1]
2 :: Two mismatches between the oligo and the target [MM2]
4 :: More than two mismatches [MM3]

0 :: The target sequence is unambiguous (all A,C,G or T) [Mute]
1 :: There are one or two ambiguous characters (eg R, Y, N, etc) in the target, likely indicating polymorphisms [AmLow]
4 :: There are three or four ambiguities in the target [AmMed]
6 :: Five or more ambiguities are present in the target. This could be a hyper-polymorphic locus, or possibly a region of poor assembly. There is a good chance that your query does not 'really' match this region [AmHi]

Homo sapiens :: A formal scientific species name [species]
# CXXC1P1 :: Pseudogenes are highlighted in pink [pseudosym]

#127 :: Genome hit count over $hc[2] times expected, based on oligo size, number of mismatches, and search space size [Expect2]
#7 :: Genome hit count over $hc[1] times expected [Expect1]
#3 :: Genome hit count within 1/$hc[0] and $hc[0] times expected [Expect0]
#2 :: Genome hit count below 1/$hc[1] times expected [Expect-1]
#1 :: Genome hit count below 1/$hc[2] times expected [Expect-2]

EOF

    return $htxt;

}

sub _help_header_row {
    my ($eh, $sname, $head) = @_;
    my @row = map { "" } (1..10);
    my @fmt = map { 'exhead' } @row;
    $row[0] = $head;
    $eh->add_row($sname, \@row, \@fmt);    
}

sub found_seqs_sheet {
    $bnch->bench_start();
    my ($eh, $allRows, $oSeqs, $olFh, $allBld, $bldLU) = @_;
    my %byOligo;
    foreach my $row (@{$allRows}) {
        next if ($row->{SKIP});
        my $oid = $row->{OligoID};
        my $dat = $byOligo{$oid} ||= {
            id   => $oid,
            seq  => $row->{OligoSeq},
            desc => $row->{OligoDesc},
            hits => {},
        };
        my $str  = $row->{Str} || 0;
        my $sdat = $dat->{hits}{$str} ||= {};
        my $hseq = $row->{HitSeq};
        my $hdat = $sdat->{$hseq} ||= {
            seq => $hseq,
            mm  => $row->{Mis},
            bld => {},
        };
        my $bld   = $row->{Build} || $unBuild;
        my $bdat  = $hdat->{bld}{$bld} ||= { };
        my $genom = $row->{GenomeFootprint};

        if (my $ll = $row->{GeneID}) {
            my $ldat = $bdat->{$ll} ||= {
                id  => $ll,
                gen => {},
            };
            if (my $ex = $row->{Exonic}) {
                $ldat->{ex} = 1 unless ($ex eq 'Intron');
            }
            if ($genom) {
                # Genomic location known
                $ldat->{gen}{$genom} = 1;
            }
            my $sym = $row->{Symbol};
            if (!$ldat->{sym} && $sym ne $ll) {
                $ldat->{sym} = $sym;
                $dat->{Syms}{uc($sym)} ||= $sym;
            }
        }
    }
    foreach my $odat (values %{$oSeqs}) {
        # Make sure any non-hit oligos are also represented
        my $oid = $odat->{id};
        $byOligo{$oid} ||= {
            id   => $oid,
            seq  => $odat->{seq},
            desc => $odat->{desc},
        };
    }
    my @os = sort keys %byOligo;
    if ($#os == -1) {
        $bnch->bench_end();
        return;
    }
    my $oname = "Found Seqs";
    my $ohead = \@os;
    $eh->sheet( -name    => $oname,
                -freeze  => [2,0],
                -cols    => \@os,
                -width   => [map { 25 } @os ]);
    $eh->add_row($oname, [map { $oSeqs->{$_}{seq} } @os],
                 [ map {'mono green' } @os]  );
    
    my @oligoTable;
    my $revcom = $args->val('revcomhtml');
    my $spanFmt = "<span class='%s'>%s</span>";
    my %legend;

    my $sty = &html_styles( $allBld );
    for my $o (0..$#os) {
        my $oid  = $os[$o];
        my $oidH = $ser->esc_xml($oid);
        my $dat  = $byOligo{$oid};
        # Populate the excel workbook
        my (%seenBuild, %seqs);
        while (my ($str, $sdat) = each %{$dat->{hits}}) {
            while (my ($hseq, $hdat) = each %{$sdat}) {
                my $mm = $hdat->{mm};
                my $exDat = $seqs{$hseq} ||= {
                    seq => $hseq,
                    num => 0,
                    mm  => $mm,
                };
                while (my ($bld, $bdat) = each %{$hdat->{bld}}) {
                    # Find number of distinct genomic loci
                    my %gLocs;
                    while (my ($ll, $ldat) = each %{$bdat}) {
                        map { $gLocs{ $_ } = 1 } keys %{$ldat->{gen}};
                    }
                    my @uLoc = keys %gLocs;
                    # If no genomic loci exist (mRNA only) then still count 1:
                    my $locNum = ($#uLoc + 1) || 1;
                    $exDat->{num} += $locNum;
                    $seenBuild{$bld}{$mm} += $locNum;
                }
                # Bubble up the summed count to the hit sequence level:
                $hdat->{num} = $exDat->{num};
            }
        }
        my @excelHits = sort {
            $b->{num} <=> $a->{num} || $a->{mm} <=> $b->{mm} || $a->{seq} cmp $b->{seq} } values %seqs;
        for my $h (0..$#excelHits) {
            my $rf    = $oligoTable[$h] ||= [ [], [] ];
            my $exDat = $excelHits[$h];
            my $mm    = $exDat->{mm};
            $rf->[0][$o] = sprintf("%s [%d]", $exDat->{seq}, $exDat->{num});
            my $fmt = "mono";
            $fmt .= " ". &{$misPick}($mm) if ($mm) ;
            $rf->[1][$o] = $fmt;
        }
        
        my $hFile = &oligo_file($oid);
        unless (open(OHTML, ">$hFile")) {
            &msg("[!!]","Failed to generate oligo HTML", $hFile, $!);
            next;
        }
        print OHTML "<html><head><title>$oidH Oligo Summary</title>\n";
        print OHTML "<link rel='shortcut icon' href='/biohtml/images/oligoDetail_16x16.png'>\n";
        print OHTML "<style>\n$sty</style>\n";
        print OHTML "</head><body>\n";

        my $fa = "";
        if (my $desc = $dat->{desc}) {
            $fa .= " <span class='desc'>$desc</desc>";
        }
        $fa .= "<br />\n";
        my $qSeq = $dat->{seq};
        $fa .= sprintf($spanFmt, 'seq', $qSeq);
        $fa .= "<br />\n";
        printf( OHTML "&gt;<span class='bigid'>%s</span>%s",
                $oidH, $fa);
        print OHTML "<h3>Oligo match sequences as extracted from the search databases</h3>\n";
        print OHTML "<i>Sequences have been reverse complemented</i><br />\n"
            if ($revcom);
        my @blds = $allBld ? @{$allBld} : sort keys %seenBuild;
        my $fullCnt = 0;
        foreach my $bld (@blds) {
            my $cnts = $seenBuild{$bld};
            next unless ($cnts);
            printf(OHTML $spanFmt, $bld, "[$bld]");
            if (my $bd = $bldLU->{$bld}) {
                if (my $tx = $bd->{taxa}) {
                    printf(OHTML $spanFmt, 'tax', " $tx");
                }
                if (my $path = $bd->{file}) {
                    printf(OHTML $spanFmt, 'file', " $path");
                }
            }
            print OHTML "<br /> &rarr; \n";
            foreach my $mm (sort { $a <=> $b } keys %{$cnts}) {
                my $mc = !$mm ? 'ext' : $mm <= 2 ? "mm$mm" : "mm3";
                printf(OHTML $spanFmt, $mc, $mm ? "mm x$mm" : 'Perfect');
                my $num = $cnts->{$mm};
                $fullCnt += $num;
                printf(OHTML " = %d Hit%s ", $num, $num == 1 ? '' : 's' );
            }
            print OHTML "<br />\n";
        }
        if ($olFh) {
            my $lnk = sprintf("<b><a href='%s.html'>%d hit%s</a></b>",
                              $oid, $fullCnt, $fullCnt == 1 ? '' : 's');
            print $olFh "&gt;<span class='normid'>$oidH</span> $lnk$fa";
        }
        $qSeq    = $su->revcom($qSeq) if ($revcom);
        my $pad  = " " x(CORE::length($qSeq));

        print OHTML "<pre>\n";
        printf(OHTML "<span class='ol'>%s</span> ",$qSeq);
        printf(OHTML "<span class='fwd'>%s &rarr;</span>\n", 
               $revcom ? 'RNA' : 'Oligo');
        my @strs = sort { $a <=> $b } keys %{$dat->{hits}};
        if ($#strs == -1) {
            printf(OHTML $spanFmt, 'nul', "No hits reported");
        }
        for my $si (0..$#strs) {
            print OHTML "\n" if ($si);
            my $str = $strs[$si];
            my ($rcs, $arr, $sty) = ($qSeq, '&rarr;', 'fwd');
            if ($str < 0) {
                $rcs = $su->revcom($rcs);
                ($arr, $sty) = ('&larr;', 'revcom');
            }
            printf(OHTML "<span class='%s'>%s</span>%s ".
                   "<span class='%s'>%s %s : Hits to %s strand</span>\n",
                   $sty, $rcs, $sty eq 'revcom' ? ($pad, 'rev') : ('', 'fwd'), 
                   $arr, $revcom ? 'Oligo' : 'RNA',
                   $str < 0 ? 'reverse' : 'forward');
            my @hits = sort { $b->{num} <=> 
                                  $a->{num} } values %{$dat->{hits}{$str}};
            foreach my $hdat (@hits) {
                my $gseq = $hdat->{seq};
                $gseq    = $su->revcom($gseq) if ($revcom);
                if (uc($gseq) eq uc($qSeq)) {
                    print OHTML "<span class='ext'>$gseq</span>";
                } else {
                    my $fa = $gseq;
                    while ($fa) {
                        if ($fa =~ /^([A-Z]+)(.*)/) {
                            $fa = $2;
                            print OHTML "<span class='mat'>$1</span>";
                        } elsif ($fa =~ /^([a-z]+)(.*)/) {
                            $fa = $2;
                            print OHTML "<span class='mis'>$1</span>";
                        } else {
                            print OHTML "<span class='err'>$fa</span>";
                            $fa = "";
                        }
                    }
                }
                my $num  = $hdat->{num};
                my $mm   = $hdat->{mm};
                printf(OHTML " [%3dx]", $num, $num == 1 ? '' : 's');
                if ($mm) {
                    my $mc = $mm <= 2 ? "mm$mm" : "mm3";
                    printf(OHTML " <span class='$mc'>mismatch x%d</span>", $mm);
                } else {
                    print OHTML  "  <i class='ext'> Perfect </i> ";
                }
                my @blds = $allBld ? @{$allBld} : sort keys %{$hdat->{bld}};
                foreach my $bld (@blds) {
                    my $bdat = $hdat->{bld}{$bld};
                    next unless ($bdat);
                    print OHTML " ";
                    printf(OHTML $spanFmt, $bld, "[$bld]");
                    my @ldats = sort { ($a->{sym} || $fakeSym) cmp ($b->{sym} || $fakeSym) || $a->{id} cmp $b->{id} } values %{$bdat};
                    foreach my $ldat (@ldats) {
                        my $ll   = $ldat->{id};
                        my $show = $ldat->{sym} || $ll;
                        printf(OHTML " <a href='../genes/%s.html'>%s</a>",
                               $ll, $show);
                        my @gHits = keys %{$ldat->{gen}};
                        if ($#gHits == -1) {
                            # No genomic locations
                            if ($ldat->{aligned}) {
                                # But the RNAs were aligned
                                # Implies that this is an exon-exon boundary
                                my $tok = sprintf
                                    ($spanFmt, 'exonexon', "Exon-Exon?");
                                print OHTML " $tok";
                                $legend{$tok} ||= "The hit may cross an exon-exon boundary: the hit is not seen in the pre-mRNA, but can be found in a genome-aligned mRNA.";
                            } else {
                                my $tok = sprintf
                                    ($spanFmt, 'rnaonly', "RNA only");
                                print OHTML " $tok";
                                $legend{$tok} ||= "The hit was only observed in a mRNA (not pre-mRNA). This is likely because genomic alignment for the locus is not available.";
                            }
                        } elsif ($#gHits > 0) {
                            $legend{"<span class='multihit'># Hits</span>"} ||= "The oligo hits this locus in two or more locations.";
                            printf(OHTML " <span class='multihit'>%d Hits</span>", $#gHits + 1);
                        }
                    }
                }
                print OHTML "\n";
            }
        }
        print OHTML "</pre>\n";
        my @leg = sort keys %legend;
        unless ($#leg == -1) {
            print OHTML "<h3>Legend</h3>\n";
            foreach my $tok (@leg) {
                print OHTML "$tok : <i>$legend{$tok}</i><br />\n";
            }
        }
        print OHTML "<p><a href='index.html'>Return to oligo index</a></p>\n";
        print OHTML "</body></html>\n";
        close OHTML;
    }
    foreach my $rf (@oligoTable) {
        $eh->add_row($oname, @{$rf});
    }
    $bnch->bench_end();
}

sub _backfill_oligo {
    my ($row) = shift;
    return {
        id   => $row->{OligoID},
        seq  => $row->{OligoSeq},
        desc => $row->{OligoDesc},
    };
}

sub query_oligo_sheet {
    my ($eh, $allRows, $allBuild, $oligos) = @_;
    my @qdats = sort { $a->{id} cmp $b->{id} } values %{$oligos};
    return if ($#qdats == -1);
    $bnch->bench_start();
    my $sname = "Query Oligos";
    my @head  = ("OligoID", "Len", "OligoSeq", "Ent1", "Ent2", "Ent3", "OligoDesc");
    my @wids  = map { $excelMeta->{$_} ? $excelMeta->{$_}[0] : 10 } @head;
    my $maxMis = 0;
    my $bfmt   = "%s+%d";
    foreach my $row (@{$allRows}) {
        my $oid = $row->{OligoID};
        my $qdat = $oligos->{$oid} ||= &_backfill_oligo( $row );
        my $mm   = $row->{Mis};
        $maxMis  = $mm if ($maxMis < $mm);
        my $tally = sprintf($bfmt, $row->{Build}, $mm);
        $qdat->{tally}{$tally}++;
        
    }
    my @fBase = map { $excelMeta->{$_} ? $excelMeta->{$_}[2] : '' } @head;

    my @hfmt = map { '' } @head;
    for my $mm (0..$maxMis) {
        push @head, "mm$mm";
        push @wids, 4;
        push @hfmt, 'cen';
    }
    my (%mmLU, @reformat);
    for my $mm (0..$maxMis) {
        foreach my $bld (@{$allBuild}) {
            my $hcol = sprintf($bfmt, $bld, $mm);
            push @head, $hcol;
            my $ind  = $mmLU{$hcol} = $#head;
            push @wids, 4;
            my $fmt  = $buildColors{$bld} ? "Rot$bld" : "RotBuild";
            push @reformat, { 
                -col => $ind, -value => $hcol, -format => $fmt };
        }
    }
    $eh->sheet( -name    => $sname,
                -freeze  => [1,2],
                -cols    => \@head,
                # -format  => \@hfmt,
                -width   => \@wids);
    # Tie in build count formats:
    $eh->sheet($sname)->set_row( 0, 72 );
    foreach my $rf (@reformat) {
        $eh->set_cell( -sheet => $sname, -row => 0, %{$rf} ); 
    }

    foreach my $qdat (@qdats) {
        my @fmt = @fBase;
        $qdat->{len} ||= CORE::length($qdat->{seq} || "") || "";
        my $qid = $qdat->{id};
        my $len = $qdat->{len};
        my @row = ($qid, $len, $qdat->{seq}, $qdat->{desc});
        if (my $url = &oligo_url($qid)) {
            $row[0] = ['url', $url, $qid];
        }
        for my $mm (0..$maxMis) {
            my $sz     = 1500000000;
            my $bp     = $len - $mm;
            my $expect = int(0.5 + 10 * $sz / (4 ** $bp)) / 10;
            push @row, $expect;
            push @fmt, 'cen';
        }
        while (my ($hcol, $ind) = each %mmLU) {
            my $num = $qdat->{tally}{$hcol} || "";
            $row[$ind] = $num;
            my $f = 'cen';
            if ($len && $hcol =~ /^(.+)\+(\d+)$/) {
                # Add a pseudocount:
                $num = ($num || 0) + 1;
                my ($bld, $mm) = ($1, $2);
                my $bp  = $len - $mm;
                my $sz  = 1500000000;
                &msg_once("[KLUDGE]", "Setting search space size to $sz for all builds");
                my $expect = 1 + ($sz / (4 ** $bp));
                my $ratio  = $num / $expect;
                my $lod4   = int(log($ratio) / ($ln4 * $lodBand));
                # &msg("$bp = $expect vs $num = $ratio = $lod4");
                if ($lod4 < -2) {
                    $lod4 = -2;
                } elsif ($lod4 > 2) {
                    $lod4 = 2;
                }
                $f = "Expect$lod4";
            }
            
            $fmt[$ind] = $f;
        }
        $eh->add_row_explicit($sname, \@row, \@fmt);
    }
    $bnch->bench_end();
}

sub dynamic_format {
    my $eh    = shift;
    my @param = @_;
    return undef if ($#param == -1);
    my $key = join("--", @param);
    unless ( $doneFmt{$key}) {
        my $name = $doneFmt{$key} = "Dynamic".++$fmtIter;
        $eh->format( -name => $name, @param );
    }
    return $doneFmt{$key};
}

sub style_sheet_from_params {
    my $styles = $usePsh->details('styles');
    my $css = "";
    foreach my $name (sort keys %{$styles}) {
        next if ($name =~ /^_/ || $name eq 'StndName');
        my $params = $styles->{$name};
        if ($name =~ /^(\.|\#)/) { 
            $name =~ s/_/ /g;
        } else {
            $name = ".$name";
        }
        my @attr;
        while (my ($css, $val) = each %{$params}) {
            next if ($css eq 'StndName');
            push @attr, sprintf("%s: %s", $css, $val);
        }
        next if ($#attr == -1);
        $css .= sprintf(" %30s { %s }\n", $name, join('; ', @attr));
    }
    return $css;
}

sub excel_formats {    
    $bnch->bench_start();
    # my ($eh, $allRows, $allBuild, $oligos) = @_;
    my $eh = shift;

    my %unknown;
    my $css2excel = {
        'background-color' => 'bg_color',
        'color'            => 'color',
        'text-align'       => 'align',
        'font-size'        => 'size',
        'font-weight'      => sub {
            my ($val, $css) = @_;
            if ($val eq 'bold') {
                return [ bold => 1 ];
            }
            $unknown{"$css : $val"}++;
            return undef;
        },
        'font-style'      => sub {
            my ($val, $css) = @_;
            if ($val eq 'italic') {
                return [ italic => 1 ];
            }
            $unknown{"$css : $val"}++;
            return undef;
        },
        'font-family'      => sub {
            my ($val, $css) = @_;
            if ($val eq 'monospace') {
                return [ font => 'Courier New' ];
            }
            $unknown{"$css : $val"}++;
            return undef;
        },
    };
    
    my $styles = $usePsh->details('styles');
    while (my ($name, $params) = each %{$styles}) {
        next if ($name =~ /^_/);
        next if ($name =~ /^\./);
        my @conf;
        while (my ($css, $val) = each %{$params}) {
            next if ($css eq 'StndName');
            next unless (defined $val);
            $val =~ s/\s.+//;
            my $ex = $css2excel->{$css};
            if ($ex) {
                if (ref($ex)) {
                    # Callback
                    if (my $c = &{$ex}( $val, $css )) {
                        push @conf, @{$c};
                    }
                    next;
                } elsif ($ex =~ /color/) {
                    my $ind = $cu->nearest_excel_color( $val );
                    push @conf, ($ex, $ind) if ($ind);
                    next;
                }
                push @conf, ($ex, $val);
            } else {
                $unknown{"$css : $val"}++;
            }
        }
        next if ($#conf == -1);
        # warn "$name = ".join(', ', @conf)."\n";
        $eh->format( -name => $name, @conf);
    }

    my @unk = sort { $unknown{$b} <=> $unknown{$a} } keys %unknown;
    &msg("Failed to map some CSS attributes to Excel:",
               map { sprintf("%d : %s", $unknown{$_}, $_ ) } @unk)
        unless ($#unk == -1);

    return;

    # Things to revisit:
 
     $eh->format( -name       => 'pseudosym',
                 -bg_color   => 'pink',
                 -align      => 'center',
                 -bold       => 0 );
 
    $eh->format( -name       => 'unknown',
                 -bg_color   => 'silver',
                 -align      => 'center',
                 -bold       => 1 );
    $eh->format( -name       => 'unkex',
                 -bg_color   => 'silver',
                 -align      => 'center',
                 -bold       => 1 );
   
    $eh->format( -name       => 'Expect0',
                 -bg_color   => 'silver',
                 -align      => 'center', );
    $eh->format( -name       => 'Expect-1',
                 -bg_color   => 'cyan',
                 -align      => 'center', );
    $eh->format( -name       => 'Expect-2',
                 -bg_color   => 'blue',
                 -color      => 'yellow',
                 -align      => 'center', );
    $eh->format( -name       => 'Expect1',
                 -bg_color   => 'pink',
                 -align      => 'center', );
    $eh->format( -name       => 'Expect2',
                 -bg_color   => 'red',
                 -color      => 'yellow',
                 -align      => 'center', );

    $eh->format( -name       => "RotBuild",
                 -align      => 'center',
                 -rotation   => -90 );


    $eh->format( -name       => 'noMatch',
                 -background => 'silver',
                 -bold       => 0);

    # Basic monospace
    $eh->format( -name       => 'mono',
                 -size       => 9,
                 -font       => 'Courier New', );
    $eh->format( -name       => 'mono green',
                 -size       => 9,
                 -background => 'green',
                 -font       => 'Courier New', );
    

     $eh->format( -name       => 'yellow',
                 -align      => 'center',
                 -background => 'yellow');

    $eh->format( -name       => 'symblock',
                 -background => 'cyan');

    $eh->format( -name       => 'orange',
                 -align      => 'center',
                 -background => 'orange');

    $eh->format( -name       => 'red',
                 -align      => 'center',
                 -background => 'red');
    

    for my $i (0..10) {
        my $gray = int(255 * $i / 10);
        my $name = sprintf("Gray%d", $i);
        my $cc = $eh->set_custom_color(40 + $i, $gray, $gray, $gray);
        my $col = ($i <= 4) ? 'white' : 'black';
        $eh->format( -name       => $name,
                     -align      => 'center',
                     -color      => $col,
                     -bg_color   => $cc );

        $eh->format( -name       => "Dec$name",
                     -align      => 'center',
                     -num_format => '0.00',
                     -color      => $col,
                     -bg_color   => $cc );
    }

    
    $bnch->bench_end();
}

sub _entropy_fmt {
    my $v = shift;
    return "" unless (defined $v);
    $v = int(0.5 + 10 * $v);
    return "DecGray$v";
}

sub find_non_genomic_entries {
    my @files = split(/[\n\r]+/, `ls -1 /gcgblast/*_pre-mRNA_*.fa`);
    my @ids;
    foreach my $file (@files) {
        
    }
}

sub START_HTML {
    return unless ($doHTML);
    my $fh = *STDOUT;
    print $fh &_html_header();
    &_create_styles();
}

sub _create_styles {
    my $conf = $usePsh->module_path( -module => $usePsh,
                                  -suffix => "styles.conf");
    my $css  = $args->val('confstyles');
    return unless ($conf && -s $conf && $css);
    return if (-s $css && (-M $css < -M $conf));
    # warn "Need to update $css";
    if (open(CSS, ">$css")) {
        print CSS &style_sheet_from_params();
        close CSS;
        chmod(0666, $css);
    } else {
        &msg("[!]","Failed to generate style sheet", $css, $!);
    }
}

sub _html_header {
    
    my $jsType = "type='text/javascript'";
    my $xtra = "";

    foreach my $url (&_parse_resource_param('stylesheet')) {
        $xtra .= sprintf("  <link href='%s' type='text/css' rel='stylesheet' media='screen' />\n", $url);
    }
    foreach my $url (&_parse_resource_param('javascript')) {
        $xtra .= sprintf
            ("  <script src='%s' type='text/javascript'></script>\n", $url);
    }
    my $rv = <<EOF;
    
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
 <head><title>Simple Oligo Search</title>
$xtra

  <link rel="shortcut icon" href="/biohtml/images/IconSOS.png" />
 </head>
 <body bgcolor="white">
EOF
 
    return $rv;
}

sub _parse_resource_param {
    my $param = shift;
    my @rv;
    my %rep = ( BETA => $isBeta ? '/biocgi/tilfordc/working' : '');
    foreach my $url ($args->each_split_val($param)) {
        next if (!$url || $url =~ /^\s*\#/);
        while ($url =~ /(\$\$([A-Z]+?)\$\$)/) {
            # Replace $$FOO$$ with $rep{FOO}
            my ($out,$key) = ($1, $2);
            my $in = $rep{$key} || "";
            $url =~ s/\Q$out\E/$in/g;
        }
        next if ($url =~ /[\'\"<>]/);
        push @rv, $url
    }
    return @rv;
}

sub HTML_INTERFACE {
    return unless ($doHTML);
    my $fh = *STDOUT;
    # &prebranch($queries);

    my $inp = "";
    foreach my $seq (sort keys %{$queries || {}}) {
        foreach my $bs (sort { uc($a->display_id()) cmp
                                   uc($b->display_id()) } @{$queries->{$seq}}) {
            $inp .= sprintf(">%s %s\n%s\n", $bs->display_id(), 
                            $bs->desc() || "", $bs->seq());
        }
    }
    $inp = $args->esc_xml($inp);
    my $strOpts = "";
    my $opt;
    foreach $opt ([-1, '-1 (Reverse)','strR'],
                  [1,'+1 (Forward)', 'strF'],
                  [0,'Both Strands','strFR']) {
        my ($val, $show, $cls) = @{$opt};
        $strOpts .= sprintf(" <option value='%s'%s", $val,
                            $val eq $geneStr ? " selected='1'" : "");
        $strOpts .= sprintf(" class='%s'", $cls) if ($cls);
        $strOpts .= sprintf(">%s</option>\n", $show);
    }
    my $modeOpts = "";
    my $chkMode  = $isAuto ? "" : $mode;
    foreach $opt (['Summary'],
                  ['Overview'],
                  ['Details'],
                  ['Genomic View'] ) {
        my ($val) = @{$opt};
        $modeOpts .= sprintf(" <option value='%s'%s>%s</option>\n", $val, 
                             $val eq $chkMode ? " selected='1'" : "", $val);
    }
    my $frmOpts = "";
    foreach my $opt (['','Auto'],
                     ['HTML'], ['Excel'], ['TSV'], ['Text'], ['JSON']) {
        my ($val, $show, $cls) = @{$opt};
        $show ||= $val;
        $frmOpts .= sprintf(" <option value='%s'%s", $val,
                            $val eq $frmReq ? " selected='1'" : "");
        $frmOpts .= sprintf(" class='%s'", $cls) if ($cls);
        $frmOpts .= sprintf(">%s</option>\n", $show);
    }
    my $mmOpts = "";
    foreach $opt ([2, '2'],
                  [1,'1'],
                  [0,'0 (Perfect Match)']) {
        my ($val, $show) = @{$opt};
        $mmOpts .= sprintf(" <option value='%s'%s>%s</option>\n", $val, 
                           $val eq $maxMiss ? " selected='1'" : "", $show);
    }
    if ($maxMiss == 3) {
        $mmOpts .= sprintf(" <option value='%s'%s>%s</option>\n", $maxMiss,
                           "selected='1'", $maxMiss);

    }

    my $extjs = &ext_setup_script( $formId );

    print $fh <<EOF;
$extjs
<form id='$formId' method='POST'>
<table class='tab'><tbody><tr><td>       
<b>Oligo requests:</b><br />
<textarea name='input' style='width:30em; font-family:monospace; height:20em; background-color:#9f9'>$inp</textarea><br />
<input style='font-size:1.5em; font-weight: bold; background-color:lime;' type='submit' value='Search' /> <span id='srchmsg'></span>
</td><td>
<b>Search Strand:</b> <select name='strand'>
$strOpts</select><br />
<b>Maximum mismatches:</b> <select name='mismatch'>
$mmOpts</select><br />
<b>Mode:</b> <select name='mode'>
$modeOpts</select>
<select name='format'>
$frmOpts</select>
<br />
<b>Limit to specific genes / RNAs:</b><br>
<textarea style='width:20em; font-family:monospace; height:8em; background-color:#9cf' name='genereq'>$keepReq</textarea><br />
<b>Exclude:</b><br>
<textarea style='width:20em; font-family:monospace; height:8em; background-color:#f66' name='exclreq'>$exclReq</textarea><br />
<a class='butt ybbutt' href='simpleOligoSearch.pl'>Clear Results / New Search</a><br />
</td></tr></tbody></table>
</form>
EOF

    if ($isBeta) {
        print $fh "<p class='beta'>This is the beta interface. It is under active development and may not always work. Unless you have a reason to use this interface, please use the <a href='/biocgi/simpleOligoSearch.pl'>released version</a> of the code</p>\n";
        # print $fh "<p style='font-size:2em;' class='alert'>23 Jan 2015 - Beta is under development and is non-functional - you <i>must</i> use the <a href='/biocgi/simpleOligoSearch.pl'>released version</a>!</p>";
    }
    
}

sub ext_setup_script {
    return "" if ($nocgi || $extDivDone);
    $extDivDone = 'extout';
    my $formId = shift;
    
    return <<EOF;
<div id='$extDivDone'></div>
<script>
function sos_fetch_resource (rsp, url, params) {
    var data;
    try {
        data = JSON.parse( rsp.responseText );
    } catch (e) {
        alert("Failed to parse JSON data from '"+url+"'");
    }
    return process_ss_json( data, null, url, params);
};
Ext.onReady(function () {
    // pgtm_basic_ajax_launch('/biohtml/temp/PgSeqHash/2015-05-06/tilfordc-2687/output/Task_1923.json', null, sos_fetch_resource);
    query_panel( '$formId' );
});
</script>
    
EOF
}

sub prebranch {
    my $foo = $args->branch(@_);
    if ($doHTML) {
        print "<pre>$foo</pre>";
    } else {
        warn $foo;
    }
}

sub contact_html {
    return <<EOF;
 <span style='color:brown; font-style: italic;'>
Please contact <a href='mailto:charles.tilford\@bms.com'>Charles Tilford</a>
if you have any questions on this tool.
</span><br />\
EOF
 
}

sub _process_gene_request {
    my ($arg, $msg) = @_;
    my ($rnas, $srcs, $gids, $subs, $stuct);
    my $pretty = "";
    my $struct;
    foreach my $req ($args->each_split_val($arg)) {
        $req = "" unless (defined $req);
        # Remove comments at end of line
        my $com = "";
        if ($req =~ /^(.*?)(\#.+)$/) {
            ($req, $com) = ($1, $2);
        }
        # Remove leading and trailing whitespace:
        $req =~ s/^\s+//; $req =~ s/\s+$//;
        unless ($req) {
            $pretty .= "$com\n" if ($com);
            next;
        }
        if ($com) {
            $com = " $com";
        } else {
            $com = "";
        }
        if ($req =~ /^([NX][MR]_|ENS[A-Z]*T\d)/) {
            $rnas ||= [];
            push @{$rnas}, $req;
            $pretty .= "$req$com\n";
            $filterNotes ||= {};
            $filterNotes->{"$msg RNA"}{$req} ||= $com;;
        } else {
            my @srcid = $srcBldSth->get_array_for_field($req);
            if ($#srcid != -1) {
                # Appears to be a build token
                $srcs ||= [];
                push @{$srcs}, @srcid;
                $pretty .= "$req$com\n";
                $filterNotes ||= {};
                $filterNotes->{"$msg Build"}{$req} ||= $com;;
                next;
            }
            my @subid = $subNameSth->get_array_for_field($req);
            if ($#subid != -1) {
                $subs ||= [];
                push @{$subs}, @subid;
                $pretty .= "$req$com\n";
                $filterNotes ||= {};
                $filterNotes->{"$msg Subject"}{$req} ||= $com;;
                next;
            }
            # Assume it is a gene.
            my $sth  = $usePsh->named_sth('Find gene by ID or alias');
            my @gids = $sth->get_array_for_field($req, $req);
            if ($#gids == -1) {
                &msg("[?]","Could not find any genes matching '$req'");
            } else {
                $gids ||= [];
                push @{$gids}, @gids;
                $pretty .= "$req$com\n";
                $filterNotes ||= {};
                $filterNotes->{"$msg Gene"}{"$req gene_id IN (".join(',',@gids).")"} ||= $com;;
            }
        }
    }
    if ($msg && $struct) {
        my @ma = ($msg);
        foreach my $k (sort keys %{$struct}) {
            push @ma, "$k:";
            push @ma, map { "  $_" } sort keys %{$struct->{$k}};
        }
        &msg("[-]", @ma);
    }
    return ($rnas, $gids, $srcs, $subs, $pretty);
}

sub note_filters {
    return if (!$filterNotes || $clean);
    return if (!$queries);
    my @u = keys %{$queries};
    return if ($#u == -1);
    my @msg;
    foreach my $filt (sort keys %{$filterNotes}) {
        push @msg, $doHTML ? "<b>$filt</b>" : $filt;
        foreach my $req (sort keys %{$filterNotes->{$filt}}) {
            my $what = $req;
            if (my $com = $filterNotes->{$filt}{$req}) {
                $what .= sprintf(" %s", $doHTML ? "<span class='note'>$com</span>" : $com);
            }
            push @msg, $what;
        }
    }
    &msg("User filters applied", @msg) unless ($args->val('norun'));
}

sub _work_dir {
    my ($s,$m,$h, $day, $mon, $year) = localtime(time);
    my $tmp = $tempDir;
    $tmp   .= sprintf("/%04d-%02d-%02d", 1900 + $year, $mon + 1, $day);
    $tmp   .= sprintf("/%s-%d", $user, $$);
    $args->assure_dir($tmp);
    return $tmp;
}

sub _temp_exclude {
    return <<EOF;
# Highly polymorphic noise:
NG_011498.1 # IL23R
NG_008875.1 # GPHN
NG_033073.1 # DOCK7
NG_008617.1 # PKD1
NG_029922.1 # HLA-DQB1
NG_029921.1 # HLA-DRB1
NG_029217.2 # HLA-A
NG_023187.1 # HLA-B
NG_029422.2 # HLA-C
NG_000007.3 # Hemoglobin B
NG_000006.1 # Hemoglobin A
EOF
    
}
