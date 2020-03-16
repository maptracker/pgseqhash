# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::PgSeqHash;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use warnings;

use JSON;

use parent qw( BMS::FriendlyDBI
               BMS::Utilities::SequenceUtilities
               BMS::Utilities::FileUtilities
               BMS::Utilities::Benchmark );

use Bio::PrimarySeq;

our $maxShownMM = 3;
our $gettingStartedProgress = 0.0001;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
    };
    bless ($self, $class); 
    $self->bench_start();
    $self->configure( @_ );
    $self->connect( @_ );
    $self->bench_end();

    return $self;
}

our $nextPrompt;

our $typePriority = {
    CodeError      => 1,
    CrypticExon    => 5,
    FullSpan       => 6,
    ExonExon       => 10,
    Exon           => 20,
    RNA            => 25,
    SpliceDonor    => 30,
    SpliceAcceptor => 30,
    SpliceJunction => 30,
    FirstExon      => 40,
    LastExon       => 50,
#    ExonInDel      => 60, # Not used
    Intron         => 70,
    Unknown        => 80,
    NoGeneInfo     => 99,
};

our $normalStrand = {
    '-1'  => -1,
    'R'   => -1,
    'REV' => -1,
    '+1'  => 1,
    '1'   => 1,
    'F'   => 1,
    'FWD' => 1,
    0     => 0,
    'FR'  => 0,
    'RF'  => 0,
};

sub DESTROY {
    my $self = shift;
    if (my $dbh = $self->{DBH}) {
        $dbh->disconnect();
    }
}

sub normalize_strand {
    my $self = shift;
    for my $i (0..$#_) {
        # Cycle through more than one argument, allows:
        # my $str = $self->normalize_strand( $args->{STR}, $args->{STRAND})
        my $req = $_[$i];
        next unless (defined $req);
        my $rv = $normalStrand->{uc($req)};
        return $rv if (defined $rv);
    }
    return undef
}

sub instance  { return shift->{INSTANCE}; }
sub dbuser    { return shift->{DBUSER}; }
sub dbpass    { return shift->{DBPASS}; }
sub task      { return shift->{TASK}; }
sub named_sth { return shift->dbh->named_sth( @_ ); }
sub prepare   { return shift->dbh->prepare( @_ ); }

sub configure {
    my $self = shift;
    $self->bench_start();
    # Try to read connection details from parameter file
    my $params = BMS::ArgumentParser->new
        ( -paramfile  => 'BMS::PgSeqHash',
          @_);
    foreach my $key ('PGPORT', 'PGHOST') {
        if (my $val = $params->val($key)) {
            $ENV{$key} = $val;
        }
    }
    $self->{PGTM}     = $params->val('pgtm');
    $self->{TASK}     = $params->val('task');
    $self->{INSTANCE} = $params->val('instance') || 'seqhash';
    $self->{DBUSER}   = $params->val('dbuser') || "";
    $self->{ADMIN}    = $params->val('admin') || '';
    $self->{ERRFILE}  = $params->val('errorfile');
    $self->{ERRMAIL}  = $params->val('adminmail');
    $self->{DBPASS}   = $params->val('dbpass');
    if (my $dbh = $self->{DBH}) {
        $dbh->close();
        delete $self->{DBH};
    }
}

sub pgtm {
    my $self = shift;
    if (defined $_[0]) {
        if ($_[0]) {
            $self->{PGTM} = $_[0];
        } else {
            delete $self->{PGTM};
        }
    }
    return $self->{PGTM};
}
sub get_srch_ids {
    my $self = shift;
    my @rv;
    my ($seq) = $self->_normalize_sequence( shift );
    if ($seq) {
        my $mm = shift || 0;
        $self->bench_start();
        my $sth = $self->dbh->named_sth('Get search_ids by sequence query');
        @rv = $sth->get_array_for_field( $seq, $mm );
        $self->bench_end();
    }
    return wantarray ? @rv : $rv[0];
}

sub set_srch_id {
    my $self = shift;
    my $req  = shift;
    my $mm   = shift || 0;
    my ($seq, $name, $desc) = $self->_normalize_sequence( $req, shift, shift );
    return 0 unless ($seq);
    my $sth = $self->dbh->named_sth('Create a new search');
    my $rv = $sth->get_single_value($seq, $mm, $name, $desc) || 0;
    return $rv;
}

sub _normalize_sequence {
    my $self = shift;
    my $req = shift;
    return ("") unless ($req);
    my ($seq, $name, $desc);
    if (ref($req)) {
        # Assume a Bio::Seq
        $seq  = $req->seq();
        $name = $req->display_id();
        $desc = $req->desc();
    } else {
        $seq  = $req;
        $name = shift;
        $desc = shift;
    }
    $seq ||= "";
    # Strip junk characters
    # Preserve dashes for the moment - do they have value?
    $seq =~ s/[^a-z\-]+//gi;
    # Translate Uracil
    $seq =~ tr/Uu/Tt/;
    # warn "$req = $seq";
    return ($seq, $name, $desc);
}

sub search_status {
    my $self = shift;
    my $srchid = shift || "";
    return 0 unless ($srchid);
    $self->bench_start();
    my $stat = shift;
    if ($stat) {
        my $sth = $self->dbh->named_sth('Set search status');
        $sth->execute($stat, $srchid);
    } else {
        my $sth = $self->dbh->named_sth('Get search status');
        $stat = $sth->get_single_value( $srchid );
    }
    $self->bench_end();
    return $stat;
}

sub clear_search {
    my $self = shift;
    $self->bench_start();
    my $srchid = shift || "";
    unless ($srchid =~ /^\d+$/) {
        $self->err("Non-integer srch_id = '$srchid'") if ($srchid);
        return ();
    }
    my $sth = $self->dbh->named_sth('Delete a search by srch_id');
    my $rv = $sth->execute( $srchid );
    $self->bench_end();
    return $rv;
}

sub search_details {
    my $self = shift;
    $self->bench_start();
    my $srchid = shift || "";
    unless ($srchid =~ /^\d+$/) {
        $self->err("Non-integer srch_id = '$srchid'") if ($srchid);
        return ();
    }
    my $sth = $self->dbh->named_sth('Get search details for srch_id');
    my $rows = $sth->selectall_arrayref( $srchid );
    $self->bench_end();
    return @{$rows->[0] || []};
}

sub bioseq_for_search_id {
    my $self = shift;
    my ($seq, $name, $descr) = $self->search_details( @_ );
    return undef unless ($seq);
    return Bio::PrimarySeq->new
            ( -seq        => $seq,
              -alphabet   => 'DNA',
              -desc       => $descr,
              -display_id => $name, );
}

sub search_details_by_name {
    my $self = shift;
    $self->bench_start();
    my $req  = shift || "";
    $req =~ s/^\s+//;
    $req =~ s/\s+$//;
    my $sth  = $self->named_sth('Get search details for name');
    $sth->execute($req, $req);
    my $found = $sth->fetchall_arrayref();
    $self->bench_end();
    return $found;   
}

sub hits_for_sequence {
    my $self = shift;
    $self->bench_start();
    my @passed = @_;
    my $pgtm   = $self->pgtm();
    unshift @passed, "-seq" if ($#passed == 0);
    my $args    = $self->parseparams( @passed );
    my $req     = $args->{SEQ} || $args->{SEQS};
    my ($seq)   = $self->_normalize_sequence( $req ); # ref($req) ? $req->bs() : $req;
    # Need to eventually make this smart
    # Either pull all possible from sourcedb table, or let user choose:
    my @wsz     = sort { $b <=> $a } (12);
    my $wsth    = $self->_get_wrd_id_sth();
    my $getHits = $self->dbh->named_sth('Get all hits for a word in a source');
    my $getSrcs = $self->dbh->named_sth('Get src_ids for word');
    my $mm      = defined $args->{MM} ? $args->{MM} : 
        defined $args->{MISMATCH} ? $args->{MISMATCH} : 0;
    
    my $rcSeq = $self->revcom( $seq );
    my $slen  = CORE::length($seq);
    my @useSz = @wsz;
    my $getWids = sub {
        # Get the word ID for an oligo
        my $seq = shift;
        my $wid = $wsth->get_single_value( $seq );
        return $wid ? ($wid) : ();
    };
    if ($slen < $wsz[0]) {
        # The query is shorter than the longest wordsize
        # We will need a trailing wildcard to recover all hits
        $self->death("Have not coded search for small sequences yet ($seq)");
    }
    my $slenMinusOne = $slen - 1;
    my $hits = {
        mm     => $mm,
        strand => {},
        query  => {
            1  => $seq,
            -1 => $rcSeq,
        },
    };
    my $promptInterval = 1;
    my @subtask = ( ["Strand",0], ["Recover",0], ["Assemble",0] );
    my $info   = "Find and store hits";
    my $status = {
        task    => $info,
        subtask => \@subtask,
    };
    $pgtm->change_task_status($status) if ($pgtm);
    my $process = $hits->{process} = {};
    my @isRCs = (0,1);
    my $numRCs = $#isRCs + 1;
    for my $irci (0..$#isRCs) {
        my $isRC  = $isRCs[$irci];
        my $str   = $isRC ? -1 : 1;
        my $niStr = $str > 0 ? '+1' : $str;
        $status->{task} = "$info : Strand $niStr";
        $nextPrompt ||= time - 1 if ($pgtm);
        $subtask[0][1] = $gettingStartedProgress;
        $subtask[0][2] = $niStr;
        $pgtm->change_task_status($status) if ($pgtm);
        foreach my $ws (@wsz) {
            $self->bench_start("Identify bait words");
            my $wsMinusOne     = $ws - 1;
            my (%sources, @allInds);
            # Find sources that are hit by the permuted words
            my @fullOligoList;
            my @offsetInds;
            for (my $i = 0; $i < $slen; $i += $ws) {
                push @offsetInds, $i;
            }
            for my $oi (0..$#offsetInds) {
                my $i    = $offsetInds[$oi];
                my $ind  = ($i + $ws < $slen) ? $i : $slen - $ws;
                my $bait = substr($isRC ? $rcSeq : $seq, $ind, $ws);
                my $pmm  = $self->permute_mismatches( $bait, $mm );
                push @allInds, $ind;
                while (my ($mmSeq, $mmDat) = each %{$pmm}) {
                    push @fullOligoList, [$ind, $mmSeq, $mmDat];
                }
            }
            my $tot = $#fullOligoList + 1;
            $self->bench_end("Identify bait words");

            # These are broken out to allow progres reporting
            $self->bench_start("Recover database sources");
            for my $mmi (0..$#fullOligoList) {
                my ($ind, $mmSeq, $mmDat) = @{$fullOligoList[$mmi]};
                foreach my $wid (&{$getWids}( $mmSeq )) {
                    my @srcs = $getSrcs->get_array_for_field($wid);
                    map { push @{$sources{$_}{$ind}}, [$wid, $mmDat] } @srcs;
                    
                    if ($nextPrompt && time >= $nextPrompt) {
                        $subtask[1][1] = &_progress_fraction
                            ( $mmi, $#fullOligoList );
                        $subtask[1][2] = $mmSeq;
                        $pgtm->change_task_status($status);
                        $nextPrompt = time + $promptInterval;
                    }
                }
            }
            $subtask[1][1] = 1;
            $subtask[1][2] = "Done";
            $subtask[2][1] = $gettingStartedProgress;
            $subtask[2][2] = "Starting...";
            $pgtm->change_task_status($status) if ($pgtm);
            
            my $expectedChunks = $#allInds;
            $self->bench_end("Recover database sources");

            # Find out which sources are ok
            my @allSources = sort { $a <=> $b } keys %sources;
            for my $asi (0..$#allSources) {
                my $srcid = $allSources[$asi];
                if ($pgtm) {
                    my $si  = $self->source_info( $srcid );
                    $subtask[2][2] = $si->{GenomeBuild} || "Unknown Build";
                    $pgtm->change_task_status($status);
                }
                my $indH  = $sources{$srcid};
                my @inds  = sort { $a <=> $b } keys %{$indH};
                # The source should be hit by every word
                unless ($#inds == $expectedChunks) {
                    $process->{"Untenable database source"}++;
                    next;
                }

                $self->bench_start("Recover subjects");
                # Ok, this source may have survivable subjects
                # Structure the subject data
                my %subjects;
                foreach my $ind (@inds) {
                    foreach my $hdat (@{$indH->{$ind}}) {
                        my $hits = $getHits->
                            get_single_value($hdat->[0], $srcid);
                        # Get the first subject ID
                        my $sid = $hits->[0];
                        my $s   = 0;
                        while ($s < $#{$hits}) {
                            $s++;
                            if (my $hpos = $hits->[$s]) {
                                # A match coordinate - adjust by the index
                                my $pos = $hpos - $ind;
                                push @{$subjects{$sid}{$pos}{$ind}}, $hdat;
                            } else {
                                # The zero sbj_id separator
                                # Get the next entry as the new sbj_id
                                $s++;
                                $sid = $hits->[$s];
                            }
                        }
                    }
                }
                $self->bench_end("Recover subjects");

                $self->bench_start("Recover positions");
                my @allsids = keys %subjects;
                my $allsidN = $#allsids + 1;
                for my $sidi (0..$#allsids) {
                    my $sid = $allsids[$sidi];
                    my $posH = $subjects{$sid};
                    while (my ($pos, $indH) = each %{$posH}) {
                        my @inds = keys %{$indH};
                        # The position should be hit by every word:
                        unless ($#inds == $expectedChunks) {
                            $process->{"Untenable subject position"}++;
                            next;
                        }
                        # Ok, we have a segment in a subject that has been
                        # covered by all word indices. Now we just need
                        # to do a final overall mismatch filter
                        # We will record all positions that are matched
                        my @matched;
                        while (my ($ind, $hdats) = each %{$indH}) {
                            foreach my $hdat (@{$hdats}) {
                                my %mmInd = map { $_ => 1 } @{$hdat->[1]};
                                for my $i (0..$wsMinusOne) {
                                    $matched[$ind + $i] ||= $mmInd{$i} ? 0 : 1;
                                }
                            }
                        }
                        unless ($#matched == $slenMinusOne) {
                            $self->msg("[?]","Failed to calculate mismatch data",
                                       "sbj_id = $sid, pos = $pos");
                        }
                        my @mmInds;
                        map { push @mmInds, $_ 
                                  unless ($matched[$_]) } (0..$#matched);
                        my $hitMM = $#mmInds + 1;
                        # If we have too many mismatches in the combined
                        # sequence, then skip
                        if ($hitMM > $mm) {
                            $process->{"Mismatch > $mm"}++;
                            next;
                        }
                        # The sequence is being kept
                        $process->{"MM$hitMM kept"}++;
                        my $phit = $hits->{strand}{$str}{$sid}{$pos} =
                            [ \@mmInds, [], $srcid ];
                        while (my ($ind, $hdats) = each %{$indH}) {
                            foreach my $hdat (@{$hdats}) {
                                push @{$phit->[1]}, [$hdat->[0], $ind];
                            }
                        }
                    }
                    if ($nextPrompt && time >= $nextPrompt) {
                        my $denom = $#allSources + 1;
                        my $frac  = ($asi + (($sidi+1) / $allsidN))/ $denom;
                        $subtask[2][1] = int(0.5 + $frac * 1000) / 1000;
                        $pgtm->change_task_status($status);
                        $nextPrompt = time + $promptInterval;
                    }
                }
                $self->bench_end("Recover positions");
                if ($pgtm) {
                    my $frac = &_progress_fraction( $asi, $#allSources);
                    $subtask[2][1] = $frac;
                    $pgtm->change_task_status($status);
                }
            }
            $subtask[2][1] = 1;
            $subtask[2][2] = "Done";
            $pgtm->change_task_status($status) if ($pgtm);
        }
    }
    $subtask[0][1] = 1;
    $subtask[0][2] = "Done";
    $pgtm->change_task_status($status) if ($pgtm);
    $self->bench_end();
    return $hits;
}

sub _progress_fraction {
    my ($ind, $len) = @_;
    return 0 if (!defined $len || $len < 0);
    $ind++;
    $len++;
    my $rv = int(0.5 + 1000 * $ind/$len) / 1000;
    return $rv || $gettingStartedProgress;
}

sub flatten_hits {
    my $self = shift;
    my ($hits, $srchid) = @_;
    unless ($hits) {
        return {
            cols => [],
            rows => [],
            flat => 1,
            type => 'null',
        };
    }
    $self->bench_start();
    my $wrdSTH  = $self->_get_word_sth();
    my $ambigLU = {};
    my $seq     = $hits->{query}{1};
    my $task    = $self->task();
    my (@rv, @got);
    my $geneStrSth = $self->dbh->named_sth
        ('Get gene strands for subject position');
    my ($bar, $nextPrompt);
    my $promptInterval = 1;
    if ($task) {
        my $bars = $task->param("progbars", []); # Will clear older bars!
        $bar  = ['Flatten data', 0];
        push @{$bars}, $bar;
        $task->write(1);
        $nextPrompt ||= time - 1;
    }
    my @allStrSubs;
    foreach my $str (sort { $b <=> $a } keys %{$hits->{strand}}) {
        foreach my $sid (sort { $a <=> $b } keys %{$hits->{strand}{$str}}) {
            push @allStrSubs, [$str, $sid];
        }
    }
    my %qChars;
    for my $ssi (0..$#allStrSubs) {
        my ($str, $sid) = @{$allStrSubs[$ssi]};
        my @chars = @{$qChars{$str} ||= [split('', uc($hits->{query}{$str}))]};
        my $cinds = $#chars;
        my $pH    = $hits->{strand}{$str}{$sid};
        while (my ($pos, $phit) = each %{$pH}) {
            $self->bench_start('Build subject sequence');
            # Build the hit words at this position into an
            # ambiguously-coded genome string
            my @ambig = map { {} } @chars;
            foreach my $wd (@{$phit->[1]}) {
                my ($wid, $ind) = @{$wd};
                my $word = $wrdSTH->get_single_value($wid);
                # $self->msg("$word = $wid [$ind]");
                my @wchar = split('', $word);
                for my $w (0..$#wchar) {
                    $ambig[ $w + $ind ]{ $wchar[$w] } = 1;
                }
            }

            # @ambig now holds every observed base at each position
            # Turn it into a subject hit string
            my $mm = 0;
            my $subSeq = "";
            for my $i (0..$#chars) {
                my $aH   = $ambig[$i] || {};
                my @c    = sort keys %{$aH};
                my $k    = join('', @c) || '-';
                my $char = $ambigLU->{$k} ||= 
                    $self->ambiguous_code( \@c );
                unless ($aH->{$chars[$i]}) {
                    # Mismatch position;
                    $char = lc($char);
                    $mm++;
                }
                $subSeq .= $char;
            }
            $subSeq = $self->revcom( $subSeq ) if ($str < 0);
            $self->bench_end('Build subject sequence');

            my $srcid = $phit->[2];
            my $end   = $pos + $cinds;
            my $fp    = $self->genomic_footprint($sid, $pos, $end, $str);
            my $genes = $self->subject_position_to_genes
                ( $sid, $pos, $end, $str );
            my %geneStrs;
            while (my ($gid, $strH) = each %{$genes || {}}) {
                map { $geneStrs{$_}++ } keys %{$strH};
            }
            my $geneStr = 
                ($geneStrs{0} || ($geneStrs{-1} && $geneStrs{1})) ? 0 : 
                $geneStrs{-1} ? -1 : $geneStrs{1} ? 1 : undef;
            #my @flat = ($subSeq, $mm, $sid, $pos, $str, 
            #            $geneStr, $srcid, $fp, \@flatGst);
            my @flat = ($subSeq, $mm, $sid, $pos, $str, 
                        $geneStr, $srcid, $fp, $genes);
            push @flat, $srchid if ($srchid);
            push @got, \@flat;
            if ($nextPrompt && time >= $nextPrompt) {
                $bar->[1] = &_progress_fraction( $ssi, $#allStrSubs);
                $bar->[2] = $self->subject_name($sid);
                $task->write(1);
                $nextPrompt = time + $promptInterval;
            }
        }
    }
    my @cols = ('SubjectSeq','MisMatch','SubjectID','SubjectPos',
                'SubjectStrand','GeneStrands','SourceID','GenomicFootprint','GeneDetails');
    push @cols, 'SearchID' if ($srchid);
    if ($task) {
        $bar->[1] = 1;
        $bar->[2] = "Done.";
        $task->write(1);
    }
    my $rv = {
        cols => \@cols,
        rows => \@got,
        flat => 1,
    };
    $self->bench_end();
    return $rv;
}

sub genomic_footprint {
    my $self = shift;
    my ($sid, $s, $e, $str) = @_;
    # $sid = subj_id for the subject
    # $s   = start coordinate of query relative to subject
    # $e   = end coordinate
    # $str = strand of query relative to subject
    return "" unless ($sid && $s);
    if (!$e) {
        $e = $s;
    } elsif ($e < $s) {
        ($s, $e) = ($e, $s);
    }
    $str ||= 1;
    my $sname = $self->subject_name($sid);
    my ($anc, @bits);
    if ($sname =~ /^(.+?):(\d+)\-(\d+)/) {
        # Contiguous genome chunk - we just need to offset the subject positions
        $anc = $1;
        my $off = $2 - 1;
        if ($s == $e) {
            push @bits, $s + $off;
        } else {
            push @bits, sprintf("%d..%d", $s + $off, $e + $off);
        }
    } elsif ($sname =~ /^(.+?):(\d+)\^(\d+)/) {
        # Junction sequence, generally exon-exon
        $anc = $1;
        my ($l, $r) = ($2, $3);
        # $l = Genome coordinate to the left of the junction
        # $r = Genome coordinate to the right
        # my @pws = $self->get_tag($sid, 'P.W');
        my @pwss = $self->get_tag($sid, 'P W S');
        if ($#pwss > 0) {
            $self->err("Multiple PWS values : ".join(',', @pwss));
            return "";
        } elsif (my $pws = $pwss[0]) {
            # The position and width of the boundary is defined
            
            # These entries are always oriented +1 relative to the
            # *RNA*.  This means that mapping to genome coordiates
            # needs to take into account the RNA-to-genome strand.
            
            if ($pws =~ /^(\d+) (\d+) (\-?1)$/) {
                my ($pl, $w, $xs) = ($1, $2, $3);
                my $isRC = $xs < 0 ? 1 : 0;
                my $lOff = $isRC ? $r + $pl : $l - $pl;
                # $pl = the SUBJECT coordinate to the left of the boundary
                # $w  = the number of bases spanning the boundary
                #       2 = a normal boundary - no RNA in between
                #           Any other value means some RNA could not be aligned
                #       5 = a boundary with 3 bases of unaligned RNA in middle
                # $xs = the strand of the gene relative to the genome
                # $isRC = 1 if the strand is -1
                # $lOff = genomic offset on left side of boundary
                if ($e <= $pl) {
                    # The segment is entirely on the "left" of the boundary
                    if ($isRC) {
                        push @bits, sprintf("%d..%d", $lOff - $e, $lOff - $s);
                    } else {
                        push @bits, sprintf("%d..%d", $lOff + $s, $lOff + $e);
                    }
                } else {
                    # At least some stuff on the "right"
                    my $pr   = $pl + $w - 1;
                    my $rOff = $isRC ? $l + $pr  : $r - $pr;
                    # $pr = subject coordinate to the right of boundary
                    # $rOff = genomic offset on the right side of boundary
                    if ($s >= $pr) {
                        # The segment is entirely on the "right" of the boundary
                        if ($isRC) {
                            push @bits, sprintf("%d..%d", $rOff -$e, $rOff -$s);
                        } else {
                            push @bits, sprintf("%d..%d", $rOff +$s, $rOff +$e);
                        }
                    } else {
                        # The alignment spans the junction, or falls inbetween
                        if ($s <= $pl) {
                            # Note the left HSP:
                            if ($isRC) {
                                my $lside = $lOff - $s;
                                push @bits, $lside == $r ? $lside : 
                                    sprintf("%d..%d", $r, $lside);
                            } else {
                                my $lside = $lOff + $s;
                                push @bits, $lside == $l ? $lside : 
                                    sprintf("%d..%d", $lside, $l);
                            }
                        }
                        if ($w > 2) {
                            # There are some bases that are unaligned and
                            # "in the middle". Make a non-standard notation
                            # for this - for example:
                            # 11732^7^19818 = There are seven bases SOMEWHERE
                            #   between 11,732 and 19,818
                            # How much is actually in the middle:
                            my $intS = ($s <= $pl) ? $pl + 1 : $s;
                            my $intE = ($e >= $pr) ? $pr - 1 : $e;
                            my $iLen = $intE - $intS +1;
                            push @bits, sprintf("%d^%d^%d", $l, $iLen, $r);
                        }
                        if ($e >= $pr) {
                            # Note the right HSP:
                            if ($isRC) {
                                my $rside = $rOff - $e;
                                push @bits, $rside == $l ? $rside : 
                                    sprintf("%d..%d", $rside, $l);
                            } else {
                                my $rside = $rOff + $e;
                                push @bits, $rside == $r ? $rside : 
                                    sprintf("%d..%d", $r, $rside);
                            }
                        }
                    }
                }
                if ($isRC) {
                    @bits = reverse @bits;
                    $str *= -1;
                }
            } else {
                $self->err("Unrecognized PWS value '$pws'");
                return "";
            }
        } else {
            $self->err("No PWS value provided for subject $sname");
           return "";
        }
    } else {
        # Treat the sequence as the footprint itself
        $anc = $sname;
        push @bits, $s == $e ? $s : "$s..$e";
    }
    my $rv = $anc . ':' . join(',', @bits);
    $rv .= '[-1]' if ($str < 0);
    return $rv;
}

sub rna_footprint {
    my $self   = shift;
    $self->bench_start();
    
    $self->bench_end();
}

sub compact_hits {
    my $self   = shift;
    my ($req) = @_;
    $self->bench_start();
    my $task = $self->task();
    my ($bar, $nextPrompt);
    my $promptInterval = 1;
    if ($task) {
        my $bars = $task->param("progbars") || $task->param("progbars", []);
        $bar  = ["Compact Results", 0];
        push @{$bars}, $bar;
        $task->write(1);
        $nextPrompt ||= time - 1;
    }
    my $data     = $self->_assure_flattened_hits( $req );
    # die $self->branch($data);
    # my $strReq = shift;
    # my $flat   = $self->flater_flat_hits( -hits => $flat, -str => $strReq );
    my $flat     = $data->{rows};
    my $flatinds = $#{$flat};
    my $Ex0col = 2;
    my $rv       = {
        cols => [qw(SearchID SourceID GeneIDsEx0 GeneIDsEx1 GeneIDsEx2 GeneIDsInt0 GeneIDsInt1 GeneIDsInt2)],
        type => 'overview',
        flat => 1,
        genes => $data->{genes},
    };
    my $rows = $rv->{rows} = [];
    if ($flatinds == -1) {
        $self->bench_end();
        return $rv;
    }
    my (@indices, @noInd);
    foreach my $colName (qw(MisMatch SearchID SourceID 
                            GeneDetails GenomicFootprint)) {
        my $ind = $self->data_index( $data, $colName );
        push @indices, $ind;
        push @noInd, $colName if ($ind < 0);
    }
    unless ($#noInd == -1) {
        $self->msg("Can not compact data, some columns missing: ".
                   join(',', @noInd));
        $self->bench_end();
        return $rv;        
    }
    my ($mmInd, $srchInd, $srcInd, $detInd, $fpInd) = @indices;
    my $sbjSTH   = $self->_get_subject_sth();
    my %struct;
    for my $ri (0..$flatinds) {
        my $row     = $flat->[$ri];
        my $mm      = $row->[$mmInd];
        my $srchid  = $row->[$srchInd] || 0;
        my $srcid   = $row->[$srcInd]  || 0;
        my $details = $row->[$detInd]  || {};
        my $fp      = $row->[$fpInd]   || "";
        
        if (!$mm) {
            $mm = 0;
        } elsif ($mm > $maxShownMM) {
            $mm = $maxShownMM;
        }
        
        my $srcTarg = $struct{$srchid}{$srcid} ||= {};
        # $self->prebranch( $data->{cols} ); $self->prebranch( $row); die;
        while (my ($gid, $gsH) = each %{$details}) {
            while (my ($geneStr, $etH) = each %{$gsH}) {
                while (my ($type, $rnaArr) = each %{$etH}) {
                    my $gTarg  = $srcTarg->{$gid}  ||= [undef, undef];
                    if ($type eq 'RNA') {
                        # These are hard to count in a normalized way
                        # We track them separately
                        my $rnaTarg = $gTarg->[1] ||= [];
                        $rnaTarg->[$mm]++;
                    } else {
                        # Should be a normalized genomic location
                        my $genTarg = $gTarg->[0] ||= {};
                        my $exType  = !$type ? 'Unk' : 
                            $type eq 'Intron' ? 'Int' : 'Ex';
                        $genTarg->{$exType}{$fp} = $mm if
                            (!defined $genTarg->{$fp} ||
                             $genTarg->{$fp} > $mm);
                    }
                }
            }
        }
        # my $exTarg  = $srcTarg->{ $exType } ||= [];
        # $exTarg->[$mm]{$gid}++;
        if ($nextPrompt && time >= $nextPrompt) {
            $bar->[1] = &_progress_fraction( $ri, $flatinds);
            $task->write(1);
            $nextPrompt = time + $promptInterval;
        }
    }
    # $self->prebranch( \%struct); die;
    # We need to now deal with the normalized information
    while (my ($srchid, $srcH) = each %struct) {
        while (my ($srcid, $srcTarg) = each %{$srcH}) {
            my %norm;
            foreach my $gid (keys %{$srcTarg}) {
                my $gTarg = $srcTarg->{$gid};
                if (my $genTarg = $gTarg->[0]) {
                    # We were able to anchor these gene hits to the genome
                    while (my ($exType, $genHits) = each %{$genTarg}) {
                        my %mms;
                        map { $mms{$_}++ } values %{$genHits};
                        while (my ($mm, $num) = each %mms) {
                            $norm{$exType}[$mm]{$gid} = $num;
                        }
                    }
                } elsif (my $rnaTarg = $gTarg->[1]) {
                    # This gene is only hit by RNAs
                    # We will simply count each mismatch class as hit once
                    foreach my $mm (0..$#{$rnaTarg}) {
                        if ($rnaTarg->[$mm]) {
                            $norm{'Ex'}[$mm]{$gid} = 1;
                        }
                    }
                } else {
                    # Huh.
                    
                    delete $srcTarg->{$gid};
                }
            }
          #  $self->prebranch( \%norm);
            my @row = ($srchid, $srcid);
            while (my ($eTyp, $mmArr) = each %norm) {
                my $offSet = $eTyp eq 'Ex' ? $Ex0col : $Ex0col + $maxShownMM;
                for my $mm (0..$#{$mmArr}) {
                    my $gidH = $norm{$eTyp}[$mm];
                    next unless ($gidH);
                    my @gids = keys %{$gidH};
                    $row[$offSet + $mm] = \@gids;
                }
            }
            push @{$rows}, \@row;
            # $srcH->{$srcid} = \%norm;
        }
    }
    if ($task) {
        $bar->[1] = 1;
        $bar->[2] = "Done.";
        $task->write(1);
    }
    return $rv;
}

sub filter_flat_hits {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $hits = $args->{HITS};
    my $rv   = $self->_assure_flattened_hits( $hits );
    my $str  = defined $args->{STR} ? $args->{STR} : $args->{STRAND};
    my $scQ  = $args->{SCORE};
    my $rnaQ = &_hash_args(  $args, 'RNA');
    my $gidQ = &_hash_args(  $args, 'GENEIDS');
    my $srcQ = &_hash_args(  $args, 'SRCIDS' );
    my $subQ = &_hash_args(  $args, 'SUBJECTS' );
    my $xRna = &_hash_args(  $args, 'NORNA' );
    my $xGid = &_hash_args(  $args, 'NOGENEIDS' );
    my $xSrc = &_hash_args(  $args, 'NOSRCIDS' );
    my $xSub = &_hash_args(  $args, 'NOSUBJECTS' );
    
    my $subInd  = $self->data_index( $rv, 'SubjectID');
    my $srcInd  = $self->data_index( $rv, 'SourceID' );
    my $detInd  = $self->data_index( $rv, 'GeneDetails' );
    my $filtDet = $str || $scQ || $rnaQ || $gidQ || $xRna || $xSub ? 1 : 0;
    
    if ($srcInd == -1 && $srcQ) {
        $self->msg("[!!]","Can not filter by Source - not in output!");
        $srcQ = undef;
    }
    if ($subInd == -1 && ($subQ || $xSub)) {
        $self->msg("[!!]","Can not filter by Subject - not in output!");
        $subQ = undef;
        $xSub = undef;
    }
    if ($detInd == -1 && $filtDet) {
        $self->msg("[!!]","Can not filter by Gene Details - not in output!");
        $filtDet = 0;
    }


    # NEED TO IMPROVE:
    # Store more detailed filter information under the {filters} hash
    
    my $filters;
    if ($srcQ) {
        $filters ||= {};
        $filters->{src_id} = join(',', keys %{$srcQ});
    }
    
    my @applied;
    push @applied, $self->standard_data_name( 'SourceID' ) if ($srcQ);
    push @applied, $self->standard_data_name( 'SubjectID' ) if ($subQ);
    if ($filtDet) {
        push @applied, $self->standard_data_name( 'GeneStrand' ) if ($str);
        push @applied, $self->standard_data_name( 'AlignmentScore' ) if ($scQ);
        push @applied, $self->standard_data_name( 'RnaAccession' ) if ($rnaQ);
        push @applied, $self->standard_data_name( 'GeneID' ) if ($gidQ);
        push @applied, $self->standard_data_name( 'ExcludedRNA' ) if ($xRna);
        push @applied, $self->standard_data_name( 'ExcludedSubjectID' )
            if ($xSub);
    }

    # Just return if no filters have been requested
    return $rv if ($#applied == -1);
 
    my $pgtm   = $self->pgtm();
    my @subtask = ( ["Filter",0] );
    my $status = {
        task    => "Filter results per user request",
        subtask => \@subtask,
    };
    my $nextPrompt;
    my $promptInterval = 1;
    if ($pgtm) {
        $pgtm->change_task_status($status);
        $nextPrompt ||= time - 1;
    }

    $rv->{filters} = \@applied;
    my $killStr = !$str ? 0 : ($str < 0) ? 1 : -1;

    my $input   = $rv->{rows};
    my $output  = $rv->{rows} = [];
    my $inds    = $#{$input};
    for my $r (0..$inds) {
        my $row = $input->[$r];
        if ($nextPrompt && time >= $nextPrompt) {
            $subtask[0][1] = &_progress_fraction( $r, $inds);
            $subtask[0][2] = "Row ".($r+1);
            $pgtm->change_task_status($status);
            $nextPrompt = time + $promptInterval;
        }
        if ($srcQ) {
            # Skip if a source whitelist is provided and row is not in it
            next unless ($srcQ->{ $row->[$srcInd] || 0});
        } elsif ($xSrc) {
            # Skip if a source blacklist provided and row is in it
            # Blacklist is ignored if whitelist provided
            next if ($srcQ->{ $row->[$srcInd] || 0});
        }
        if ($subQ) {
            # Skip if subject whitelist, and row not in it
            next unless ( $subQ->{ $row->[$subInd] || 0});
        } elsif ($xSub) {
            # Skip if subject blacklist and row is in it:
            # Blacklist is ignored if whitelist provided
            next if ( $xSub->{ $row->[$subInd] || 0});
        }
        unless ($filtDet) {
            # No details, keep the row as-is
            push @{$output}, $row;
            next;
        }
        my $details = $row->[$detInd];
        # Lack of details means we need to reject:
        next unless ($details);
        my $kept;
        while (my ($gid, $strH) = each %{$details}) {
            my $isOk;
            if ($gidQ) {
                # If a gene_id whitelist is provided, skip if not listed:
                next unless ($gidQ->{$gid});
            } elsif ($xGid) {
                # Skip if listed on a gene_id blacklist:
                # Blacklist is ignored if whitelist provided
                next if ($xGid->{$gid});
            }
            if ($str) {
                # Strand filter provided
                if ($strH->{$str}) {
                    # This gene is on the correct strand
                    # Eliminate any hits on the other strand
                    delete $strH->{$killStr};
                } else {
                    # Wrong strand. Ignore hit.
                    next;
                }
            }
            if ($rnaQ || $scQ || $xRna) {
                # We are filtering by RNA name and/or score
                my $keptRNA;
                $isOk = 0;
                while (my ($kstr, $etH) = each %{$strH}) {
                    while (my ($et, $rnaDat) = each %{$etH}) {
                        foreach my $rd (@{$rnaDat}) {
                            my $accV = $rd->[0] || "";
                            my $accU = $accV; $accU =~ s/\.\d+$//;
                            if ($rnaQ) {
                                # Skip if we have requested an RNA accession
                                # and this row does not match
                                next unless ($rnaQ->{$accV} || $rnaQ->{$accU});
                            } elsif ($xRna) {
                                # RNA blacklist provided
                                # Blacklist is ignored if whitelist provided
                                next if ($xRna->{$accV} || $xRna->{$accU});
                            }
                            if ($scQ) {
                                # Skip if a minimum score is requested and
                                # the row falls below
                                my $sc = $rd->[3] || 0;
                                next if ($sc < $scQ);
                            }
                            # The row survived the filters
                            $keptRNA ||= {};
                            push @{$keptRNA->{$kstr}{$et}}, $rd;
                        }
                    }
                }
                if ($keptRNA) {
                    $details->{$gid} = $keptRNA;
                    $isOk = 1;
                }
            }
            if ($gidQ && !$isOk) {
                $isOk = $gidQ->{$gid} ? 1 : 0;
            }
            next if (defined $isOk && !$isOk);
            # All filters passed (Gene and RNA considered independently)
            $kept ||= {};
            $kept->{$gid} = $details->{$gid};
        }
        if ($kept) {
            # At least one gene/RNA passed the filter(s)
            push @{$output}, $row;
            $row->[$detInd] = $kept;
        }
    }
    # warn "$#{$rv} => $#filtered    ($str || $rnaQ  || $gidQ || $scQ);";
    if ($pgtm) {
        $subtask[0][1] = 1;
        $subtask[0][2] = "Done";
        $pgtm->change_task_status($status);
    }
    return $rv;
}

sub source_id_callback {
    my $self    = shift;
    my ($dataReq, $colReq) = @_;
    my $data    = $self->_assure_flattened_hits( $dataReq );
    my $srcCol  = 'SourceID';
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $colReq, $srcCol, [qw(SourcePath GenomeBuild Species)],
          sub {
              my $col = shift;
              if ($col eq 'SourcePath') {
                  return sub { return shift->{Path} };
              } elsif ($col eq 'GenomeBuild') {
                  return sub { return shift->{GenomeBuild} };
              } elsif ($col eq 'Species') {
                  return sub { return shift->{Species} };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $rows     = $data->{rows};
    
    my $sbjSTH   = $self->_get_subject_sth();
    foreach my $row (@{$rows}) {
        my $srcid = $row->[$srcInd];
        next unless ($srcid);
        my $info  = $self->source_info( $srcid );
        while (my ($ind, $cb) = each %{$needed}) {
            $row->[$ind] = &{$cb}( $info );
        }
    }
    $self->bench_end();
    return;
}

sub overview_count_callback {
    my $self    = shift;
    my ($dataReq, $colReq) = @_;
    my $data    = $self->_assure_flattened_hits( $dataReq );
    my ($needed) = $self->_find_needed_expansion_columns
        ( $data, $colReq, undef, [qw(Ex0 Ex1 Ex2 Int0 Int1 Int2)],
          sub {
              my ($col, $data) = @_;
              my $srcCol = 'GeneIDs'.$col;
              my $ind    = $self->data_index( $data, $srcCol );
              if ($ind != -1) {
                  return sub {
                      my $row = shift;
                      if (my $arr = $row->[$ind]) {
                          return $#{$arr} + 1;
                      } else {
                          return "";
                      }
                  };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $rows     = $data->{rows};
    foreach my $row (@{$rows}) {
        while (my ($ind, $cb) = each %{$needed}) {
            $row->[$ind] = &{$cb}( $row );
        }
    }
    $self->bench_end();
    return;
}

sub subject_id_callback {
    my $self    = shift;
    my ($dataReq, $colReq) = @_;
    my $data    = $self->_assure_flattened_hits( $dataReq );
    my $srcCol  = 'SubjectID';
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $colReq, $srcCol, [qw(SubjectAccession SubjectLength)],
          sub {
              my $col = shift;
              if ($col eq 'SubjectAccession') {
                  return sub { return shift->[0] };
              } elsif ($col eq 'SubjectLength') {
                  return sub { return shift->[2] };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $rows     = $data->{rows};
    my $sbjSTH   = $self->_get_subject_sth();
    foreach my $row (@{$rows}) {
        my $sid   = $row->[$srcInd];
        $sbjSTH->execute( $sid );
        my $sRows = $sbjSTH->fetchall_arrayref();
        if ($#{$sRows} == 0) {
            # Unique DB entry found
            while (my ($ind, $cb) = each %{$needed}) {
                $row->[$ind] = &{$cb}( $sRows );
            }
        } else {
            $self->msg_once("[??]","Failed to find $srcCol = $sid");
        }
    }
    $self->bench_end();
    return;
}

sub genomic_footprint_callback {
    my $self    = shift;
    my ($dataReq, $colReq) = @_;
    my $data    = $self->_assure_flattened_hits( $dataReq );
    my $srcCol  = 'GenomicFootprint';
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $colReq, $srcCol, [qw(GenomeAccession GenomeStrand GenomeStart GenomeEnd)],
          sub {
              my $col = shift;
              if ($col eq 'GenomeAccession') {
                  return sub { return shift->[0] };
              } elsif ($col eq 'GenomeStrand') {
                  return sub { return shift->[1] };
              } elsif ($col eq 'GenomeStart') {
                  return sub { return shift->[2] };
              } elsif ($col eq 'GenomeEnd') {
                  return sub { return shift->[3] };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $rows     = $data->{rows};
    foreach my $row (@{$rows}) {
        my $fp   = $row->[$srcInd];
        next unless ($fp);
        my @gDat;
        # Parse the footprint to extract coordinate information
        if ($fp =~ /^([^:]+):(.+)$/) {
            $gDat[0] = $1; # Accession
            my $ses  = $2; # Footprint coordinates
            if ($ses =~ /(.+)\[([\+\-]?1)\]$/) {
                # Strand token on end: [1], [+1], [-1]
                $ses = $1;
                $gDat[1] = 0 + $2;
            } else {
                $gDat[1] = 1;
            }
            my @bits;
            foreach my $pos (split(',', $ses)) {
                if ($pos =~ /^(\d+)$/) {
                    # Single base
                    push @bits, $pos;
                } elsif ($pos =~ /^(\d+)(\.\.|\^)(\d+)$/) {
                    push @bits, ($1, $3);
                }
            }
            @bits = sort {$a <=> $b} @bits;
            $gDat[2] = $bits[0];
            $gDat[3] = $bits[-1];
        }
        while (my ($ind, $cb) = each %{$needed}) {
            $row->[$ind] = &{$cb}( \@gDat );
        }
    }
    $self->bench_end();
    return;
}

sub ambiguous_count_callback {
    my $self = shift;
    my ($dataReq, $colReq) = @_;
    my $data = $self->_assure_flattened_hits( $dataReq );
    my $srcCol  = 'SubjectSeq';
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $colReq, $srcCol, [qw(AmbiguousCount)], sub {
              my $col = shift;
              if ($col eq 'AmbiguousCount') {
                  return sub {
                      my $seq = shift;
                      if ($seq) {
                          $seq = uc($seq);
                          $seq =~ s/[^KMRYSWBVHDXN]//g;
                          return CORE::length($seq);
                      }
                      return undef;
                  };
              }
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $rows     = $data->{rows};
    foreach my $row (@{$rows}) {
        while (my ($ind, $cb) = each %{$needed}) {
            $row->[$ind] = &{$cb}( $row->[$srcInd] );
        }
    }
    $self->bench_end();
    return;
    
}

sub gene_details_callback {
    my $self    = shift;
    my ($dataReq, $colReq) = @_;
    my $data    = $self->_assure_flattened_hits( $dataReq );
    my $srcCol  = 'GeneDetails';
    my ($needed, $srcInd) = $self->_find_needed_expansion_columns
        ( $data, $colReq, $srcCol, [qw(HitType SubType GeneStrand GeneAccession GeneSymbol HumanSymbol HumanAccession HumanScore GeneDescription RnaFootprint ExonNumber RnaCount RnaStart RnaEnd Notes ExonIntron GeneID)],
          sub {
              my $col = shift;
              if ($col eq 'HitType') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{typ};
                  };
              } elsif ($col eq 'SubType') {
                  # Pulls out CDS, UTR5, etc
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      my $typ = $rna->{typ};
                      if ($#{$typ} != -1 && $typ->[0] =~ / (.+)/) {
                          # warn "$typ->[0] = $1\n";
                          return $1;
                      } else {
                          return "";
                      }
                  };
              } elsif ($col eq 'GeneStrand') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{str};
                  };
              } elsif ($col eq 'GeneAccession') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{acc};
                  };
              } elsif ($col eq 'GeneSymbol') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{sym};
                  };
              } elsif ($col eq 'GeneID') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{geneid};
                  };
              } elsif ($col eq 'HumanSymbol') {
                  return sub {
                      my ($self, $info) = @_;
                      # die $self->branch( -ref => $info, -skipkey => [qw(rna eth)]);
                      return $info->{humansym};
                  };
              } elsif ($col eq 'HumanAccession') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{humanacc};
                  };
              } elsif ($col eq 'HumanScore') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{humanscore};
                  };
              } elsif ($col eq 'GeneDescription') {
                  return sub {
                      my ($self, $info) = @_;
                      return $info->{desc};
                  };
              } elsif ($col eq 'RnaAccession') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{acc};
                  };
              } elsif ($col eq 'RnaFootprint') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{fp};
                  };
              } elsif ($col eq 'ExonNumber') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{ex};
                  };
              } elsif ($col eq 'RnaCount') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $#{$rna->{acc} || []} + 1;
                  };
              } elsif ($col eq 'RnaStart') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{s};
                  };
              } elsif ($col eq 'RnaEnd') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{e};
                  };
              } elsif ($col eq 'Notes') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      return $rna->{note};
                  };
              } elsif ($col eq 'ExonIntron') {
                  return sub {
                      my ($self, $info) = @_;
                      my $rna = $self->_rna_details( $info );
                      my $typ = lc($rna->{typ}[0] || "");
                      $typ =~ s/ .+//;
                      my $exn = $rna->{ex}[0];
                      #die $self->branch($rna);
                      my $rv;
                      if (!$typ || !defined $exn || $exn eq '') {
                          return "";
                      } elsif ($typ eq 'intron') {
                          $rv = sprintf("intron.%dto%d", $exn, $exn+1);
                      } elsif ($typ eq 'exon' || 
                               $typ eq 'firstexon' || $typ eq 'lastexon') {
                          $rv = sprintf("exon.%d", $exn);
                      } elsif ($typ eq 'exonexon') {
                          $rv = sprintf("exon.%d|exon.%d", $exn, $exn + 1);
                      } elsif ($typ =~ /^splice(\S+)/) {
                          my $side = $1;
                          if ($side eq 'donor') {
                              $rv = sprintf("exon%d|intron.%dto%d",
                                             $exn, $exn, $exn+1);
                          } elsif ($side eq 'acceptor') {
                              $rv = sprintf("intron.%dto%d|exon%d",
                                             $exn-1, $exn, $exn);
                          } else {
                              $rv = sprintf("exon%d|intron.?", $exn);
                          }
                      } elsif ($typ eq 'crypticexon') {
                          $rv = sprintf("exon.%d?", $exn);
                      }
                      $rv ||= sprintf("exon?.%d?", $exn);
                      # warn "$typ [$exn] = $rv\n";
                  };
              }
              return undef;
          });
    # If we already have all columns just return:
    return unless ($needed);
    
    # Iterate over the rows and apply callbacks
    $self->bench_start();
    my $pgtm   = $self->pgtm();
    my @subtask = ( ["Gene",0] );
    my $status = {
        task    => "Add gene details",
        subtask => \@subtask,
    };
    my $nextPrompt;
    my $promptInterval = 1;
    if ($pgtm) {
        $pgtm->change_task_status($status);
        $nextPrompt ||= time - 1;
    }
    my $input   = $data->{rows};
    my $output  = $data->{rows} = [];
    my $inds    = $#{$input};
    for my $ri (0..$inds) {
        my $row     = $input->[$ri];
        my $details = $row->[$srcInd] || {};
        my @gids    = keys %{$details};
        foreach my $gid (@gids) {
            my @gStrs = keys %{$details->{$gid}};
            foreach my $geneStr (@gStrs) {
                my $etH = $details->{$gid}{$geneStr};
                my $geneRow = [ @{$row} ];
                # Make a surgical slice of the details hash:
                $geneRow->[$srcInd] = { $gid => { $geneStr => $etH } };
                push @{$output}, $geneRow;
                my $info = {
                    str    => $geneStr,
                    eth    => $etH,
                    geneid => $gid,
                };
                if ($gid) {
                    my ($acc, $sym, $desc, $tax, $ali, $hacc, $hsym, $hsc) =
                        $self->gene_info( $gid );
                    if ($acc) {
                        $info->{acc}    = $acc;
                        $info->{sym}    = $sym;
                        $info->{desc}   = $desc;
                        # Why no taxa or aliases?
                        $info->{humanacc}   = $hacc;
                        $info->{humansym}   = $hsym;
                        $info->{humanscore} = $hsc;
                    } else {
                        $self->msg_once("[??]","Failed to find gene_id = $gid");
                    }
                }
                while (my ($ind, $cb) = each %{$needed}) {
                    $geneRow->[$ind] = &{$cb}( $self, $info );
                }
                if ($nextPrompt && time >= $nextPrompt) {
                    $subtask[0][1] = &_progress_fraction( $ri, $inds);
                    $subtask[0][2] = $info->{sym} || 'NoSymbol';
                    $pgtm->change_task_status($status);
                    $nextPrompt = time + $promptInterval;
                }
            }
        }
        # die $self->branch($row);
    }
    if ($pgtm) {
        $subtask[0][1] = 1;
        $subtask[0][2] = "Done";
        $pgtm->change_task_status($status);
    }
}

sub _rna_details {
    my $self = shift;
    my $info = shift;
    unless ($info->{rna}) {
        my $etH    = $info->{eth};
        my @types  = $self->_sorted_types( [ keys %{$etH} ] );
        my (@accs, @rfps, @starts, @ends, @exons, %notes);
        foreach my $type (@types) {
            # Sort the entries by exon number
            my @rdats;
            foreach my $rdat (@{$etH->{$type}}) {
                # Sort numerically, accomodating intron notation of 17^18 :
                my $ex = $type eq 'RNA' ? -1 :
                    ($rdat->[1] && $rdat->[1] =~ /^(\d+)/) ? $1 : 0;
                push @rdats, [$rdat, $ex];
            }
            foreach my $rdata (sort { $a->[1] <=> $b->[1] } @rdats) {
                my ($rdat, $ex) = @{$rdata};
                my $acc = $rdat->[0] || "";
                push @accs, $acc || '?';
                my $rfp = $acc;
                my ($s, $e);
                if ($acc) {
                    if (my $ft = $rdat->[2]) {
                        $rfp .= ":$ft";
                        my @se = sort { $a <=> $b } split(/[^\d]+/, $ft);
                        $s = $se[0];
                        $e = $se[-1];
                    }
                    $rfp .= '[-1]' if ($info->{str} < 0);
                }
                push @rfps, $rfp || '?';
                push @exons, $ex < 1 ? '' : $ex;
                push @starts, defined $s ? $s : '?';
                push @ends, $e || '?';
                if (my $note = $rdat->[4]) {
                    $notes{$note} = 1;
                }
            }
        }
        $info->{rna} = {
            acc  => \@accs,
            fp   => \@rfps,
            ex   => \@exons,
            s    => \@starts,
            e    => \@ends,
            typ  => \@types,
            note => join('; ', sort keys %notes ),
        };
    }
    return $info->{rna};
}

sub _find_needed_expansion_columns {
    my $self = shift;
    $self->bench_start();
    my ($data, $request, $srcCol, $ac, $logic, $depends) = @_;
    my $srcInd;
    if ($srcCol) {
        $srcInd  = $self->data_index( $data, $srcCol );
        if ($srcInd == -1) {
            $self->msg("[?]","Can not expand $srcCol - source column absent");
            $self->bench_end();
            return ();
        }
    }
    my @allCols = map { $self->standard_data_name($_) } @{$ac};
    # warn join(" + ", @{$request});
    $request  ||= \@allCols;
    my %wanted  = map { $self->standard_data_name($_) => 1 } @{$request};
    if ($depends) {
        foreach my $col (@{$request}) {
            if (my $parents = $depends->{$col}) {
                map { $wanted{$self->standard_data_name($_)} = 1 } @{$parents};
            }
        }
    }
    my ($needed, %clear);
    foreach my $col (@allCols) {
        next unless ($wanted{$col});
        # Do we already have it?
        my $ind = $self->data_index( $data, $col );
        unless ($ind == -1) {
            # The column is already in the header
            unless ($data->{need}{$col}) {
                next;
            }
            # There was a flag that we should still calculate the column
            $clear{$col}++;
        }
          
        # Not there, need to add it
        if (my $cb = &{$logic}( $col, $data )) {
            $needed ||= {};
            my $newInd = $self->data_index( $data, $col, 'add' );
            $needed->{$newInd} = $cb;
        } else {
            $srcCol ||= '?';
            $self->err("No callback defined to expand $srcCol to $col");
        }
    }
    foreach my $done (keys %clear) {
        delete $data->{need}{$done};
        # $self->msg("[+]","Ready to add $done");
    }
    $self->bench_end();
    return ($needed, $srcInd);
}

sub _subject_base_class_hash {
    my $self = shift;
    unless ($self->{SBCH}) {
        my $sbch = $self->{SBCH} = {};
        # Upper case bases are normal, unclassed
        map { $sbch->{$_} = '' } qw(A C G T);
        # Upper case ambiguities are special:
        map { $sbch->{$_} = 'mmamb' } split('', 'KMRYSWBVHDXN-');
        # Lower case nucleotide letters are mismatched:
        map { $sbch->{$_} = 'mmseq' } split('', 'acgtkmryswbvhdxn');
        # Gap
        $sbch->{'-'} = 'mmgap';
    }
    return $self->{SBCH};
}

sub data_to_sencha {
    my $self = shift;
    my @passed = @_;
    unshift @passed, "-data" if ($#passed == 0);
    my $args = $self->parseparams( @passed );
    my $data = $args->{DATA};
    my $rv = { type => 'null', rows => [] };
    if ($data) {
        foreach my $k (keys %{$data}) {
            next if ($k eq 'rows' || $k eq 'cols');
            $rv->{$k} = $data->{$k};
        }
        my @cols;
        foreach my $stnd (@{$data->{cols} || []}) {
            my $ext = $self->get_metadata('columns', $stnd, 'extfield');
            push @cols, $ext || $stnd;
        }
        my $inds = $args->{INDS} || [0..$#cols];
        $rv->{cols}  = [ map { $cols[$_] } @{$inds} ];
        $rv->{names} = [ map { $data->{cols}[$_] } @{$inds} ];
        my $reformat = $args->{REFORMAT} || {};
        foreach my $row (@{$data->{rows} || []}) {
            my %hash;
            foreach my $c (@{$inds}) {
                my $val = $row->[$c];
                if (my $cb = $reformat->{$c}) {
                    $val = &{$cb}( $self, $val );
                }
                $hash{$cols[$c]} = $val;
            }
            push @{$rv->{rows}}, \%hash;
        }
    }
    return $rv;
}

sub _assure_flattened_hits {
    my $self = shift;
    my $hitReq = shift;
    return $self->flatten_hits() unless ($hitReq);
    $self->bench_start();
    my $flat;
    if (my $reqRef = ref($hitReq)) {
        if ($reqRef eq 'HASH') {
            if (exists $hitReq->{flat}) {
                # Already flattened
                $flat = $hitReq;
            } else {
                $flat = $self->flatten_hits( $hitReq );
            }
        }
    }
    unless ($flat) {
        $self->msg("[!]","Not sure how to flatten hit request '$hitReq'");
        $flat = $self->flatten_hits();
    }
    $self->bench_end();
    return $flat;
}

sub subject_position_to_genes {
    my $self = shift;
    $self->bench_start();
    my ($sid, $s, $e, $hstr) = @_;
    # $sid  = subj_id of subject
    # $s    = query start relative to subject
    # $e    = query end
    # $hstr = query strand relative to subject
    
    my $gsth = $self->dbh->named_sth
        ('Get gene locations for subject position');
    $gsth->execute( $sid, $e, $s );
    my $grows = $gsth->fetchall_arrayref();
    my %genes;
    for my $i (0..$#{$grows}) {
        my ($locid, $gid, $gstr, $json) = @{$grows->[$i]};
        # warn join(" , ", map { defined $_ ? $_ : 'UNDEF' }  @{$grows->[$i]});
        # $gid  = gene_id for the gene
        # $gstr = gene strand relative to subject
        # $json = JSON text describing RNA arangement
        my $types = $self->exons_to_hit_type( $s, $e, $gstr, $json, $sid );
        unless ($types) {
            # This should not happen! gene_loc should only have recovered
            # hits that had at least one overlap!
            my ($acc, $sym, $desc) = $self->gene_info( $gid );
            my $sname = $self->subject_name($sid);
            $sym = $sym ? " [$sym]" : "";
            $self->err("Failed to find overlap for [$s..$e] on sbj_id = $sid",
                       $sname, "gene_id = $gid $acc$sym\n");
            next;
        }
        my $str   = $gstr * $hstr;
        # $str = strand of gene relative to query
        while (my ($type, $rnaArr) = each %{$types}) {
            push @{$genes{$gid}{$str}{$type}}, @{$rnaArr};
        }
    }
    # We will keep the data structured
    $self->bench_end();
    return \%genes;

    
    my @rv;
    foreach my $gid (sort { $a <=> $b } keys %genes) {
        foreach my $str (sort { $b <=> $a } keys %{$genes{$gid}}) {
            my %ut = map { $_ => 1 } @{$genes{$gid}{$str}};
            my @types = $self->_sorted_types( [ keys %ut ] );
            my @etids = map { $self->exon_type_id( $_ ) } @types;
            push @rv, [$gid, $str, @etids];
        }
    }
    $self->bench_end();
    return \@rv;
}

sub _sorted_types {
    my $self = shift;
    my $arr  = shift || [];
    my @rv;
    foreach my $t (@{$arr}) {
        my $tp = $typePriority->{$t};
        unless ($tp) {
            if ($t =~ /^(\S+) /) {
                # "Exon CDS" or "Exon UTR5|CDS"
                my $base = $1;
                $typePriority->{$t} = $typePriority->{$base};
            }
            $tp = $typePriority->{$t} ||= 999;
        }
        push @rv, $t;
    }
    return sort { $typePriority->{$a} <=> $typePriority->{$b} } @rv;
}

sub exon_type_id {
    my $self = shift;
    my $req  = shift || "";
    unless (defined $self->{EXONTYPEIDS}{$req}) {
        my $sth  = $self->dbh->named_sth('Get/Make exon type id');
        my $etid = $sth->get_single_value($req) || 0;
        $self->{EXONTYPEIDS}{$req} = $etid;
    }
    return $self->{EXONTYPEIDS}{$req};
}

sub et_id_to_type {
    my $self = shift;
    my $req  = shift || 0;
    unless ($self->{ETIDNAMES}{$req}) {
        my $sth  = $self->dbh->named_sth('Get exon type name for id');
        my $etid = $sth->get_single_value($req) || "Unknown";
        $self->{ETIDNAMES}{$req} = $etid;
    }
    return $self->{ETIDNAMES}{$req};
}

sub _decode_json_details {
    my $self = shift;
    my $json = shift;
    return ( { NoGeneInfo  => [ ["", "", "", 0, 
                                 "No exon structure provided for gene"] ] },
             "") unless ($json);
    my $details;
    eval {
        $details  = decode_json( $json );
    };
    return ($details, undef) if ($details);
    return ( { CodeError => [ ["", "", "", 0, 
                               "Malformed gene exon structure data"] ] },
             "Malformed gene data structure");
}

sub exons_to_hit_type {
    my $self = shift;
    $self->bench_start();
    my ($s, $e, $gstr, $json, $sid) = @_;
    # $s    = query start relative to subject
    # $e    = query end
    # $gstr = Strand of the gene on the subject
    # $json = JSON text describing RNA arangement
    # $sid  = subj_id of subject sequence
    my ($details, $err) = $self->_decode_json_details( $json );
    if (defined $err) {
        $self->msg_once("[!]", $err, "'$json'") if ($err);
        $self->bench_end();
        return $details;
    }

    my $isNotGenomic = $details->{type};
    my $isRC = ($gstr < 0) ? 1 : 0;
    my $types;
    foreach my $rdat (@{$details->{rna}}) {
        my ($racc, $sc, $exons, $enum, $notes, $cds) = @{$rdat};
        unless ($e > $exons->[0][0] && $s < $exons->[-1][1]) {
            # The [s,e] range does not cover this transcript
            next;
        }
        # Our [s,e] range is overlapping this transcript
        for my $ei (0..$#{$exons}) {
            my ($l, $r, $rs, $re) = @{$exons->[$ei]};
            # $l   = Left flanking position of exon
            # $r   = Right flanking position of exon
            #        So an exon from gDNA 100..200 -> $l = 99, $r = 201
            # $rs  = RNA start coordintate (RNA coordinate space)
            # $re  = RNA end coordinate
            #        So $rs = 451 = This exon starts at RNA base 451
            if ($e <= $l) {
                # Our query end position falls before this exon begins
                # We have explored past our query - we do not need to
                # consider any other exons
                if ($ei) {
                    # This is not the first exon = so we must be in an intron
                    my $rfp;
                    if ($rs) {
                        # We can calculate an RNA footprint, though it will
                        # just be "between here and there"
                        my $prior = $exons->[$ei-1];
                        if ($isRC) {
                            $rfp = sprintf("%d^%d", $re, $prior->[2]);
                        } else {
                            $rfp = sprintf("%d^%d", $prior->[3], $rs);
                        }
                    }
                    # Get the intron number:
                    my $inum = $isRC ? $#{$exons} + 1 - $ei : $ei;
                    $types ||= {};
                    push @{$types->{Intron}}, [$racc, $inum, $rfp, $sc];
                }
                last;
            }
            if ($s >= $r) {
                # Our query start position is still after this exon
                # We have not yet found an exon on or past the query,
                # keep looking
                next;
            }

            # We are overlapping this exon
            
            # Calculate the exon number:
            my $inum = $isRC ? $#{$exons} + 1 - $ei : $ei + 1;
            # Calculate the offset from the LEFT of the subject exon (0 based):
            my $off = $s - $l - 1;
            # ---------[$l                        $r]-----------
            #                  |$s..$e|
            #          <-$off->
            
            # The length of our query segment, minus one
            my $lenM1 = $e - $s;
            my ($type, $rfp);
            if ($s > $l && $e < $r) {
                # Our [s,e] range is FULLY contained
                if ($isNotGenomic) {
                    # This is a segment of RNA not aligned to genome (1)
                    # Or a raw RNA (2)
                    $type = ($isNotGenomic eq 'ExEx') ? 'CrypticExon' : 'RNA';
                } else {
                    # Standard exon
                    $type = 'Exon';
                }
                $rfp = &_rna_footchunk( $rs, $re, $off, $lenM1, $isRC) if ($rs);
            } else {
                # Partially overlapping - at least some of [s,e]
                # is NOT in the exon boundary. Tease out what it is...
                if ($isNotGenomic) {
                    if ($isNotGenomic eq 'ExEx') {
                        # Special Exon-Exon boundary sequence (spliced RNA)
                        $type = 'ExonExon';
                        # Here, the RNA positions are also left/right values
                        # So we should 'shrink' the $rs/$re coordinates
                        $rfp = &_rna_footchunk( $rs + 1, $re - 1, $off,
                                                $lenM1, $isRC) if ($rs);
                        # Need to extract the exon number
                        $inum = $enum || 0;
                    } else {
                        # If we are here, this really should be a transcript,
                        # which should be a monolithic piece of RNA.
                        # The strand should ALWAYS be +1 in such cases, but
                        # we will treat it as varaible to be safe.
                        my ($trimS, $trimE) = ($s, $e);
                        if ($s <= $l && $e >= $r) {
                            # We overhang BOTH sides
                            # Do not expect this to happen often, but plausible
                            # with miRNAs
                            $type = 'FullSpan';
                            ($trimS, $trimE) = ($l + 1, $r - 1);
                        } else {
                            my $side = 1; # 1 = on right side of subject
                            if ($s <= $l) {
                                $side = -1; # -1 = on left side of subject
                                # Need to flush the query to the left side:
                                $trimS = $l + 1;
                            } elsif ($e >= $r) {
                                # We should always end up here if the prior
                                # conditional fails
                                # Flush the right query coordinate:
                                $trimE = $r - 1;
                            }
                            
                            $side   *= ($gstr || 1); # Now -1 = RNA start
                            if ($side < 0) {
                                # The 'start' of the RNA'
                                $type = 'FirstExon';
                            } else {
                                $type = 'LastExon';
                            }
                            if ($rs) {
                                # Recalculate the exon-centric offset:
                                my $trOff   = $trimS - $l - 1;
                                my $trLenM1 = $trimE - $trimS;
                                # The RNA footprint WILL NOT COVER the
                                # entire query - there should be an overhang
                                $rfp = &_rna_footchunk
                                    ( $rs, $re, $trOff, $trLenM1, $isRC);
                            }
                        }
                    }
                } elsif ($s <= $l && $e >= $r) {
                    # We are hanging off BOTH sides of the exon
                    $type = 'FullSpan';
                    if ($rs) {
                        my $diff = $l - $s + 1;
                        $rfp = &_rna_footchunk($rs, $re, $off + $diff,
                                               $r - $l - 1, $isRC);
                    }
                } elsif ($s <= $l) {
                    # Hangs off left side of the exon
                    if ($ei) {
                        $type = $isRC ? 'SpliceDonor' : 'SpliceAcceptor';

                        # $type = 'SpliceJunction';
                        # We are on two exon numbers now. We will join
                        # with a decimal to still allow sorting. A
                        # leading zero will indicate the side that is
                        # intronic.
                        
                        #$inum = $isRC ?
                        #    sprintf("%d.0%d", $inum, $inum ) :
                        #    sprintf("0%d.%d", $inum - 1, $inum) ;
                            
                        
                    } else {
                        # Left-most exon
                        $type = $isRC ? 'LastExon' : 'FirstExon';
                    }
                    if ($rs) {
                        my $diff = $l - $s + 12;
                        $rfp = &_rna_footchunk($rs, $re, $off + $diff,
                                               $lenM1 - $diff, $isRC);
                    }
                } elsif ($e >= $r) {
                    # Hangs off right side
                    # Special case handle if this is the last exon
                    if ($ei == $#{$exons}) {
                        # Right-most exon
                        $type = $isRC ? 'FirstExon' : 'LastExon';
                    } else {
                        $type = $isRC ? 'SpliceAcceptor' : 'SpliceDonor';
                        
                        #$type = 'SpliceJunction';
                        # We are on two exon numbers now. We will join
                        # with a decimal to still allow sorting. A
                        # leading zero will indicate the side that is
                        # intronic.
                        
                        #$inum = $isRC ?
                        #    sprintf("0%d.%d", $inum - 1, $inum) :
                        #    sprintf("%d.0%d", $inum, $inum);
                             
                    }
                    if ($rs) {
                        my $diff = $e - $r + 1;
                        $rfp = &_rna_footchunk($rs, $re, $off,
                                               $lenM1 - $diff, $isRC);
                    }
                } else {
                    # Uh. Should not be here.
                    $type = "Unknown";
                }
            }
            if ($cds && $type && $type ne 'Intron') {
                # The ORF start/end range is defined
                my ($cs, $ce) = @{$cds};
                my @cbits;
                if ($s <= $ce && $e >= $cs) {
                    # The CDS is covered by at least one base
                    push @cbits, "CDS";
                }
                if ($s < $cs) {
                    # The right UTR is covered
                    if ($isRC) {
                        push @cbits, "UTR3";
                    } else {
                        unshift @cbits, "UTR5";
                    }
                }
                if ($e > $ce) {
                    # The left UTR is covered
                    if ($isRC) {
                        unshift @cbits, "UTR5";
                    } else {
                        push @cbits, "UTR3";
                    }
                }
                $type .= ' ' .join('|', @cbits);
            }
            $types ||= {};
            push @{$types->{$type}}, [$racc, $inum, $rfp, $sc, $notes];
            last;
        }
    }
    unless ($types) {
        # Huh. We really should have found SOMETHING.
        # There is one valid exception - a genomic chunk that contains
        # MORE THAN ONE COPY of a gene. Example:
        
        # 14.Rnor_5:5556762-6520637.FR  uncharacterized LOC103690086
        #  ["XR_595692.1",100,[[0,887,476,1361],
        #    [2752,2839,390,475],[3661,4051,1,389]]],
        #  ["XR_595692.1",100,[[71719,72606,476,1361],
        #    [74471,74558,390,475],[75380,75770,1,389]]] ]}
        
        # Two non-overlapping copies of a transcript, both with perfect (100%)
        # alignments to neighboring segments on Rnor_5 Chr14. They are pulled
        # into the same segment because of one or more overlapping other genes

        my %sides;
        foreach my $rdat (@{$details->{rna}}) {
            my $exons = $rdat->[2];
            if ($e <= $exons->[0][0]) {
                # Our query is to the left of this RNA
                $sides{L} = 1;
            } elsif ($s >= $exons->[-1][1]) {
                # Query to the right
                $sides{R} = 1;
            } else {
                # THIS SHOULD NOT HAPPEN
                # If we are neither to the left nor right, we should have
                # caught that above.
                $types = { CodeError => [ ["", "", "", 0, "Unexpected result in exons_to_hit_type"] ] };
                last;
            }
        }
        if ($sides{L} && $sides{R}) {
            # Yup, there it is: we are both to the left AND to the right
            # of this particular gene. So we are falling between two copies
            # I do not think this will generally be 'notable', so we will
            # simply set the $types array as an empty hash structure, which
            # will suppress errors later on
            $types ||= {};
        } elsif (!$types) {
            # Huh. Assuming we are running this method from 
            # subject_position_to_genes(), we should have found SOMETHING

            # Leave $types null as a flag
        }
    }
    $self->bench_end();
    return $types;
}

sub _rna_footchunk {
    my ($rs, $re, $off, $lenM1, $isRC) = @_;
    # $rs    = low RNA coordinate for exon
    # $re    = high RNA coordinate for exon
    # $off   = left query offset of the hit with regard to the exon
    # $lenM1 = length of query, minus 1
    if ($isRC) {
        # The exon is reversed
        my $rpos = $re - $off; 
        return sprintf("%d..%d", $rpos - $lenM1, $rpos);
    } else {
        my $rpos = $rs + $off; # This is the START coordinate
        return sprintf("%d..%d", $rpos, $rpos + $lenM1);
    }
}

sub _array_args {
    my ($args, $param) = @_;
    if (my $val = $args->{$param}) {
        if (ref($val)) {
            return $#{$val} == -1 ? undef : $val;
        }
        return [$val];
    }
    return undef;
}

sub _hash_args {
    my $arr = &_array_args( @_ );
    return undef unless ($arr);
    return { map { $_ => 1 } @{$arr} };
}

sub _search_hit_where_clause {
    my $self = shift;
    $self->bench_start();
    
    my $args = $self->parseparams( @_ );
    my (@where, @binds, @what, $filter);
    # $needsFilter indicates if finer-grained filtering will be needed on
    # the recovered rows
    my @needsFilter;
    if (my $srchid = $args->{SRCHID} || $args->{SEARCHID}) {
        push @where, 'sh.srch_id = ?';
        push @binds, $srchid;
        my ($seq, $name, $descr, $date, $mm) = $self->search_details( $srchid );
        push @what, "Seq: $seq";
        $filter ||= [];
        push @{$filter}, 'SearchID';
    }
    my $mm = $args->{MM}; $mm = $args->{MISMATCH} unless (defined $mm);
    if (defined $mm) {
        push @where, 'sh.mm <= ?';
        push @binds, $mm;
        push @what, "MM <= $mm";
        $filter ||= [];
        push @{$filter}, 'MisMatch';
    }

    if (my $str = $self->normalize_strand
        ($args->{GENESTRAND} || $args->{STRAND} || $args->{STR})) {
        push @where, "sh.genestrand != ?";
        push @binds, $str == -1 ? 1 : 1;
        push @needsFilter, 'strand';
        $filter ||= [];
        push @{$filter}, 'GeneStrand';
    }


    # The positive selective requests will be OR'ed into @geneBits
    my (@geneBits, @geneBinds);
    if (my $rnaQ = &_array_args( $args, 'RNA' )) {
        foreach my $rq (@{$rnaQ}) {
            push @geneBits, "upper(sh.details) LIKE ?";
            push @geneBinds, uc("%${rq}%");
        }
        push @what, "RNA = ".join('/', @{$rnaQ});
        push @needsFilter, 'rnaid';
        $filter ||= [];
        push @{$filter}, 'RnaAccession';
    }
    my $subQ   = &_array_args( $args, 'SUBJECTS' );
    my $noSubQ = &_array_args( $args, 'NOSUBJECTS' );
    if ($subQ) {
        $filter ||= [];
        push @{$filter}, 'SubjectID';
    }
    if (my $geneQ = &_array_args( $args, 'GENEIDS' )) {
        my $findSbj = $self->dbh->named_sth('Get subjects for gene_id');
        foreach my $gq (@{$geneQ}) {
            next unless ($gq =~ /^\d+$/);
            my @sids = $findSbj->get_array_for_field( $gq );
            # Only consider if at least one subject has this gene on it:
            next if ($#sids == -1);
            push @geneBits, "sh.details LIKE ?";
            push @geneBinds, uc("%\"${gq}\"%");
            $subQ ||= [];
            push @{$subQ}, @sids;
        }
        push @what, "gene_id = ".join('/', @{$geneQ});
        $filter ||= [];
        push @{$filter}, 'GeneID';
        push @needsFilter, 'geneid';
        # $args->msg_once("[TODO]","Extend gene_id filters to subj_id for extra speed");
    }
    if ($subQ) {
        push @where, sprintf("sh.sbj_id IN (%s)", join
                             (',', map {'?'} @{$subQ}));
        push @binds, sort { $a <=> $b } @{$subQ};
        push @what, "subj_id = ".join('/', @{$subQ});
        $noSubQ = undef;
    }
    if ($noSubQ) {
        push @where, sprintf("sh.sbj_id NOT IN (%s)", join
                             (',', map {'?'} @{$noSubQ}));
        push @binds, @{$noSubQ};
        push @what, "subj_id != ".join('/', @{$noSubQ});
        $filter ||= [];
        push @{$filter}, 'ExcludedSubjectID';
    }
    if (my $srcQ = &_array_args( $args, 'SRCIDS' )) {
        push @geneBits, sprintf("sh.src_id IN (%s)", join
                                (',', map {'?'} @{$srcQ}));
        push @geneBinds, @{$srcQ};
        push @what, "src_id = ".join('/', @{$srcQ});
        $filter ||= [];
        push @{$filter}, 'SourceID';
    }
    unless ($#geneBits == -1) {
        push @where, '('.join(' OR ', @geneBits).')';
        push @binds, @geneBinds;
    }
    my $wc = join(' AND ', @where);
    return ($wc || "", \@binds, \@what, join(' ', @needsFilter), $filter);
}

sub recover_search_hits {
    my $self = shift;
    $self->bench_start();
    
    my $args = $self->parseparams( @_ );
    my $pgtm   = $self->pgtm();
    my ($where, $binds, $what) = $self->_search_hit_where_clause( @_ );
    # "sh.genestrand, sh.src_id, sh.footprint, sh.details, sh.srch_id
    my @defCol = qw(SubjectSeq MisMatch SubjectID SubjectPos SubjectStrand 
                    GeneStrands SourceID GenomicFootprint GeneDetails SearchID);
    if (my $creq = $args->{COLS}) {
        @defCol = ref($creq) ? @{$creq} : split(/[\t\s\,]/, $creq);
    }
    my (@cols, @dbCols);
    foreach my $col (@defCol) {
        my $cname = $self->standard_data_name( $col );
        unless ($cname) {
            $args->msg_once("[?]", "Unknown column request '$cname'");
            next;
        }
        my $dbcol = $self->get_metadata('columns', $cname, 'dbcol');
        unless ($dbcol) {
            $self->msg_once("[?]", "No database column assigned to $cname");
            next;
        }
        push @cols, $cname;
        # Normalize to table search_hit:
        $dbcol =~ s/^.+\.//;
        push @dbCols, "sh.$dbcol";
    }

    my $rv = {
        cols => \@cols,
        rows => [],
        flat => 1,
        type => 'detailed',
    };
      
    unless ($where) {
        $self->bench_end();
        $pgtm->change_task_status({error => "No parameters provided for recovering old results"}) if ($pgtm);
        return $rv;
    }
    my @subtask = ( ["Query DB", 0] );
    my $status = {
        task    => "Recover hits from database",
        subtask => \@subtask,
    };
    $pgtm->change_task_status($status) if ($pgtm);
    my $sql = "SELECT ".join(', ', @dbCols)." FROM search_hit sh WHERE $where";
    my $sth = $self->dbh->prepare( -name => "Recover old search data",
                                   -sql => $sql,
                                   -limit => $args->{LIMIT} || 0);
    
    $self->msg("[SQL]", $sth->show_and_explain( $binds))
        if ($args->{DUMPSQL});
    my $rows   = $sth->selectall_arrayref( @{$binds} );
    my $detInd = $self->data_index( $rv, 'GeneDetails' );
    # $self->prebranch($rv);
    if ($detInd != -1) {
        # Normalize the JSON details to a Perl structure
        for my $i (0..$#{$rows}) {
            my $row = $rows->[$i];
            my $details;
            if (my $json = $row->[$detInd]) {
                eval {
                    $details = decode_json( $json );
                };
                unless ($details) {
                    $self->msg_once("Failed to parse details JSON", "'$json'");
                }
            }
            $row->[$detInd] = $details || {};
        }
    }
    $subtask[0][1] = 1;
    $subtask[0][2] = "Done";
    $pgtm->change_task_status($status) if ($pgtm);
    $rv->{rows} = $rows;
    # $self->prebranch($rv); die;
    $self->bench_end();
    return $rv;
}

sub search_summary {
    my $self = shift;
    my ($where, $binds, $what, $doFilter, $filters) =
        $self->_search_hit_where_clause( @_ );
    my %struct;
    my $maxMM = 0;
    my $rv = {
        flat => 1,
        type => 'summary',
        filters => $filters,
    };
    my $args = $self->parseparams( @_ );
    if (!$doFilter || $doFilter eq 'strand') {
        # We can get what we need by a simple SQL query
        my $sql =
            "SELECT srch_id, src_id, mm, count(mm) FROM search_hit sh".
            " WHERE $where GROUP BY srch_id, src_id, mm".
            " ORDER BY srch_id, src_id, mm ASC";
        my $sth = $self->dbh->prepare( -name => "Summarize search hits",
                                       -sql  => $sql );
        $self->msg("[SQL]", $sth->show_and_explain( $binds))
            if ($args->{DUMPSQL});
        my $rows = $sth->selectall_arrayref( @{$binds} );
        foreach my $row (@{$rows}) {
            my ($srch_id, $src_id, $mm, $num) = @{$row};
            $mm = $maxShownMM if ($mm > $maxShownMM);
            $struct{$srch_id}{$src_id}{$mm} += $num;
        }
    } else {
        # We will need to get all rows, then filter them in detail
        my $data = $self->recover_search_hits
            ( @_, -cols => [qw(SearchID SourceID MisMatch SubjectID)], );
        $self->filter_flat_hits( -hits => $data, @_ );
        my $subInd  = $self->data_index( $data, 'SubjectID');
        my $srcInd  = $self->data_index( $data, 'SourceID' );
        my $mmInd   = $self->data_index( $data, 'MisMatch' );
        my $srchInd = $self->data_index( $data, 'SearchID' );
        my $rows    = $data->{rows};
        foreach my $row (@{$rows}) {
            my ($srch_id, $src_id, $mm, $subj) =
                ($row->[$srchInd], $row->[$srcInd], $row->[$mmInd], 
                 $row->[$subInd]);
            $mm = $maxShownMM if ($mm > $maxShownMM);
            # We need to build this as a hash to be consistent with
            # direct DB recovery. The direct SQL counts unique subjects,
            # but it is possible for one subject to have multiple hits
            $struct{$srch_id}{$src_id}{$mm}{$subj} = 1;
        }
        $self->err("Working on filtered summary. Where clause:", $where);
    }
    # Now reduce the data to a table
    my $tab = $rv->{rows} = [];
    foreach my $srchid (sort {$a <=> $b } keys %struct) {
        foreach my $srcid (sort {$a <=> $b } keys %{$struct{$srchid}}) {
            my @row = ($srchid, $srcid);
            push @{$tab}, \@row;
            foreach my $mm (keys %{$struct{$srchid}{$srcid}}) {
                $maxMM = $mm if ($maxMM < $mm);
                my $num = $struct{$srchid}{$srcid}{$mm};
                if (ref($num)) {
                    my @u = keys %{$num};
                    $num = $#u + 1;
                }
                $row[$mm + 2] = $num;
            }
        }
    }
    my @cols = qw(SearchID SourceID);
    push @cols, map { "MM$_" } (0..$maxMM);
    $rv->{cols} = \@cols;
    # $self->prebranch($rv); die;
    return $rv;
}

sub store_search_hits {
    my $self = shift;
    my ($data) = @_;
    return unless ($data);
    $self->bench_start();
    my $dbh   = $self->dbh;
    # Determine if the data are already in the right order for bulk loading:
    my @dbCol  = map { "sh.".lc($_) } $dbh->column_order( 'search_hit' );
    my @inDat  = map { $self->data_index( $data, $_ ) } @dbCol;
    my $expect = join(' ', (0..$#dbCol));
    my $actual = join(' ', @inDat);
    my $reMap;
    unless ($expect eq $actual) {
        # The data block is not in the same order as the database. We will
        # have to rearrange the column order for loading
        my %colMap;
        my $hasNull = 0;
        my $maxInd  = $#dbCol;
        for my $in (0..$maxInd) {
            # my $dbColName = $dbCol[$in];
            my $inInd = $inDat[$in];
            if ($inInd == -1) {
                $hasNull++;
            } else {
                $colMap{$inInd} = $in;
            }
        }
        $reMap = sub {
            my $inRow = shift;
            my @outRow;
            # Pad out to maximum length
            $outRow[$maxInd] = undef;
            while (my ($inInd, $outInd) = each %colMap) {
                $outRow[$outInd] = $inRow->[$inInd];
            }
            return \@outRow;
        };
        # die $self->branch( {dbCol => \@dbCol, inDat => \@inDat, colMap => \%colMap});
    }
    my %dbColH  = map { ($dbCol[$_] || '') => $_ } (0..$#dbCol);
    my $detInd  = $dbColH{'sh.details'};
    unless ($detInd) {
        $self->err("[?]","Failed to find gene details column for data");
        $self->bench_end();
        return -1;
    }
    my $srchInd = $self->data_index( $data, "SearchID" );
    unless (defined $srchInd) {
        $self->err("[?]","Failed to find search ID for data");
        $self->bench_end();
        return -1;
    }
    my %seen;
    map { $seen{ $_->[$srchInd] ||0 }++ } @{$data->{rows}};

   
    # We will load in blocks to allow task tickets to be updated
    my $task  = $self->task();
    my $blk   = 10000;
    my $hlen  = $#{$data->{rows}};
    
    my $bar;
    if ($task) {
        my $bars = $task->param("progbars") || $task->param("progbars", []);
        $bar  = ['Record results', 0];
        push @{$bars}, $bar;
        $task->write(1);
    }
    my $delSTH = $dbh->named_sth('Clear search hits');
    $dbh->begin_work;
    foreach my $srchid (keys %seen) {
        $delSTH->execute( $srchid );
    }
    for (my $i = 0; $i <= $hlen; $i += $blk) {
        my $e = $i + $blk - 1;
        $e = $hlen if ($e > $hlen);
        my @toLoad;
        for my $j ($i..$e) {
            my $row = $reMap ? &{$reMap}($data->{rows}[$j]) : $data->{rows}[$j];
            if (ref($row->[$detInd])){
                # GeneDetails are represented as a JSON structure
                # Serialize
                $row->[$detInd] = encode_json( $row->[$detInd] );
            }
            push @toLoad, $row;
        }
        $self->dbh->insert_array( 'search_hit', \@toLoad );
        if ($task) {
            $bar->[1] = &_progress_fraction( $e, $hlen );
            $task->write(1);
        }
    }
    $dbh->commit;
    if ($task) {
        $bar->[0] = "Database updated.";
        $bar->[1] = 1;
        $task->write(1);
    }
    $self->bench_end();
    return $hlen + 1;
}

=head2 _careful_mismatch

Counting mismatches is complicated by ambiguities in the source
database - it is possible that a specific location on the genome gets
'hit twice', with differing alignment scores for each hit. For
example

   Query     AAATAAAATAAA
   Genome    AAAYAAAAYAAA

   Word1     AAAcAAAAcAAA   MM = 2
   Word2     AAATAAAAcAAA   MM = 1
   Word3     AAAcAAAATAAA   MM = 1
   Word4     AAATAAAATAAA   MM = 1

We also have uneven coverage. For example a 12bp wordsize with a 15bp oligo

=cut

sub _careful_mismatch {
    my $self = shift;
    my ($fd, $mm, $ws, $pos) = @_;
    # If the crude measure is ok, do not bother checking deeper:
    return 0 if ($fd->[0] <= $mm);
    $self->bench_start();
    $ws--;
    # We will track all positions (indices) that have succesfully matched:
    my $matchedPos = $fd->[2] ||= [];
    foreach my $wd (@{$fd->[1]}) {
        # If we have already parsed this entry ignore it
        next if ($wd->[3]);
        $wd->[3] = 1;
        my $ind  = $wd->[1];
        # Note the non-matching positions
        my %notMatched = map { $_ => 1 } @{$wd->[2]};
        for my $i (0..$ws) {
            $matchedPos->[$ind + $i] ||= $notMatched{$i} ? 0 : 1;
        }
    }
    my $newMM = 0;
    map { $newMM++ unless ($matchedPos->[$_]) } (0..$#{$matchedPos});
    $fd->[0] = $newMM;
    $self->bench_end();
    return $newMM > $mm ? $newMM : 0;
}

sub permute_mismatches {
    # ~ 30ms for 12bp + 2 mismatch
    my $self = shift;
    my ($base, $mm) = @_;
    my %rv;
    return \%rv unless ($base);
    $self->bench_start();
    # Begin by expanding any ambiguities in the initial base sequence
    # Length is zero indexed:
    my $slen  = CORE::length($base) - 1;
    my $prior = [];
    foreach my $seq ( $self->expand_ambiguous_sequence( uc($base) ) ) {
        # Prior tracks the sequences that will undergo next round of mutation
        # [The sequence, and the last position it was permuted at]:
        push @{$prior}, [ $seq, 0 ];
        # RV is a hash keyed to sequences, pointing to an array of
        # zero-indexed mismatch locations
        $rv{$seq} = [  ];
    }
    if ($mm) {
        # Note which characters are considred mismatches for each position
        # There may be a more elegant way of doing this...
        my %all = map { $_ => 1 } qw(A C G T);
        my @mismatchChars  = map { { %all } } (0..$slen);
        foreach my $seq ( keys %rv) {
            my @chars = split('', $seq);
            for my $i (0..$#mismatchChars) {
                delete $mismatchChars[$i]{ $chars[$i] };
            }
        }
        # Normalize the hashes to arrays
        for my $i (0..$#mismatchChars) {
            $mismatchChars[$i] = [ sort keys %{$mismatchChars[$i]} ];
        }
        # Now start building all sequences with up to $mm mismatches
        for my $mismatch (1..$mm) {
            my @current;
            for my $p (0..$#{$prior}) {
                my ($pSeq, $lp) = @{$prior->[$p]};
                my @pmm = @{$rv{$pSeq}};
                for my $i ($lp..$slen) {
                    my $mmChars = $mismatchChars[$i];
                    foreach my $base (@{$mmChars}) {
                        my $mutSeq = $pSeq;
                        substr($mutSeq, $i, 1) = $base;
                        # I do not think we will revisit mutations with the
                        # current logic. But put in a safety catch in case.
                        next if ($rv{$mutSeq});
                        $rv{$mutSeq} = [ @pmm, $i];
                        push @current, [$mutSeq, $i + 1];
                    }
                }
            }
            $prior = \@current;
        }
    }
    $self->bench_end();
    return \%rv;
}

sub _permuted_mismatches_to_text {
    my $self = shift;
    my $mmh  = shift;
    my @seqs = sort { $mmh->{$a}[0] <=> $mmh->{$b}[0] ||
                          $a cmp $b } keys %{$mmh || {}};

    my $text = "";
    my %cnt;
    foreach my $seq (@seqs) {
        my ($mm, $mmLoc) = @{$mmh->{$seq}};
        $cnt{$mm}++;
        my $show = $seq;
        foreach my $i (@{$mmLoc}) {
            substr($show, $i, 1) = lc(substr($show, $i, 1));
        }
        $text .= sprintf("  %s %2d\n", $show, $mm);
    }
    my $head = sprintf("%d permutations:\n", $#seqs + 1);
    foreach my $mm (sort {$a <=> $b} keys %cnt) {
        $head .= sprintf("   %dmm = %d\n", $mm, $cnt{$mm});
    }
    return $head . $text;
}

*dbh = \&connect;
*dbi = \&connect;
sub connect {
    my $self = shift;
    unless ($self->{DBH}) {
        $self->bench_start();
        # Try to read connection details from parameter file
        my $intName  = $self->instance();
        my $dbType   = "dbi:Pg:service=$intName";

        my $dbh      = BMS::FriendlyDBI->connect
            ($dbType, undef, undef, {
                RaiseError  => 0,
                PrintError  => 0,
                LongReadLen => 100000,
                RowCacheSize => 1,
                AutoCommit  => 1, },
             -errorfile => $self->{ERRFILE},
             -adminmail => $self->{ERRMAIL}, );
        $dbh->schema( $self->schema() );

        $self->{DBH} = $dbh;
        $self->_set_statement_data();
        $self->bench_end();
    }
    return $self->{DBH};
}

sub disconnect {
    my $self = shift;
    if (my $dbh = $self->{DBH}) {
        $dbh->disconnect();
    }
}

sub src_id {
    my $self = shift;
    my ($path, $ws) = @_;
    return 0 unless ($path && $ws);
    my $sth = $self->dbh->named_sth
        ('Get or create a src_id for a DB source and wordsize');
    return $sth->get_single_value( $path, $ws );
}


sub sbj_id {
    my $self = shift;
    my ($name, $srcid, $len) = @_;
    if (!$name) {
        return 0;
    } elsif ($srcid) {
        my $sth = $self->dbh->named_sth
            ('Get a sbj_id for subject - create if absent');
        return $sth->get_single_value( $name, $srcid, $len );
    } else {
        # Just query
        my $sth = $self->dbh->named_sth
            ('Get a sbj_id for subject - no creation or source');
        my @found = $sth->get_array_for_field( $name );
        return $#found == 0 ? $found[0] : 0;
    }
}

sub subject_name {
    my $self = shift;
    my $id  = shift || "";
    unless ($id =~ /^\d+$/) {
        return $id;
    }
    my $sth = $self->_get_subject_sth();
    $sth->execute( $id );
    my $sRows = $sth->fetchall_arrayref();
    return $sRows->[0][0] || "";
}

sub gene_info {
    my $self = shift;
    my @rv;
    if (my $gid = shift) {
        $self->bench_start();
        my $locSTH = $self->dbh->named_sth('Get Gene Information');
        $locSTH->execute($gid);
        my $gRows = $locSTH->fetchall_arrayref();
        @rv = @{$gRows->[0] || []};
        # Deal with ill-advised BMS symbol tokens:
        $rv[1] =~ s/[\*\?]$// if ($rv[1]);
        $self->bench_end()
    }
    return wantarray ? @rv : \@rv;
}

sub source_info {
    my $self = shift;
    my $srcid = shift || 0;
    my $wsreq = shift;
    my $key   = $wsreq ? "$srcid\t$wsreq" : $srcid;
    unless ($self->{SRCINFO}{$key}) {
        my $sth;
        if ($wsreq) {
            # Path and wordsize
            $sth = $self->_get_source_by_path_sth();
            $sth->execute( $srcid, $wsreq );
        } else {
            # src_id request
            $sth = $self->_get_source_sth();
            $sth->execute( $srcid );
        }
        $sth->execute( $srcid );
        my $sRows = $sth->fetchall_arrayref();
        my ($sid, $path, $ws, $num, $chr, $mod) = @{$sRows->[0] || []};
        $ws ||= 0;
        my $info = 
            $self->{SRCINFO}{$srcid} = 
            $self->{SRCINFO}{"$path\t$ws"} = {
            src_id      => $sid,
            Path        => $path,
            WordSize    => $ws,
            Entries     => $num,
            Characters  => $chr,
            ModDate     => $mod,
        };
        my $tags = $self->tagvals( $srcid );
        foreach my $tv (@{$tags}) {
            my ($t, $v) = @{$tv};
            $info->{$t} ||= $v;
        }
        unless ($info->{Species}) {
            my $hack = $path;
            $hack =~ s/.+\///;
            if ($hack =~ /^(.+)_refseq/i) {
                my $species = lc($1);
                $species =~ s/_/ /g;
                substr($species, 0, 1) = uc(substr($species, 0, 1));
                $info->{Species} = $species;
            }
        }
    }
    return $self->{SRCINFO}{$key};
}

sub word_for_id {
    my $self = shift;
    my $id   = shift || "";
    unless ($id =~ /^\d+$/) {
        return $id;
    }
    my $sth = $self->_get_word_sth();
    $sth->execute( $id );
    return $sth->get_single_value($id) || "";
}

sub prior_search_ids {
    my $self = shift;
    my ($seq) = $self->_normalize_sequence( shift );
    my @rv;
    if ($seq) {
        my $chk = $self->dbh->named_sth('Get searches by sequence query');
        $chk->execute($seq);
        @rv = $chk->get_array_for_field();
    }
    return wantarray ? @rv : \@rv;
}

sub get_tag {
    my $self = shift;
    my @rv;
    if (my $objID = shift) {
        if (my $tag = shift) {
            my $sth = $self->dbh->named_sth('Get values for a specific tag');
            @rv = $sth->get_array_for_field($objID, $tag);
        }
    }
    return wantarray ? @rv : $rv[0];
}

sub set_tag {
    my $self = shift;
    if (my $objID = shift) {
        if (my $tag = shift) {
            my $val = shift;
            if (defined $val) {
                my $sth = $self->dbh->named_sth('Set a tag/val pair');
                $sth->execute( $objID, $tag, $val );
            }
        }
    }
}

sub clear_tag {
    my $self = shift;
    if (my $objID = shift) {
        if (my $tag = shift) {
            my $sth = $self->dbh->named_sth('Clear a tag for object');
            $sth->execute( $objID, $tag );
        }
    }
}

sub gene_id {
    my $self = shift;
    my $acc = shift;
    return 0 unless ($acc);
    my $sth = $self->dbh->named_sth('Get gene_id for accession');
    $sth->execute( $acc );
    return $sth->get_single_value($acc) || 0;
}

sub set_gene_loc_details {
    my $self = shift;
    my ($sid, $gid, $l, $r, $str, $sc, $details) = @_;
    return unless ($sid && $gid);
    return unless (defined $l && defined $r);
    my $sth = $self->dbh->named_sth('Set Gene Location details');
    $sth->execute( $sid, $gid, $l, $r, $str, $sc, $details );
}

sub set_gene_meta {
    my $self = shift;
    my $args   = $self->parseparams( @_ );
    my $gid;
    if ($gid = $args->{GID}) {
        # directly passed
    } elsif (my $acc = $args->{ACC}) {
        unless ($gid = $self->gene_id( $acc )) {
            $self->msg_once("[!]","Failed to get gene_id for '$acc'");
            return;
        }
    } else {
        $self->msg_once("[!]","Failed to set_gene_meta()",
                        "Provide either -id or -acc");
        return;
    }
    if (my $sym = $args->{SYM}) {
        my $sth = $self->dbh->named_sth("Set Gene Symbol");
        $sth->execute($sym, $gid);
    }
    if (my $desc = $args->{DESC} || $args->{DESCR}) {
        my $sth = $self->dbh->named_sth("Set Gene Description");
        $sth->execute($desc, $gid);
    }
}

sub tagvals {
    my $self = shift;
    my $rv = [];
    if (my $id = shift) {
        my $sth = $self->dbh->named_sth("Get tag/value pairs");
        $sth->execute($id);
        $rv = $sth->fetchall_arrayref();
    }
    return $rv;
}


sub _get_wrd_id_sth {
    return shift->dbh->named_sth
        ('Get a wrd_id from a word - no creation');
}

sub _set_word_sth {
    return shift->dbh->named_sth
        ('Get a wrd_id from a word - create if absent');
}

sub _get_word_sth {
    return shift->dbh->named_sth('Get a word from a wrd_id');
}

sub _get_subject_sth {
    return shift->dbh->named_sth('Get subject information for sbj_id');
}

sub _get_source_sth {
    return shift->dbh->named_sth('Get source information for src_id');
}

sub _get_source_by_path_sth {
    return shift->dbh->named_sth
        ('Get source information for path and wordsize');
}

sub _set_srchit_sth {
    return shift->dbh->named_sth('Set a hit for a word');
}


sub set_hits {
    my $self = shift;
    $self->death("NEED TO IMPLEMENT");
}

sub schema {
    my $tables = {

        sourcedb => {
            name => 'sourcedb',
            com  => "A fasta flat file representing a searchable DB",
            index => {
                src_path_ind   => {
                    cols => [ 'path', 'wordsize' ],
                    unique => 1,
                },
            },
            pkey  => 'src_id',
            cols  => [ ['src_id', 'integer', 
                        'Primary key ID representing this object', {
                            SEQUENCE => 'global_seq',
                        } ],
                       ['path', 'text', 
                        'The path to the file'],
                       ['num', 'integer', 
                        'The number of subject entries in the database'],
                       ['chars', 'integer', 
                        'The number of subject nucleotides'],
                       ['wordsize', 'integer', 
                        'The size of the hash words used in the DB'],
                       ['mod_date', 'varchar(25)', 
                        'The time the file was last modified'],
                       ],
        },
        
        subject => {
            name => 'subject',
            com  => "A sequence entry from one of the source databases",
            index => {
                sbj_primary   => {
                    cols => [ 'name', 'src_id' ],
                    unique => 1,
                },
                sbj_by_src   => {
                    cols => [ 'src_id' ],
                },
                sbj_uc_name   => {
                    cols => [ 'upper(name)', ],
                },  

            },
            fkey  => {
                src_id  => 'sourcedb.src_id',
            },
            pkey  => 'sbj_id',
            cols  => [ ['sbj_id', 'integer', 
                        'Primary key ID representing this object', {
                            SEQUENCE => 'global_seq',
                        } ],
                       ['src_id', 'integer', 
                        'Fkey pointing to the source database'],
                       ['name', 'text', 
                        'The name / display_id of the sequence'],
                       ['len', 'integer', 
                        'The length of the subject'],
                       ],
        },

        search_seq => {
            name => 'search_seq',
            com  => "A sequence used to search the database",
            index => {
                srchseq_seq   => {
                    cols => [ 'upper(seq)', ],
                },
                srchseq_name   => {
                    cols => [ 'upper(name)', ],
                },
            },
            pkey  => 'srch_id',
            cols  => [ ['srch_id', 'integer', 
                        'Primary key ID representing this object', {
                            SEQUENCE => 'global_seq',
                        } ],
                       ['seq', 'text', 
                        'The sequence of the search sequence'],
                       ['name', 'text', 
                        'The accession of the search sequence'],
                       ['descr', 'text', 
                        'The description of the search sequence'],
                       ['daterun', 'date', 
                        'The time that the search was run'],
                       ['mm', 'int', 
                        'Maximum number of allowed mismatches'],
                       ['status', 'integer', 
                        '-1 indicates completion. Other values indicate PID running the search'],
                       ],
        },

        search_hit => {
            name => 'search_hit',
            com  => "A stored hit for a searched oligo",
            index => {
                sh_srchid => { cols => [ 'srch_id', 'mm' ] },
                sh_sbjid  => { cols => [ 'sbj_id' ] },
                sh_srcid  => { cols => [ 'src_id' ] },
            },
            fkey  => {
                srch_id  => 'search_seq.srch_id ON DELETE CASCADE',
                src_id   => 'sourcedb.src_id',
                sbj_id   => 'subject.sbj_id',
            },
            cols  => [ ['srch_id', 'integer', 
                        'FKEY pointing back to search_seq' ],
                       ['subseq', 'text', 
                        'The aligned sequence from the subject'],
                       ['mm', 'int', 
                        'Number of mismatches'],
                       ['sbj_id', 'int', 
                        'FKEY pointing to the hit subject'],
                       ['pos', 'int',
                        'The starting position of the alignment'],
                       ['strand', 'int', 
                        'The strand of the alignment'],
                       ['genestrand', 'int', 
                        'The strand of overlapping genes. 0 is both'],
                       ['src_id', 'int', 
                        'FKEY pointing to the hit source'],
                       ['footprint', 'text',
                        'Hit footprint in absolute genomic coordinates'],
                       ['gene_str_type', 'integer[]', 
                        '1D array of genes, using -2 as separators. gene_id, strand, and one or more exon_type ids'],
                       ['details', 'text', 
                        'A JSON string describing genome-, gene- and rna-level details of the hit'],
                       ],
        },

        gene => {
            name => 'gene',
            com  => "Details for a gene",
            index => {
                gene_acc_ind   => {
                    cols => [ 'upper(acc)', ],
                    unique => 1,
                },
                gene_alias_ind => {
                    cols => ['aliases'],
                    using => 'gin',
                },
            },
            pkey  => 'gene_id',
            cols  => [ ['gene_id', 'integer', 
                        'Primary key referencing the gene', {
                            SEQUENCE => 'global_seq',
                        } ],
                       ['acc', 'varchar(50)', 
                        'Gene accession'],
                       ['sym', 'varchar(50)', 
                        'Main gene symbol' ],
                       ['taxa', 'varchar(50)', 
                        'Scientific species name'],
                       ['descr', 'text',
                        'Description for the gene'],
                       ['aliases', 'varchar(50)[]',
                        'Synonyms and alternative symbols'],
                       ['hs_gene', 'varchar(50)',
                        'The human orthologue'],
                       ['hs_sym', 'varchar(50)',
                        'Human orthologue gene symbol (best available)'],
                       ['hs_score', 'numeric',
                        'Value between 0-100, a percentage score relating the similarity to the human orthologue'],
                       ],
        },

        geneloc => {
            name => 'geneloc',
            com  => "Location of a gene within a subject",
            index => {
                geneloc_sbj_ind   => {
                    cols => [ 'sbj_id', ],
                },
                gene_primary_ind   => {
                    cols => [ 'gene_id', 'lft', 'sbj_id' ],
                    unique => 1,
                },
            },
            fkey  => {
                sbj_id   => 'subject.sbj_id',
                gene_id  => 'gene.gene_id',
            },
            pkey  => 'geneloc_id',
            cols  => [ ['geneloc_id', 'integer', 
                        'Primary key referencing the gene at this location', {
                            SEQUENCE => 'global_seq',
                        } ],
                       ['sbj_id', 'int', 
                        'FKEY pointing to the subject the gene resides on'],
                       ['gene_id', 'int', 
                        'FKEY pointing to the gene entry' ],
                       ['lft', 'int', 
                        'The subject coordinate to the left of the gene'],
                       ['rgt', 'int',
                        'The subject coordinate to the right of the gene'],
                       ['strand', 'int', 
                        'The strand of the gene on the subject'],
                       ['score', 'real',
                        'An alignment score, generally a percent match'],
                       ['details', 'text', 
                        'A JSON string holding exon-level information for individual transcripts'],
                       ],
        },

        word => {
            name => 'word',
            com  => "Normalization for sequence words",
            index => {
                word_primary   => {
                    cols => [ 'word', ],
                    unique => 1,
                },
            },
            pkey  => 'wrd_id',
            cols  => [ ['wrd_id', 'integer', 
                        'Primary key ID representing the word', {
                            SEQUENCE => 'word_seq',
                        } ],
                       ['word', 'varchar(100)', 
                        'The sequence being hit' ],
                       ],
        },

        wordsrc => {
            name => 'wordsrc',
            com  => "All locations of a specific word in a source",
            index => {
                wsrc_primary   => {
                    cols => [ 'wrd_id', 'src_id' ],
                    unique => 1,
                },
            },
            fkey  => {
                # Keep these informal for speed?
                #wrd_id  => 'word.wrd_id',
                #sbj_id  => 'subject.sbj_id',
            },
            cols  => [ ['wrd_id', 'integer', 
                        'Fkey pointing to the sequence being hit' ],
                       ['src_id', 'integer', 
                        'Fkey pointing to the subject sequence'],
                       ['num', 'integer', 
                        'Number of distinct hits. Excludes exon junction locations'],
                       ['hits', 'integer[]', 
                        '1D array of hits. A sbj_id is followed by the coordinates found in that subject. Multiple subjects are separated by zeros'],
                       ],
        },
        
        tagval => {
            name  => 'tagval',
            com   => 'Tag-value pairs assigned to arbitrary object',
            index => {
                tagval_primary   => {
                    cols => [ qw(obj_id) ],
                },
                tagval_bytag   => {
                    cols => [ qw(tag obj_id val) ],
                    unique => 1,
                },
                tagval_byval   => {
                    cols => [ qw(val obj_id) ],
                },
            },
            cols  => [['obj_id', 'integer',
                       'The PKEY for the object being tagged', ],
                      ['tag', 'text',
                       'The name of the tag' ],
                      ['val', 'text',
                       'The value assigned to the tag' ],
                      ],
        },

        exon_type => {
            name  => 'exon_type',
            com   => 'Simple table to allow exon types to be integerized',
            index => {
                et_primary   => {
                    cols => [ qw(exontype) ],
                },
            },
            pkey  => 'et_id',
            cols  => [['et_id', 'integer',
                       'The PKEY for the exon type being referenced',  {
                           SEQUENCE => 'global_seq',
                       }],
                      ['exontype', 'text',
                       'The human-readable exon type name' ],
                      ],
        },

        quiet_tagval_write => {
            name  => 'quiet_tagval_write',
            com   => 'Write a tag/value pair to the DB without complaining',
            db    => 'postgres',
            args  => [ oid    => 'INT',
                       tagtxt => 'TEXT',
                       valtxt => 'TEXT' ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
BEGIN
  PERFORM 0 FROM tagval WHERE obj_id = oid AND tag = tagtxt AND val = valtxt;
IF NOT FOUND THEN
  -- It is still possible for a violation to occur here, but we do not care
  -- The unique constraint prevents duplicate rows in DB, and this function
  -- minimizes complaints in the error log
  INSERT INTO tagval (obj_id, tag, val) VALUES (oid, tagtxt, valtxt);
  RETURN 1;
END IF;
  -- The value already exists in database, leave happy
  RETURN 0;
END;
"
},

        quiet_geneloc_details => {
            name  => 'quiet_geneloc_details',
            com   => 'Write data for a gene assignment',
            db    => 'postgres',
            args  => [ sid    => 'INT',
                       gid    => 'INT',
                       l      => 'INT',
                       r      => 'INT',
                       str    => 'INT',
                       sc     => 'real',
                       det    => 'text' ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
BEGIN
  PERFORM 0 FROM geneloc WHERE sbj_id = sid AND gene_id = gid AND lft = l;
IF NOT FOUND THEN
  -- Make sure the word / source / subject row is in the database
  INSERT INTO geneloc (sbj_id, gene_id, lft) VALUES (sid, gid, l);
END IF;
  -- Now update the other data
  UPDATE geneloc SET rgt = r, strand = str, score = sc, details = det
   WHERE sbj_id = sid AND gene_id = gid AND lft = l;
  RETURN 0;
END;
"
},

        quiet_gene_id => {
            name  => 'quiet_gene_id',
            com   => 'Make or get a gene entry',
            db    => 'postgres',
            args  => [ name   => 'varchar', ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
DECLARE rv integer;
BEGIN
  SELECT gene_id INTO rv FROM gene 
   WHERE upper(acc) = upper(name);
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this name and source
  BEGIN
    INSERT INTO gene (acc) VALUES (name);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the id that should have been auto-generated
  SELECT gene_id INTO rv FROM gene 
   WHERE upper(acc) = upper(name);
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate gene_id';
  RETURN 0;
END;
"
},

        quiet_exontype_id => {
            name  => 'quiet_exontype_id',
            com   => 'Make or get a exon type entry',
            db    => 'postgres',
            args  => [ val   => 'text', ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
DECLARE rv integer;
BEGIN
  SELECT et_id INTO rv FROM exon_type 
   WHERE exontype = val;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this text value
  BEGIN
    INSERT INTO exon_type (exontype) VALUES (val);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the id that should have been auto-generated
  SELECT et_id INTO rv FROM exon_type 
   WHERE exontype = val;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate et_id';
  RETURN 0;
END;
"
},

        quiet_wordsrc_write => {
            name  => 'quiet_wordsrc_write',
            com   => 'Write hits for a word in a database',
            db    => 'postgres',
            args  => [ wid    => 'INT',
                       sid    => 'INT',
                       count  => 'INT',
                       poses  => 'INT[]' ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
BEGIN
  PERFORM 0 FROM wordsrc WHERE wrd_id = wid AND src_id = sid;
IF NOT FOUND THEN
  -- Make sure the word / source / subject row is in the database
  INSERT INTO wordsrc (wrd_id, src_id) VALUES (wid, sid);
END IF;
  -- Now update the coordinates
  UPDATE wordsrc SET num = count, hits = poses
   WHERE wrd_id = wid AND src_id = sid;
  RETURN 0;
END;
"
},

        quiet_subject_write => {
            name  => 'quiet_subject_write',
            com   => 'Write a subject entry to the database',
            db    => 'postgres',
            args  => [ srcid    => 'INT',
                       nametxt  => 'TEXT' ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
BEGIN
  PERFORM 0 FROM subject WHERE src_id = srcid AND name = nametxt;
IF NOT FOUND THEN
  -- It is still possible for a violation to occur here, but we do not care
  -- The unique constraint prevents duplicate rows in DB, and this function
  -- minimizes complaints in the error log
  INSERT INTO subject (src_id, name) VALUES (srcid, nametxt);
  RETURN 1;
END IF;
  -- The value already exists in database, leave happy
  RETURN 0;
END;
"
},

        quiet_word_write => {
            name  => 'quiet_word_write',
            com   => 'Fetch or create a pkey entry for a word',
            db    => 'postgres',
            args  => [ seq    => 'varchar(100)' ],
            retval    => 'INT',
            language  => 'plpgsql',
            function  =>
"
DECLARE rv integer;
BEGIN
  SELECT wrd_id INTO rv FROM word WHERE word = seq;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this name and source
  BEGIN
    INSERT INTO word (word) VALUES (seq);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the pop_id that should have been auto-generated
  SELECT wrd_id INTO rv FROM word WHERE word = seq;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate wrd_id';
  RETURN 0;
END;
"
},

    source_db_id => {
        name  => 'source_db_id',
        com   => 'Fetch or create a pkey entry for a source',
        db    => 'postgres',
        args  => [ fp => 'TEXT',
                   ws => 'INT' ],
        retval    => 'INT',
        language  => 'plpgsql',
        function  =>
"
DECLARE rv integer;
BEGIN
  SELECT src_id INTO rv FROM sourcedb WHERE path = fp AND wordsize = ws;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this name and source
  BEGIN
    INSERT INTO sourcedb (path, wordsize) VALUES (fp, ws);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the pop_id that should have been auto-generated
  SELECT src_id INTO rv FROM sourcedb WHERE path = fp AND wordsize = ws;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate src_id';
  RETURN 0;
END;
"

},
        
    create_search_entry => {
        name  => 'create_search_entry',
        com   => 'Create a new search entry',
        db    => 'postgres',
        args  => [ dna  => 'TEXT',
                   mis  => 'INT',
                   id   => 'TEXT',
                   d    => 'TEXT' ],
        retval    => 'INT',
        language  => 'plpgsql',
        function  =>
"
DECLARE rv integer;
DECLARE n date;
BEGIN
  SELECT now INTO n FROM now();
  INSERT INTO search_seq (seq, name, descr, mm, daterun, status)
  VALUES (dna, id, d, mis, n, 0);
  SELECT srch_id INTO rv FROM search_seq WHERE seq = dna AND daterun = n;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate srch_id';
  RETURN 0;
END;
"
},

    subject_id => {
        name  => 'subject_id',
        com   => 'Fetch or create a pkey entry for a subject',
        db    => 'postgres',
        args  => [ nm  => 'TEXT',
                   sid => 'INT',
                   l => 'INT'],
        retval    => 'INT',
        language  => 'plpgsql',
        function  =>
"
DECLARE rv integer;
BEGIN
  SELECT sbj_id INTO rv FROM subject WHERE name = nm AND src_id = sid;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this name and source
  BEGIN
    INSERT INTO subject (name, src_id, len) VALUES (nm, sid, l);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the pop_id that should have been auto-generated
  SELECT sbj_id INTO rv FROM subject WHERE name = nm AND src_id = sid;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate src_id';
  RETURN 0;
END;
"
},
    
    v_gl => {
        name  => 'v_gl',
        com   => 'Hman-readable geneloc entries',
        view  =>
"
SELECT g.acc, g.sym, g.taxa, s.name, gl.lft, gl.rgt, gl.strand, gl.score,
       gl.geneloc_id, gl.sbj_id, gl.gene_id
  FROM gene g, geneloc gl, subject s
 WHERE g.gene_id = gl.gene_id
   AND s.sbj_id = gl.sbj_id
 ORDER BY g.sym, s.name
"
},

    v_fasta => {
        name  => 'v_fasta',
        com   => 'Report search_seq as fasta format',
        view  =>
"
SELECT '>' || name || COALESCE(' ' || descr,'') || '\n' || seq AS fasta,
  daterun, srch_id, mm, status 
  FROM search_seq ORDER BY daterun DESC;
"
},



    v_func => {
        name  => 'v_func',
        com   => 'Summarizes functions',
        db    => 'postgres',
        view  =>
"
SELECT proname || '(\n  ' || array_to_string(proargnames, ',\n  ') || ')' AS Function, 
       prosrc AS Source
  FROM pg_catalog.pg_namespace n
  JOIN pg_catalog.pg_proc p
    ON pronamespace = n.oid
 WHERE nspname = 'public';
"
},

    v_tab => {
        name  => 'v_tab',
        com   => 'Summarizes activities on tables',
        db    => 'postgres',
        view  =>
"
SELECT relname,  seq_scan, idx_scan, 
       n_tup_ins AS Inserts, n_tup_upd AS Updates, n_tup_del AS Deletes,
       to_char(last_analyze, 'YYYY Mon DD') AS Analyzed,
       to_char(last_vacuum, 'YYYY Mon DD') AS Vacuumed
  FROM pg_stat_all_tables where schemaname = 'public'
 ORDER BY relname
"
},

v_ind => {
    name  => 'v_ind',
    com   => 'Summarizes size and location of indices',
    db    => 'postgres',
    view  =>
"
 SELECT c.relname AS Index, tc.relname AS Table, 
        ts.spcname AS tablespace,
  c.relpages::double precision / 1000::double precision AS kilopages,
  floor(c.reltuples / 1000::double precision) AS kilotuples,
  pg_size_pretty(pg_total_relation_size(c.oid)) AS disk
   FROM pg_class tc, pg_namespace ns, pg_index ix, pg_class c
   LEFT OUTER JOIN pg_tablespace ts
     ON (c.reltablespace = ts.oid)
  WHERE ns.oid = c.relnamespace 
    AND ns.nspname NOT IN ('pg_catalog', 'information_schema', 'pg_toast')
    AND tc.oid     = ix.indrelid
    AND c.oid      = ix.indexrelid
  ORDER BY c.reltuples DESC
"
},

v_xid => {
    name  => 'v_xid',
    com   => 'Report on the transaction IDs for each table',
    db    => 'postgres',
    requires => [ 'queries' ],
    view  =>
"
 SELECT c.relname, ts.spcname AS tablespace,
  c.relpages::double precision / 1000::double precision AS kilopages,
  floor(c.reltuples / 1000::double precision) AS kilotuples,
  age(c.relfrozenxid)::double precision /1000000::double precision AS mega_xid,
  pg_size_pretty(pg_total_relation_size(c.relname::text)) AS disk
   FROM pg_namespace ns, pg_class c
   LEFT OUTER JOIN pg_tablespace ts
     ON (c.reltablespace = ts.oid)
  WHERE ns.oid = c.relnamespace
    AND ns.nspname NOT IN ('pg_catalog', 'information_schema', 'pg_toast')
    AND c.relkind  = 'r'
  ORDER BY c.reltuples DESC
"
},

v_wait => {
    name  => 'v_wait',
    com   => 'Find queries that are not immediately returning',
    db    => 'postgres',
    requires => [ 'queries' ],
    view  =>
"
SELECT count(queries.current_query) AS count,
       floor(100::double precision * (avg(date_part('minutes'::text, queries.query_age) * 60::double precision + date_part('seconds'::text, queries.query_age)) / 60::double precision)) / 100::double precision AS minutes,
        queries.current_query
   FROM queries
  GROUP BY queries.current_query
  ORDER BY floor(100::double precision * (avg(date_part('minutes'::text, queries.query_age) * 60::double precision + date_part('seconds'::text, queries.query_age)) / 60::double precision)) / 100::double precision DESC;
"
},

    v_size => {
        name  => 'v_size',
        com   => 'Show size of installed postgres databases',
        db    => 'postgres',
        view  =>
"
SELECT datid, datname, 
       pg_size_pretty(pg_database_size(datname)) AS size_on_disk
  FROM pg_stat_database
 ORDER BY pg_database_size(datname) DESC;
"
},

    queries => {
        name  => 'queries',
        com   => 'Shows Postgres SQL statements currently running for ALL databases',
        db    => 'postgres',
        view  =>
"
SELECT pg_stat_activity.datname, pg_stat_activity.usename,
       date_trunc('second', now() - pg_stat_activity.query_start) AS query_age,
       date_trunc('second', now() - pg_stat_activity.backend_start) AS backend_age, btrim(pg_stat_activity.query) AS current_query
   FROM pg_stat_activity
  WHERE pg_stat_activity.query <> '<IDLE>'::text
    AND pg_stat_activity.state = 'active'
  ORDER BY date_trunc('second', now() - pg_stat_activity.query_start),
           date_trunc('second', now() - pg_stat_activity.backend_start)
"
},

};
    
    return $tables;
}

sub prebranch {
    my $self = shift;
    printf("<pre>%s</pre>", $self->branch( @_ ));
}

sub user_col_prefix { return 'USER-'; }

sub standard_data_name {
    my $self = shift;
    my $req  = shift || "";
    my $colDet = $self->details('columns');
    if (my $stnd = $colDet->{__NAMELOOKUP}{lc($req)}) {
        return $stnd;
    }
    my $metaPrfx = $self->user_col_prefix();
    if ($req =~ /^\Q$metaPrfx\E/) {
        # This is a user-defined column, keep as-is
        return $req;
    } else {
        $self->msg_once("[?]","Unknown column '$req' requested");
        # warn "Stack Trace:";
        # die $self->branch($colDet);
        # warn "Foo";
        return "";
    }
}

sub get_metadata {
    my $self = shift;
    my ($type, $cname, $tag) = @_;
    return undef unless ($type && $cname);
    my $colDet = $self->details($type);
    my $rv = {};
    if (my $cName  = $self->standard_data_name( $cname )) {
        $rv = $colDet->{$cName} || {};
    } else {
        $self->msg_once("[?]","Unrecognized '$type' request for '$cname'");
    }
    if ($tag) {
        $rv = $rv->{lc($tag)};
    }
    return $rv;
}

sub data_index {
    my $self = shift;
    my ($data, $colReq, $addIfNeeded) = @_;
    return -1 unless ($data && $colReq);
    my $stnd = $self->standard_data_name( $colReq );
    return -1 unless ($stnd);
    my $lu = $data->{colH} ||= { map { 
        $self->standard_data_name($data->{cols}[$_]) => $_ + 1 }
                                 (0..$#{$data->{cols}})
    };
    # die $self->branch($data);
    if (my $num = $lu->{$stnd}) {
        return $num -1;
    } elsif ($addIfNeeded) {
        # Add the column to the end
        push @{$data->{cols}}, $stnd;
        my $sz = $#{$data->{cols}};
        $lu->{$stnd} = $sz + 1;
        $data->{need}{$stnd} = 1;
        # $self->msg("[+]","Column $colReq ($stnd) added at $sz");
        return $sz;
    } else {
        return -1;
    }
}

sub _clear_data_index_lookup {
    my $self = shift;
    my $data = shift;
    delete $data->{colH} if ($data);
}

sub custom_column_order {
    my $self = shift;
    my ($data, $req, $more) = @_;
    return ([], [], []) unless ($data);
    my $cols = $data->{cols};
    my $num  = $#{$cols};
    my (@names, @srcInds, %expanders);
    # If no special request, then just use all columns in data set:
    $req ||= $cols;
    my (%alreadyThere, @newColumn);
    foreach my $name (@{$req}) {
        my $stnd = $self->standard_data_name( $name );
        my $ind  = $self->data_index( $data, $name, 'add' );
        if ($ind == -1) {
            # This is not a recognized column name
            push @{$cols}, undef;
            $ind = $#{$cols};
        }
        $names[$ind] = $name;
        push @srcInds, $ind;
        $alreadyThere{$stnd}++;
    }
    # $self->prebranch({cols => $cols, request => $req, names => \@names, at => \%alreadyThere});
    if ($more) {
        # In addition to the default/requested columns, some are being
        # added to the end. This is useful if you want to just add a few
        # columns but otherwise do not want to disrupt the default
        foreach my $name (@{$more}) {
            my $stnd = $self->standard_data_name( $name );
            my $ind  = $self->data_index( $data, $name, 'add' );
            if ($ind == -1) {
                # Not a recognized column name
                push @{$cols}, undef;
                $ind = $#{$cols};
            } elsif ($alreadyThere{$stnd}++) {
                # Already present in output
                next;
            }
            $names[$ind] = $name;
            push @srcInds, $ind;
        }
    }
    
    # Check if we need to run any expander methods. New columns will
    # have been flagged with a 'need' tag
    # die $self->branch($data);
    while (my ($cname, $needFlag) = each %{$data->{need} || {}}) {
        next unless ($needFlag);
        # warn $cname;
        my $cb = $self->get_metadata('columns', $cname, 'expandCB');
        next unless ($cb);
        my $targ = $expanders{$cb.""} ||= {
            cb   => $cb,
            cols => [],
        };
        push @{$targ->{cols}}, $cname;
    }
    # For generation of CSS styles and Excel format names, some columns
    # specify a shorter alias for the column name. Set those as a parallel
    # name array
    my $clsNames = $data->{classNames} = [];
    my $colDet = $self->details('columns');
    foreach my $ind (@srcInds) {
        my $cname = $cols->[$ind];
        my $clsNm = $self->get_metadata('columns', $cname, 'classPrefix');
        $clsNames->[$ind] = $clsNm || $cname;
    }
    # warn $self->branch(\%expanders);
    return (\@names, \@srcInds, [values %expanders]);
}

sub json_encoded_columns {
    my $self = shift;
    my ($data) = @_;
    my %jsonInds;
    return \%jsonInds if (!$data);
    my @cols = @{$data->{cols}};
    for my $ind (0..$#cols) {
        if (my $cb = $self->get_metadata('columns', $cols[$ind], 'isRef')) {
            $jsonInds{$ind} = $cb;
        }
    }
    return \%jsonInds;
}

sub column_format_callbacks {
    my $self = shift;
    my ($data, $fmt, $colNames) = @_;
    # Build a hash keyed to column index and pointing to a callback method
    my %fmtInds;
    return \%fmtInds unless ($data && $fmt);
    my @cols = @{$data->{cols} || []};
    $fmt = lc($fmt);
    my $cbk = $fmt.'cb';
    for my $ind (0..$#cols) {
        my $cname = $cols[$ind];
        if (my $cb = $self->get_metadata( 'columns', $cname, $cbk )) {
            # There is a callback for this column in this format
            $fmtInds{$ind} = $cb;
            if (my $newName = $self->get_metadata
                ('columns', $cname, $fmt.'cbcol')
                || $self->get_metadata('columns', $cname, 'cbcol')) {
                # We should also rename the column if a new column name
                # is specified for the callback
                $colNames->[$ind] = $newName if ($colNames);
            }
        }
    }
    return \%fmtInds;
}

sub column_prefixed_value_callback {
    # Callback used to generate style or format names from controlled
    # vocabulary columns like strand or mismatch
    my $self = shift;
    my ($val, $data, $col) = @_;
    # The format name is the column name concatenated with the value
    return "" unless (defined $val);
    my $colNames = $data->{classNames} || $data->{cols} || [];
    my $cn = $colNames->[$col];
    return $cn ? $cn . $val : "";
}

sub _first_array_value_callback {
    # Callback used to generate style or format names from controlled
    # vocabulary columns like HitType that exist as an array.
    my $self = shift;
    my ($val, $data, $col, $row, $src) = @_;
    # If the cell is populated, return the first value
    $src  ||= $data->{rows};
    my $raw = $src->[$row][$col];
    return $raw ? $raw->[0] : "";
}

sub _gene_id_list_to_text_callback {
    my $self  = shift;
    my $valIn = shift;
    my $locSTH = $self->dbh->named_sth('Get Gene Information');
    my (@syms, @ids);
    foreach my $gid (@{$valIn || []}) {
        $locSTH->execute($gid);
        my $gRows = $locSTH->fetchall_arrayref();
        my ($acc, $sym) = @{$gRows->[0] || []};
        if ($sym) {
            if ($sym =~ /^LOC\d+$/) {
                # Really just the accession
                push @ids, $sym;
            } else {
                $sym =~ s/[\*\?]$//;
                push @syms, $sym;
            }
        } else {
            push @ids, $acc || 'UNK';
        }
    }
    return join(', ', (sort @syms), (sort @ids)) || "";
}

sub _join_array_comma_callback {
    my $self  = shift;
    my $valIn = shift;
    my @list;
    foreach my $val (@{$valIn}) {
        push @list, defined $val ? $val : '?';
    }
    return join (', ', @list) || "";
}

sub _serialize_json_callback {
    my $self  = shift;
    my $valIn = shift;
    return encode_json($valIn || "");
}


sub details {
    my $self = shift;
    my $type = shift;
    my $det   = $self->{DETAILS}{$type || ""};
    unless ($det && $det->{__FILEREAD}) {
        $det = $self->_get_or_create_details( $det, $type );
        my $colFile = $self->module_path( -module => $self,
                                          -suffix => "$type.conf");
        $self->read_conf_file( $colFile, $det );
        $det->{__FILEREAD} = 1;
    }
    return $det;
}

sub _get_or_create_details {
    my $self = shift;
    my ($obj, $type) = @_;
    if (!$obj) {
        $type ||= 'NotSpecified';
        return $self->{DETAILS}{$type} ||= {
            __TYPE       => $type,
            __NAMELOOKUP => {},
        };
    } elsif (!ref($obj)) {
        # Assume this is a type name
        return $self->_get_or_create_details( undef, $obj );
    }
    # Already instantiated
    return $obj;
}

sub read_conf_file {
    my $self = shift;
    my ($file, $detReq) = @_;
    my $det = $self->_get_or_create_details( $detReq );
    $det->{__SOURCE} = $file;
    unless ($file && -s $file) {
        $self->msg("[?]","Failed to read configuration file",
                   $file || "Unspecified file!", "Does not exist");
        return $det;
    }
    unless (open(FILE,"<$file")) {
        $det->{__ERROR} = $!;
        $self->msg("[?]","Failed to read configuration file", $file, $!);
        return $det;
    }
    my $colName = "";
    my %gather;
    while (<FILE>) {
        s/[\n\r]+$//;
        # Skip blank lines and those starting with #
        next if (/^\s*$/ || /^\#/);
        if (/^\s*(\S+)\s*$/) {
            # Treat a single word as a column / parameter name
            $colName = $1;
            next;
        }
        if (/^\s+(\S+):\s+(.+?)\s*?$/) {
            $gather{$colName}{lc($1)} = $2;
        } else {
            $self->msg("[?]","Oddly formatted configuration line '$_'");
        }
    }
    close FILE;
    while (my ($cn, $keyH) = each %gather) {
        while (my ($key, $val) = each %{$keyH}) {
            $self->reset_detail( $cn, $key, $val, $det, $keyH );
        }
    }
    return $det;
}

our $aliasKeys = { map { lc($_) => 1 } qw(dbcol alias)};

sub _subdetail_object {
    my $self = shift;
    my ($det, $subkey) = @_;
    my $rv = $det->{$subkey};
    unless ($rv) {
        my $skey = 'StndName';
        $rv = $det->{$subkey} ||= {
            $skey => $subkey,
        };
        $self->_set_detail_alias( $det, $subkey, $skey );
        # warn sprintf("Created %s : %s", $det->{__TYPE}, $subkey);
    }
    return $rv;
}

sub reset_detail {
    my $self = shift;
    my ($colName, $k, $v, $detReq, $otherKeys) = @_;
    return unless ($colName && $k && defined $v);
    $k = lc($k);
    my $det  = $self->_get_or_create_details( $detReq );
    my $targ = $self->_subdetail_object( $det, $colName );
    $otherKeys ||= {};
    
    # Special handle / safety check some keys:
    if ($k eq 'dbcol') {
        $v = lc($v);
        if ($targ->{$k} && $targ->{$k} ne $v) {
        } else {
            $targ->{$k} = $v;
        }
    } elsif ($k eq 'isref') {
        # The isRef tag flags columns that are normally complex objects
        # The value is a callback that returns an empty structure,
        # used when a null value is needed for the column
        $v = uc($v);
        if ($v eq 'ARRAY') {
            $targ->{$k} = sub { return []; }
        } elsif ($v eq 'HASH') {
            $targ->{$k} = sub { return {}; }
        } else {
            $self->msg("[?]","Column $targ->{StndName} has unrecognized isRef value of '$v'");
        }
    } elsif ($k eq 'alias') {
        # If commas are seen, split on them, othewise split on whitespace
        my $splitter = ($v =~ /,/) ? '\s*,\s*' : '[\s\t]+';
        $v = lc($v);
        $v =~ s/^\s+//;
        $v =~ s/\s+$//;
        foreach my $ali ( split(/$splitter/, $v)) {
            push @{$targ->{$k}}, $ali if ($ali);
        }
    } elsif ($k =~ /^(\S+)cb$/) {
        # Callback method for manipulating values, generally for output
        # but also for expanding one column into others
        my $what = $1;
        my $r = ref($v);
        my $cb;
        if ($r eq 'CODE') {
            # A code reference is being passed directly
            $cb = $v;
        } elsif ($v eq '$self') {
            # Return the value itself
            $cb = sub { my $self = shift; return shift; }
        } elsif ($cb = $self->can($v)) {
            # method is being defined by name
            unless ($v =~ /callback$/) {
                # Insist on only method names with 'callback' at the end
                # This is to prevent random calling of potentially destructive
                # or security related methods, like dbpass
                $self->msg("[!]","Callback for '$what' on column $targ->{StndName} requested for illegal method '$v' - method name must end with 'callback'");
                $cb = undef;
            }
        } elsif ($v =~ /^BIN\s+(.+)/) {
            # Request to return binned values based on numeric comparison
            my $binReq  = $1;
            my $nullVal = $otherKeys->{$what.'cbnull'} || "";
            my @bins;
            foreach my $br (split(/[\s,]+/, $binReq)) {
                next if (!defined $br || (!$br && $br ne '0'));
                my $bin = $br;
                my $val = $bin;
                if ($bin =~ /^(.+?):(\S+)/) {
                    ($bin, $val) = ($1, $2);
                }
                if ($bin =~ /^\-?(\d+|\d+\.\d+)$/) {
                    push @bins, [$bin, $val];
                } else {
                    $self->msg("[?]","Non-numeric bin value '$bin' on '$what' callback function for column $targ->{StndName}");
                }
            }
            if ($#bins != -1) {
                @bins = sort { $a->[0] <=> $b->[0] } @bins;
                $cb = sub {
                    my $self = shift;
                    my $val = shift;
                    return $nullVal if (!defined $val);
                    return $nullVal unless ($val =~ /^\-?(\d+|\d+\.\d+)$/);
                    for (my $i = 0; $i < $#bins; $i++) {
                        return $bins[$i][1] if ($val <= $bins[$i][0]);
                    }
                    return $bins[-1][1];
                };
            }
        } elsif ($v =~ /^\s*\'([^\']*)\'\s*$/ ||
                 $v =~ /^\s*\'([^\']*)\'\s*$/) {
            # A quoted value is treated as a literal, make a stub callback
            my $val = $1;
            $cb = sub { return $val; };
            # warn "$k = $val";
        }
        if ($cb) {
            # We have a code reference for the callback
            if ($what eq 'classname') {
                # BMS::ExcelHelper does not care what you name formats
                # But CSS does. Do some regexp to assure compliance
                # http://www.w3.org/TR/CSS21/grammar.html#scanner
                # https://stackoverflow.com/questions/448981/what-characters-are-valid-in-css-class-selectors
                $targ->{$k} = sub {
                    my $classname = &{$cb}( @_ ) || "";
                    # Allow only letters, numbers, underscores or dashes:
                    $classname =~ s/[^a-zA-Z0-9\-_]+//g;
                    # No numbers or dashes at the front. Technically a
                    # dash is sometimes ok but would require extra parsing
                    # so I am just going to ignore that possibility
                    $classname =~ s/^[0-9\-]//;
                    return $classname;
                };
            } else {
                $targ->{$k} = $cb;
            }
        } else {
            $self->msg("[?]","Failed to find '$what' callback function for column $targ->{StndName}");
        }
    } else {
        $targ->{$k} = $v;
    }
    $self->_set_detail_alias( $det, $colName, $k ) if ($aliasKeys->{$k});
}

sub _set_detail_alias {
    my $self = shift;
    my ($det, $k, $sk) = @_;
    my $src  = $self->_subdetail_object( $det, $k );
    my $stnd = $src->{StndName};
    my $vals = $src->{$sk};
    return unless ($vals);
    $vals  = [ $vals ] unless (ref($vals));
    my $lu = $det->{__NAMELOOKUP};
    $self->err("") unless ($lu);
    foreach my $ali (@{$vals}) {
        next unless ($ali);
        $ali = lc($ali);
        my $prior = $lu->{$ali};
        if (defined $prior) {
            if ($prior ne $stnd && $prior ne "") {
                # This alias is assigned to multiple stnds - not allowed
                $self->msg("[?]", "Alias $ali was pointing to both $prior and $stnd");
                $lu->{$ali} = "";
            }
        } else {
            $lu->{$ali} = $stnd;
            # warn "Alias $ali = $stnd\n";
        }
    }
}

sub _set_statement_data {
    my $self = shift;
    my $sthTxt = <<EOF;

# Get or create a src_id for a DB source and wordsize 2
SELECT source_db_id(?,?)

# Get a sbj_id for subject - create if absent 3
SELECT subject_id(?,?,?)

# Get a sbj_id for subject - no creation 3
SELECT sbj_id FROM subject WHERE name = ? AND src_id = ?

# Get a sbj_id for subject - no creation or source 3
SELECT sbj_id FROM subject WHERE name = ?

# Get a word from a wrd_id 3
SELECT word FROM word WHERE wrd_id = ?

# Get a wrd_id from a word - create if absent 3
SELECT quiet_word_write(?)

# Get a wrd_id from a word - no creation 3
SELECT wrd_id FROM word WHERE word = upper(?)


# Get src_ids for word 2
SELECT src_id FROM wordsrc WHERE wrd_id = ?

# Get all hits for a word 2
SELECT src_id, num, hits FROM wordsrc WHERE wrd_id = ?

# Get all hits for a word in a source
SELECT hits FROM wordsrc WHERE wrd_id = ? AND src_id = ?

# Set a hit for a word 2
SELECT quiet_wordsrc_write(?,?,?,?)

# Get subject information for sbj_id 1
SELECT name, src_id, len FROM subject WHERE sbj_id = ?

# Get sbj_id for name 1
SELECT sbj_id FROM subject WHERE upper(name) = upper(?)

# Get source information for src_id 1
SELECT src_id, path, wordsize, num, chars, mod_date FROM sourcedb
 WHERE src_id = ?

# Get source information for path and wordsize 1
SELECT src_id, path, wordsize, num, chars, mod_date FROM sourcedb
 WHERE path = ? AND wordsize = ?

# Set Gene Location details 2
SELECT quiet_geneloc_details(?,?,?,?,?,?,?)

# Get gene locations for subject position 2
SELECT geneloc_id, gene_id, strand, details FROM geneloc 
 WHERE sbj_id = ? AND lft < ? AND rgt > ?

# Get gene strands for subject position 3
SELECT DISTINCT strand FROM geneloc WHERE sbj_id = ? AND lft < ? AND rgt > ?

# Get gene_id for accession 2
SELECT quiet_gene_id(?)

# Set Gene Symbol 2
UPDATE gene SET sym = substring(? from 1 for 50) WHERE gene_id = ?

# Set Gene Taxa 2
UPDATE gene SET taxa = substring(? from 1 for 50) WHERE gene_id = ?

# Set Gene Description 2
UPDATE gene SET descr = ? WHERE gene_id = ?

# Set Gene Aliases 2
UPDATE gene SET aliases = ? WHERE gene_id = ?

# Get Gene Information 3
SELECT acc, sym, descr, taxa, aliases, hs_gene, hs_sym, hs_score
  FROM gene WHERE gene_id = ?

# Set Human Gene
UPDATE gene SET hs_gene = ? WHERE gene_id = ?
# Set Human Symbol
UPDATE gene SET hs_sym = ? WHERE gene_id = ?
# Set Human Score
UPDATE gene SET hs_score = ? WHERE gene_id = ?


# Get subjects for gene_id
SELECT sbj_id FROM geneloc WHERE gene_id = ?

# Get search_ids by sequence query 1
SELECT srch_id FROM search_seq
 WHERE upper(seq) = upper(?) AND mm >= ?
ORDER BY daterun desc;

# Get search details for srch_id 1
SELECT seq, name, descr, daterun, mm FROM search_seq WHERE srch_id = ?

# Delete a search by srch_id 1
DELETE FROM search_seq WHERE srch_id = ?

# Clear search hits 1
DELETE FROM search_hit WHERE srch_id = ?

# Get search details for name 1
SELECT seq, name, descr, daterun, mm FROM search_seq
 WHERE upper(name) = upper(?)
UNION
SELECT ss.seq, ss.name, ss.descr, ss.daterun, ss.mm
 FROM search_seq ss, tagval tv
WHERE upper(tv.val) = upper(?)
  AND tv.tag = 'Other ID'
  AND tv.obj_id = ss.srch_id


# Set search status 2
UPDATE search_seq SET status = ? WHERE srch_id = ?

# Get search status 2
SELECT status FROM search_seq WHERE srch_id = ?

# Create a new search 1
SELECT create_search_entry(?,?,?,?)

# Summarize hits for srch_id 1
SELECT src_id, mm, count(mm) FROM search_hit
 WHERE srch_id = ? AND mm <= ? AND genestrand != ?
 GROUP BY src_id, mm
 ORDER BY src_id, mm ASC

# Store a search hit 2
INSERT INTO search_hit (srch_id, subseq, mm, sbj_id, 
                        pos, strand, genestrand, footprint)
 VALUES (?,?,?,?,?,?,?)

# Set a tag/val pair 2 ignore=duplicate key value
SELECT quiet_tagval_write(?,?,?)

# Get tag/value pairs 2
SELECT tag, val FROM tagval WHERE obj_id = ?

# Clear a tag for object 2
DELETE FROM tagval WHERE obj_id = ? AND tag = ?

# Get values for a specific tag 2
SELECT val FROM tagval WHERE obj_id = ? AND tag = ?

# Get/Make exon type id 3
SELECT quiet_exontype_id( ? )

# Get exon type name for id 3
SELECT exontype FROM exon_type WHERE et_id = ?

# Find gene by ID or alias 1
SELECT gene_id FROM gene WHERE acc = ? OR ? = ANY (aliases)

# Get source ID for build 1
SELECT tv.obj_id FROM tagval tv, sourcedb s
 WHERE tag = 'GenomeBuild'
   AND upper(tv.val) = upper(?)
   AND tv.obj_id = s.src_id


EOF

    $self->dbh->parse_sth_text_block($sthTxt);
}
