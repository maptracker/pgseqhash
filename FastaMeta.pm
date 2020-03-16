# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use DB_File;

use BMS::ErrorInterceptor;
use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::Benchmark;
use Bio::DB::Fasta;

use vars qw(@ISA);
@ISA      = qw(BMS::ErrorInterceptor
               BMS::Utilities::SequenceUtilities
               BMS::Utilities::Benchmark);

our $strandLookup = {
    'F'  => 1,
    'R'  => -1,
    'FR' => 0,
    'RF' => 0,
    '0'  => 0,
    '-1' => -1,
    '+1' => 1,
    '1'  => 1,
    ''   => 1,
};

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
    };
    bless ($self, $class);
    my $args = $self->parseparams( @_ );
    if (my $fa = $args->{FASTA}) {
        $self->death("-fasta '$fa' does not exist") unless (-s $fa);
        my $m;
        if ($fa =~ /(.+)\.dbfile$/) {
            $m = $fa;
            $fa = $1;
        } else {
            $m = "$fa.dbfile";
        }
        $self->{FASTA} = $fa;
        if (-s $m) {
            $self->{META} = $m;
            my $tieHash = $self->{TIEHASH} = {};
            tie(%{$tieHash}, 'DB_File', $m, O_CREAT|O_RDWR, 0644, $DB_BTREE ) ||
                $self->death("Failed to tie DBFILE hash", $m, $!);
        } else {
            $self->msg("[!]", "No metadata file available for fasta", $fa);
        }
        $self->{DBFASTA} = Bio::DB::Fasta->new( $fa );
        $self->debug->skip_key([qw(DBFASTA TIEHASH)], 'global');
    } else {
        $self->death("Please provide -fasta (path to relevant fasta file)");
    }
    return $self;
}

*dbfile = \&meta;
*tie    = \&tiehash;
sub db      { return shift; }
sub fasta   { return shift->{FASTA}; }
sub meta    { return shift->{META}; }
sub dbfasta { return shift->{DBFASTA}; }
sub tiehash { return shift->{TIEHASH}; }

sub data {
    my $self = shift;
    unless ($self->{DATA}) {
        $self->bench_start();
        my $tie  = $self->tiehash();
        my $json = $tie->{ GLOBAL_DATA };
        $self->{DATA} = $json ? decode_json( $json ) : {};
        $self->bench_end();
    }
    return $self->{DATA};
}

sub val {
    my $self = shift;
    my $key  = shift;
    my $rv   = "";
    if ($key) {
        my $data = $self->data();
        # All keys should be uppercase and underscored, eg FASTA_FILE:
        $key     =~ s/\s+/_/g;
        $rv      = $data->{uc($key)};
        $rv      = "" unless (defined $rv);
    }
    return $rv;
}

sub numeric_strand {
    my $self = shift;
    my $str  = shift;
    return 0 unless (defined $str);
    return $strandLookup->{uc($str)} || 0;
}

sub fetch {
    my $self = shift;
    my $req   = shift;
    return undef unless ($req);
    $self->bench_start();
    my $class;
    my @try = ($req);
    if ($req =~ /(\d+)\^(\d+)/) {
        my ($l,$r) = ($1, $2);
        if ($r == $l + 1) {
            $class = 'BMS::FastaMeta::GenomeInsertion';
        } else {
            $class = 'BMS::FastaMeta::SpliceJunction';
        }
    } elsif ($req =~ /^(LOC|ENS[A-Z]*G)/i) {
        $class = 'BMS::FastaMeta::Gene';        
    } elsif ($req =~ /^([NX][MR]|ENS[A-Z]*T)/i) {
        $class = 'BMS::FastaMeta::RNA';
        if ($req =~ /(.+)\.(\d+)$/) {
            # Include unversioned ID
            push @try, $1;
        }
    } else {
        $class = 'BMS::FastaMeta::gDNA';        
    }
    foreach my $id (@try) {
        my $obj  = BMS::FastaMeta::Common->new
            ( -id => $id, -db => $self, -class => $class );
        next unless ($obj->json());
        $self->bench_end();
        return $obj;
    }
    $self->bench_end();
    return undef;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Common;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use Scalar::Util qw(weaken);
use parent qw( BMS::Utilities::SequenceUtilities
               BMS::Utilities::Benchmark
               BMS::ErrorInterceptor );

use JSON;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
    };
    bless ($self, $class);
    my $args = $self->parseparams( @_ );
    my $id   = $args->{ID};
    $id      = "" unless (defined $id);
    $self->{ID} = $id;
    if (my $db = $args->{DB}) {
        weaken( $self->{DB} = $db);
    } else {
        $self->death("Can not create object for '$id' without specifying -db",
                     "Please provide the base FastaMeta object");
    }
    if (my $creq = $args->{CLASS}) {
        bless($self, $creq);
    }
    return $self;
}

sub db      { return shift->{DB}; }
sub id      { return shift->{ID}; }
sub fastaid { return shift->{ID}; }
sub type    { return 'Unknown'; }
*dbfile = \&meta;
*tie    = \&tiehash;
sub fasta   { return shift->{DB}->fasta(); }
sub meta    { return shift->{DB}->meta(); }
sub dbfasta { return shift->{DB}->dbfasta(); }
sub tiehash { return shift->{DB}->tiehash(); }

sub acc  {
    my $self = shift;
    my $rv   = $self->id();
    if (my $vers = $self->val('v')) { $rv .= ".$vers"; }
    return $rv;
}

sub subseq_genomic {
    my $self = shift;
    my ($s, $e) = @_;
    if ($self->can('offset')) {
        if (my $off = $self->offset()) {
            $s -= $off;
            $e -= $off;
        }
    }
    return $self->subseq($s, $e);
}

sub subseq {
    my $self   = shift;
    $self->bench_start();
    my ($s, $e) = @_;
    my $id = $self->fastaid();
    my $ss = $self->dbfasta()->seq($id, $s, $e) || "";
    # $self->msg("[DEBUG]",sprintf("%s:%d..%d = %s", $id, $s, $e, $ss));
    $self->bench_end();
    return $ss;
}

sub json {
    my $self = shift;
    unless (defined $self->{JSON}) {
        my $tie = $self->tiehash();
        $self->{JSON} = $tie->{ $self->id() };
    }
    return $self->{JSON};
}

sub data {
    my $self = shift;
    unless ($self->{DATA}) {
        $self->bench_start();
        my $json = $self->json();
        $self->{DATA} = $json ? decode_json( $json ) : {};
        $self->bench_end();
    }
    return $self->{DATA};
}

sub val {
    my $self = shift;
    my $key  = shift;
    my $rv   = "";
    if ($key) {
        my $data = $self->data();
        $rv      = $data->{$key};
        $rv      = "" unless (defined $rv);
    }
    return $rv;
}

sub to_one_line {
    my $self = shift;
    my $rv   = sprintf("%s (%s)", $self->id(), $self->type());
    return $rv;
}

sub to_text {
    my $self = shift;
    return $self->to_one_line() . "\n";
}

sub supporting_json {
    my $self = shift;
    my %rv = ( $self->id() => $self->json);
    foreach my $method ('each_rna', 'each_gene', 'each_gdna') {
        if (my $cb = $self->can($method)) {
            foreach my $obj ( &{$cb}($self) ) {
                $rv{ $obj->id() } ||= $obj->json();
            }
        }
    }
    return \%rv;
}

sub _sublocs {
    my $self = shift;
    my $loc  = shift;
    return () unless ($loc);
    my @sl = $loc->can('sub_Location') ? $loc->sub_Location() : ($loc);
    return @sl;
}

sub location_text {
    my $self = shift;
    my @sls  = $self->_sublocs( @_ );
    return "" if ($#sls == -1);
    my @bits;
    foreach my $loc (@sls) {
        my ($s, $e) = ($loc->start(), $loc->end() );
        if ($s == $e) {
            push @bits, $s;
        } elsif ($loc->location_type() eq 'BETWEEN') {
            
        }
        
    }
    return join(',', @bits) || "";
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::BioObject;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::Common
               BMS::FastaMeta::Aligned);

*description = \&desc;
sub type   { return 'BioObject'; }
*symbol = \&sym;
sub sym    { return shift->val('sym'); }
sub desc   { return shift->val('desc'); }
sub taxa   { return shift->val('taxa'); }

sub to_one_line {
    my $self = shift;
    my $opts = shift || "";
    my $rv   = $self->id();
    if (my $sym = $self->sym()) {
        $rv .= " [$sym]";
    }
    # $rv .= sprintf(" (%s)", $self->type());
    if (my $fp = $self->footprint()) {
        $rv .= " $fp->[0]" unless ($opts =~ /nochr/);
        $rv .= sprintf(" %s..%s %s", map { $_ || '?' } @{$fp}[1..3]);
    }
    if (my $sc   = $self->score()) {
        $rv     .= sprintf(" %.1f%%", $sc);
    }
    if (my $desc = $self->desc()) {
        $rv .= " $desc";
    }
    return $rv;
}

sub to_text {
    my $self = shift;
    my $rv   = $self->to_one_line()."\n";
    return $rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Anchored;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sub anchor_id { return shift->chrBuild(); }

sub _loc_object {
    my $self = shift;
    my ($s, $e, $str, $ins) = @_;
    my $loc = BMS::FastaMeta::FMLoc->new
        ( -strand    => $str, 
          -ins_count => $ins,
          -start     => $s,
          -end       => $e );
    $loc->db_source( $self, $s, $e, $str, $ins );
    return $loc;
}

sub _finalize_locations {
    my $self = shift;
    my ($locs, $str) = @_;
    return () if ($#{$locs} == -1);
    my @rv;
    if ($#{$locs} != 0) {
        # Two or more locations, need to make a split
        # $locs = [ sort { $a->start() <=> $b->start() } @{$locs} ];
        # $locs = [ reverse @{$locs} ] if ($str < 0);
        my $split = BMS::FastaMeta::FMSplit->new( -splittype => 'join' );
        my @order = $str < 0 ?
            sort { $b->start() <=> $a->start() } @{$locs} :
            sort { $a->start() <=> $b->start() } @{$locs};
        foreach my $loc (@order) {
            $split->add_sub_Location( $loc );
        }
        $locs = [ $split ];
    }
    # $self->msg("[DEBUG]","_finalize_locations : ".$self->type().":".$self->id()." Str$str");
    my $aid = $self->anchor_id();
    map { $_->seq_id( $aid ) } @{$locs};
    return wantarray ? @{$locs} : $locs;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::FMSplit;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
use Bio::Location::Split;

@ISA      = qw(Bio::Location::Split);


sub new {
    my ($proto, @args) = @_;
    my $self   = Bio::Location::Split->new( @args );
    bless( $self, 'BMS::FastaMeta::FMSplit');
    return $self;
}

sub ignore {
    my $self = shift;
    my @slocs = $self->sub_Location();
    return $#slocs == -1 ? 1 : 0;
}

sub to_BMSstring {
    my $self = shift;
    my @bits;
    foreach my $loc ($self->sub_Location()) {
        push @bits, $loc->to_BMSstring();
    }
    @bits = reverse @bits if ($self->strand() < 0);
    my $rv = join(',', @bits);
    if (my $id = $self->seq_id()) {
        $rv = "$id:$rv";
    }
    if ($self->strand() < 0) {
        $rv .= '[-1]';
    }
    return $rv;
}

sub add_sub_Location {
    my $self = shift;
    my $loc  = shift;
    $self->SUPER::add_sub_Location($loc) unless ($loc->ignore());
    return $loc->is_subloc(1);
}

sub seq_id {
    my $self = shift;
    if (defined $_[0]) {
        $self->{SEQ_ID} = $_[0];
    }
    return $self->{SEQ_ID} || undef;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::FMLoc;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
use Bio::Location::Simple;

@ISA      = qw(Bio::Location::Simple);

sub new {
    my ($proto, @args) = @_;
    my $self   = Bio::Location::Simple->new( @args );
    bless( $self, 'BMS::FastaMeta::FMLoc');
    my ($ins) = $self->_rearrange([qw(INS_COUNT)],@args);
    if ($ins) {
        # $self->location_type('^');
        $self->ins_count( $ins );
    }
    return $self;
}

sub ignore { return shift->{IGNORE}; }

sub is_subloc {
    my $self = shift;
    if (defined $_[0]) {
        $self->{ISSUBLOC} = $_[0];
    }
    return $self->{ISSUBLOC} || 0;
}

sub ins_count {
    my $self = shift;
    if (defined $_[0]) {
        $self->{INSCOUNT} = $_[0];
    }
    return $self->{INSCOUNT} || 0;
}

sub db_source {
    my $self = shift;
    if ($_[0]) {
        $self->{DBSOURCE} = [ @_ ];
    }
    return @{$self->{DBSOURCE} || []};
}

sub to_BMSstring { 
    my ($self) = @_;
    my ($s, $e) = ($self->start, $self->end);
    my $rv = "";
    if ($s == $e) {
        $rv = $s;
    } elsif (my $ins = $self->ins_count()) {
        $rv = sprintf("%d^%d^%d", $s, $ins, $e);
    } else {
        $rv = sprintf("%d..%d", $s, $e);
    }
    unless ($self->is_subloc()) {
        if (my $id = $self->seq_id()) {
            $rv = "$id:$rv";
        }
        if ($self->strand() < 0) {
            $rv .= '[-1]';
        }
    }
    return $rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Aligned;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
use Bio::Location::Simple;
use Bio::Location::Fuzzy;

@ISA      = qw(BMS::FastaMeta::Common BMS::FastaMeta::Anchored);

sub fullchr { return shift->{FULLCHR}; }
sub se      { return shift->{SE}; }
sub str     { return shift->{STR}; }
sub score   { return shift->{SC}; }
sub offset  { return shift->{OFF}; }
sub _parse_fullchr {
    my $self = shift;
    my $fc   = $self->fullchr() || "";
    if ($fc =~ /^([^\.]+)\.([^\.]+)\.([^\.]+)\.(.+)$/) {
        my $tax = $1;
        $self->{GTYPE}     = $2;
        $self->{SHORTCHR}  = $3;
        $self->{BUILD}     = $4;
        $tax =~ s/_/ /g;
        substr($tax, 0, 1) = uc(substr($tax, 0, 1));
        $self->{BLDTAX}    = $tax;
    } else {
        $self->{GTYPE}     = "";
        $self->{SHORTCHR}  = "";
        $self->{BUILD}     = $self->val('build');
        $self->{BLDTAX}    = "";
    }
}

sub chr {
    my $self = shift;
    $self->_parse_fullchr() unless (defined $self->{SHORTCHR});
    return $self->{SHORTCHR};
}

*bld = \&build;
sub build {
    my $self = shift;
    unless (defined $self->{BUILD}) {
        $self->_parse_fullchr();
    }
    return $self->{BUILD};
}

sub chrBuild {
    my $self = shift;
    my ($chr, $bld) = ($self->chr(), $self->bld());
    return ($chr && $bld) ? "$chr.$bld" : $chr ? "Chr$chr" : "?.?";
}

sub footprint {
    my $self = shift;
    my $rv;
    if (my $fullchr = $self->fullchr()) {
        my $se = $self->se();
        if (!$se) {
            $se = [];
        } elsif (my $off = $self->offset()) {
            map { $_ -= $off } @{$se};
        }
        $rv = [ $fullchr, $se->[0], $se->[1], $self->str() ];
    }
    return $rv;
}

sub _donate_align_data {
    my $self = shift;
    my $targ = shift;
    foreach my $key (qw(SC STR FULLCHR OFF SE)) {
        $targ->{$key} = $self->{$key};
    }
    if (my $oneHsp = $self->{HSP}) {
        $targ->{HSP} = $oneHsp;
    }
    if ($targ->can('strand')) {
        $targ->strand($self->{STR});
    }
}

sub each_alignment {
    my $self = shift;
    if ($self->{HSP}) {
        # This object is already aligned to something
        return ($self);
    }
    my $hsps = $self->{HSPS};
    return () if (!$hsps || $#{$hsps} == -1);
    $self->bench_start();
    my $fullchr = $self->fullchr();
    my $off     = $self->offset();
    my $db      = $self->db();
    my @rv;
    foreach my $hsp (@{$hsps}) {
        my $rid = $hsp->{rna};
        my $rna = $db->fetch( $rid );
        next unless ($rna);
        my $str     = $hsp->{str};
        $rna->{SC}  = $hsp->{sc};
        $rna->{STR} = $str;
        $rna->{FULLCHR} = $fullchr;
        $rna->{OFF} = $off;
        $rna->{HSP} = $hsp->{hsp};
        $rna->{SE}  = $str < 0 ? 
            [$hsp->[-1][0], $hsp->[0][1]] : [$hsp->[0][0], $hsp->[-1][1]];
        push @rv, $rna;
    }
    $self->bench_end();
    return @rv;
}

sub each_genome_alignment_text {
    my $self = shift;
    my @rv;
    foreach my $aln ($self->each_alignment()) {
        my $at = "";
        if ($aln->chr()) {
            $at = $aln->chrBuild();
        } elsif ($at = $aln->fullchr()) {
        } else {
            $at = "UnkChr";
        }
        $at .= sprintf(":%d-%d %s", map { $_ || '?' } (@{$aln->{SE}}, $aln->{STR}));
        if (my $sc = $aln->{SC}) {
            $at .= sprintf(" %.1f%%", $sc);
        }
        my $hNum = $#{$aln->{HSP} || []} + 1;
        $at .= sprintf(" %d HSP%s", $hNum, $hNum == 0 ? '' : 's');
        push @rv, $at;
    }
    return @rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Gene;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::BioObject);

sub type { return 'Gene'; }
# sub strand { return 

sub each_gene { return (shift); }

sub each_gdna {
    my $self = shift;
    my $data = $self->data();
    my @rv;
    foreach my $gid (@{$data->{gdna} || []}) {
        my $gdna = $self->db->fetch( $gid );
        push @rv, $gdna if ($gdna);
    }
    return @rv;
}

sub each_rna {
    my $self = shift;
    my $data = $self->data();
    my $db   = $self->db();
    my @rv;
    if (my $hsps = $self->{HSPS}) {
        foreach my $hd (@{$hsps}) {
            my $rna = $db->fetch( $hd->{rna} );
            next unless ($rna);
            $self->_donate_align_data( $rna );
            my $hsp = $rna->{HSP} = $hd->{hsp};
            $rna->{SE}  = $hd->{str} < 0 ? 
                [$hsp->[-1][0], $hsp->[0][1]] : [$hsp->[0][0], $hsp->[-1][1]];
            push @rv, $rna;
        }
       # die $self->branch(\@rv);
    } else {
        foreach my $rid (@{$data->{rna} || []}) {
            my $rna = $db->fetch( $rid );
            push @rv, $rna if ($rna);
        }
    }
    # Sort by variant number
    my @sorter;
    foreach my $rna (@rv) {
        my $var = $rna->variant() || "";
        my ($m, $l, $r) = ($var, "", "");
        if ($m =~ /^([^\d]+)$/          || # Totally non-numeric : XYZ
            $m =~ /^([^\d]+)(.*?)$/     || # (X)(4)
            $m =~ /^(.*?[^\d])(\d.*?)$/    # (X)(4) , (427p)(1)
            ) {     # 
            ($l, $m) = ($1, $2 || "");
        }
        if ($m =~ /^(.*?)([^\d]+)$/) {
            ($m, $r) = ($1 || "", $2);
        }
        if ($m && $m =~ /[^\d]/) {
            $l = $l.$m;
            $m = "";
        }
        $m ||= 999;
        my $accNum = $self->id();
        $accNum    =~ s/\.\d+$//;
        $accNum    =~ s/[^\d]+//g;
        my $s = sprintf("%3s%03d%3s%015d", $l, $m, $r, $accNum);
        # warn sprintf("%20s %10s = '%s'\n", $self->id(), $var, $s);
        push @sorter, [$s, $rna];
    }
    @rv = map { $_->[1] } sort {$a->[0] cmp $b->[0] } @sorter;
    return @rv;
}

sub to_text {
    my $self = shift;
    my $rv = $self->to_one_line()."\n";
    foreach my $rna ($self->each_rna()) {
        $rv .= "  ".$rna->acc();
        if (my $v = $rna->variant()) {
            $rv .= " variant $v";
        }
        $rv .= "\n";
        foreach my $gdna ($rna->each_gdna()) {
            foreach my $at ($gdna->each_genome_alignment_text()) {
                $rv .= "    $at\n";
            }
        }
    }
    return $rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::RNA;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::BioObject BMS::FastaMeta::Anchored);

sub type { return 'RNA'; }

sub each_rna { return (shift); }

*len = \&length;
sub length {
    my $self   = shift;
    my ($s, $e) = @_;
    my $id = $self->fastaid();
    my $l  = $self->dbfasta()->length($id) || 0;
    return $l;
}

sub strand {
    my $self = shift;
    if (defined $_[0]) {
        $self->{STRAND} = $strandLookup->{$_[0]} || 0;
    }
    return $self->{STRAND} || 0;
}

*vers = \&version;
sub version  { return shift->val('v'); }

*idv = \&versioned_id;
sub anchor_id { return shift->versioned_id(); }
sub fastaid   { return shift->versioned_id(); }
sub versioned_id {
    my $self = shift;
    my $rv = $self->id();
    if (my $v = $self->version()) {
        $rv .= ".$v";
    }
    return $rv;
}


sub anchor_coordinates { return shift->rna_coordinates( @_ ); }
sub rna_coordinates {
    my $self = shift;
    my ($s, $e, $str) = @_; # 1-indexed character position of query
    $str    = $strandLookup->{$str || ''} || 1;
    my $loc = $self->_loc_object( $s, $e, $str);
    my @rv  = $self->_finalize_locations( [$loc], $str );
    return @rv;
}

sub each_overlapping_rna {
    # Return itself
    my $self = shift;
    my $hit = shift;
    return () unless ($hit);
    $self->exon_class( $self->type() );
    # The strand is the strand of the RNA on the parent
    # Since the parent is self-referential, the strand is 1
    $self->strand( 1 );
    return $self;
}

sub variant {
    my $self = shift;
    unless (defined $self->{VARIANT}) {
        my $d = $self->desc();
        my $iso = "";
        if ($d =~ /transcript variant ([\dA-Z]+)/i) {
            $iso = $1;
        }
        $self->{VARIANT} = $iso;
    }
    return $self->{VARIANT};
}

sub exon_class {
    my $self = shift;
    if (defined $_[0]) {
        $self->{EXON_CLASS} = $_[0];
    }
    return $self->{EXON_CLASS} || 0;
}

sub each_gdna {
    my $self = shift;
    my $rv   = $self->{GDNAS};
    return @{$rv} if ($rv);
    my $gene = $self->gene();
    $rv      = $self->{GDNAS} = [];
    return @{$rv} unless ($gene);
    my $id   = $self->id();
    my $db   = $self->db();
    # Find all the gDNA locations for the RNA's gene
    foreach my $gdna ($gene->each_gdna()) {
        # Just because a gene aligns to a stretch of gDNA does not mean
        # that all RNAs do. Make sure this RNA is actually at that location
        # Also, one RNA could be aligned multiple times to a location
        foreach my $rna ($gdna->each_rna()) {
            # Skip if it is not our RNA
            next if ($rna->id() ne $id);
            # print "\n\nRNA each_gdna(rna)\n".$self->branch($rna)."\n";
            # Make a fresh copy of the gDNA object:
            my $copy = $db->fetch( $gdna->id() );
            next unless ($copy);
            # Assign alignment data for this particular RNA
            $rna->_donate_align_data( $copy );
            push @{$rv}, $copy;
            # die $self->branch($copy);
        }
    }
    return @{$rv};
}

sub gene {
    my $self = shift;
    my $data = $self->data();
    my $gene = $self->db->fetch( $data->{gene} );
    return $gene;
}

sub each_gene {
    my $self = shift;
    my $gene = $self->gene();
    return $gene ? ($gene) : ();
}

sub to_one_line {
    my $self = shift;
    my $rv   = sprintf("%s (%s)", $self->versioned_id(), $self->type());
    if (my $v = $self->variant()) {
        $rv .= " Variant $v";
    }
    return $rv;
}

sub to_text {
    my $self = shift;
    my $rv = $self->to_one_line()."\n";
    if (my $gene = $self->gene()) {
        $rv .= "  ".$gene->to_one_line()."\n";
    } else {
        $rv .= "  /NO GENE/\n";
    }
    foreach my $gdna ($self->each_gdna()) {
        foreach my $at ($gdna->each_genome_alignment_text()) {
            $rv .= "    $at\n";
        }
    }
    return $rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::gDNA;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::Common
               BMS::FastaMeta::Aligned);

sub type { return 'gDNA'; }
sub fullchr  { return shift->val('chr'); }

sub ses {
    my $self = shift;
    unless ($self->{SES}) {
        my @ses;
        my $id = $self->id();
        if ($id =~ /(\d+)\-(\d+)\.([FR]+)$/) {
            @ses = ($1, $2, $3);
        } else {
            $self->msg_once("Unable to parse StartEndStrand", $self->id());
            die $self->branch($self);
        }
        $self->{SES} = \@ses;
    }
    return @{$self->{SES}};
}

sub length {
    my $self = shift;
    unless ($self->{LENGTH}) {
        my ($s,$e) = $self->ses();
        $self->{LENGTH} = $s ? $e - $s + 1 : 1
    }
    return $self->{LENGTH};
}

sub offset {
    my $self = shift;
    unless (defined $self->{OFF}) {
        my ($off) = $self->ses();
        $self->{OFF} = $off ? $off - 1 : 0;
    }
    return $self->{OFF};
}

sub each_gdna { return ( shift ); }

sub each_gene {
    my $self = shift;
    unless ($self->{GENES}) {
        $self->bench_start();
        my $data    = $self->data();
        my $fullchr = $self->fullchr();
        my $off     = $self->offset();
        my $db      = $self->db();
        my @rv;
        foreach my $gd (@{$data->{genes} || []}) {
            my $gobj         = $db->fetch( $gd->{id} );
            next unless ($gobj);
            $gobj->{SC}      = $gd->{sc};
            $gobj->{STR}     = $strandLookup->{ $gd->{str} || 0} || 0;
            $gobj->{FULLCHR} = $fullchr;
            $gobj->{SE}      = [ $gd->{s}, $gd->{e} ];
            $gobj->{OFF}     = $off;
            my $hsps = $gobj->{HSPS} = $gd->{hsps};
            # Normalize strand letters
            foreach my $hsp (@{$hsps || []}) {
                $hsp->{str} = $strandLookup->{ $hsp->{str} || 0} || 0;
            }
            push @rv, $gobj;
        }
        $self->{GENES} = \@rv;
        $self->bench_end();
    }
    return @{$self->{GENES}};
}

*citgc = \&char_index_to_genome_coordinate;
sub char_index_to_genome_coordinate {
    my $self = shift;
    my $ind  = shift;
    return undef if (!defined $ind);
    return $ind + $self->offset() + 1;
}

*anchor_coordinates = \&genome_coordinates;
sub genome_coordinates {
    my $self = shift;
    my ($s, $e, $str) = @_;
    my @rv;
    if (defined $s && defined $e) {
        $str = $strandLookup->{$str || ''} || 1;
        my $off = $self->offset();
        push @rv, $self->_loc_object( $off + $s, $off + $e, $str);
    }
    @rv = $self->_finalize_locations( \@rv, $str );
    return @rv;
}

sub to_text {
    my $self = shift;
    # die $self->json();
    my @genes = $self->each_gene();
    my $rv = sprintf("%s : %d Gene%s\n", $self->to_one_line(), $#genes + 1,
                     $#genes == 0 ? '' : 's');
    foreach my $gene (@genes) {
        $rv .= "  ".$gene->to_one_line( 'nochr')."\n";
    }
    return $rv;
}

sub each_rna {
    my $self = shift;
    unless ($self->{RNAS}) {
        $self->bench_start();
        my @rv;
        foreach my $gene ($self->each_gene()) {
            foreach my $rna ( $gene->each_rna() ) {
                push @rv, $rna;
                # $gene->_donate_align_data( $rna );
                $rna->strand( $rna->{STR} );
                # die $self->branch($rna);
            }
        }
        $self->{RNAS} = \@rv;
        $self->bench_end();
    }
    return @{$self->{RNAS}};
}

sub each_overlapping_rna {
    my $self = shift;
    $self->bench_start();
    my $hit     = shift;
    my $keepStr = $strandLookup->{ shift || 0} || 0;
    my $off     = $self->offset();
    my ($s, $e) = ($hit->start() + $off, $hit->end() + $off);
    my $qStr    = $hit->strand();
    my $pad     = $hit->pad();
    my @rv;
    # my $dbg = sprintf("Qry:[%d..%d]\n", $s, $e);
    foreach my $rna ($self->each_rna()) {
        my $rStr = $rna->strand();
        next if ($keepStr && $rStr * $qStr != $keepStr);
        my @hsps = sort { $a->[0] <=> $b->[0] } @{$rna->{HSP} || []};
        next if ($#hsps == -1);
        # Skip any locations that do not completely overlap:
        next if (($s + $pad) < $hsps[0][0] || ($e - $pad) > $hsps[-1][1]);
        # Find out if it is in an exon or not
        my $exCls = "Intron";

        my $hsp;
        for my $h (0..$#hsps) {
            $hsp = $hsps[$h];
            # HSP is to the right of hit:
            next if ($hsp->[1] < $s);
            # HSP is to the left of hit - nothing left to do:
            last if ($hsp->[0] > $e);
            # We overlap - how?
            if ($s >= $hsp->[0] && $e <= $hsp->[1]) {
                # Fully contained
                $exCls = "Exon";
            } else {
                # Partial overlap with the exon
                if ($h == 0 && $s < $hsp->[0]) {
                    # Dangling over the 'left' end
                    $exCls = $rStr < 0 ? "LastExon" : "FirstExon";
                } elsif ($h == $#hsps && $e > $hsp->[-1]) {
                    # Dangling over the 'right' end
                    $exCls = $rStr < 0 ? "FirstExon" : "LastExon";
                } else {
                    $exCls = "Splice";
                }
            }
            last;
        }
        # $dbg .= sprintf("  %s [%d..%d] %s\n%s", $exCls, $hsp->[0], $hsp->[1], $rna->id(), join('', map { "    $_->[0]..$_->[1]\n" } @hsps));
        $rna->exon_class($exCls);
        push @rv, $rna;
    }
    # warn "$dbg\n"; warn "\n";
    $self->bench_end();
    return @rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Interruption;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::Common
               BMS::FastaMeta::Aligned);

our $intronOffset   = 5;
sub type    { return 'Interruption'; }

sub fullchr {
    my $self = shift;
    unless (defined $self->{FULLCHR}) {
        my $par = $self->db->fetch( $self->val('par') );
        $self->{FULLCHR} = $par ? $par->fullchr() || "" : "";
    }
    return $self->{FULLCHR};
}

sub each_rna {
    my $self = shift;
    unless ($self->{RNAS}) {
        $self->bench_start();
        my $exons = $self->val('ex');
        my @rnas;
        foreach my $id (sort keys %{$exons}) {
            my $rna = $self->db->fetch( $id );
            next unless ($rna);
            $rna->strand( 1 );
            $rna->exon_class( $self->exon_class() );
            my $ex1 = $exons->{$id};
            push @rnas, $rna;
        }
        $self->{RNAS} = \@rnas;
        $self->bench_end();
    }
    return @{$self->{RNAS}};
#    die $self->branch($exons);
}

sub each_gene {
    my $self = shift;
    unless ($self->{GENES}) {
        my %genes;
        foreach my $rna ($self->each_rna()) {
            if (my $gene = $rna->gene()) {
                $genes{$gene->id()} ||= $gene;
            }
        }
        $self->{GENES} = [ values %genes ];
    }
    return @{$self->{GENES}};
}

sub each_overlapping_rna {
    my $self = shift;
    my $hit     = shift;
    my $keepStr = $strandLookup->{ shift || 0} || 0;
    my @rv;
    foreach my $rna ($self->each_rna()) {
        next if ($keepStr && $rna->strand() * $hit->strand() != $keepStr);
        push @rv, $rna;
    }
    return @rv;
}

sub lrspw {
    my $self = shift;
    unless ($self->{LRSPW}) {
        # L = Genome coordinate just to left (RNA perspective) of exon
        # R = Genome coordinate to right of exon
        # S = Strand of RNA on genome
        #     NOT CORRECT - S always equals 'F'
        # P = 1-indexed coordinate of just-to-left character in fasta record 
        # W = Width of Left-Right boundary in fasta file
        #     2 = normal = RNA is nicely aligned to just exons
        #     7 = there are 5 additional RNA bases between the exons
        my @lrspw;
        my $id = $self->id();
        if ($id =~ /(\d+)\^(\d+)(?:\.\d+)?\.([FR]+)$/) {
            @lrspw = ($1, $2, $strandLookup->{$3});
            my ($p, $w) = ($self->val('pos'), $self->val('w') || 2 );
            if ($p && $w) {
                push @lrspw, ($p, $w);
            } else {
                $self->msg_once("Unable to parse PositionWidth", $self->id());
            }
        } else {
            $self->death("Unable to parse LeftRightStrand", $self->id());
        }
        $self->{LRSPW} = \@lrspw;
    }
    return @{$self->{LRSPW}};
}

sub anchor_coordinates { return shift->genome_coordinates( @_ ); }
sub genome_coordinates {
    my $self = shift;
    $self->bench_start();
    my ($s, $e, $str) = @_; # 1-indexed character position of query
    $str = $strandLookup->{$str || ''} || 1;
    my ($l, $r, $sStr, $p, $w) = $self->lrspw();
    # $str *= $sStr;
    my $lOff = $l - $p;
    my @rv;
    # warn "($s, $e, $str) : ($l, $r, $sStr, $p, $w)\n";
    if (!defined $s || !defined $e) {
        # Eh. Can not do anything
    } elsif ($e <= $p) {
        # The coordinate does not span the splice junction, left side only
        push @rv, $self->_loc_object( $lOff + $s, $lOff + $e, $str);
        # We IGNORE this result. Hit should be caught by pre-mRNA entries
        $rv[-1]{IGNORE} = 'Redundant to pre-mRNA';
    } else {
        my $pw   = $p + $w - 1; # 1-index char pos of start of left exon
        my $rOff = $r - $pw;
        if ($s >= $pw) {
            # The coordinate does not span the splice junction, right side only
            push @rv, $self->_loc_object( $rOff + $s, $rOff + $e, $str);
            # We IGNORE this result. Hit should be caught by pre-mRNA entries
            $rv[-1]{IGNORE} = 'Redundant to pre-mRNA';
        } elsif ($s > $p && $e < $pw) {
            # Coordinates are fully between the two exons
            # Normally this would be caught by the pre-mRNA entries
            # In this case it implies some unaligned sequence in the intron
            my $start = $lOff + $intronOffset;
            push @rv, $self->_loc_object( $l, $r, $str, $e - $s + 1 );
        } else {
            # At least part of the span is crossing the exon-exon boundary
            if ($s <= $p) {
                # There is a segment on the left
                push @rv, $self->_loc_object( $lOff + $s, $lOff + $p, $str);
            }
            if ($w > 2) {
                # At least part of the span is in between the two exons
                # We need to find how much
                my $loc = $self->_interrupted_location( $s, $e, $str );
                push @rv, $loc;
            }
            if ($e >= $pw) {
                # There is a segment on the right
                push @rv, $self->_loc_object( $rOff + $pw, $rOff + $e , $str);
            }
        }
    }
    @rv = $self->_finalize_locations( \@rv, $str );
    $self->bench_end();
    return wantarray ? @rv : \@rv;
}

sub _interrupted_count {
    # Count the number of characters falling between genome blocks
    my $self = shift;
    my ($s, $e) = @_; # 1-indexed character position of query
    my ($l, $r, $sStr, $p, $w) = $self->lrspw();
    my $pw    = $p + $w - 1; # 1-index char pos of start of left exon
    my $is    = $s <= $p  ? $p  + 1 : $s;
    my $ie    = $e >= $pw ? $pw - 1 : $e;
    return $ie - $is + 1;
}

sub _interrupted_location {
    my $self = shift;
    my ($s, $e, $str) = @_; # 1-indexed character position of query
    my $len = $self->_interrupted_count( $s, $e );
    return () unless ($len);
    my ($l, $r, $sStr, $p, $w) = $self->lrspw();
    return $self->_loc_object( $l, $r, $str, $len);
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::SpliceJunction;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::Interruption );

sub type       { return 'SpliceJunction'; }
sub exon_class { return 'ExonExon'; }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::GenomeInsertion;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use vars qw(@ISA);
@ISA      = qw(BMS::FastaMeta::Interruption );

sub type       { return 'GenomeInsertion'; }
sub exon_class { return 'ExonInDel'; }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FastaMeta::Hit;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use BMS::ErrorInterceptor;
use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::Benchmark;

use vars qw(@ISA);
@ISA      = qw(BMS::ErrorInterceptor
               BMS::Utilities::SequenceUtilities
               BMS::Utilities::Benchmark);
sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
    };
    bless ($self, $class);
    $self->bench_start();
    my $args = $self->parseparams( @_ );
    $self->query_seq( $args->{QUERYSEQ} );
    $self->query_id( $args->{QUERYID} );
    $self->query_desc( $args->{QUERYDESC} );
    $self->subject( $args->{SUBJECT} );
    $self->start( $args->{START} );
    $self->end( $args->{END} );
    $self->strand( $args->{STRAND} );
    $self->{PAD} = $args->{PAD} || 0;
    $self->bench_end();
    return $self;
}

# Allows for slight overhang on edges of RNAs
sub pad { return shift->{PAD}; }

sub anchor_locations {
    my $self = shift;
    unless ($self->{ANC_LOC}) {
        $self->bench_start();
        if (my $subObj = $self->subject()) {
            my ($s, $e, $str) =($self->start(), $self->end(), $self->strand());
            # warn $subObj->to_one_line()."\n";
            my @locs = $subObj->anchor_coordinates($s, $e, $str);
            if ($#locs != 0) {
                my $id = $subObj->id();
                $self->msg("[!]", ($#locs == -1 ? "Failed to get genome location" : "Multiple genome locations")." for $id $s..$e [$str]");
            }
            $self->{ANC_LOC} = \@locs;
        }
        $self->bench_end();
    }
    return @{$self->{ANC_LOC} || []};
}

sub each_overlapping_rna {
    my $self = shift;
    if (my $subObj = $self->subject()) {
        return $subObj->each_overlapping_rna($self, @_);
    }
    return ();
}

sub row_data {
    my $self = shift;
    my $gTxt    = $self->genomic_footprint();
    return () if ($gTxt =~ /^IGNORE/);
    $self->bench_start();
    my $keepStr = $strandLookup->{ $_[0] || 0} || 0;
    my $qid     = $self->query_id();
    my $qdesc   = $self->query_desc();
    my $qseq    = $self->query_seq();
    my ($s, $e) = ($self->start(), $self->end());
    my $qstr    = $self->strand();
    my $subObj  = $self->subject();
    my $sid     = $subObj->fastaid();
    my $bld     = $subObj->build();
    my $sseq    = $self->subject_seq();
    my $mm      = $self->mismatch();

    my %byGene;
    foreach my $rna ($self->each_overlapping_rna( $keepStr)) {
        # warn "  ".$rna->id()."\n";
        my $gene = $rna->gene();
        next unless ($gene);
        my $gid  = $gene->id();
        my $desc = $gene->desc();
        # We expect all the RNAs to be the same strand
        # Inverted duplications could result in overlapping assignments
        # to both strands:
        my $rstr = $rna->strand() * $qstr;
        my $gdat = $byGene{$gene->id()}{$rstr} ||= {
            rnas => {},
            excl => {},
            row  => {
                OligoID         => $qid,
                GeneID          => $gid,
                HitID           => $sid,
                Build           => $bld,
                Symbol          => $gene->sym(),
                Start           => $s,
                End             => $e,
                Str             => $rstr,
                GenomeFootprint => $gTxt,
                OligoSeq        => $qseq,
                OligoDesc       => $qdesc,
                GeneDesc        => $desc,
                Mis             => $mm,
                HitSeq          => $sseq,
            },
        };
        # RNA specific information structured in hashes:
        $gdat->{rnas}{ $rna->id() }++;
        $gdat->{excl}{ $rna->exon_class() }++;
    }
    my @rv;
    foreach my $strDat (values %byGene) {
        foreach my $gdat (values %{$strDat}) {
            my $row        = $gdat->{row};
            # Join RNA specific information into single string
            $row->{RnaID}  = join(' | ', sort keys %{$gdat->{rnas}});
            $row->{Exonic} = join(' | ', sort keys %{$gdat->{excl}});
            push @rv, $row;
        }
    }
    $self->bench_end();
    return wantarray ? @rv : \@rv;
}

sub genomic_footprint {
    my $self = shift;
    my @gLocs   = $self->anchor_locations();
    if ($#gLocs == -1) {
        return "IGNORE : Unanchored";
    } elsif (my $ign = $gLocs[0]->ignore()) {
        return "IGNORE : $ign";
    }
    return join(' / ', map { $_->to_BMSstring() } @gLocs);
}

sub query_seq {
    my $self = shift;
    if ( $_[0]) {
        $self->{QUERY_SEQ} = $_[0];
    }
    return $self->{QUERY_SEQ} || "";
}

sub subject {
    my $self = shift;
    if ($_[0]) {
        $self->{SUBOBJ} = $_[0];
    }
    return $self->{SUBOBJ} || undef;
}

sub _seq_details {
    my $self = shift;
    unless ($self->{SEQDETAILS}) {
        $self->bench_start();
        my $seq   = "";
        my $qseq  = $self->query_seq();
        my @qchar = split('', uc($qseq));
        my $str   = $self->strand();
        my @mm;
        if ( my $subObj  = $self->subject()) {
            my ($s, $e) = ($self->start, $self->end);
            my $prePad = "";
            if ($s < 1) {
                # There is a left side overhang on the sequence
                # We will represent it with gap characters
                $prePad = '-' x (1 - $s);
                $s = 1;
            }
            $seq = $subObj->subseq( $s, $e );
            $seq = $prePad . $seq if ($prePad);
            if (my $pad = $#qchar + 1 - CORE::length($seq)) {
                # Not enough characters in the sequence
                # This could be a right side overhang, or something bad
                # has happened with sequence recovery. Bio::DB::Fasta
                # sometimes has corrupted indices.
                $seq .= ('-' x $pad);
            }
            if ($str < 0) {
                $seq = $subObj->revcom( $seq );
            }
            my @schar = split('', uc($seq));
            for my $i (0..$#qchar) {
                my $qc = $qchar[$i];
                my $sc = $schar[$i];
                push @mm, $i if (!$sc || $sc ne $qc);
            }
        } else {
            $self->msg("[!!]", "Failed to find subject!");
            @mm = (0..$#qchar);
        }
        # Note mismatches in lowercase
        my $misSeq = uc($seq);
        map { substr($misSeq, $_, 1) = lc(substr($misSeq, $_, 1)) } @mm;
        $self->{SEQDETAILS} = [ $misSeq, \@mm, $seq ];
        $self->bench_end();
    }
    return @{$self->{SEQDETAILS}};
}

sub subject_seq {
    my $self = shift;
    my ($seq) = $self->_seq_details();
    return $seq;
}

sub mismatch {
    my $self = shift;
    my ($seq, $mm) = $self->_seq_details();
    return wantarray ? @{$mm} : $#{$mm} + 1;
}

sub start {
    my $self = shift;
    if (defined $_[0]) {
        $self->{START} = $_[0];
    }
    return $self->{START};
}

sub end {
    my $self = shift;
    if (defined $_[0]) {
        $self->{END} = $_[0];
    }
    return $self->{END};
}

sub strand {
    my $self = shift;
    if (defined $_[0]) {
        $self->{STRAND} = $strandLookup->{$_[0]};
    }
    return $self->{STRAND} || 0;
}

sub query_id {
    my $self = shift;
    if ($_[0]) {
        $self->{QUERY_ID} = $_[0];
    }
    return $self->{QUERY_ID} || "";
}

sub query_desc {
    my $self = shift;
    if (defined $_[0]) {
        $self->{QUERY_DESC} = $_[0];
    }
    return $self->{QUERY_DESC} || "";
}

sub foo {
    my $self = shift;
    if (defined $_[0]) {
        $self->{FOO} = $_[0];
    }
    return $self->{FOO} || "";
}

