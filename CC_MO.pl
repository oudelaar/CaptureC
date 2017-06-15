#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use String::Approx "amatch";

# This script performs the analysis of a Capture-C experiment. It takes as input the aligned sam file, the oligo file (CCanalyser format), and a file with the
# coordinates of all the restriction fragments in the genome of the restriction enzyme used for digestion (format chr:start-stop\n).
# For more info on file formats and how to generate these, see http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/UserManualforCaptureCanalysis.pdf.

# The script analyses the sam file and filters out:
# 1. Duplicate reads. For flashed reads duplicate removal is based on a unique start and stop of all mapped fragments; for non-flashed reads only the start of the
#    two separate reads is considered and needs to be unique. The script has an option for more stringent analysis, in which reads with only 1 (-stringent 1) or
#    2 (-stringent 2) differences are also removed.
# 2. Proximity exclusion. All mapped restriction fragments that fall within a 1kb window from any of the viewpoints used will be removed.
# The remaining read fragments are considered reporters and are mapped to the corresponding restriction fragments. 
# The script contains a 'globin option' to combine the separate globin tracks into a combined track.

# The script is dependent on the ucsctools module. It is also necessary to specify a folder containing a file with the size of the genome ($bigwig_folder).

# Output:
# 1. Report with basic statistics: total aligned reads, duplicates, proximity exclusion, reporters in cis and trans
# 2. Oligo-specific txt file with reporter fragments in the format "read_name \t frag1-coord_frag2-coord_..."  
# 3. Oligo-specific wig file
# 4. Oligo-specific gff file

# The script takes 5-10 min to run when using the default duplicate filtering. Stringent duplicate filtering is very slow and will take overnight for samples
# sequenced on a miseq (sam file up to 3 Gb). This is probably only required when low cell numbers are used. 

# Example of run command:
# module load ucsctools
# perl CC_MO.pl -sam your_input_sam_file.sam -o your_input_oligo_file.txt -r mm9_dpnII_coordinates.txt -pf /public/oudelaar/ -name short_sample_name &

# contact: oudelaar@well.ox.ac.uk

#specify
my $public_folder = "/public/oudelaar/";
my $public_url = "sara.molbiol.ox.ac.uk/public/oudelaar/";
my $bigwig_folder = "/t1-data/user/config/bigwig/";
my $email = 'oudelaar@well.ox.ac.uk';

my $genome = "mm9";
my $globin = 0;
my $stringent = 0;

&GetOptions
(
	"sam=s" =>\ my $sam_file, 	    # -sam          sam file
    "o=s" => \my $oligo_file,       # -o            oligo file (CCanalyser format)           
    "r=s" =>\ my $dig_genome,       # -r            file with restriction coordinates in genome 
    "name=s" =>\ my $input_name,    # -name         name of the experiment, can be whatever you like, will be used for names of output files
    "genome=s" =>\ $genome,         # -genome       genome, eg mm9 or hg18
    "pf=s" =>\ $public_folder,      # -pf		    your public folder (e.g. /public/username/)
    "pu=s" =>\ $public_url,         # -pu           your public url (e.g. sara.molbiol.ox.ac.uk/public/username/)
    "globin=i" =>\ $globin,         # -globin       1 = combine Hba-1 and Hba-2; 2 = also combine Hbb-b1 and Hbb-b2
    "stringent=i" =>\ $stringent,   # -stringent    default / 0 = normal duplicate removal, 1 = stringent duplicate removal (1 letter difference),
                                    #               2 = superstringent duplicate removal (2 letters difference) 
);

#open filehandles
open (SAMFH, $sam_file) or die "can't open sam file $sam_file";
open (OLIFH, $oligo_file) or die "can't open oligo file $oligo_file";
open (GENFH, $dig_genome) or die "can't open genome file with restriction coordinates $dig_genome";

#generate output directory and files
my $name = "undef";
if ($stringent == 0) {
    $name = "$input_name\_CC";
}
if ($stringent == 1) {
    $name = "$input_name\_CC_str1";
}
if ($stringent == 2) {
    $name = "$input_name\_CC_str2";
}

my $path = "undef"; 
if ($sam_file =~ /(.*\/)(\V++)/) {
    $path = $1;
    }
my $dir = "$path/$name";
unless (-d $dir) {
    mkdir $dir;
}

my $report = "$dir/$name\_report.txt";
open (REP, ">$report") or die "can't open output report file $report";

#store oligo coordinates in hash
my %oligos;
while (my $line = <OLIFH>) {
    chomp $line;
    my ($id, $chr1, $start1, $stop1, $chr2, $start2, $stop2, $SNP_loc, $SNP_base) = split (/\t/, $line);
    $oligos{$id}{"chr"}=$chr1;
    $oligos{$id}{"start"}=$start1;
    $oligos{$id}{"stop"}=$stop1;
    $oligos{$id}{"start_prox"}=$start2;
    $oligos{$id}{"stop_prox"}=$stop2;
}

#store restriction fragment coordinates in hash and sort in ascending order, to use for binary search
my %RE_hash;
while (my $line = <GENFH>) {
    chomp $line;
    my ($chr, $start, $stop) = split (/\W/, $line);
    push @{$RE_hash{$chr}}, $start;
}

my @chr_index = keys %RE_hash;
foreach my $chr (@chr_index) {
    @{$RE_hash{$chr}} = sort {$a <=> $b} @{$RE_hash{$chr}};
}

#store read info (name, pair-mate, fragment number [after in silico digest]) and coordinates in hash, and check if the read contains a capture
my %data_hash; my %counter; 
while (my $line = <SAMFH>) {
    chomp $line;
    my ($name, $flag, $chr, $start, $map_qual, $cigar, $mate_ref, $mate_start, $mate_insert, $seq, $seq_qual, $opt1, $opt2, $opt3) = split (/\t/, $line);
    my $read_name = "undefined"; my $stop = 0;
    if ($chr =~ /chr.*/ and $name =~ /(.*):PE(\d++):(\d++)$/) {
        $read_name = $1; my $mate = $2; my $frag_nr= $3;
        $chr =~ /chr(.*)/; $chr = $1;
        $data_hash{$read_name}{$mate}{$frag_nr}{"chr"} = $chr;
        $data_hash{$read_name}{$mate}{$frag_nr}{"start"} = $start;
    if ($cigar =~/(\d++)(.*)/) {
        $stop = $1 + $start;                                                        #first number in cigar string ($1) is alignment length; if case of indels, the part of the read until first indel 
        $data_hash{$read_name}{$mate}{$frag_nr}{"stop"}= $stop;                     #is used to map the read to the restriction fragment and to generate the coordinate string for duplicate removal
        }
    foreach my $id (keys %oligos) {                                                 #check if read contains a capture
        if ($chr eq $oligos{$id}{"chr"} and $start >= $oligos{$id}{"start"} - 1 and $stop <= $oligos{$id}{"stop"} + 1) {
            $data_hash{$read_name}{"label"} = $id;
            }
        }
    }
}

#generate coord strings for duplicate removal
foreach my $read_name (keys %data_hash) {
    for (my $mate = 1; $mate < 3; $mate++) { 
        foreach my $frag_nr (sort {$data_hash{$read_name}{$mate}{$a} <=> $data_hash{$read_name}{$mate}{$b}} keys %{$data_hash{$read_name}{$mate}}) {
            my $chr = $data_hash{$read_name}{$mate}{$frag_nr}{"chr"};
            my $start = $data_hash{$read_name}{$mate}{$frag_nr}{"start"};
            my $stop = $data_hash{$read_name}{$mate}{$frag_nr}{"stop"};
            if (exists($data_hash{$read_name}{"2"})) {              #non-flashed reads: only use start of the two reads for duplicate removal
                if (exists($data_hash{$read_name}{$mate}{"0"})) {   #check if first part of the read (frag0) is mapped, if not use second (frag1) or third (frag2) part
                    if ($frag_nr == 0) {
                        $data_hash{$read_name}{"coord"} .= $chr.":".$start."_";                 
                    }
                }   
                elsif (exists($data_hash{$read_name}{$mate}{"1"})) {
                    if ($frag_nr == 1) {
                        $data_hash{$read_name}{"coord"} .= $chr.":".$start."_";
                    }
                }
                elsif (exists($data_hash{$read_name}{$mate}{"2"})) {
                    if ($frag_nr == 2) {
                        $data_hash{$read_name}{"coord"} .= $chr.":".$start."_";
                    }
                }
            }
            else {                                                  #flashed reads: use start and stop of all mapped fragments for duplicate removal
                $data_hash{$read_name}{"coord"} .= $chr.":".$start."-".$stop."_";
            }
        }
    }
}

#remove duplicates
my %coord_hash;             
foreach my $key (keys %data_hash) {
    $counter{"1Total aligned reads"}++;
    if (exists $coord_hash{$data_hash{$key}{"coord"}}) {
        delete $data_hash{$key};
    }
    else {
    $coord_hash{$data_hash{$key}{"coord"}} = 1;
    }
}

if ($stringent == 1) {                                      #remove duplicates with one letter difference in coordinate string
    foreach my $key (keys %data_hash) {
        my $coordinate = $data_hash{$key}{"coord"};
        LABEL:
        foreach my $stored_coord (keys %coord_hash) {
            if (length($coordinate) == length($stored_coord) and $coordinate ne $stored_coord and amatch ($coordinate, ['1'], $stored_coord)) {
                delete $data_hash{$key};
                last LABEL;
            }
        }
    }
}

if ($stringent == 2) {                                      #remove duplicates with two letters difference in coordinate string
    foreach my $key (keys %data_hash) {
        my $coordinate = $data_hash{$key}{"coord"};
        LABEL:
        foreach my $stored_coord (keys %coord_hash) {
            if (length($coordinate) == length($stored_coord) and $coordinate ne $stored_coord and amatch ($coordinate, ['2'], $stored_coord)) {
                delete $data_hash{$key};
                last LABEL;
            }
        }
    }
}
    
$counter{"2Unique aligned reads"} = scalar keys %data_hash;
$counter{"3Duplicated aligned reads"} = $counter{"1Total aligned reads"} - $counter{"2Unique aligned reads"};

#map read to restriction fragment using the binary_search subroutine and delete fragments that fall within the proximity exclusion
my %frag_hash;
my %frag_counter;
my $chr = 0;
my $frags_counted = 0;
foreach my $read_name (keys %data_hash) {
    $frags_counted = 0;
    for (my $mate = 1; $mate < 3; $mate++) {
        foreach my $frag (sort {$data_hash{$read_name}{$mate}{$a} <=> $data_hash{$read_name}{$mate}{$b}} keys %{$data_hash{$read_name}{$mate}}) {   
            $frags_counted++;
            foreach my $oligo (keys %oligos) {
            $chr = $oligos{$oligo}{"chr"};
            
            #capture fragments
            if (exists ($data_hash{$read_name}{$mate}{$frag}{"chr"}) and $data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $data_hash{$read_name}{$mate}{$frag}{"start"} >= $oligos{$oligo}{"start"} - 1 and $data_hash{$read_name}{$mate}{$frag}{"stop"} <= $oligos{$oligo}{"stop"} + 1) {
                $data_hash{$read_name}{$mate}{$frag}{"label"} = "capture $oligo";
                $counter{"4Unique reads with capture $oligo"} ++; 
            }
            
            #proximity exclusion                                                                                               
            else {
                if (exists ($data_hash{$read_name}{$mate}{$frag}{"stop"})) {
                    my ($start_frag, $end_frag) = binary_search(\@{$RE_hash{$chr}}, $data_hash{$read_name}{$mate}{$frag}{"start"}, $data_hash{$read_name}{$mate}{$frag}{"stop"}, \%counter);                    #map read to restriction fragment
                    unless ($start_frag =~ /error/) {
                        if ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag <= $oligos{$oligo}{"start_prox"} - 1 and $end_frag >= $oligos{$oligo}{"start_prox"} - 1) {    #proximity exclusion if end of fragment is in proximity zone 
                            $counter{"5aProximity excluded fragments total"}++;
                            $counter{"5bProximity excluded fragments $oligo"}++;
                            delete $data_hash{$read_name}{$mate}{$frag};
                            }
                        elsif ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag >= $oligos{$oligo}{"start_prox"} - 1 and $end_frag <= $oligos{$oligo}{"stop_prox"} + 1) {  #proximity exclusion if entire fragment is in proximity zone
                            $counter{"5aProximity excluded fragments total"}++;
                            $counter{"5bProximity excluded fragments $oligo"}++;
                            delete $data_hash{$read_name}{$mate}{$frag};
                            }
                        elsif ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"} and $start_frag <= $oligos{$oligo}{"stop_prox"} - 1 and $end_frag >= $oligos{$oligo}{"stop_prox"} + 1) {   #proximity exclusion if start of fragment is in proximity zone
                            $counter{"5aProximity excluded fragments total"}++;
                            $counter{"5bProximity excluded fragments $oligo"}++;
                            delete $data_hash{$read_name}{$mate}{$frag};
                            }
                        
                        #reporter fragments (remaining fragments)
                        else {
                            if (exists($data_hash{$read_name}{"label"}) and $data_hash{$read_name}{"label"} eq $oligo) {
                                unless (exists ($data_hash{$read_name}{$mate}{$frag}{"label"}) and $data_hash{$read_name}{$mate}{$frag}{"label"} =~ /capture/) {
                                    $data_hash{$read_name}{$mate}{$frag}{"label"} = "reporter $oligo";
                                    $data_hash{$read_name}{$mate}{$frag}{"mapped_coord"} = "$data_hash{$read_name}{$mate}{$frag}{\"chr\"}:$start_frag-$end_frag";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    $frag_counter{"$frags_counted"}++;
}

#print Dumper \%data_hash;

#generate %frag_hash containing filtered reporters needed for output
foreach my $read_name (keys %data_hash) {
    for (my $mate = 1; $mate < 3; $mate++) {                                                                                                        
        foreach my $frag (sort {$data_hash{$read_name}{$mate}{$a} <=> $data_hash{$read_name}{$mate}{$b}} keys %{$data_hash{$read_name}{$mate}}) {   
            foreach my $oligo (keys %oligos) {
                if (exists ($data_hash{$read_name}{$mate}{$frag}{"label"}) and $data_hash{$read_name}{$mate}{$frag}{"label"} eq "reporter $oligo") {
                    
                    unless (grep( /$data_hash{$read_name}{$mate}{$frag}{"mapped_coord"}/, @{$frag_hash{$oligo}{$read_name}} ) ) {                  #to make sure that a restriction fragment can only be reported once for each capture read
                        push (@{$frag_hash{$oligo}{$read_name}}, "$data_hash{$read_name}{$mate}{$frag}{\"mapped_coord\"}");
                        $counter{"6aReporter fragments total"}++;
                        $counter{"6bReporter fragments $oligo"}++;
                        if ($data_hash{$read_name}{$mate}{$frag}{"chr"} eq $oligos{$oligo}{"chr"}) {
                            $counter{"6cReporter fragments $oligo cis"}++;
                        }
                        else {
                            $counter{"6dReporter fragments $oligo trans"}++;
                        }
                    }
                }
            }
        }
    }
}
                
#combine globins
if ($globin > 0) {                                                                          #if globin is 1 or 2, combine Hba
    $oligos{"Hba_comb"}{"chr"} = 11;
    foreach my $oligo (keys %frag_hash) {
        unless ($oligo eq "Hba_comb") {
            foreach my $read (keys %{$frag_hash{$oligo}}) {
                if ($oligo =~ /Hba/) {
                    push (@{$frag_hash{"Hba_comb"}{$read}}, @{$frag_hash{$oligo}{$read}});
                }
            }
        }
    } 
}

if ($globin == 2) {                                                                          #if globin is 2, also combine Hbb
    $oligos{"Hbb_comb"}{"chr"} = 7;
    foreach my $oligo (keys %frag_hash) {
        unless ($oligo eq "Hbb_comb") {
            foreach my $read (keys %{$frag_hash{$oligo}}) {
                if ($oligo =~ /Hbb/) {
                    push (@{$frag_hash{"Hbb_comb"}{$read}}, @{$frag_hash{$oligo}{$read}});
                }
            }
        }
    } 
}

#print basic statistics to report file
print REP "Counters:\n";
foreach my $key (sort keys %counter) {
    printf REP "%-8s %s\n", $key, $counter{$key};
}

#print REP "Number of fragments per read:\n";
#foreach my $key (sort keys %frag_counter) {
#    printf REP "%-8s %s\n", $key, $frag_counter{$key};
#}

print REP "\n\nTracks:\n";

#generate oligo-specific fragment output file 
foreach my $oligo (keys %frag_hash) {
    my $frag_output = "$dir/$name\_$oligo\_fragments.txt";
    open (FRAG_OUT, ">$frag_output") or die "can't open fragment output file $frag_output";
    foreach my $read (keys %{$frag_hash{$oligo}}) {
        print FRAG_OUT "$read\t";
        my $length = scalar @{$frag_hash{$oligo}{$read}};
        for (my $i = 0; $i < $length; $i++ ) {
            my ($chr, $start, $stop) = split(/\W/, ${frag_hash{$oligo}{$read}}[$i]);            
            print FRAG_OUT "$chr:$start-$stop\_";
        }
        print FRAG_OUT "\n";   
    }
}

#generate gff hash 
my %gff_hash;
foreach my $oligo (keys %frag_hash) {
    foreach my $read (keys %{$frag_hash{$oligo}}) {
        my $length = @{$frag_hash{$oligo}{$read}};
        for(my $i = 0; $i < $length; $i++) {
            my $frag = ${frag_hash{$oligo}{$read}}[$i];
            ${$gff_hash{$oligo}{$frag}}++;            
        }
    }
}

#generate oligo-specific gff output file 
foreach my $oligo (keys %gff_hash) {
    my $gff_output = "$dir/$name\_$oligo\.gff";
    open (GFF_OUT, ">$gff_output") or die "can't open fragment output file $gff_output";
    foreach my $coord (sort keys %{$gff_hash{$oligo}}) {
        my ($chr, $dpn_start, $dpn_stop) = split(/\W/, $coord);
        $dpn_start = $dpn_start + 2;
        $dpn_stop = $dpn_stop + 2;
        unless ($chr =~ /M|Y/) {                                                                            
            print GFF_OUT "chr$chr\tco-cap\t$name\t$dpn_start\t$dpn_stop\t${$gff_hash{$oligo}{$coord}}\t+\t0\t.\n";
        }
    }
}

#generate tracks
frag_to_wig(\%frag_hash,"");

unless (open(TRACKHUBA, ">$public_folder/$name\_hub.txt")){print "Cannot open file $public_folder/$name\_tracks.txt\n";exit;}
unless (open(TRACKHUBB, ">$public_folder/$name\_genomes.txt")){print "Cannot open file $public_folder/$name\_tracks.txt\n";exit;}
unless (open(TRACKHUBC, ">$public_folder/$name\_tracks.txt")){print "Cannot open file $public_folder/$name\_tracks.txt\n";exit;}

print TRACKHUBA "hub $name
shortLabel $name
longLabel $name\_CaptureC
genomesFile http://$public_url/$name\_genomes.txt
email $email";

print TRACKHUBB "genome $genome
trackDb http://sara.molbiol.ox.ac.uk$public_folder/$name\_tracks.txt";

foreach my $oligo (keys %oligos) {
    print TRACKHUBC
    "track $name\_$oligo
    type bigWig
    longLabel CC_$name\_$oligo
    shortLabel $name\_$oligo
    bigDataUrl http://$public_url/$name\_$oligo.bw
    visibility hide
    priority 200
    color 0,0,0
    autoScale on
    alwaysZero on

"
}

print REP "\n\nThe track data hub can be found at:
http://$public_url/$name\_hub.txt
Paste this link in the UCSC genome browser.\n";

#subroutines

#binary search: plots digested readas to restriction fragments by performing a binary search of a sorted array that contains the restriction fragments
sub binary_search {
    my ($chr_array, $start, $stop, $counter_hash) = @_;
    my $mid_value = ($start + $stop)/2;
    my $first_pos = 0;
    my $last_pos = scalar @$chr_array - 1; 
    my $counter =0;
    if (($mid_value < $$chr_array[$first_pos]) or ($mid_value > $$chr_array[$last_pos])) {
        $$counter_hash{"Binary search error: search outside range of restriction enzyme coords"}++;
        return ("error", "error")}
    for (my $i=0; $i<99; $i++) {
        my $mid_search = int(($first_pos + $last_pos) / 2);
        if ($$chr_array[$mid_search] > $$chr_array[$mid_search+1]) {
            $$counter_hash{"Binary search error: restriction enzyme array coordinates not in ascending order"}++;
            return ("error", "error")}
        if (($$chr_array[$mid_search] <= $mid_value) and ($$chr_array[$mid_search+1] > $mid_value)) {    #maps the mid point of the read to a fragment
            if (($$chr_array[$mid_search] <= $start+2) and ($$chr_array[$mid_search+1] >= $stop-2)) {    #checks the whole read is on the fragment +/-2 to allow for the dpnII overlaps
                return ($$chr_array[$mid_search], $$chr_array[$mid_search+1]-1)
                }
            else {
                $$counter_hash{"Binary search error: fragment overlaps multiple restriction sites"}++;
                return ("error", "error");
                }
        }       
        elsif ($$chr_array[$mid_search] > $mid_value) {
            $last_pos = $mid_search-1;
            }    
        elsif ($$chr_array[$mid_search] < $mid_value) {
            $first_pos = $mid_search+1;
            }
        else {
            $$counter_hash{"Binary search error: end of loop reached"}++;
            }
    }
    $$counter_hash{"Binary search error: couldn't map read to fragments"}++;
    return ("error", "error")
}

#frag_to_wig: outputs restriction fragments in data hash (format %hash{$oligo}{$read}{$1:1000-2000}) to wig format  
sub frag_to_wig {                   
    my ($hashref, $file_name) = @_; 
    foreach my $oligo (keys %$hashref) {
        foreach my $read (keys %{$$hashref{$oligo}}) {
            my $length = scalar @{$$hashref{$oligo}{$read}};
            for (my $i = 0; $i < $length; $i++ ) {
                my ($chr, $start, $stop) = split(/\W/, $$hashref{$oligo}{$read}[$i]);
                #$chr =~ s/chr//gi;
                if ($chr eq $oligos{$oligo}{"chr"}) {                                   #only plot in cis to make script run faster
                    my @range = ($start..$stop);
                    for(my $j = 0; $j < $#range; $j++) {
                        ${$coord_hash{$oligo}{$chr}{$range[$j]}}++;
                    }
                }
            }
        }
        my $tracks_out = "$name\_$oligo$file_name";
        open (TRACKS_OUT, ">$dir/$tracks_out.wig");
        foreach my $chr (sort {$coord_hash{$oligo}{$a} <=> $coord_hash{$oligo}{$b}} keys %{$coord_hash{$oligo}}) {
            print TRACKS_OUT "variableStep  chrom=chr$chr\n";
            foreach my $coord (sort {$a <=> $b} keys %{$coord_hash{$oligo}{$chr}}) {
                my $count = ${$coord_hash{$oligo}{$chr}{$coord}};
                print TRACKS_OUT "$coord\t$count\n";
            }
        }
        system ("wigToBigWig -clip $dir/$tracks_out.wig $bigwig_folder/$genome\_sizes.txt $dir/$tracks_out.bw") == 0 or die "couldn't bigwig files\n";
        system ("mv $dir/$tracks_out.bw $public_folder") == 0 or die "couldn't move files\n";		
        system ("chmod 755 $public_folder/$tracks_out.bw") == 0 or die "couldn't chmod files\n";   
        print REP "track type=bigWig name=\"$tracks_out\" description=\"CaptureC_$oligo $tracks_out\" bigDataUrl=http://$public_url/$tracks_out.bw\n";
    }
}

    

