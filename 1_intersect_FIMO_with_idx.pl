use strict;
use warnings;
use Getopt::Std;
use File::Basename;


#reading all the commans line options
my %opt;
getopt('irud',\%opt);

&help_message if (defined $opt{h});

####### Variable declaration ######
my $up_dist = 20; my $down_dist = 20; 

# reading the directory that has the index files in  it
my $idxdir = $opt{'i'};
$idxdir = check_dir($idxdir);
# reading the directory that has the reference files in  it
my $refdir = $opt{'r'};
$refdir = check_dir($refdir);

if (defined $opt{u}){
    $up_dist = $opt{'u'}
}

if (defined $opt{d}){
    $down_dist = $opt{'d'}
}



opendir IDX, $idxdir || die "Cannot open the directory";

while((my $filename = readdir(IDX))){
    # checking for .idx / .tab file in the folder
    my $gff;
    next if(check_extension($filename) < 1 || check_extension($filename) > 3);
    if((check_extension($filename) == 2) || (check_extension($filename) == 3)){
        $gff = 0
    }
    elsif(check_extension($filename) == 1){
        $gff = 1
    }
    my $base_name = return_basename($filename);
    
    #creating  a tmp folder to store tmp files
    unless(-d $idxdir."tmp/"){
    system("mkdir ".$idxdir."tmp/");
    }
    
    # creating an output folder for each index file present in the folder.
    unless(-d $idxdir."output_intersect/"){
    system("mkdir ".$idxdir."output_intersect/");
    }
    print "\nINFO: Index file detected as input! Abort if incorrect\n";
    print "\nINFO: Input tag/index file = ".$filename;
    open IN, $idxdir.$filename || die "File not found";
    # open up the index/tag file and divide it into fwd and rev tag GFF files
    #print "\nINFO: Splitting the tags file into fwd and rev gff files\n";
    open OUT1,">".$idxdir."tmp/tags.gff" || die "File not found";
    #open OUT2,">".$idxdir."tmp/rev.tmp" || die "File not found";
    if($gff != 1){
        while(<IN>){
        chomp($_);
        next if(($_ =~ /^#/) || ($_ =~ /^chrom/));
        my @cols = split(/\t/,$_);
        if($cols[2] > 0){
            print OUT1 $cols[0]."\t.\t.\t".$cols[1]."\t".($cols[1]+1)."\t".$cols[2]."\t+\t.\t.\n";
        }
        if($cols[3] > 0){
            
            print OUT1 $cols[0]."\t.\t.\t".$cols[1]."\t".($cols[1]+1)."\t".$cols[3]."\t-\t.\t.\n";
        }
     }
    }
    
    close(IN);close(OUT1);
    # read the index file. Moving on to the refernece directory.
    opendir REF, $refdir || die "Cannot open the directory";
    while( (my $refname = readdir(REF))){ # Splitting each reference file into sense and anti-sense files.
        next if(check_extension($refname) != 1);
        my $ref_basename = return_basename($refname);
        print "INFO: Reference file = ".$refname."\nINFO: Splitting it into sense and antisesne files\n";
        open IN, $refdir.$refname || die "File not found";
        open OUT2,">".$idxdir."tmp/ref.gff" || die "File not found";
        #open OUT2,">".$idxdir."tmp/antisense.tmp" || die "File not found"; 
        while(<IN>){
            chomp($_);
            next if(($_ =~ /^#/) || ($_ =~ /^chrom/));
            my @cols = split(/\t/,$_);
            # First create sesnse and antisense ref feature file and add/subtract the up/down-stream distance from the start co-ordinate.
            if($cols[6] eq "+"){
                if($cols[3] - $up_dist <= 0){
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t1\t".($cols[3]+$down_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3]."-".$cols[4]."\t".$cols[8]."\n";
                }
                else{
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t".($cols[3]-$up_dist)."\t".($cols[3]+$down_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3]."-".$cols[4]."\t".$cols[8]."\n";
                }
            
            }
            else{
                if($cols[4] - $down_dist <= 0){
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t1\t".($cols[4]+$up_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3]."-".$cols[4]."\t".$cols[8]."\n";
                }
                else{
                    print OUT2 $cols[0]."\t".$cols[1]."\t".$cols[2]."\t".($cols[4]-$down_dist)."\t".($cols[4]+$up_dist)."\t".$cols[5]."\t".$cols[6]."\t".$cols[3]."-".$cols[4]."\t".$cols[8]."\n";
                }
            
            }
         # completed splitting the reference into sense and anti-sense gff files   
        }
        close(IN);close(OUT2);
    # Running intersectBed on the fwd/sense, rev/sense, fwd/antisense and rev/antisense and then joining fwd/sense with rev/antisense and rev/sense with fwd/antisense.
    # Feature with no tags would have "0" as overlap and features with tags would have >0 as overlap.
    #print "/usr/local/bin/closestBed -a ".$idxdir."fwd.gff -b ".$refdir."sense.gff  >".$idxdir."output/fwd_sense.txt\n";
    print "INFO: Running intersectBed on the generated files\n";
    system("/usr/local/bin/intersectBed -wao  -a ".$idxdir."tmp/ref.gff -b ".$idxdir."tmp/tags.gff  >".$idxdir."output_intersect/intersect_output.txt");
    
  } # end of the loop over reference directory.
    print "INFO: Your output is present in ".$idxdir."output_intersect/\n";
    system("rm -fr ".$idxdir."tmp/");
      
}
 


###################### Subroutine declarations ###################

sub print_hash{
    
    my($hash1,$OUT2) = @_;
    my @order;
    foreach my $value (sort {$$hash1{$b} <=> $$hash1{$a} }keys %$hash1)
    {
     print $OUT2 $value."\n";
     my @tmp = split(/\t/,$value);
     push(@order,shift(@tmp));
     #print $OUT2 $$hash1{$value}."\t".$value."\n";
    }
    return(@order);
}

sub print_hash_by_order{
    my($hash,$order,$OUT2) = @_;
    foreach my $id (@$order){
        print $OUT2 $id."\t".$$hash{$id}."\n";
        
    }
    
}
sub print_array{
    
	my ($array_to_print,$OUT2) = @_;
	foreach my $val (@$array_to_print)
	{
		print $OUT2 "\t".$val;
	}
	print  $OUT2 "\n";
}


sub print_header{
    my $OUT2 = shift;
    
    print $OUT2 "Uniqe ID\tGWEIGHT\t";
    for(my $i = (-1)*$up_dist;$i<=$down_dist;$i++){
        print $OUT2 $i."\t";
    }
    print $OUT2 "\nEWEIGHT\t\t";
    
    for(my $i = (-1)*$up_dist;$i<=$down_dist;$i++){
        print $OUT2 "1\t";
    }
    print $OUT2 "\n";
    
}


sub check_dir{
    my $dir = shift;
    my $tmp = substr($dir,-1,1);
    if($tmp eq "/"){
        return($dir);
    }
    else{
        return($dir."/");
    }
}


sub return_basename{
    my $filename = shift;
    my @suffixes = (".gff",".idx",".tab",".txt",".tmp",".tab");
    my $basename = basename($filename, @suffixes);
    return $basename;
    
}
sub check_extension{
    my $filename = shift;
    my @suffixes = (".gff",".idx",".tab",".txt");
    my ($fname,$path,$suffix) = fileparse($filename,@suffixes);
    if($suffix eq ".gff"){
        return 1;
    }
    elsif($suffix eq ".idx"){
        return 2;
    }
    elsif($suffix eq ".tab"){
        return 3;
    }
    elsif($suffix eq ".txt"){
        return 4;
    }
    else{
        return 0;
    }
}

sub help_message {
  print qq{
    Program: intersect_FIMO_with_idx.pl (Calculate distance from reference feature and generate CDT files)
    Contact: Rohit Reja <rzr142\@psu.edu>
    Usage:   intersect_FIMO_with_idx.pl -i <path_to_index_file> -r <path_to_REF_feature_files>

    Options: -i <path1>     path to the folder with index files in .idx or .tab format 
             -r <path2>     path to the folder with reference feature files in gff format.
             -u <number>    Upstream distance to go, defualt = 20bp.
             -d <number>    Downstream distance to go, default = 20bp.
             -h             help

    Example:
      perl intersect_FIMO_with_idx.pl -i  /usr/local/folder -r /folder1/ref/
      
  
  };
  exit;
}







