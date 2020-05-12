#!/usr/bin/perl -w

###############################################
#
# This CGI perl script parses primer design
# job requests v2.1 
# Damian Yap, PhD
# 21 Sep 2013 & 1 Oct 2013
#
##############################################

use DateTime;

 $|=1;            # Flush immediately.

print "Content-Type: text/plain\n\n"; 

print "\n=============================================\n"; 

# get total arg passed to this script
my $total = $#ARGV + 1;
 
# get script name
my $scriptname = $0;
 
print "Total args passed to $scriptname : $total\n";

if ($total lt "5" ) {
    print "\n Usage:";
    print "\n perl $scriptname Project={Project} Sample={Sample} Ref={hg18_or_hg19} Miseq={Miseq_run_in_bp} p3settings=v2\n\n\n";
    exit;
    }

#
# Now the variable $total has your input data.
# Create your associative array.
#
    foreach $pair (@ARGV) {
        if ($pair =~ /(.*)=(.*)/) {  # found key=value;
        ($key,$value) = ($1,$2);     # get key, value.
        $value =~ s/\+/_/g;  # substitute spaces for + signs.
        $value =~ s/%(..)/pack('c',hex($1))/eg;
        $value =~ s/\|/_/g;  # substitute | for _ signs.
        $value =~ s/\ /_/g;  # substitute [space] for _ signs.
        $inputs{$key} = $value;   # Create Associative Array.
        }
    }
 
print "Project: $inputs{'Project'}\n";
print "Sample: $inputs{'Sample'}\n";
print "Genome Ref: $inputs{'Ref'}\n";
print "MiSeq run length: $inputs{'Miseq'}\n";
print "Primer3 Settings: $inputs{'p3settings'}\n\n\n";
 

#############################################################################
# This section appends all the requests to out log file for tracking purposes
 
print "Writing to /home/dyap/PrimerDesign/v1.5_requests\n"; 

open (MYARchIVE,">> /home/dyap/v1.5_requests") || die "Cannot open Archive";
print MYARchIVE "|";
print MYARchIVE "$inputs{'Project'}|";
print MYARchIVE "$inputs{'Sample'}|";
print MYARchIVE "$inputs{'Ref'}|";
print MYARchIVE "$inputs{'Miseq'}|";
print MYARchIVE "$inputs{'p3settings'}|";
print MYARchIVE "\n";
close MYARchIVE;

######################################################################
# This section writes the list to the file to be used for primer design

$file = $inputs{'Sample'};
$dir = $inputs{'Project'};
$sample = $inputs{'Sample'};
$Miseq = $inputs{'Miseq'};

$p3set = $inputs{'p3settings'};
if ($p3set eq "v0") {
    print "\n Using v0 settings"; 
    }
if ($p3set eq "v0-1") {
    print "\n Using v0-1 settings"; 
    }
if ($p3set eq "v0-2") {
    print "\n Using v0-2 settings"; 
    }
if ($p3set eq "v2") {
    print "\n Using v2 settings - default + iterative"; 
    }


print "\n-----------------------------------------------------------\n"; 
print "View the output online in another browser window when done:\n";
print "http://godel.cluster.bccrc.ca/workflow/primer3check/\n"; 
print "\n-----------------------------------------------------------\n"; 

print "\nWriting list of positions to ${'file'} in ${'dir'}/${'sample'}\n"; 

system("mkdir  /home/dyap/Projects/PrimerDesign/$dir");
system("mkdir  /home/dyap/Projects/PrimerDesign/$dir/positions");

$pos = $inputs{'Ref'};
print "\nThis is the value of hg = ${'pos'}\n";
print "\nThis is the value of sample = ${'sample'}\n";
print "\nThis is the value of dir = ${'dir'}\n";
print "\nThis is the value of file = ${'file'}\n";

if ($pos eq "NULL") {
    	print "\nFor privacy issues, please do not submit any patient, private or confidential data using this service. ";
	print "\nExiting now....";
    	exit;
    	}

if ($pos eq "hg19") {
	print "File to use /home/dyap/Projects/PrimerDesign/$dir/positions/$file-hg19";
	print "Please save your file in this location and press ENTER";
        $ans = <>;
	system ("export sample=${'sample'}; export Project=${'dir'}; export name=${'file'};/home/dyap/Scripts/v1.5_pipeline/v1.5_format_webinput.sh");

	}

if ($pos eq "hg18") {
    	print "\nConverting hg18->hg19, please double check the data when done. ";

	print "File to use /home/dyap/Projects/PrimerDesign/$dir/positions/$file-hg18";
	print "Please save your file in this location and press ENTER";
	$ans = <>;
	system ("export sample=${'sample'}; export Project=${'dir'}; export name=${'file'};/home/dyap/Scripts/v1.5_pipeline/v1.5_liftover.sh");
    }


print "\n=============================================\n"; 

print "Calling Rscript...\n"; 

print("Rscript /home/dyap/Scripts/v1.5_pipeline/v1.5_GetSeq.R --no-save --no-restore --args ${'dir'}/${'sample'}/${'file'}\n");
system("Rscript /home/dyap/Scripts/v1.5_pipeline/v1.5_GetSeq.R --no-save --no-restore --args ${'dir'}/${'sample'}/${'file'}");

print "\n=============================================================\n"; 

print "\nBack to perl script here\n"; 

print "=============================================================\n"; 
print "\n\nDesigning primers...\n (This takes a while please wait):\n"; 

system ("export sample=${'sample'}; export p3set=${'p3set'}; export hg=${'pos'};export Project=${'dir'}; export name=${'file'}; export Miseq=${'Miseq'};/home/dyap/Scripts/v1.5_pipeline/v1.5_primer3pipeline.sh");

print "=============================================================\n"; 

print "Primer Design job queued. Please check this address for results when done.\n"; 
print "http://godel.cluster.bccrc.ca/workflow/primer3check/\n"; 

print "=============================================================\n"; 

exit;


