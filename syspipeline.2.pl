#!/usr/bin/perl

#Author:David Selassie Opoku
#Date started: July 15, 2011
#syspipeline.pl -- This program receives a list of fastq files and process them
#                  using the program BWA and SAMTOOLS to process a list of sorted
#                  bam files for use in later processes.
#Updates: July 18, 2011
#         1)Implemented code to use Perl module Getopt::Std to accept options and values from command line for system call.
#         2)Included print statements with "Command:" delimiter before each system call.

use strict;
use warnings;
use Getopt::Std;#Perl module to get program option values from command line.
 
#Declare hash for bwa aln options.
my %bwaoptions=();
getopts("u:t:n:o:e:", \%bwaoptions);#populations "%bwaoptions hash" with bwa aln options and values.

#This subroutine is called by the main subroutine and checks the command line arguments for valid inputs.This is done as follows:
#                1)Checksthat all "bwa aln" command options and values have been entered
#                2)Checks that FASTQ FILE LIST has been entered and IS NOT EMPTY.
#                3)Checks that REFERENCE GENOME File has been entered and IS NOT EMPTY.
#                4)Checks that TOO MANY ARGUMENTS are not entered in command line.

sub VerifyInputs{
    print "Verifying Inputs...\n";
    #### -u T (user set parameters)
    if ($bwaoptions{u} eq "T"){
    	if(scalar(keys %bwaoptions) !=5){
        	die "exec -u T -t 10 -n 4 -o 1 -e 2 fastq_list_file ref_fa \n";
	}
	else{
		die "Must specify a value for -t\n" unless defined $bwaoptions{t};
                die "Must specify a value for -n\n" unless defined $bwaoptions{n};
                die "Must specify a value for -o\n" unless defined $bwaoptions{o};
                die "Must specify a value for -e\n" unless defined $bwaoptions{e};
	}	
     }
     ### -u F (program set parameters
     elsif ($bwaoptions{u} eq "F"){
	if(scalar(keys %bwaoptions) !=2){
		 die "exec -u F -t 10 fastq_list_file ref_fa \n";
	}
	else{
		die "Must specify a value for -t\n" unless defined $bwaoptions{t};
	}
     }
     ### no -u provided Or -u X
     else{
		print "exec -u T -t 10 -n 4 -o 1 -e 2 fastq_list_file ref_fa \n";
		print "OR\n";
		die "exec -u F -t 10 fastq_list_file ref_fa \n";
	}	    
    
    if($#ARGV <1){
	die "$0 requires a FILE WITH LIST OF FASTQ FILES and a REFERENCE FASTA FILE\n";
    }
    elsif($#ARGV == 1){
	die "The fastq list file  is EMPTY\n" if -z $ARGV[0];
	die "The reference fasta is EMPTY\n" if -z $ARGV[1];
    }
    else {   
	 die "$0 cannot take MORE THAN 2 arguments!\n";
    }
    print "$0 Inputs Verified...\n";
}


sub GetOptionValues{
	my $readLength =shift;
	if( $readLength < 50){
		$bwaoptions{n} = 1;
		$bwaoptions{o} = 0;
		$bwaoptions{e} = 0;
	}
        elsif( $readLength < 80){
                $bwaoptions{n} = 2;
                $bwaoptions{o} = 0;
                $bwaoptions{e} = 0;
        }
	elsif( $readLength < 100){
		$bwaoptions{n} = 3;
                $bwaoptions{o} = 0;
                $bwaoptions{e} = 0;
	}
	elsif( $readLength < 125){
		$bwaoptions{n} = 4;
                $bwaoptions{o} = 1;
                $bwaoptions{e} = 2;
	}
        else {
                $bwaoptions{n} = 6;
                $bwaoptions{o} = 2;
                $bwaoptions{e} = 2;
        }
}


#The subroutine, SystemCommands opens individual fastq files and carries out the following system commands:
#                1)bwa aln
#                2)bwa samse
#                3)C++ pipeline to select unique,non-redundant reads based on a start position key.
#                4)C++ pipeline to select unique,non-redundant reads based on an end position key.
#                5)samtools view
#                6)samtools sort
#                7)rm command to remove all transitionary files leaving the SAI and SORTED BAM files.
sub SystemCommands{
    open(FASTQLIST, $ARGV[0]) || die "Can't open $ARGV[0]\n";
    my $refFASTA = $ARGV[1];#declaring scalar for REFERENCE GENOME FILE(fasta format).
    my $count=1;
    while(<FASTQLIST>){
	my $FASTQ =$_;
	chomp($FASTQ);
        my $readLength;
        $readLength= `head -2 $FASTQ | tail -1 | wc -c`;
        $readLength--;
	if($bwaoptions{u} ne "T"){
 		&GetOptionValues($readLength);
	}
	# The following scalars are declared for the t,n,o & e options of the "bwa aln" command
	my $T=$bwaoptions{t};#scalar for -t option value
        my $N=$bwaoptions{n};#scalar for -n option value
        my $O=$bwaoptions{o};#scalar for -o option value
        my $E=$bwaoptions{e};#scalar for -e option value
	

	my $prefix="n".$N."."."o".$O."."."e".$E.".";#create prefix for output file.
	my $OUTPUT;  
        if ($FASTQ =~m|^(.*[/\\])(.*)\.\w+|)#grab file name from pathway using regular expression. $2 scalar stores file name.
        {
         $OUTPUT = $prefix.$2;
        }
        elsif($FASTQ =~ m|(.*)\.\w+|) #no pathway, has suffix
        {
         $OUTPUT = $prefix.$1;
        }
        elsif($FASTQ =~ m|^(.*[/\\])(.*)|)  #pathway, no suffix
        {
         $OUTPUT = $prefix.$2;
        }
        else
        {
         $OUTPUT = $prefix.$FASTQ;
        }

	print "Processing FASTQ FILE $count:$FASTQ\n";


	my $SAItemp=$OUTPUT.".sai";#create output name for SAI file from "bwa aln" command.
	print "Command:bwa aln -t $T -n $N -o $O -e $E -f $SAItemp $refFASTA $FASTQ\n";
	system("bwa aln -t $T -n $N -o $O -e $E -f $SAItemp $refFASTA $FASTQ"); #bwa alignment system call
	
	
	my $SAMtemp = $OUTPUT.".sam";#create output name for temporary SAM file from "bwa samse" command.
	print "Command:bwa samse -f $SAMtemp $refFASTA $SAItemp $FASTQ\n";
	system("bwa samse -f $SAMtemp $refFASTA $SAItemp $FASTQ");#bwa samse system call


	my $UNQNRsam = $OUTPUT.".uniq.nonRed.sam";#create ouput name for unique,non-redundant SAM file from two C++ pipelines.
	print "Command:./rmRedunSam2 $SAMtemp tempCpp.sam\n";
	system("./rmRedunSam2 $SAMtemp tempCpp.sam");#system call for C++ pipeline to select unique, non-redundant reads based on start position key.
	print "Command:./rmRedunSam3 tempCpp.sam $UNQNRsam $N\n";
	system("./rmRedunSam3 tempCpp.sam $UNQNRsam $N");#system call for C++ pipeline to select unique, non-redundant reads based on end position key.


	my $UNQNRbam = $OUTPUT.".uniq.nonRed.bam";#create output name for unique,non-redudant BAM file from "samtools view" command.
	print "Command:samtools view -u -S -o $UNQNRbam $UNQNRsam\n";
	system("samtools view -u -S -o $UNQNRbam $UNQNRsam");#samtools view system call.


	my $UNQNRsorted = $OUTPUT.".uniq.nonRed.sorted";#create output name for sorted BAM file from "samtools sort" command.
	print "Command:samtools sort $UNQNRbam $UNQNRsorted\n";
	system("samtools sort $UNQNRbam $UNQNRsorted");#samtools sort system call.
        
        print "$FASTQ\t$readLength(bp)\t-n $N -o $O -e $E\t";
 
	my $lineNum = `cat $FASTQ | wc -l`;
	chomp($lineNum);
	print int($lineNum)/4, "\t"; 
	
        my $comm = `awk '\$2 == 0 || \$2 == 16' $SAMtemp | wc -l`;
	chomp($comm);
        print $comm,"\t";
        print $comm/$lineNum * 400, "%\t";

        $comm= `grep  "	XT:A:U	" $SAMtemp | wc -l`;
        chomp($comm);
        print $comm,"\t";
        print $comm/$lineNum * 400, "%\t";
        
        $comm= `grep  "	XT:A:U	" $UNQNRsam | wc -l`;
        chomp($comm);
        print $comm,"\t";      
        print $comm/$lineNum * 400, "%\n";
	print "Command:rm $SAMtemp tempCpp.sam $UNQNRsam $UNQNRbam\n";
	system("rm $SAMtemp tempCpp.sam $UNQNRsam $UNQNRbam");#removing intermediary files to leave 
	$count++;
	print "\n\n";
    }
}
 
sub main{
    &VerifyInputs;
    &SystemCommands;
}
    
&main;
