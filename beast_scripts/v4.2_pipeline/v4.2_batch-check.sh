#!/bin/sh

# Name of Project
Project="cellmix"
# MiSeq reads
reads=200
# Project Directory

dir="/home/dyap/Projects/PrimerDesign"
p3dir=$dir"/"$Project"/primer3"
inpath="/home/dyap/dyap_temp/vcf"

# Input the names of the samples, they are processed in series
# To do - rewrite to make them run in parallel

# 2 Dec 2014 (after using old shared positions)
input=$inpath"/L2-H-shared-masked_check.vcf"
bname=`basename $input | awk -F"." '{print $1}'`
sample=$bname

out=$inpath"/"$bname

##################################
###   Change these variables    ##
echo -e "SA_chr_pos"'\t'"REF"'\t'"ALT"'\t'"Prob"'\t'"Geno"'\t'"Maskseq" > $out.tmp
vcf-query $out.vcf -f 'SA040:SA036\_chr%CHROM\_%POS\t%REF\t%ALT\t%INFO/PR\t%INFO/GT\t%INFO/FS\n' >> $out.tmp
##################################

# This are the passing conditions, please change as necessary
# GT (ie $5) must be 0/1 ie het

# This picks out those that PASS the conditions
awk -F"\t" '{ if ($5=="0/1") print $0}' $out.tmp | tail -n +1 > $out.txt

infile=$out".txt"
echo $infile

primfile="/home/dyap/Projects/PrimerDesign/cellmix/Shared_Primers"
outfile=$inpath"/"$sample"_primer3_output.txt.cat" 

        for i in `cat $infile | awk -F"\t" '{print $1}' `
                do

        ####################################################
        #   This parser must be changed for varying input  #
        ####################################################


                ID=`grep -m1 $i $infile | awk -F"\t" '{print $1}'`
                Ref=`grep -m1 $i $infile | awk -F"\t" '{print $2}'`
                Alt=`grep -m1 $i $infile | awk -F"\t" '{print $3}'`
                Pos=301 ## by definition

		LEFT=`grep -m1 $i $primfile | awk -F"\t" '{print $2}'`
		RIGHT=`grep -m1 $i $primfile | awk -F"\t" '{print $3}'`

                ####################################################################
                # clean up to make unrecognised characters and spaces N in the Seq #
                ####################################################################

#               Seq=`grep -m1 $i $infile | awk -F"\t" '{print $6}' | tr -d '[' | tr -d ']' | tr -c 'ATCGatcg' 'N' `
                Seq=`grep -m1 $i $infile | awk -F"\t" '{print $6}' | tr -d '['  | sed 's/\/.*\]//g'`


                ################################################################################
                # For matching of context, we get the 5' 10 bp of the SNV or nreakpoint or end #
                ################################################################################

                llim=`echo "$Pos - 5" | bc`
                rlim=`echo "$Pos + 5" | bc`
#               echo $llim, $rlim
                range=`echo $llim"-"$rlim`
#               echo $range

                if [[ "$Pos" -eq "target" ]];
                        then continue
                        else

                #################################################
                # This is the context matching for viewing only #
                #################################################

                cxtSeq=`echo $Seq | cut -c$range`
                fi

                seqLength=`echo $Seq | wc -c`

        if [[ "$seqLength" -lt "$Miseq" ]];
                # Not enough design space, next.
                then continue;
                else
        {
        echo "SEQUENCE_ID="$ID;
        echo "SEQUENCE_TEMPLATE="$Seq;

        ###################################################################
        # The Context, ALT, REF alleles are carried through as P3_COMMENT #
        ###################################################################

        echo "P3_COMMENT="$cxtSeq"@"REF_"$Ref":"VB_"$Alt":"


        #############################################
        # Primer3 Parameters specific to the expt   #
        # Note see ~/Projects/PrimerDesign/settings #
        #     for general settings for Primer3      #
        #############################################

        # <start>,<length>
        plim=`echo "$Pos - 10" | bc`
        echo "SEQUENCE_TARGET="$plim",20";

        # Specifies the optimal primer sizes
        echo "PRIMER_OPT_SIZE=22";
        echo "PRIMER_MIN_SIZE=18";
        echo "PRIMER_MAX_SIZE=26";

        # Specifies the number of unknown bases in any primer
        # For N Masking to prevent primer overlap
        # This must be uncommented 
        echo "PRIMER_NUM_NS_ACCEPTED=0";

        # Specifies the product range (or ranges)
        # Must be within range of $reads
        # dlim1=`echo "0.9 * $reads" | bc` # For automated calculation etc
        echo "PRIMER_PRODUCT_SIZE_RANGE=170-195 145-220";

        echo "PRIMER_GC_CLAMP=2";
        echo "PRIMER_NUM_RETURN=5";
        echo "PRIMER_EXPLAIN_FLAG=1";
	echo "PRIMER_LEFT_EXPLAIN=ok 1";
	echo "PRIMER_RIGHT_EXPLAIN=ok 1";
	echo "PRIMER_PAIR_EXPLAIN=ok 1";
	echo "PRIMER_LEFT_NUM_RETURNED=1";
	echo "PRIMER_RIGHT_NUM_RETURNED=1";
	echo "PRIMER_INTERNAL_NUM_RETURNED=0";
	echo "PRIMER_PAIR_NUM_RETURNED=1";
	echo "PRIMER_PAIR_0_PENALTY=N/A";
	echo "PRIMER_LEFT_0_PENALTY=N/A";
	echo "PRIMER_RIGHT_0_PENALTY=N/A";
	echo "PRIMER_LEFT_0_SEQUENCE="$LEFT;
	echo "PRIMER_RIGHT_0_SEQUENCE="$RIGHT;

        echo "=";
        } >> $outfile
        fi
                done


#-----------------------------DO NOT EDIT ANYTHING BELOW THIS LINE -------------------------------------

	 echo `date` \r >> ~/Projects/PrimerDesign/jobs.log ; echo $sample "started processing." >> ~/Projects/PrimerDesign/jobs.log

	if [ -d $dir"/"$Project ]; then
			mkdir $p3dir
		else
			mkdir $dir"/"$Project
			mkdir $p3dir
	fi

# copies the final version of the primer design process to the primer3 project directory
	cp $inpath"/"$sample"_primer3_output.txt.cat" $p3dir"/"$sample"_primer3_output.txt"

	~/Scripts/v4.2_pipeline/v4.2_display_pipeline.sh $Project $sample $reads
	echo `date` \r >> ~/Projects/PrimerDesign/jobs.log; echo $sample "completed using v4.2 ." >> ~/Projects/PrimerDesign/jobs.log

exit;

