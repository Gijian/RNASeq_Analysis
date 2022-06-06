#!/bin/bash

if [ ! -d "./fastp" ] ## Create folder if it does not exist to store summary output from fastp
	then
		mkdir ./fastp
fi

if [ ! -d "./hisat_summary" ] ## Create folder if it does not exist to store summary output from hisat alignment
	then
		mkdir ./hisat_summary
fi

if [ ! -d "./bam" ] ## Create folder if it does not exist to store output BAM files of all files processed
	then
		mkdir ./bam
fi

## This loop processes all fq.gz files sharing directory with this script
## File format given by NGS output is [name]_1.fq.gz (read1) and [name]_2.fq.gz (read2)
## This will also only loop based on # of read1 file as each sample contain read1 and read2 file
## Rather than piping all the commands in 1 line, I split it up as it can take a while per sample and can periodically error out on my computer.  This would allow me to some what resume.

for i in `ls -1 *_1.fq.gz | sed 's/\_1.fq.gz//'` ## Looping based on read1 file but truncating the name from [name]_1.fq.gz to just [name]
do
	echo $i ## Verbal parameter to tell me what sample it is processing
	date +%r ## Print starting date/time
  
  ## Check to see if the BAM file for this sample already exist.  If so, it would skip. This would allow for script to resume of something happen during multiple file processing.
  
	if [ -s ./bam2/$i\.bam ]  
	then
		echo "BAM alignment exist"
	else
	  printf "Running Fastp $nam at $(date +%r)\n" ## Print out what it is doing
  
  ## I find fastp is a convenient tools for trimming and cleaning up poor quality bases.
  ## The output would create a (read1) temp1.fq.gz and (read2) temp2.fq.gz which would be used for the next step in alignment.  These files would be overwritten and deleted at the end.
	
    fastp -w $(nproc) --detect_adapter_for_pe -q 30 -i $i\_1.fq.gz -I $i\_2.fq.gz -o temp1.fq.gz -O temp2.fq.gz -h /fastp/$i\_fastp.html 
	  printf "Fastp Completed\n\n"
  
  ## Alignment with HISAT2 on mouse MM10 genome reference store on my computer which was downloaded at http://daehwankimlab.github.io/hisat2/download/
  ## The process will return summary that would be placed in the hisat_summary folder created.  The temp SAM file is large and will be overwritten after each sample and at the end of run.
  ## Added --tmo because I only have interest in transcript region for DGE analysis

    printf "Running HISAT2 $nam at $(date +%r)\n" ## Print out time of HISAT2 starting alignment 
    hisat2 -p $(nproc) -t -x /mnt/c/ref/Mouse/grcm38_tran/genome_tran -1 temp1.fq.gz -2 temp2.fq.gz -S temp.sam --tmo --summary-file ./hisat_summary/$i\_hisat.txt 
	  printf "SAM alignment completed\n\n"
	  printf "Converting to BAM at $(date +%r)\n"
    
  ## Samtools sort and conversion to BAM is under a while loop because for some reason it randomly error out at different time point and I cannot find the solution.
  ## -f 0x2 is used to keep reads that are paired.  Can be ommited.  
  
		while ! samtools view -@ $(nproc) -h -q 30 -f 0x2 temp.sam | samtools sort -@ $(nproc) -o ./bam/$i\.bam -
		do
			echo "Samtools error...as always...repeating..."
			sleep 10
		done
	  printf "Completed Alignment at $(date +%r)\n\n"  ## Final time output to contrast with the starting time to see how long it took to process
	fi
done

## Remove all the temp files once the processing is completed
rm temp1.fq.gz
rm temp2.fq.gz
rm temp.sam
