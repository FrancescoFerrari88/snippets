#merge bam files in series
while read line; do read sampleA sampleB outbam rest <<< "$line"; echo samtools merge $outbam $sampleA $sampleB; done < <(cat mESC_pasted.txt)
