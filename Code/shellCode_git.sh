# RNAseq process
trim_galore -j 20 -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired "$FQ1" "$FQ2" -o "$output_folder"

hisat2 -p 20 -x "$genome_index" -1 "${output_folder}/${Al_fq1}" -2 "${output_folder}/${Al_fq2}" -S "${output_folder}/${SampleName}.sam"

samtools view -@ 20 -bS "${output_folder}/${SampleName}.sam" > "${output_folder}/${SampleName}.bam"

samtools sort -@ 20 -o "${output_folder}/${SampleName}_sorted.bam" "${output_folder}/${SampleName}.bam"

samtools index -@ 20 "${output_folder}/${SampleName}_sorted.bam"

featureCounts -p -t exon -g gene_id -T 40 -a "$gtf_file" -o "${output_folder}/${SampleName}_counts.txt" "${output_folder}/${SampleName}_sorted.bam"	
