# Identifying_redundancy
Small script designed to take a single genome assembly and identify contigs which show redundancy with larger contigs through BLAST alignment <br/>
This is useful for regions that are repetitive and/or high copy number which can be assembled several times in different contigs <br/>
For example this often happens when the mitochondrial genome is assembled <br/>

ALTERNATIVE FOUND AT A LATER DATE AND EASIER::<br/>
Just install funannotate and then run the function funannotate clean <br/>
Provides an output and some verbatim telling you how many contigs were removed and the identity and coverage <br/> 
Below I use the threshold of 99% covered and 99% basepair identity <br/>
Much quicker and easier (wish I knew of this before) <br/>

        mamba create -n funannotate -c bioconda funannotate <br/>
        funannotate clean -i genome.fa -o genome.cleaned.fa -p 99 -c 99 

For this you need only two tools installed and in path, BLAST and samtools

For identifying redundancy we can use 4 parameters:
1. The lower limit of identity that should be considered relevant for the within genome comparison = identity
2. The amount of the contig that should be aligned i.e how much should be contained within another contig = overlap
3. The distance between alignments that should be closed to combined alignments = gap
4. The largest size of contigs that should be considered potentially redundant ( i.e we do not expect large contigs to be contained within even larger ones, plus these would be very obvious in dotplots) = contigsize

Therefore only alignments with an identity greater than *identity*% are kept <br/>
Single contiguous non-redudant alignments are calculated after filling gaps <= *gap* <br/>
Single alignments which cover > *overlap*% of the contig are identified <br/>
If this percentage cover is greater than that of the corresponding contig (i.e the contig is smaller than that which is had aligned to), then it is considered the redudndant contig <br/>


Below I will use a 99% identity and overlap wih a 10bp gap and only look for redudancy in contigs smaller than 500kb <br/>
Note if wanting to go smaller than 99%, for either overlap of identity, do not add decimal points <br/>
For example: An overlap of at greater than 99.9% would be overlap="999" as this is translated to .999 <br/>

    genome="genome.fa"
    overlap="99"
    identity="99"
    gap="10"
    contigsize="500000"
    
    #just generate a prefix for the file naming by removing any fasta suffix
    genome2=$( echo $genome | awk -F ".fa" '{print $1}'  )
    #generate a file output region which be be filled with intermediate files such as the BLAST alignments etc
    mkdir ${genome2}_redundancy
    #remove the summary output file incase the same output directory is being used again
    rm ${genome2}_redundancy/OFINTEREST.${identity}_${overlap}.tsv
    #index the genome in order to the get the contig sizes and names easily
    samtools faidx ${genome}
    #take the contigs smaller than the contigsize value and blast each one seperately against the whole genome
    cat ${genome}.fai | awk -v contigsize="$contigsize" '{if($2 < contigsize) print $0}' | while read line
    do
      contig=$( echo $line | awk '{print $1}' )
      size=$( echo $line | awk '{print $2}' )
      samtools faidx $genome $contig > ${genome2}_redundancy/temp.fa
      makeblastdb -in ${genome2}_redundancy/temp.fa -input_type fasta -parse_seqids -dbtype nucl
      blastn -db ${genome2}_redundancy/temp.fa  -query $genome -outfmt 6 | awk -v contig="$contig" '{if($1 != contig) print}' >  ${genome2}_redundancy/${contig}.mtBLAST_raw.tsv
      #take only alignments greater than the identity value, rearrange inverted alignments, sort them and then generate non-redundant single alignments that are joint if there is a distance between them <= the gap value
      cat ${genome2}_redundancy/${contig}.mtBLAST_raw.tsv | awk -v identity=".$identity" '{if($3 > identity) print $0 }' | awk '{if($9 < $10) {print $1"\t"$9"\t"$10} else {print $1"\t"$10"\t"$9}}' | sort -k1,1 -k2,2n | bedtools merge -d ${gap} > ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.bed
      #remove file incase the same output directory is being used again
      rm ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.proportion.bed
      #get the size of the corresponding contig in order to calulate which contig is covered more by the alignment
      cat ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.bed | cut -f1 | sort -u | while read contig2
      do
        size2=$( cat ${genome}.fai | awk -v contig2="$contig2" '{if($1 == contig2) print $2}' )
        cat ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.bed | awk -v contig2="$contig2" -v size2="$size2" -v contig="$contig" -v size="$size" '{if($1 == contig2) print $2"\t"$3"\t"($3-$2)+1"\t"contig"\t"size"\t"(($3-$2)+1)/size"\t"contig2"\t"size2"\t"(($3-$2)+1)/size2}' >> ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.proportion.bed
      done
      #take all alignments of interest, in that they are greater than the overlap value and that the contig in the 4th column is covered more by the alignment than the other contig
      cat ${genome2}_redundancy/${contig}.mtBLAST_raw.merged2k.proportion.bed  | awk -v contig="$contig" -v size="$size" -v overlap=".$overlap" '{if($6 > overlap && $6 > $9)print $0}'  >> ${genome2}_redundancy/OFINTEREST.${identity}_${overlap}.tsv
    done
    #get a list of the contigs to remove
    list=$( cat ${genome2}_redundancy/OFINTEREST.${identity}_${overlap}.tsv | cut -f4 | sort -u )
    #take the opposite for those to keep
    listkeep=$(  cut -f1 ${genome}.fai | grep -v "${list}" )
    #extract only those contigs to keep
    samtools faidx ${genome} ${listkeep} > ${genome2}.nonredundant.fa

The output is the \*nonredundant.fa file which has the redudant contigs removed <br/>
And also a directorry containing all the preliminary files <br/>
One preliminary file 'OFINTEREST.\*.tsv' is the alignments considered of interest to the removal process, i.e all those that matched the final criteria <br/>
The rest are files per contig considered or redundancy such as the raw blast match filtered for only matches above the identity parameter
