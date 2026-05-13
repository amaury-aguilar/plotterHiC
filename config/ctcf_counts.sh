

cat <(awk '{OFS="\t"; print $1,$2-500,$3+500,$4,$5,$6,$2,$3,"192,91,56",1,$3-$2,0,"plus"}' Ctcf.ATACseq.Ctcf_WT_KO.motifs_bed6.plus.bed) <(awk '{OFS="\t"; print $1,$2-500,$3+500,$4,$5,$6,$2,$3,"181,250,92",1,$3-$2,0,"minus"}' Ctcf.ATACseq.Ctcf_WT_KO.motifs_bed6.minus.bed) > Ctcf.ATACseq.Ctcf_WT_KO.Ctcf_motifs.plus_minus.1kb_extended.bed
 
