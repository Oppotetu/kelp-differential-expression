
diamond blastx --db ./phaeophyceae_db.dmnd -q ./light1_extracted.fasta \
-o matches.tsv --outfmt 6 qseqid sseqid bitscore evalue pident salltitles \
--ultra-sensitive 

