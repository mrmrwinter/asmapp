install.packages("circlize")

library(circlize)

asm = "../concat_mrl20k_q15l15k_pepper/outputs/redundans/scaffolds.reduced.fasta"

random_bed = generateRandomBed()
chordDiagram(random_bed)
