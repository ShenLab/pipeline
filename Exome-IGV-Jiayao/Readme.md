Step 1: Prepare:
	Input:
		1. variant csv file with header 'Chrom', 'Pos', 'SampleID'
		2. Pedigree file with 'SampleID', 'Father', 'Mother'
		3. BamList contains bam locations.
	Function:
		link variant with sample and parents, link samples with their bam name.
		Generate script to run bamout bams.


Step 2: Bamout:
	Input:
		script produced by previous step. More, contain a Target Inverval and bam.
	Function:
		Generate bamout of a Region-Sample pair.
		Link bamout with variant

Step 3: IGV:
	Input:
		Variant list with format: chrom:pos    bam,bamout,fa_bam,fa_bamout,mo_bam,mo_bamout
		Bamlist with bam locations, should be original bams concatnate with bamout bams.
		