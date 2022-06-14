import os
import pybedtools



def clean_gff(gff_in, gff_out, f_tag):
	# Clean up GFF3 files so they work with pybedtools
	with open(gff_in, 'r') as gff_i:
		with open(gff_out, 'w') as gff_o:
			data = gff_i.readlines()
			clean_data = []
			for line in data:
				tab_cnt = len(line.split('\t'))
				feat_str = '\t' + f_tag + '\t'
				if feat_str in line: # only add feature lines
					clean_data.append(line)
			gff_o.write(''.join(clean_data))



def compare_gffs(orf_gff_file, feat_gff_file, tag='UNKNOWN'):
	""" This function is meant to compare predicted ORF GFFs
		to other predicted feature GFFs, i.e., rRNA or tRNA.
		Produces a GFF that contains overlaps and saves in 
		original ORF directory.
	"""

	# Get sample ID and ORF save path
	sv_path = os.path.dirname(feat_gff_file)
	sample_id = os.path.basename(orf_gff_file).rsplit('.', 2)[0]

	# Clean input GFF3 files
	orf_clean_file = os.path.join(sv_path, sample_id + '.annot.cds.gff')
	feat_clean_file = os.path.join(sv_path, sample_id + '.' + tag + '.clean.gff')
	clean_gff(orf_gff_file, orf_clean_file, 'CDS')
	clean_gff(feat_gff_file, feat_clean_file, tag)
	
	# Load GFF3s
	orf_gff = pybedtools.BedTool(orf_clean_file)
	feat_gff = pybedtools.BedTool(feat_clean_file)

	# Find overlaps between ORFs and Feature GFF
	orf_overlaps = orf_gff.intersect(feat_gff, u=True)
	# Save overlaps for each comparison
	orf_overlaps.saveas(os.path.join(sv_path, sample_id + '.cds.' + tag + '.overlaps.gff'))

	# Filter out overlapping ORFs
	#orf_rrna_filtered = orf_gff.subtract(rrna_gff, A=True)
	#orf_trna_filtered = orf_rrna_filtered.subtract(trna_gff, A=True)
	#orf_trna_filtered.saveas(os.path.join(sv_path, sample_id + '.cds.no_overlaps.gff'))

	# ORF overlap stats
	with open(os.path.join(sv_path, sample_id + '.' + tag + '.feature.log'), 'w') as l_out:
		l_out.write("Raw predicted ORFs: " + str(orf_gff.count()) + '\n')
		l_out.write("Predicted " + tag + ": " + str(feat_gff.count()) + '\n')
		l_out.write("ORFs that overlap with " + tag + ": " + str(orf_overlaps.count()) + '\n')

