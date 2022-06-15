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



def compare_gffs(diag_path, orf_gff_file, feat_gff_file, tag='UNKNOWN'):
	"""
	This function is meant to compare predicted ORF GFFs
	to other predicted feature GFFs, i.e., rRNA or tRNA.
	Produces a GFF that contains overlaps and saves in 
	diagnostics directory.
	"""

	# Get sample ID and ORF/Feature save path
	orf_base_file = os.path.basename(orf_gff_file)
	feat_base_file = os.path.basename(feat_gff_file)
	sample_id = orf_base_file.rsplit('.', 2)[0]

	# Clean input GFF3 files
	orf_tmp_file = os.path.join(diag_path, orf_base_file)
	feat_tmp_file = os.path.join(diag_path, feat_base_file)
	clean_gff(orf_gff_file, orf_tmp_file, 'CDS')
	clean_gff(feat_gff_file, feat_tmp_file, tag)
	
	# Load GFF3s
	orf_gff = pybedtools.BedTool(orf_tmp_file)
	feat_gff = pybedtools.BedTool(feat_tmp_file)

	# Find overlaps between ORFs and Feature GFF
	orf_overlaps = orf_gff.intersect(feat_gff, u=True)
	# Save overlaps for each comparison
	orf_overlaps.saveas(os.path.join(diag_path, sample_id + '.ORF_' + tag + '_overlaps.gff'))

	# ORF overlap stats
	with open(os.path.join(diag_path, sample_id + '.' + tag + '_overlap.log'), 'w') as l_out:
		l_out.write("Raw predicted ORFs: " + str(orf_gff.count()) + '\n')
		l_out.write("Predicted " + tag + ": " + str(feat_gff.count()) + '\n')
		l_out.write("ORFs that overlap with " + tag + ": " + str(orf_overlaps.count()) + '\n')

	return orf_tmp_file, [tag, feat_tmp_file]


def subtract_gffs(diag_path, orf_gff_file, feat_file_dict):
	"""
	This function takes the ORF GFF in the diagnostics
	dir and a dictionary containing features, i.e. tRNAs
	and rRNAs, and subtracts the features from the ORFs.
	"""

	# Get sample ID and ORF save path
	orf_base_file = os.path.basename(orf_gff_file)
	sample_id = orf_base_file.rsplit('.', 2)[0]
	
	# Load ORF GFF3
	orf_gff = pybedtools.BedTool(orf_gff_file)

	for tag in feat_file_dict:
		feat_gff_file = feat_file_dict[tag]
		feat_gff = pybedtools.BedTool(feat_gff_file)
		# Filter out overlapping ORFs
		orf_gff = orf_gff.subtract(feat_gff, A=True)

	orf_gff.saveas(os.path.join(diag_path, sample_id + '.ORFs_minus_features.gff'))