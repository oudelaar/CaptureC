# CaptureC

# This script performs the analysis of a Capture-C experiment. It takes as input the aligned sam file, the oligo file (CCanalyser format), and a file with the coordinates of all the restriction fragments in the genome of the restriction enzyme used for digestion (format chr:start-stop\n).
# For more info on file formats and how to generate these, see http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/UserManualforCaptureCanalysis.pdf.

# The script analyses the sam file and filters out:
# 1. Duplicate reads. For flashed reads duplicate removal is based on a unique start and stop of all mapped fragments; for non-flashed reads only the start of the two separate reads is considered and needs to be unique. The script has an option for more stringent analysis, in which reads with only 1 (-stringent 1) or 2 (-stringent 2) differences are also removed.
# 2. Proximity exclusion. All mapped restriction fragments that fall within a 1kb window from any of the viewpoints used will be removed.
# The remaining read fragments are considered reporters and are mapped to the corresponding restriction fragments. 
# The script contains a 'globin option' to combine the separate globin tracks into a combined track.

# The script is dependent on the ucsctools module. It is also necessary to specify a folder containing a file with the size of the genome ($bigwig_folder).

# Output:
# 1. Report with basic statistics: total aligned reads, duplicates, proximity exclusion, reporters in cis and trans
# 2. Oligo-specific txt file with reporter fragments in the format "read_name \t frag1-coord_frag2-coord_..."  
# 3. Oligo-specific wig file
# 4. Oligo-specific gff file

# The script takes 5-10 min to run when using the default duplicate filtering. Stringent duplicate filtering is very slow and will take overnight for samples sequenced on a miseq (sam file up to 3 Gb). This is probably only required when low cell numbers are used. 

# Example of run command:
# module load ucsctools
# perl CC_MO.pl -sam your_input_sam_file.sam -o your_input_oligo_file.txt -r mm9_dpnII_coordinates.txt -pf /public/oudelaar/ -name short_sample_name &

# contact: oudelaar@well.ox.ac.uk
