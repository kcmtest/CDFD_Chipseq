set -e
set -u
set -o pipefail -o errexit -o nounset


#for peak in *bed
for peak in *narrowPeak
do
	file_name=$(basename $peak .narrowPeak)
	bedtools intersect -a "$peak" -b /run/media/punit/data2/PRIMARY_CELL_ATAC_SEQ/test/DISEASE_DIFFBIND/hg38-blacklist.v2.bed -v \
	> ${file_name}_flt.bed
done
