### SCRIPT — STARLING RADseq DATA ###

#################################### PRE-PROCESSING ####################################

# download data from server
wget -q -c -O 9131_3270_67666_CB3C6ANXX_NHI1_ATCACG_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1432436124&refid=363594"
wget -q -c -O 9131_3270_67667_CB3C6ANXX_NHI2_CGATGT_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=564768908&refid=363595"
wget -q -c -O 9131_3270_67668_CB3C6ANXX_NHI3_TTAGGC_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1844585339&refid=363596"
wget -q -c -O 9131_3270_67669_CB3C6ANXX_NHI4_TGACCA_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1426440080&refid=363597"
wget -q -c -O 9131_3270_67670_CB3C6ANXX_NHI5_ACAGTG_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=352663514&refid=363598"
wget -q -c -O 9131_3270_67671_CB3C6ANXX_NHI6_GCCAAT_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=770739288&refid=363599"
wget -q -c -O 9131_3270_67672_CB3C6ANXX_NHI7_CAGATC_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1601033762&refid=363600"
wget -q -c -O 9131_3270_67673_CB3C6ANXX_NHI8_ACTTGA_R1.fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=527257196&refid=363601"

# unzip files
gunzip *.fastq.gz

# quality check using FastQC

fastqc 9131_3270_67666_CB3C6ANXX_NHI1_ATCACG_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67667_CB3C6ANXX_NHI2_CGATGT_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67668_CB3C6ANXX_NHI3_TTAGGC_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67669_CB3C6ANXX_NHI4_TGACCA_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67670_CB3C6ANXX_NHI5_ACAGTG_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67671_CB3C6ANXX_NHI6_GCCAAT_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67672_CB3C6ANXX_NHI7_CAGATC_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC
fastqc 9131_3270_67673_CB3C6ANXX_NHI8_ACTTGA_R1.fastq --outdir=home/nrh44/EUST_RADseq_2018/FastQC

# All files were of acceptable quality and passed filters based on expectations for RAD-Seq Data
# quality score across all bases ~ 35 --- all ok!

# Initial Filtering and processing using FASTX-Toolkit

# Trim 3' end of all reads to a length of 97 bp (FASTX Trimmer):

fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67666_CB3C6ANXX_NHI1_ATCACG_R1.fastq -o NH_I1.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67667_CB3C6ANXX_NHI2_CGATGT_R1.fastq -o NH_I2.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67668_CB3C6ANXX_NHI3_TTAGGC_R1.fastq -o NH_I3.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67669_CB3C6ANXX_NHI4_TGACCA_R1.fastq -o NH_I4.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67670_CB3C6ANXX_NHI5_ACAGTG_R1.fastq -o NH_I5.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67671_CB3C6ANXX_NHI6_GCCAAT_R1.fastq -o NH_I6.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67672_CB3C6ANXX_NHI7_CAGATC_R1.fastq -o NH_I7.fastq
fastx_trimmer -f 1 -l 97 -Q33 -i 9131_3270_67673_CB3C6ANXX_NHI8_ACTTGA_R1.fastq -o NH_I8.fastq


# Two-step filtering process:

# 1. Eliminate 100% of sequences with Phred quality scores below 10 

fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I1.fastq -o NH_I1filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I2.fastq -o NH_I2filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I3.fastq -o NH_I3filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I4.fastq -o NH_I4filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I5.fastq -o NH_I5filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I6.fastq -o NH_I6filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I7.fastq -o NH_I7filter_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i NH_I8.fastq -o NH_I8filter_1.fastq

# 2. then 5% of reads with Phred quality scores below 20

fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I1filter_1.fastq -o NH_I1filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I2filter_1.fastq -o NH_I2filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I3filter_1.fastq -o NH_I3filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I4filter_1.fastq -o NH_I4filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I5filter_1.fastq -o NH_I5filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I6filter_1.fastq -o NH_I6filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I7filter_1.fastq -o NH_I7filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i NH_I8filter_1.fastq -o NH_I8filter_2.fastq

# Demultiplex sequences using process_radtags (STACKS):

mkdir L1_raw
mkdir L2_raw
mkdir L3_raw
mkdir L4_raw
mkdir L5_raw
mkdir L6_raw
mkdir L7_raw
mkdir L8_raw

mv NH_I1filter_2.fastq L1_raw
mv NH_I2filter_2.fastq L2_raw
mv NH_I3filter_2.fastq L3_raw
mv NH_I4filter_2.fastq L4_raw
mv NH_I5filter_2.fastq L5_raw
mv NH_I6filter_2.fastq L6_raw
mv NH_I7filter_2.fastq L7_raw
mv NH_I8filter_2.fastq L8_raw

mkdir demultifilter

# run process_radtags
# upload input text files with barcodes (for -b)

/programs/stacks/bin/process_radtags -p /workdir/nrh44/L1_raw -b /workdir/nrh44/EUST_RAD_index1.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L2_raw -b /workdir/nrh44/EUST_RAD_index2.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L3_raw -b /workdir/nrh44/EUST_RAD_index3.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L4_raw -b /workdir/nrh44/EUST_RAD_index4.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L5_raw -b /workdir/nrh44/EUST_RAD_index5.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L6_raw -b /workdir/nrh44/EUST_RAD_index6.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L7_raw -b /workdir/nrh44/EUST_RAD_index7.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -p /workdir/nrh44/L8_raw -b /workdir/nrh44/EUST_RAD_index8.txt -o /workdir/nrh44/demultifilter -r -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina

# trim all sequences at their 3' end to the length of the shortest sequences

mkdir demultiplexed_trimmed

fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0504.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0504.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0509.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0509.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0107.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0107.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0402.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0402.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0505.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0505.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0203.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0203.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0403.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0403.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0110.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0110.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0107.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0107.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0203.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0203.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0402.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0402.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0303.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0303.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0506.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0506.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0402.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0402.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0106.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0106.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0107.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0107.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0501.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0501.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0108.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0108.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0108.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0108.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0405.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0405.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0507.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0507.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0205.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0205.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0106.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0106.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0109.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0109.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0503.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0503.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0106.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0106.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0201.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0201.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0207.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0207.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0106.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0106.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0302.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0302.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0206.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0206.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0402.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0402.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0304.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0304.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NM0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NM0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0309.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0309.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0501.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0501.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0406.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0406.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0402.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0402.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0409.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0409.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0201.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0201.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0201.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0201.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0504.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0504.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0110.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0110.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0305.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0305.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0504.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0504.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0203.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0203.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0310.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0310.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0304.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0304.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0201.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0201.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0204.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0204.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0403.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0403.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0307.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0307.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0407.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0407.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0109.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0109.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0406.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0406.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0407.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0407.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEE301.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEE301.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0301.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0301.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0301.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0301.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0404.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0404.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0502.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0502.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0202.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0202.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0503.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0503.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0408.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0408.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0203.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0203.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0305.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0305.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0107.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0107.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0205.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0205.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0201.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0201.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0207.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0207.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0404.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0404.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0403.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0403.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0508.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0508.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0404.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0404.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0405.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0405.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0106.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0106.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0507.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0507.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_MO0108.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_MO0108.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_KS0502.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_KS0502.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0308.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0308.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0209.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0209.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0406.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0406.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0405.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0405.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEN401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEN401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0208.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0208.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEW504.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEW504.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0202.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0202.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0202.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0202.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0303.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0303.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0209.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0209.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0410.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0410.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEE303.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEE303.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0403.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0403.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0103.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0103.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0202.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0202.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NC0302.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NC0302.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEN403.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEN403.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0204.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0204.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEN404.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEN404.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NV0306.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NV0306.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEW503.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEW503.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WA0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WA0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0205.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0205.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEW501.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEW501.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0104.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0104.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEE302.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEE302.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0206.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0206.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CO0102.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CO0102.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IL0401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IL0401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0401.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0401.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEE304.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEE304.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0502.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0502.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0204.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0204.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0501.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0501.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0506.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0506.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_ID0503.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_ID0503.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NEW502.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NEW502.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0109.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0109.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_CA0105.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_CA0105.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0107.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0107.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0208.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0208.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NH0505.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NH0505.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0207.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0207.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_NY0404.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_NY0404.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_IA0210.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_IA0210.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0108.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0108.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_AZ0101.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0101.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_TX0210.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_TX0210.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0208.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0208.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0206.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0206.fq
fastx_trimmer -Q33 -l90 -i /workdir/nrh44/demultifilter/EUST_WI0203.fq -o /workdir/nrh44/demultiplexed_trimmed/EUST_WI0203.fq

# confirm that sequences are 90 bp in length

wc -L demultiplexed_trimmed/*.fq

#################################### ASSEMBLY ####################################

# index the reference genome, -f means fasta input, EUST is the indexed file
# makes 6 output files with .bt2 as output

bowtie2-build -f /home/nrh44/EUSTref.fna EUST

# make folder for mapped files

mkdir /workdir/nrh44/mapped

# align every read to reference genome
# used end-to-end option (doesn't clip mismatches at the end)

bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0102.fq -S /workdir/nrh44/mapped/EUST_NM0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0105.fq -S /workdir/nrh44/mapped/EUST_NM0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0504.fq -S /workdir/nrh44/mapped/EUST_KS0504.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0509.fq -S /workdir/nrh44/mapped/EUST_KS0509.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0107.fq -S /workdir/nrh44/mapped/EUST_NM0107.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0103.fq -S /workdir/nrh44/mapped/EUST_MO0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0402.fq -S /workdir/nrh44/mapped/EUST_IL0402.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0505.fq -S /workdir/nrh44/mapped/EUST_KS0505.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0104.fq -S /workdir/nrh44/mapped/EUST_MO0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0203.fq -S /workdir/nrh44/mapped/EUST_NC0203.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0403.fq -S /workdir/nrh44/mapped/EUST_IL0403.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0110.fq -S /workdir/nrh44/mapped/EUST_AZ0110.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0107.fq -S /workdir/nrh44/mapped/EUST_MO0107.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0203.fq -S /workdir/nrh44/mapped/EUST_IA0203.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0402.fq -S /workdir/nrh44/mapped/EUST_ID0402.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0303.fq -S /workdir/nrh44/mapped/EUST_NV0303.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0506.fq -S /workdir/nrh44/mapped/EUST_KS0506.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0102.fq -S /workdir/nrh44/mapped/EUST_WA0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0402.fq -S /workdir/nrh44/mapped/EUST_CA0402.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0106.fq -S /workdir/nrh44/mapped/EUST_NM0106.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0107.fq -S /workdir/nrh44/mapped/EUST_CO0107.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0501.fq -S /workdir/nrh44/mapped/EUST_KS0501.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0108.fq -S /workdir/nrh44/mapped/EUST_CO0108.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0108.fq -S /workdir/nrh44/mapped/EUST_WA0108.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0405.fq -S /workdir/nrh44/mapped/EUST_IL0405.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0103.fq -S /workdir/nrh44/mapped/EUST_NC0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0507.fq -S /workdir/nrh44/mapped/EUST_KS0507.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0205.fq -S /workdir/nrh44/mapped/EUST_TX0205.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0103.fq -S /workdir/nrh44/mapped/EUST_WA0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0104.fq -S /workdir/nrh44/mapped/EUST_NM0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0106.fq -S /workdir/nrh44/mapped/EUST_CO0106.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0109.fq -S /workdir/nrh44/mapped/EUST_MO0109.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0503.fq -S /workdir/nrh44/mapped/EUST_NH0503.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0106.fq -S /workdir/nrh44/mapped/EUST_WA0106.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0201.fq -S /workdir/nrh44/mapped/EUST_IA0201.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0101.fq -S /workdir/nrh44/mapped/EUST_NM0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0207.fq -S /workdir/nrh44/mapped/EUST_TX0207.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0106.fq -S /workdir/nrh44/mapped/EUST_AZ0106.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0302.fq -S /workdir/nrh44/mapped/EUST_NV0302.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0206.fq -S /workdir/nrh44/mapped/EUST_TX0206.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0402.fq -S /workdir/nrh44/mapped/EUST_NH0402.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0304.fq -S /workdir/nrh44/mapped/EUST_NV0304.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0105.fq -S /workdir/nrh44/mapped/EUST_MO0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NM0103.fq -S /workdir/nrh44/mapped/EUST_NM0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0309.fq -S /workdir/nrh44/mapped/EUST_NV0309.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0501.fq -S /workdir/nrh44/mapped/EUST_ID0501.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0406.fq -S /workdir/nrh44/mapped/EUST_IL0406.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0402.fq -S /workdir/nrh44/mapped/EUST_NY0402.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0409.fq -S /workdir/nrh44/mapped/EUST_NY0409.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0201.fq -S /workdir/nrh44/mapped/EUST_TX0201.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0201.fq -S /workdir/nrh44/mapped/EUST_NC0201.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0504.fq -S /workdir/nrh44/mapped/EUST_ID0504.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0110.fq -S /workdir/nrh44/mapped/EUST_WA0110.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0101.fq -S /workdir/nrh44/mapped/EUST_CO0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0305.fq -S /workdir/nrh44/mapped/EUST_NC0305.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0504.fq -S /workdir/nrh44/mapped/EUST_NH0504.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0203.fq -S /workdir/nrh44/mapped/EUST_TX0203.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0310.fq -S /workdir/nrh44/mapped/EUST_NV0310.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0102.fq -S /workdir/nrh44/mapped/EUST_MO0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0101.fq -S /workdir/nrh44/mapped/EUST_MO0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0304.fq -S /workdir/nrh44/mapped/EUST_NC0304.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0201.fq -S /workdir/nrh44/mapped/EUST_CA0201.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0204.fq -S /workdir/nrh44/mapped/EUST_IA0204.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0403.fq -S /workdir/nrh44/mapped/EUST_ID0403.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0307.fq -S /workdir/nrh44/mapped/EUST_NV0307.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0407.fq -S /workdir/nrh44/mapped/EUST_NY0407.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0109.fq -S /workdir/nrh44/mapped/EUST_WA0109.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0406.fq -S /workdir/nrh44/mapped/EUST_NY0406.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0407.fq -S /workdir/nrh44/mapped/EUST_IL0407.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEE301.fq -S /workdir/nrh44/mapped/EUST_NEE301.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0301.fq -S /workdir/nrh44/mapped/EUST_NV0301.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0104.fq -S /workdir/nrh44/mapped/EUST_CA0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0301.fq -S /workdir/nrh44/mapped/EUST_NC0301.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0404.fq -S /workdir/nrh44/mapped/EUST_NH0404.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0502.fq -S /workdir/nrh44/mapped/EUST_ID0502.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0202.fq -S /workdir/nrh44/mapped/EUST_TX0202.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0503.fq -S /workdir/nrh44/mapped/EUST_KS0503.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0401.fq -S /workdir/nrh44/mapped/EUST_NH0401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0408.fq -S /workdir/nrh44/mapped/EUST_NY0408.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0203.fq -S /workdir/nrh44/mapped/EUST_CA0203.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0305.fq -S /workdir/nrh44/mapped/EUST_NV0305.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0107.fq -S /workdir/nrh44/mapped/EUST_WA0107.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0101.fq -S /workdir/nrh44/mapped/EUST_CA0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0401.fq -S /workdir/nrh44/mapped/EUST_NY0401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0205.fq -S /workdir/nrh44/mapped/EUST_IA0205.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0201.fq -S /workdir/nrh44/mapped/EUST_WI0201.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0207.fq -S /workdir/nrh44/mapped/EUST_IA0207.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0404.fq -S /workdir/nrh44/mapped/EUST_ID0404.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0403.fq -S /workdir/nrh44/mapped/EUST_CA0403.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0508.fq -S /workdir/nrh44/mapped/EUST_KS0508.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0101.fq -S /workdir/nrh44/mapped/EUST_WA0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0105.fq -S /workdir/nrh44/mapped/EUST_WA0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0404.fq -S /workdir/nrh44/mapped/EUST_IL0404.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0105.fq -S /workdir/nrh44/mapped/EUST_CO0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0405.fq -S /workdir/nrh44/mapped/EUST_NH0405.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0106.fq -S /workdir/nrh44/mapped/EUST_MO0106.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0103.fq -S /workdir/nrh44/mapped/EUST_AZ0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0507.fq -S /workdir/nrh44/mapped/EUST_NH0507.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_MO0108.fq -S /workdir/nrh44/mapped/EUST_MO0108.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_KS0502.fq -S /workdir/nrh44/mapped/EUST_KS0502.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0308.fq -S /workdir/nrh44/mapped/EUST_NV0308.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0209.fq -S /workdir/nrh44/mapped/EUST_IA0209.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0104.fq -S /workdir/nrh44/mapped/EUST_NC0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0102.fq -S /workdir/nrh44/mapped/EUST_NC0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0102.fq -S /workdir/nrh44/mapped/EUST_AZ0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0406.fq -S /workdir/nrh44/mapped/EUST_NH0406.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0401.fq -S /workdir/nrh44/mapped/EUST_CA0401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0405.fq -S /workdir/nrh44/mapped/EUST_NY0405.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEN401.fq -S /workdir/nrh44/mapped/EUST_NEN401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0208.fq -S /workdir/nrh44/mapped/EUST_TX0208.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0104.fq -S /workdir/nrh44/mapped/EUST_CO0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEW504.fq -S /workdir/nrh44/mapped/EUST_NEW504.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0202.fq -S /workdir/nrh44/mapped/EUST_WI0202.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0202.fq -S /workdir/nrh44/mapped/EUST_CA0202.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0101.fq -S /workdir/nrh44/mapped/EUST_NC0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0303.fq -S /workdir/nrh44/mapped/EUST_NC0303.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0209.fq -S /workdir/nrh44/mapped/EUST_TX0209.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0103.fq -S /workdir/nrh44/mapped/EUST_CO0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0410.fq -S /workdir/nrh44/mapped/EUST_NY0410.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0105.fq -S /workdir/nrh44/mapped/EUST_AZ0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEE303.fq -S /workdir/nrh44/mapped/EUST_NEE303.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0403.fq -S /workdir/nrh44/mapped/EUST_NY0403.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0103.fq -S /workdir/nrh44/mapped/EUST_CA0103.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0202.fq -S /workdir/nrh44/mapped/EUST_IA0202.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NC0302.fq -S /workdir/nrh44/mapped/EUST_NC0302.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEN403.fq -S /workdir/nrh44/mapped/EUST_NEN403.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0204.fq -S /workdir/nrh44/mapped/EUST_TX0204.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEN404.fq -S /workdir/nrh44/mapped/EUST_NEN404.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NV0306.fq -S /workdir/nrh44/mapped/EUST_NV0306.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEW503.fq -S /workdir/nrh44/mapped/EUST_NEW503.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WA0104.fq -S /workdir/nrh44/mapped/EUST_WA0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0205.fq -S /workdir/nrh44/mapped/EUST_WI0205.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEW501.fq -S /workdir/nrh44/mapped/EUST_NEW501.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0104.fq -S /workdir/nrh44/mapped/EUST_AZ0104.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEE302.fq -S /workdir/nrh44/mapped/EUST_NEE302.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0206.fq -S /workdir/nrh44/mapped/EUST_IA0206.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CO0102.fq -S /workdir/nrh44/mapped/EUST_CO0102.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IL0401.fq -S /workdir/nrh44/mapped/EUST_IL0401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0401.fq -S /workdir/nrh44/mapped/EUST_ID0401.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEE304.fq -S /workdir/nrh44/mapped/EUST_NEE304.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0502.fq -S /workdir/nrh44/mapped/EUST_NH0502.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0204.fq -S /workdir/nrh44/mapped/EUST_WI0204.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0501.fq -S /workdir/nrh44/mapped/EUST_CA0501.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0506.fq -S /workdir/nrh44/mapped/EUST_NH0506.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_ID0503.fq -S /workdir/nrh44/mapped/EUST_ID0503.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NEW502.fq -S /workdir/nrh44/mapped/EUST_NEW502.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0109.fq -S /workdir/nrh44/mapped/EUST_AZ0109.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_CA0105.fq -S /workdir/nrh44/mapped/EUST_CA0105.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0107.fq -S /workdir/nrh44/mapped/EUST_AZ0107.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0208.fq -S /workdir/nrh44/mapped/EUST_IA0208.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NH0505.fq -S /workdir/nrh44/mapped/EUST_NH0505.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0207.fq -S /workdir/nrh44/mapped/EUST_WI0207.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_NY0404.fq -S /workdir/nrh44/mapped/EUST_NY0404.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_IA0210.fq -S /workdir/nrh44/mapped/EUST_IA0210.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0108.fq -S /workdir/nrh44/mapped/EUST_AZ0108.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_AZ0101.fq -S /workdir/nrh44/mapped/EUST_AZ0101.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_TX0210.fq -S /workdir/nrh44/mapped/EUST_TX0210.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0208.fq -S /workdir/nrh44/mapped/EUST_WI0208.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0206.fq -S /workdir/nrh44/mapped/EUST_WI0206.sam
bowtie2 --phred33 --sensitive -x EUST -U /workdir/nrh44/demultiplexed_trimmed/EUST_WI0203.fq -S /workdir/nrh44/mapped/EUST_WI0203.sam

# manually remove samples with less than 15 mb of data, since poor quality data doesn't run well through STACKS
# all individuals had over 40 mb of data, still 159 individuals

# ref_map below runs the core pipeline 
# step 1: pstacks to build loci from reference and call SNPs
# step 2: cstacks assembles catalog
# step 3: sstacks matches to that catalog


# -m is minimum number of identical, raw reads to create a stack

mkdir /workdir/nrh44/eustrad
ref_map.pl -m 10 -n 5 -T 8 -B eust_radtags -S -b 1 -D "eust_radtags" -i 1 -o /workdir/nrh44/eustrad \
-s /workdir/nrh44/mapped/EUST_NM0102.sam \
-s /workdir/nrh44/mapped/EUST_NM0105.sam \
-s /workdir/nrh44/mapped/EUST_KS0504.sam \
-s /workdir/nrh44/mapped/EUST_KS0509.sam \
-s /workdir/nrh44/mapped/EUST_NM0107.sam \
-s /workdir/nrh44/mapped/EUST_MO0103.sam \
-s /workdir/nrh44/mapped/EUST_IL0402.sam \
-s /workdir/nrh44/mapped/EUST_KS0505.sam \
-s /workdir/nrh44/mapped/EUST_MO0104.sam \
-s /workdir/nrh44/mapped/EUST_NC0203.sam \
-s /workdir/nrh44/mapped/EUST_IL0403.sam \
-s /workdir/nrh44/mapped/EUST_AZ0110.sam \
-s /workdir/nrh44/mapped/EUST_MO0107.sam \
-s /workdir/nrh44/mapped/EUST_IA0203.sam \
-s /workdir/nrh44/mapped/EUST_ID0402.sam \
-s /workdir/nrh44/mapped/EUST_NV0303.sam \
-s /workdir/nrh44/mapped/EUST_KS0506.sam \
-s /workdir/nrh44/mapped/EUST_WA0102.sam \
-s /workdir/nrh44/mapped/EUST_CA0402.sam \
-s /workdir/nrh44/mapped/EUST_NM0106.sam \
-s /workdir/nrh44/mapped/EUST_CO0107.sam \
-s /workdir/nrh44/mapped/EUST_KS0501.sam \
-s /workdir/nrh44/mapped/EUST_CO0108.sam \
-s /workdir/nrh44/mapped/EUST_WA0108.sam \
-s /workdir/nrh44/mapped/EUST_IL0405.sam \
-s /workdir/nrh44/mapped/EUST_NC0103.sam \
-s /workdir/nrh44/mapped/EUST_KS0507.sam \
-s /workdir/nrh44/mapped/EUST_TX0205.sam \
-s /workdir/nrh44/mapped/EUST_WA0103.sam \
-s /workdir/nrh44/mapped/EUST_NM0104.sam \
-s /workdir/nrh44/mapped/EUST_CO0106.sam \
-s /workdir/nrh44/mapped/EUST_MO0109.sam \
-s /workdir/nrh44/mapped/EUST_NH0503.sam \
-s /workdir/nrh44/mapped/EUST_WA0106.sam \
-s /workdir/nrh44/mapped/EUST_IA0201.sam \
-s /workdir/nrh44/mapped/EUST_NM0101.sam \
-s /workdir/nrh44/mapped/EUST_TX0207.sam \
-s /workdir/nrh44/mapped/EUST_AZ0106.sam \
-s /workdir/nrh44/mapped/EUST_NV0302.sam \
-s /workdir/nrh44/mapped/EUST_TX0206.sam \
-s /workdir/nrh44/mapped/EUST_NH0402.sam \
-s /workdir/nrh44/mapped/EUST_NV0304.sam \
-s /workdir/nrh44/mapped/EUST_MO0105.sam \
-s /workdir/nrh44/mapped/EUST_NM0103.sam \
-s /workdir/nrh44/mapped/EUST_NV0309.sam \
-s /workdir/nrh44/mapped/EUST_ID0501.sam \
-s /workdir/nrh44/mapped/EUST_IL0406.sam \
-s /workdir/nrh44/mapped/EUST_NY0402.sam \
-s /workdir/nrh44/mapped/EUST_NY0409.sam \
-s /workdir/nrh44/mapped/EUST_TX0201.sam \
-s /workdir/nrh44/mapped/EUST_NC0201.sam \
-s /workdir/nrh44/mapped/EUST_ID0504.sam \
-s /workdir/nrh44/mapped/EUST_WA0110.sam \
-s /workdir/nrh44/mapped/EUST_CO0101.sam \
-s /workdir/nrh44/mapped/EUST_NC0305.sam \
-s /workdir/nrh44/mapped/EUST_NH0504.sam \
-s /workdir/nrh44/mapped/EUST_TX0203.sam \
-s /workdir/nrh44/mapped/EUST_NV0310.sam \
-s /workdir/nrh44/mapped/EUST_MO0102.sam \
-s /workdir/nrh44/mapped/EUST_MO0101.sam \
-s /workdir/nrh44/mapped/EUST_NC0304.sam \
-s /workdir/nrh44/mapped/EUST_CA0201.sam \
-s /workdir/nrh44/mapped/EUST_IA0204.sam \
-s /workdir/nrh44/mapped/EUST_ID0403.sam \
-s /workdir/nrh44/mapped/EUST_NV0307.sam \
-s /workdir/nrh44/mapped/EUST_NY0407.sam \
-s /workdir/nrh44/mapped/EUST_WA0109.sam \
-s /workdir/nrh44/mapped/EUST_NY0406.sam \
-s /workdir/nrh44/mapped/EUST_IL0407.sam \
-s /workdir/nrh44/mapped/EUST_NEE301.sam \
-s /workdir/nrh44/mapped/EUST_NV0301.sam \
-s /workdir/nrh44/mapped/EUST_CA0104.sam \
-s /workdir/nrh44/mapped/EUST_NC0301.sam \
-s /workdir/nrh44/mapped/EUST_NH0404.sam \
-s /workdir/nrh44/mapped/EUST_ID0502.sam \
-s /workdir/nrh44/mapped/EUST_TX0202.sam \
-s /workdir/nrh44/mapped/EUST_KS0503.sam \
-s /workdir/nrh44/mapped/EUST_NH0401.sam \
-s /workdir/nrh44/mapped/EUST_NY0408.sam \
-s /workdir/nrh44/mapped/EUST_CA0203.sam \
-s /workdir/nrh44/mapped/EUST_NV0305.sam \
-s /workdir/nrh44/mapped/EUST_WA0107.sam \
-s /workdir/nrh44/mapped/EUST_CA0101.sam \
-s /workdir/nrh44/mapped/EUST_NY0401.sam \
-s /workdir/nrh44/mapped/EUST_IA0205.sam \
-s /workdir/nrh44/mapped/EUST_WI0201.sam \
-s /workdir/nrh44/mapped/EUST_IA0207.sam \
-s /workdir/nrh44/mapped/EUST_ID0404.sam \
-s /workdir/nrh44/mapped/EUST_CA0403.sam \
-s /workdir/nrh44/mapped/EUST_KS0508.sam \
-s /workdir/nrh44/mapped/EUST_WA0101.sam \
-s /workdir/nrh44/mapped/EUST_WA0105.sam \
-s /workdir/nrh44/mapped/EUST_IL0404.sam \
-s /workdir/nrh44/mapped/EUST_CO0105.sam \
-s /workdir/nrh44/mapped/EUST_NH0405.sam \
-s /workdir/nrh44/mapped/EUST_MO0106.sam \
-s /workdir/nrh44/mapped/EUST_AZ0103.sam \
-s /workdir/nrh44/mapped/EUST_NH0507.sam \
-s /workdir/nrh44/mapped/EUST_MO0108.sam \
-s /workdir/nrh44/mapped/EUST_KS0502.sam \
-s /workdir/nrh44/mapped/EUST_NV0308.sam \
-s /workdir/nrh44/mapped/EUST_IA0209.sam \
-s /workdir/nrh44/mapped/EUST_NC0104.sam \
-s /workdir/nrh44/mapped/EUST_NC0102.sam \
-s /workdir/nrh44/mapped/EUST_AZ0102.sam \
-s /workdir/nrh44/mapped/EUST_NH0406.sam \
-s /workdir/nrh44/mapped/EUST_CA0401.sam \
-s /workdir/nrh44/mapped/EUST_NY0405.sam \
-s /workdir/nrh44/mapped/EUST_NEN401.sam \
-s /workdir/nrh44/mapped/EUST_TX0208.sam \
-s /workdir/nrh44/mapped/EUST_CO0104.sam \
-s /workdir/nrh44/mapped/EUST_NEW504.sam \
-s /workdir/nrh44/mapped/EUST_WI0202.sam \
-s /workdir/nrh44/mapped/EUST_CA0202.sam \
-s /workdir/nrh44/mapped/EUST_NC0101.sam \
-s /workdir/nrh44/mapped/EUST_NC0303.sam \
-s /workdir/nrh44/mapped/EUST_TX0209.sam \
-s /workdir/nrh44/mapped/EUST_CO0103.sam \
-s /workdir/nrh44/mapped/EUST_NY0410.sam \
-s /workdir/nrh44/mapped/EUST_AZ0105.sam \
-s /workdir/nrh44/mapped/EUST_NEE303.sam \
-s /workdir/nrh44/mapped/EUST_NY0403.sam \
-s /workdir/nrh44/mapped/EUST_CA0103.sam \
-s /workdir/nrh44/mapped/EUST_IA0202.sam \
-s /workdir/nrh44/mapped/EUST_NC0302.sam \
-s /workdir/nrh44/mapped/EUST_NEN403.sam \
-s /workdir/nrh44/mapped/EUST_TX0204.sam \
-s /workdir/nrh44/mapped/EUST_NEN404.sam \
-s /workdir/nrh44/mapped/EUST_NV0306.sam \
-s /workdir/nrh44/mapped/EUST_NEW503.sam \
-s /workdir/nrh44/mapped/EUST_WA0104.sam \
-s /workdir/nrh44/mapped/EUST_WI0205.sam \
-s /workdir/nrh44/mapped/EUST_NEW501.sam \
-s /workdir/nrh44/mapped/EUST_AZ0104.sam \
-s /workdir/nrh44/mapped/EUST_NEE302.sam \
-s /workdir/nrh44/mapped/EUST_IA0206.sam \
-s /workdir/nrh44/mapped/EUST_CO0102.sam \
-s /workdir/nrh44/mapped/EUST_IL0401.sam \
-s /workdir/nrh44/mapped/EUST_ID0401.sam \
-s /workdir/nrh44/mapped/EUST_NEE304.sam \
-s /workdir/nrh44/mapped/EUST_NH0502.sam \
-s /workdir/nrh44/mapped/EUST_WI0204.sam \
-s /workdir/nrh44/mapped/EUST_CA0501.sam \
-s /workdir/nrh44/mapped/EUST_NH0506.sam \
-s /workdir/nrh44/mapped/EUST_ID0503.sam \
-s /workdir/nrh44/mapped/EUST_NEW502.sam \
-s /workdir/nrh44/mapped/EUST_AZ0109.sam \
-s /workdir/nrh44/mapped/EUST_CA0105.sam \
-s /workdir/nrh44/mapped/EUST_AZ0107.sam \
-s /workdir/nrh44/mapped/EUST_IA0208.sam \
-s /workdir/nrh44/mapped/EUST_NH0505.sam \
-s /workdir/nrh44/mapped/EUST_WI0207.sam \
-s /workdir/nrh44/mapped/EUST_NY0404.sam \
-s /workdir/nrh44/mapped/EUST_IA0210.sam \
-s /workdir/nrh44/mapped/EUST_AZ0108.sam \
-s /workdir/nrh44/mapped/EUST_AZ0101.sam \
-s /workdir/nrh44/mapped/EUST_TX0210.sam \
-s /workdir/nrh44/mapped/EUST_WI0208.sam \
-s /workdir/nrh44/mapped/EUST_WI0206.sam \
-s /workdir/nrh44/mapped/EUST_WI0203.sam 

# now you've called SNPs, but you want to correct for quality
# so, we'll run rxstacks after running the core pipeline in ref_map.pl

# RxStacks makes corrections to genotype and haplotype calls in individual samples
# corrections include SNP model corrections, Log likelihood filtering, confounded locus filter, haplotype pruning

# --lnl_lim is minimum log likelihood requirement to keep a catalog locus
# --lnl_lim -100.0 sets the likelihood limit to -100 for the catalog locus to be kept. Lower (more negative) than this are thrown away

# Pre filter = 19679 loci, Post filter = ?? loci


mkdir /workdir/nrh44/rxstacks/
/programs/stacks/bin/rxstacks -b 1 -P /workdir/nrh44/eustrad -o /workdir/nrh44/rxstacks/ --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.1 --lnl_lim -300.0 --lnl_dist -t 24 --verbose
# t is threads, adjust for machine you're using

# once rxstacks is run to filter genotypes, cstacks and sstacks should be run again

# cstacks: builds a catalog from a set of samples. Creates a set of consensus loci.

/programs/stacks/bin/cstacks -b 1 -o /workdir/nrh44/rxstacks -p 15 -n 5 \
-s /workdir/nrh44/rxstacks/EUST_NM0102 \
-s /workdir/nrh44/rxstacks/EUST_NM0105 \
-s /workdir/nrh44/rxstacks/EUST_KS0504 \
-s /workdir/nrh44/rxstacks/EUST_KS0509 \
-s /workdir/nrh44/rxstacks/EUST_NM0107 \
-s /workdir/nrh44/rxstacks/EUST_MO0103 \
-s /workdir/nrh44/rxstacks/EUST_IL0402 \
-s /workdir/nrh44/rxstacks/EUST_KS0505 \
-s /workdir/nrh44/rxstacks/EUST_MO0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0203 \
-s /workdir/nrh44/rxstacks/EUST_IL0403 \
-s /workdir/nrh44/rxstacks/EUST_AZ0110 \
-s /workdir/nrh44/rxstacks/EUST_MO0107 \
-s /workdir/nrh44/rxstacks/EUST_IA0203 \
-s /workdir/nrh44/rxstacks/EUST_ID0402 \
-s /workdir/nrh44/rxstacks/EUST_NV0303 \
-s /workdir/nrh44/rxstacks/EUST_KS0506 \
-s /workdir/nrh44/rxstacks/EUST_WA0102 \
-s /workdir/nrh44/rxstacks/EUST_CA0402 \
-s /workdir/nrh44/rxstacks/EUST_NM0106 \
-s /workdir/nrh44/rxstacks/EUST_CO0107 \
-s /workdir/nrh44/rxstacks/EUST_KS0501 \
-s /workdir/nrh44/rxstacks/EUST_CO0108 \
-s /workdir/nrh44/rxstacks/EUST_WA0108 \
-s /workdir/nrh44/rxstacks/EUST_IL0405 \
-s /workdir/nrh44/rxstacks/EUST_NC0103 \
-s /workdir/nrh44/rxstacks/EUST_KS0507 \
-s /workdir/nrh44/rxstacks/EUST_TX0205 \
-s /workdir/nrh44/rxstacks/EUST_WA0103 \
-s /workdir/nrh44/rxstacks/EUST_NM0104 \
-s /workdir/nrh44/rxstacks/EUST_CO0106 \
-s /workdir/nrh44/rxstacks/EUST_MO0109 \
-s /workdir/nrh44/rxstacks/EUST_NH0503 \
-s /workdir/nrh44/rxstacks/EUST_WA0106 \
-s /workdir/nrh44/rxstacks/EUST_IA0201 \
-s /workdir/nrh44/rxstacks/EUST_NM0101 \
-s /workdir/nrh44/rxstacks/EUST_TX0207 \
-s /workdir/nrh44/rxstacks/EUST_AZ0106 \
-s /workdir/nrh44/rxstacks/EUST_NV0302 \
-s /workdir/nrh44/rxstacks/EUST_TX0206 \
-s /workdir/nrh44/rxstacks/EUST_NH0402 \
-s /workdir/nrh44/rxstacks/EUST_NV0304 \
-s /workdir/nrh44/rxstacks/EUST_MO0105 \
-s /workdir/nrh44/rxstacks/EUST_NM0103 \
-s /workdir/nrh44/rxstacks/EUST_NV0309 \
-s /workdir/nrh44/rxstacks/EUST_ID0501 \
-s /workdir/nrh44/rxstacks/EUST_IL0406 \
-s /workdir/nrh44/rxstacks/EUST_NY0402 \
-s /workdir/nrh44/rxstacks/EUST_NY0409 \
-s /workdir/nrh44/rxstacks/EUST_TX0201 \
-s /workdir/nrh44/rxstacks/EUST_NC0201 \
-s /workdir/nrh44/rxstacks/EUST_ID0504 \
-s /workdir/nrh44/rxstacks/EUST_WA0110 \
-s /workdir/nrh44/rxstacks/EUST_CO0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0305 \
-s /workdir/nrh44/rxstacks/EUST_NH0504 \
-s /workdir/nrh44/rxstacks/EUST_TX0203 \
-s /workdir/nrh44/rxstacks/EUST_NV0310 \
-s /workdir/nrh44/rxstacks/EUST_MO0102 \
-s /workdir/nrh44/rxstacks/EUST_MO0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0304 \
-s /workdir/nrh44/rxstacks/EUST_CA0201 \
-s /workdir/nrh44/rxstacks/EUST_IA0204 \
-s /workdir/nrh44/rxstacks/EUST_ID0403 \
-s /workdir/nrh44/rxstacks/EUST_NV0307 \
-s /workdir/nrh44/rxstacks/EUST_NY0407 \
-s /workdir/nrh44/rxstacks/EUST_WA0109 \
-s /workdir/nrh44/rxstacks/EUST_NY0406 \
-s /workdir/nrh44/rxstacks/EUST_IL0407 \
-s /workdir/nrh44/rxstacks/EUST_NEE301 \
-s /workdir/nrh44/rxstacks/EUST_NV0301 \
-s /workdir/nrh44/rxstacks/EUST_CA0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0301 \
-s /workdir/nrh44/rxstacks/EUST_NH0404 \
-s /workdir/nrh44/rxstacks/EUST_ID0502 \
-s /workdir/nrh44/rxstacks/EUST_TX0202 \
-s /workdir/nrh44/rxstacks/EUST_KS0503 \
-s /workdir/nrh44/rxstacks/EUST_NH0401 \
-s /workdir/nrh44/rxstacks/EUST_NY0408 \
-s /workdir/nrh44/rxstacks/EUST_CA0203 \
-s /workdir/nrh44/rxstacks/EUST_NV0305 \
-s /workdir/nrh44/rxstacks/EUST_WA0107 \
-s /workdir/nrh44/rxstacks/EUST_CA0101 \
-s /workdir/nrh44/rxstacks/EUST_NY0401 \
-s /workdir/nrh44/rxstacks/EUST_IA0205 \
-s /workdir/nrh44/rxstacks/EUST_WI0201 \
-s /workdir/nrh44/rxstacks/EUST_IA0207 \
-s /workdir/nrh44/rxstacks/EUST_ID0404 \
-s /workdir/nrh44/rxstacks/EUST_CA0403 \
-s /workdir/nrh44/rxstacks/EUST_KS0508 \
-s /workdir/nrh44/rxstacks/EUST_WA0101 \
-s /workdir/nrh44/rxstacks/EUST_WA0105 \
-s /workdir/nrh44/rxstacks/EUST_IL0404 \
-s /workdir/nrh44/rxstacks/EUST_CO0105 \
-s /workdir/nrh44/rxstacks/EUST_NH0405 \
-s /workdir/nrh44/rxstacks/EUST_MO0106 \
-s /workdir/nrh44/rxstacks/EUST_AZ0103 \
-s /workdir/nrh44/rxstacks/EUST_NH0507 \
-s /workdir/nrh44/rxstacks/EUST_MO0108 \
-s /workdir/nrh44/rxstacks/EUST_KS0502 \
-s /workdir/nrh44/rxstacks/EUST_NV0308 \
-s /workdir/nrh44/rxstacks/EUST_IA0209 \
-s /workdir/nrh44/rxstacks/EUST_NC0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0102 \
-s /workdir/nrh44/rxstacks/EUST_AZ0102 \
-s /workdir/nrh44/rxstacks/EUST_NH0406 \
-s /workdir/nrh44/rxstacks/EUST_CA0401 \
-s /workdir/nrh44/rxstacks/EUST_NY0405 \
-s /workdir/nrh44/rxstacks/EUST_NEN401 \
-s /workdir/nrh44/rxstacks/EUST_TX0208 \
-s /workdir/nrh44/rxstacks/EUST_CO0104 \
-s /workdir/nrh44/rxstacks/EUST_NEW504 \
-s /workdir/nrh44/rxstacks/EUST_WI0202 \
-s /workdir/nrh44/rxstacks/EUST_CA0202 \
-s /workdir/nrh44/rxstacks/EUST_NC0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0303 \
-s /workdir/nrh44/rxstacks/EUST_TX0209 \
-s /workdir/nrh44/rxstacks/EUST_CO0103 \
-s /workdir/nrh44/rxstacks/EUST_NY0410 \
-s /workdir/nrh44/rxstacks/EUST_AZ0105 \
-s /workdir/nrh44/rxstacks/EUST_NEE303 \
-s /workdir/nrh44/rxstacks/EUST_NY0403 \
-s /workdir/nrh44/rxstacks/EUST_CA0103 \
-s /workdir/nrh44/rxstacks/EUST_IA0202 \
-s /workdir/nrh44/rxstacks/EUST_NC0302 \
-s /workdir/nrh44/rxstacks/EUST_NEN403 \
-s /workdir/nrh44/rxstacks/EUST_TX0204 \
-s /workdir/nrh44/rxstacks/EUST_NEN404 \
-s /workdir/nrh44/rxstacks/EUST_NV0306 \
-s /workdir/nrh44/rxstacks/EUST_NEW503 \
-s /workdir/nrh44/rxstacks/EUST_WA0104 \
-s /workdir/nrh44/rxstacks/EUST_WI0205 \
-s /workdir/nrh44/rxstacks/EUST_NEW501 \
-s /workdir/nrh44/rxstacks/EUST_AZ0104 \
-s /workdir/nrh44/rxstacks/EUST_NEE302 \
-s /workdir/nrh44/rxstacks/EUST_IA0206 \
-s /workdir/nrh44/rxstacks/EUST_CO0102 \
-s /workdir/nrh44/rxstacks/EUST_IL0401 \
-s /workdir/nrh44/rxstacks/EUST_ID0401 \
-s /workdir/nrh44/rxstacks/EUST_NEE304 \
-s /workdir/nrh44/rxstacks/EUST_NH0502 \
-s /workdir/nrh44/rxstacks/EUST_WI0204 \
-s /workdir/nrh44/rxstacks/EUST_CA0501 \
-s /workdir/nrh44/rxstacks/EUST_NH0506 \
-s /workdir/nrh44/rxstacks/EUST_ID0503 \
-s /workdir/nrh44/rxstacks/EUST_NEW502 \
-s /workdir/nrh44/rxstacks/EUST_AZ0109 \
-s /workdir/nrh44/rxstacks/EUST_CA0105 \
-s /workdir/nrh44/rxstacks/EUST_AZ0107 \
-s /workdir/nrh44/rxstacks/EUST_IA0208 \
-s /workdir/nrh44/rxstacks/EUST_NH0505 \
-s /workdir/nrh44/rxstacks/EUST_WI0207 \
-s /workdir/nrh44/rxstacks/EUST_NY0404 \
-s /workdir/nrh44/rxstacks/EUST_IA0210 \
-s /workdir/nrh44/rxstacks/EUST_AZ0108 \
-s /workdir/nrh44/rxstacks/EUST_AZ0101 \
-s /workdir/nrh44/rxstacks/EUST_TX0210 \
-s /workdir/nrh44/rxstacks/EUST_WI0208 \
-s /workdir/nrh44/rxstacks/EUST_WI0206 \
-s /workdir/nrh44/rxstacks/EUST_WI0203 

# run sstacks: sets of stacks constructed by ustacks can be searched agains the catalog produced by cstacks. 

/programs/stacks/bin/sstacks -b 1 -c /workdir/nrh44/rxstacks/batch_1 -p 15 -o /workdir/nrh44/rxstacks \
-s /workdir/nrh44/rxstacks/EUST_NM0102 \
-s /workdir/nrh44/rxstacks/EUST_NM0105 \
-s /workdir/nrh44/rxstacks/EUST_KS0504 \
-s /workdir/nrh44/rxstacks/EUST_KS0509 \
-s /workdir/nrh44/rxstacks/EUST_NM0107 \
-s /workdir/nrh44/rxstacks/EUST_MO0103 \
-s /workdir/nrh44/rxstacks/EUST_IL0402 \
-s /workdir/nrh44/rxstacks/EUST_KS0505 \
-s /workdir/nrh44/rxstacks/EUST_MO0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0203 \
-s /workdir/nrh44/rxstacks/EUST_IL0403 \
-s /workdir/nrh44/rxstacks/EUST_AZ0110 \
-s /workdir/nrh44/rxstacks/EUST_MO0107 \
-s /workdir/nrh44/rxstacks/EUST_IA0203 \
-s /workdir/nrh44/rxstacks/EUST_ID0402 \
-s /workdir/nrh44/rxstacks/EUST_NV0303 \
-s /workdir/nrh44/rxstacks/EUST_KS0506 \
-s /workdir/nrh44/rxstacks/EUST_WA0102 \
-s /workdir/nrh44/rxstacks/EUST_CA0402 \
-s /workdir/nrh44/rxstacks/EUST_NM0106 \
-s /workdir/nrh44/rxstacks/EUST_CO0107 \
-s /workdir/nrh44/rxstacks/EUST_KS0501 \
-s /workdir/nrh44/rxstacks/EUST_CO0108 \
-s /workdir/nrh44/rxstacks/EUST_WA0108 \
-s /workdir/nrh44/rxstacks/EUST_IL0405 \
-s /workdir/nrh44/rxstacks/EUST_NC0103 \
-s /workdir/nrh44/rxstacks/EUST_KS0507 \
-s /workdir/nrh44/rxstacks/EUST_TX0205 \
-s /workdir/nrh44/rxstacks/EUST_WA0103 \
-s /workdir/nrh44/rxstacks/EUST_NM0104 \
-s /workdir/nrh44/rxstacks/EUST_CO0106 \
-s /workdir/nrh44/rxstacks/EUST_MO0109 \
-s /workdir/nrh44/rxstacks/EUST_NH0503 \
-s /workdir/nrh44/rxstacks/EUST_WA0106 \
-s /workdir/nrh44/rxstacks/EUST_IA0201 \
-s /workdir/nrh44/rxstacks/EUST_NM0101 \
-s /workdir/nrh44/rxstacks/EUST_TX0207 \
-s /workdir/nrh44/rxstacks/EUST_AZ0106 \
-s /workdir/nrh44/rxstacks/EUST_NV0302 \
-s /workdir/nrh44/rxstacks/EUST_TX0206 \
-s /workdir/nrh44/rxstacks/EUST_NH0402 \
-s /workdir/nrh44/rxstacks/EUST_NV0304 \
-s /workdir/nrh44/rxstacks/EUST_MO0105 \
-s /workdir/nrh44/rxstacks/EUST_NM0103 \
-s /workdir/nrh44/rxstacks/EUST_NV0309 \
-s /workdir/nrh44/rxstacks/EUST_ID0501 \
-s /workdir/nrh44/rxstacks/EUST_IL0406 \
-s /workdir/nrh44/rxstacks/EUST_NY0402 \
-s /workdir/nrh44/rxstacks/EUST_NY0409 \
-s /workdir/nrh44/rxstacks/EUST_TX0201 \
-s /workdir/nrh44/rxstacks/EUST_NC0201 \
-s /workdir/nrh44/rxstacks/EUST_ID0504 \
-s /workdir/nrh44/rxstacks/EUST_WA0110 \
-s /workdir/nrh44/rxstacks/EUST_CO0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0305 \
-s /workdir/nrh44/rxstacks/EUST_NH0504 \
-s /workdir/nrh44/rxstacks/EUST_TX0203 \
-s /workdir/nrh44/rxstacks/EUST_NV0310 \
-s /workdir/nrh44/rxstacks/EUST_MO0102 \
-s /workdir/nrh44/rxstacks/EUST_MO0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0304 \
-s /workdir/nrh44/rxstacks/EUST_CA0201 \
-s /workdir/nrh44/rxstacks/EUST_IA0204 \
-s /workdir/nrh44/rxstacks/EUST_ID0403 \
-s /workdir/nrh44/rxstacks/EUST_NV0307 \
-s /workdir/nrh44/rxstacks/EUST_NY0407 \
-s /workdir/nrh44/rxstacks/EUST_WA0109 \
-s /workdir/nrh44/rxstacks/EUST_NY0406 \
-s /workdir/nrh44/rxstacks/EUST_IL0407 \
-s /workdir/nrh44/rxstacks/EUST_NEE301 \
-s /workdir/nrh44/rxstacks/EUST_NV0301 \
-s /workdir/nrh44/rxstacks/EUST_CA0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0301 \
-s /workdir/nrh44/rxstacks/EUST_NH0404 \
-s /workdir/nrh44/rxstacks/EUST_ID0502 \
-s /workdir/nrh44/rxstacks/EUST_TX0202 \
-s /workdir/nrh44/rxstacks/EUST_KS0503 \
-s /workdir/nrh44/rxstacks/EUST_NH0401 \
-s /workdir/nrh44/rxstacks/EUST_NY0408 \
-s /workdir/nrh44/rxstacks/EUST_CA0203 \
-s /workdir/nrh44/rxstacks/EUST_NV0305 \
-s /workdir/nrh44/rxstacks/EUST_WA0107 \
-s /workdir/nrh44/rxstacks/EUST_CA0101 \
-s /workdir/nrh44/rxstacks/EUST_NY0401 \
-s /workdir/nrh44/rxstacks/EUST_IA0205 \
-s /workdir/nrh44/rxstacks/EUST_WI0201 \
-s /workdir/nrh44/rxstacks/EUST_IA0207 \
-s /workdir/nrh44/rxstacks/EUST_ID0404 \
-s /workdir/nrh44/rxstacks/EUST_CA0403 \
-s /workdir/nrh44/rxstacks/EUST_KS0508 \
-s /workdir/nrh44/rxstacks/EUST_WA0101 \
-s /workdir/nrh44/rxstacks/EUST_WA0105 \
-s /workdir/nrh44/rxstacks/EUST_IL0404 \
-s /workdir/nrh44/rxstacks/EUST_CO0105 \
-s /workdir/nrh44/rxstacks/EUST_NH0405 \
-s /workdir/nrh44/rxstacks/EUST_MO0106 \
-s /workdir/nrh44/rxstacks/EUST_AZ0103 \
-s /workdir/nrh44/rxstacks/EUST_NH0507 \
-s /workdir/nrh44/rxstacks/EUST_MO0108 \
-s /workdir/nrh44/rxstacks/EUST_KS0502 \
-s /workdir/nrh44/rxstacks/EUST_NV0308 \
-s /workdir/nrh44/rxstacks/EUST_IA0209 \
-s /workdir/nrh44/rxstacks/EUST_NC0104 \
-s /workdir/nrh44/rxstacks/EUST_NC0102 \
-s /workdir/nrh44/rxstacks/EUST_AZ0102 \
-s /workdir/nrh44/rxstacks/EUST_NH0406 \
-s /workdir/nrh44/rxstacks/EUST_CA0401 \
-s /workdir/nrh44/rxstacks/EUST_NY0405 \
-s /workdir/nrh44/rxstacks/EUST_NEN401 \
-s /workdir/nrh44/rxstacks/EUST_TX0208 \
-s /workdir/nrh44/rxstacks/EUST_CO0104 \
-s /workdir/nrh44/rxstacks/EUST_NEW504 \
-s /workdir/nrh44/rxstacks/EUST_WI0202 \
-s /workdir/nrh44/rxstacks/EUST_CA0202 \
-s /workdir/nrh44/rxstacks/EUST_NC0101 \
-s /workdir/nrh44/rxstacks/EUST_NC0303 \
-s /workdir/nrh44/rxstacks/EUST_TX0209 \
-s /workdir/nrh44/rxstacks/EUST_CO0103 \
-s /workdir/nrh44/rxstacks/EUST_NY0410 \
-s /workdir/nrh44/rxstacks/EUST_AZ0105 \
-s /workdir/nrh44/rxstacks/EUST_NEE303 \
-s /workdir/nrh44/rxstacks/EUST_NY0403 \
-s /workdir/nrh44/rxstacks/EUST_CA0103 \
-s /workdir/nrh44/rxstacks/EUST_IA0202 \
-s /workdir/nrh44/rxstacks/EUST_NC0302 \
-s /workdir/nrh44/rxstacks/EUST_NEN403 \
-s /workdir/nrh44/rxstacks/EUST_TX0204 \
-s /workdir/nrh44/rxstacks/EUST_NEN404 \
-s /workdir/nrh44/rxstacks/EUST_NV0306 \
-s /workdir/nrh44/rxstacks/EUST_NEW503 \
-s /workdir/nrh44/rxstacks/EUST_WA0104 \
-s /workdir/nrh44/rxstacks/EUST_WI0205 \
-s /workdir/nrh44/rxstacks/EUST_NEW501 \
-s /workdir/nrh44/rxstacks/EUST_AZ0104 \
-s /workdir/nrh44/rxstacks/EUST_NEE302 \
-s /workdir/nrh44/rxstacks/EUST_IA0206 \
-s /workdir/nrh44/rxstacks/EUST_CO0102 \
-s /workdir/nrh44/rxstacks/EUST_IL0401 \
-s /workdir/nrh44/rxstacks/EUST_ID0401 \
-s /workdir/nrh44/rxstacks/EUST_NEE304 \
-s /workdir/nrh44/rxstacks/EUST_NH0502 \
-s /workdir/nrh44/rxstacks/EUST_WI0204 \
-s /workdir/nrh44/rxstacks/EUST_CA0501 \
-s /workdir/nrh44/rxstacks/EUST_NH0506 \
-s /workdir/nrh44/rxstacks/EUST_ID0503 \
-s /workdir/nrh44/rxstacks/EUST_NEW502 \
-s /workdir/nrh44/rxstacks/EUST_AZ0109 \
-s /workdir/nrh44/rxstacks/EUST_CA0105 \
-s /workdir/nrh44/rxstacks/EUST_AZ0107 \
-s /workdir/nrh44/rxstacks/EUST_IA0208 \
-s /workdir/nrh44/rxstacks/EUST_NH0505 \
-s /workdir/nrh44/rxstacks/EUST_WI0207 \
-s /workdir/nrh44/rxstacks/EUST_NY0404 \
-s /workdir/nrh44/rxstacks/EUST_IA0210 \
-s /workdir/nrh44/rxstacks/EUST_AZ0108 \
-s /workdir/nrh44/rxstacks/EUST_AZ0101 \
-s /workdir/nrh44/rxstacks/EUST_TX0210 \
-s /workdir/nrh44/rxstacks/EUST_WI0208 \
-s /workdir/nrh44/rxstacks/EUST_WI0206 \
-s /workdir/nrh44/rxstacks/EUST_WI0203 
# this is now in rxstacks2 directory in SNPcalling

### run populations ###
# pop map includes all 17 putative populations

# make populations.txt file (individual name *tab* string; example EUST_WI0203	wi), where string varies based on populations you're calling  

# symbolic links to input files needed (since only using these for information, not changing directly)
ln -s /home/nrh44/EUST_RADseq_2018/Filtering/SNPcalling/rxstacks2/* ./
ln -s /home/nrh44/EUST_RADseq_2018/Populations/populations.txt ./
ln -s /home/nrh44/EUST_RADseq_2018/Populations/populationsF.txt ./
ln -s /home/nrh44/EUST_RADseq_2018/Populations/populations1.txt ./


# note that we want all SNPs *and* one SNP because different downstream analyses use different data
# example: STRUCTURE wants SNPs that aren't in tight linkage, so use one SNP file

# run populations - pick all SNPs per locus - one population
# -r --- percentage of individuals in a population (80%)
# -p --- minimum number of populations a locus must be present in (1) 
# -m --- minimum stack depth (10)
# -t --- threads to run program on

# all SNPs for Bayescan

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/tsv -M /workdir/nrh44/populationsF.txt -r 0.8 -p 1 -m 10 -t 15 --structure --vcf --fasta_strict

mv batch_1.sumstats.tsv EUSTallSNPr8.sumstats.tsv
mv batch_1.haplotypes.tsv EUSTallSNPr8.haplotypes.tsv
mv batch_1.strict.fa EUSTallSNPr8.strict.fa
mv batch_1.structure.tsv EUSTallSNPr8.structure.tsv
mv batch_1.vcf EUSTallSNPr8.vcf
mv batch_1.hapstats.tsv EUSTallSNPr8.hapstats.tsv
mv batch_1.sumstats_summary.tsv EUSTallSNPr8.sumstats_summary.tsv

# checked for missing data in R
# removed TX0202 from populations.txt so that not overly penalizing for missing data
# re-ran code with updated populations*.txt file


### checking coverage again
# get depth information
cp *.bt2 /workdir/nrh44
cp *.fq /workdir/nrh44
vcftools --vcf EUSTallSNPr8.vcf --depth -c > EUSTallSNPr8_depth_summary.txt
samtools depth EUST_AZ0101cov.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
7.518

# get coverage
# map and sort
bowtie2 -x EUST --no-unal -U EUST_AZ0101.fq -S EUST_AZ0101cov.sam -p 12
samtools sort EUST_AZ0101cov.sam -o EUST_AZ0101cov.bam
# create samtools index
samtools index EUST_AZ0101cov.bam
# get total number of bases covered at MIN_COVERAGE_DEPTH or higher 
MIN_COVERAGE_DEPTH=5
samtools mpileup EUST_AZ0101cov.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l

3569540

# get length of reference genome
bowtie2-inspect -s EUST | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}'

1036755994

coverage = (3569540*27)/1036755994 = 0.9% 


bowtie2 -x EUST --no-unal -U EUST_NY0408.fq -S EUST_NY0408cov.sam -p 12
samtools sort EUST_NY0408cov.sam -o EUST_NY0408cov.bam
samtools index EUST_NY0408cov.bam
samtools mpileup EUST_NY0408cov.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l

MIN 5
1009803

bowtie2 -x EUST --no-unal -U EUST_NM0104.fq -S EUST_NM0104cov.sam -p 12
samtools sort EUST_NM0104cov.sam -o EUST_NM0104cov.bam
samtools index EUST_NM0104cov.bam
samtools mpileup EUST_NM0104cov.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$4>=X' | wc -l

MIN 5
1014210

# check for relatedness
vcftools --vcf EUSTallSNPr8.vcf --relatedness
# IA0203 and IA0204 are related: 0.954346
# remove IA0204 from populations.txt

# re-run all SNPS and overwrite EUSTallSNPr8 with problem individuals removed
# this has no MAF filter!!!!

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations1.txt -r 0.8 -p 1 -m 10 -t 15 --structure --vcf --fasta_strict

mv batch_1.sumstats.tsv EUSTallSNPr8.sumstats.tsv
mv batch_1.haplotypes.tsv EUSTallSNPr8.haplotypes.tsv
mv batch_1.strict.fa EUSTallSNPr8.strict.fa
mv batch_1.structure.tsv EUSTallSNPr8.structure.tsv
mv batch_1.vcf EUSTallSNPr8.vcf
mv batch_1.hapstats.tsv EUSTallSNPr8.hapstats.tsv
mv batch_1.sumstats_summary.tsv EUSTallSNPr8.sumstats_summary.tsv

cp EUSTallSNPr8* /home/nrh44/EUST_RADseq_2018/Populations/forBayescanr8allSNP &

# re-run with MAF 0.1, all pops
# use this for EEMS, other analyses!

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations.txt -r 0.8 -p 1 -m 10 --min_maf 0.1 -t 15 --structure --vcf --fasta_strict

mv batch_1.sumstats.tsv EUSTallSNPr8maf01p.sumstats.tsv
mv batch_1.haplotypes.tsv EUSTallSNPr8maf01p.haplotypes.tsv
mv batch_1.strict.fa EUSTallSNPr8maf01p.strict.fa
mv batch_1.structure.tsv EUSTallSNPr8maf01p.structure.tsv
mv batch_1.vcf EUSTallSNPr8maf01p.vcf
mv batch_1.hapstats.tsv EUSTallSNPr8maf01p.hapstats.tsv
mv batch_1.sumstats_summary.tsv EUSTallSNPr8maf01p.sumstats_summary.tsv
mv batch_1.populations.log EUSTallSNPr8maf01p.populations.log 

cp EUSTallSNPr8maf01p* /home/nrh44/EUST_RADseq_2018/Populations/forBayescanr8allSNPmaf01p &

# one SNP for STRUCTURE - pick the first SNP per locus
# -f gives p-value correction to Fst score 
# popmap with one putative population
# test for sensitivity in allele frequency for STRUCTURE, PCA

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations1.txt -f p_value --p_value_cutoff 0.05 -r 0.8 -p 1 -m 5 --min_maf 0.3 -t 15 --structure --vcf --fstats --write_single_snp &

mv batch_1.sumstats.tsv EUSToneSNPr8maf03.sumstats.tsv
mv batch_1.sumstats_summary.tsv EUSToneSNPr8maf03.sumstats_summary.tsv
mv batch_1.phistats.tsv EUSToneSNPr8maf03.phistats.tsv
mv batch_1.hapstats.tsv EUSToneSNPr8maf03.hapstats.tsv
mv batch_1.haplotypes.tsv EUSToneSNPr8maf03.haplotypes.tsv
mv batch_1.strict.fa EUSToneSNPr8maf03.strict.fa
mv batch_1.structure.tsv EUSToneSNPr8maf03.structure.tsv
mv batch_1.vcf EUSToneSNPr8maf03.vcf
mv batch_1.populations.log EUSToneSNPr8maf03.populations.log 

cp EUSToneSNPr8maf03* /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF03

# one SNP for STRUCTURE with populations specified in pop map, maf cut-off 0.3

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations.txt -f p_value --p_value_cutoff 0.05 -r 0.8 -p 1 -m 5 --min_maf 0.3 -t 15 --structure --vcf --fstats --write_single_snp &

mv batch_1.sumstats.tsv EUSToneSNPr8maf03p.sumstats.tsv
mv batch_1.sumstats_summary.tsv EUSToneSNPr8maf03p.sumstats_summary.tsv
mv batch_1.phistats.tsv EUSToneSNPr8maf03p.phistats.tsv
mv batch_1.hapstats.tsv EUSToneSNPr8maf03p.hapstats.tsv
mv batch_1.haplotypes.tsv EUSToneSNPr8maf03p.haplotypes.tsv
mv batch_1.structure.tsv EUSToneSNPr8maf03p.structure.tsv
mv batch_1.vcf EUSToneSNPr8maf03p.vcf
mv batch_1.populations.log EUSToneSNPr8maf03p.populations.log 

cp batch_1.phistats_*.tsv /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF03p
cp batch_1.fst_*.tsv /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF03p
cp EUSToneSNP* /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF03p

# with pop map, lower MAF cut-off of 0.1

/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations.txt -f p_value --p_value_cutoff 0.05 -r 0.8 -p 1 -m 5 --min_maf 0.1 -t 15 --structure --vcf --fstats --write_single_snp &

mv batch_1.sumstats.tsv EUSToneSNPr8maf01p.sumstats.tsv
mv batch_1.sumstats_summary.tsv EUSToneSNPr8maf01p.sumstats_summary.tsv
mv batch_1.phistats.tsv EUSToneSNPr8maf01p.phistats.tsv
mv batch_1.hapstats.tsv EUSToneSNPr8maf01p.hapstats.tsv
mv batch_1.haplotypes.tsv EUSToneSNPr8maf01p.haplotypes.tsv
mv batch_1.structure.tsv EUSToneSNPr8maf01p.structure.tsv
mv batch_1.vcf EUSToneSNPr8maf01p.vcf
mv batch_1.populations.log EUSToneSNPr8maf01p.populations.log 

cp batch_1.phistats_*.tsv /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF01p
cp batch_1.fst_*.tsv /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF01p
cp EUSToneSNP* /home/nrh44/EUST_RADseq_2018/Populations/forStructureMAF01p

# do rare alleles drive patterns?

vcftools --vcf EUSTallSNPr8.vcf --max-maf 0.3 --min-alleles 2 --max-alleles 2 --recode --out EUSTallSNPr8_rare.vcf
vcftools --vcf EUSTallSNPr8.vcf --max-maf 0.1 --min-alleles 2 --max-alleles 2 --recode --out EUSTallSNPr8_rarer.vcf

# only non-missing sites for RDA

vcftools --vcf EUSTallSNPr8.vcf --max-missing 1 --recode --out EUSTallSNPr8.nomiss

vcftools --vcf EUSTallSNPr8.nomiss.vcf.recode.vcf --012 --out EUSTallSNPr8.nomiss

# input for GPhoCS
vcftools --vcf EUSTallSNPr8.vcf --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --out EUSTallSNPr8_strict.vcf

# bumped up r to 0.95, make sure all one population
/programs/stacks-1.44/bin/populations -b 1 -P /workdir/nrh44/ -M /workdir/nrh44/populations1.txt -r 0.95 -p 1 -m 10 -t 15 --structure --vcf --fasta_strict &> gphocsfasta.log &

mv batch_1.sumstats.tsv EUSTallSNPr95.sumstats.tsv
mv batch_1.sumstats_summary.tsv EUSTallSNPr95.sumstats_summary.tsv
mv batch_1.strict.fa EUSTallSNPr95.strict.fa
mv batch_1.hapstats.tsv EUSTallSNPr95.hapstats.tsv
mv batch_1.haplotypes.tsv EUSTallSNPr95.haplotypes.tsv
mv batch_1.structure.tsv EUSTallSNPr95.structure.tsv
mv batch_1.vcf EUSTallSNPr95.vcf
mv batch_1.populations.log EUSTallSNPr95.populations.log 

# send .strict.fa to Leo for conversion

### POP GEN SUMMARY STATS

vcftools --vcf EUSTallSNPr8.vcf --freq --out EUSTallSNPr8 &


vcftools --vcf EUSToneSNPr8maf01p.vcf --freq --out EUSToneSNPr8maf01p &> RADfreq.log &
vcftools --vcf EUSToneSNPr8maf01p.vcf --het --out EUSToneSNPr8maf01p &> RADhet.log &
vcftools --vcf EUSToneSNPr8maf01p.vcf --hardy --out EUSToneSNPr8maf01p &> RADhardy.log &
vcftools --vcf EUSToneSNPr8maf01p.vcf --TajimaD 25000 --out EUSToneSNPr8maf01p25kb &> RADtajimad.log &
vcftools --vcf EUSToneSNPr8maf01p.vcf --window-pi 25000 --out EUSToneSNPr8maf01p25kb &> RADpi25kb.log &

vcftools --vcf EUSTallSNPr8maf01p.vcf --freq --out EUSTallSNPr8maf01p &> RADfreq.all.log &
vcftools --vcf EUSTallSNPr8maf01p.vcf --het --out EUSTallSNPr8maf01p &> RADhet.all.log &
vcftools --vcf EUSTallSNPr8maf01p.vcf --hardy --out EUSTallSNPr8maf01p &> RADhardy.all.log &
vcftools --vcf EUSTallSNPr8maf01p.vcf --TajimaD 25000 --out EUSTallSNPr8maf01p25kb &> RADtajimad.all.log &
vcftools --vcf EUSTallSNPr8maf01p.vcf --window-pi 25000 --out EUSTallSNPr8maf01p25kb &> RADpi25kb.all.log &

vcftools --vcf EUSTallSNPr8maf01p.vcf --fst-window-size 25000 --weir-fst-pop popeast.txt --weir-fst-pop popwest.txt --out EUSTallSNPr8maf01pfst25kb &> RADfst25kb.all.log

vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop CA.txt --out az_v_ca
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop CO.txt --out az_v_co
# run PCA in R



cp EUSTallSNPr8maf01p* /home/nrh44/EUST_RADseq_2018/popgen_stats

# convert to all plink formats possible
vcftools --vcf EUSTallSNPr8maf01p.vcf --out EUSTallSNPr8maf01p.plink --plink
plink --file EUSTallSNPr8maf01p.plink --out EUSTallSNPr8maf01p --recodeA
plink --map EUSTallSNPr8maf01p.plink.map --ped EUSTallSNPr8maf01p.plink.ped --make-bed --out EUSTallSNPr8maf01p.plink



vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop IA.txt --out az_v_ia
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop ID.txt --out az_v_id
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop IL.txt --out az_v_il
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop KS.txt --out az_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop MO.txt --out az_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NC.txt --out az_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NE.txt --out az_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NH.txt --out az_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NM.txt --out az_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NV.txt --out az_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop NY.txt --out az_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop TX.txt --out az_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop WA.txt --out az_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop AZ.txt --weir-fst-pop WI.txt --out az_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop CO.txt --out ca_v_co
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop IA.txt --out ca_v_ia
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop ID.txt --out ca_v_id
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop IL.txt --out ca_v_il
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop KS.txt --out ca_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop MO.txt --out ca_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NC.txt --out ca_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NE.txt --out ca_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NH.txt --out ca_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NM.txt --out ca_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NV.txt --out ca_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop NY.txt --out ca_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop TX.txt --out ca_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop WA.txt --out ca_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CA.txt --weir-fst-pop WI.txt --out ca_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop IA.txt --out co_v_ia
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop ID.txt --out co_v_id
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop IL.txt --out co_v_il
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop KS.txt --out co_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop MO.txt --out co_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NC.txt --out co_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NE.txt --out co_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NH.txt --out co_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NM.txt --out co_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NV.txt --out co_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop NY.txt --out co_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop TX.txt --out co_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop WA.txt --out co_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop CO.txt --weir-fst-pop WI.txt --out co_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop ID.txt --out ia_v_id
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop IL.txt --out ia_v_il
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop KS.txt --out ia_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop MO.txt --out ia_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NC.txt --out ia_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NE.txt --out ia_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NH.txt --out ia_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NM.txt --out ia_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NV.txt --out ia_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop NY.txt --out ia_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop TX.txt --out ia_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop WA.txt --out ia_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IA.txt --weir-fst-pop WI.txt --out ia_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop IL.txt --out id_v_il
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop KS.txt --out id_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop MO.txt --out id_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NC.txt --out id_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NE.txt --out id_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NH.txt --out id_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NM.txt --out id_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NV.txt --out id_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop NY.txt --out id_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop TX.txt --out id_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop WA.txt --out id_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop ID.txt --weir-fst-pop WI.txt --out id_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop KS.txt --out il_v_ks
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop MO.txt --out il_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NC.txt --out il_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NE.txt --out il_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NH.txt --out il_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NM.txt --out il_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NV.txt --out il_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop NY.txt --out il_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop TX.txt --out il_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop WA.txt --out il_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop IL.txt --weir-fst-pop WI.txt --out il_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop MO.txt --out ks_v_mo
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NC.txt --out ks_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NE.txt --out ks_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NH.txt --out ks_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NM.txt --out ks_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NV.txt --out ks_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop NY.txt --out ks_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop TX.txt --out ks_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop WA.txt --out ks_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop KS.txt --weir-fst-pop WI.txt --out ks_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NC.txt --out mo_v_nc
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NE.txt --out mo_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NH.txt --out mo_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NM.txt --out mo_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NV.txt --out mo_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop NY.txt --out mo_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop TX.txt --out mo_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop WA.txt --out mo_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop MO.txt --weir-fst-pop WI.txt --out mo_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop NE.txt --out nc_v_ne
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop NH.txt --out nc_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop NM.txt --out nc_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop NV.txt --out nc_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop NY.txt --out nc_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop TX.txt --out nc_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop WA.txt --out nc_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NC.txt --weir-fst-pop WI.txt --out nc_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop NH.txt --out ne_v_nh
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop NM.txt --out ne_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop NV.txt --out ne_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop NY.txt --out ne_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop TX.txt --out ne_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop WA.txt --out ne_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NE.txt --weir-fst-pop WI.txt --out ne_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop NM.txt --out nh_v_nm
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop NV.txt --out nh_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop NY.txt --out nh_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop TX.txt --out nh_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop WA.txt --out nh_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NH.txt --weir-fst-pop WI.txt --out nh_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NM.txt --weir-fst-pop NV.txt --out nm_v_nv
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NM.txt --weir-fst-pop NY.txt --out nm_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NM.txt --weir-fst-pop TX.txt --out nm_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NM.txt --weir-fst-pop WA.txt --out nm_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NM.txt --weir-fst-pop WI.txt --out nm_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NV.txt --weir-fst-pop NY.txt --out nv_v_ny
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NV.txt --weir-fst-pop TX.txt --out nv_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NV.txt --weir-fst-pop WA.txt --out nv_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NV.txt --weir-fst-pop WI.txt --out nv_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NY.txt --weir-fst-pop TX.txt --out ny_v_tx
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NY.txt --weir-fst-pop WA.txt --out ny_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop NY.txt --weir-fst-pop WI.txt --out ny_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop TX.txt --weir-fst-pop WA.txt --out tx_v_wa
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop TX.txt --weir-fst-pop WI.txt --out tx_v_wi
vcftools --vcf EUSTallSNPr8maf01p.vcf --weir-fst-pop WA.txt --weir-fst-pop WI.txt --out wa_v_wi


# remove symlinks 
find -type l -delete

### STRUCTURE

# use admixture model, correlated allele frequencies model due to expected shared ancestry
# at least 20 iterations
# create three .txt files to specify data (infile) and parameters (mainparams & extraparams)
#in directory with <infile>, mainparams.txt, extraparams.txt
structure >& structurelog &

# run to infer lambda first
# run with MAF filter
# then include estimated lambda in extraparams file
#line is exactly what is written below
#define LAMBDA <integer>

# zip all of _f output files into one .zip folder and upload to http://taylor0.biology.ucla.edu/structureHarvester/ 
# import into structure harvester

# use CLUMPP to identify clusters
export PATH=/programs/CLUMPP_Linux64.1.1.2:$PATH
CLUMPP paramfile

# use distruct 
export PATH=/programs/distruct1.1:$PATH
distructLinux1.1 drawparams

### fineRADstructure

# convert input file (rawSTACKSoutput_haplotypes.tsv) to correct format
# use Stacks2fineRAD.py script provided (need to edit first line for mac computers)
python Stacks2fineRAD.py -i EUSToneSNPr8maf01p.haplotypes.tsv -n 10 -m 30

mv EUSToneSNPr8maf01p.haplotypes.tsv.fineRADpainter.lociFilt.samples30%missFilt.txt EUSToneSNPr8maf01p.fineRADpainter.txt

# to prepare, run
export PATH=/programs/fineRADstructure/bin:$PATH

# run RADpainter to calculate coancestry matrix
RADpainter paint EUSToneSNPr8maf01p.fineRADpainter.txt >& RADpainterlog &
# assign individuals to populations
finestructure -x 100000 -y 100000 -z 1000 EUSToneSNPr8maf01p.fineRADpainter_chunks.out EUSToneSNPr8maf01p.fineRADpainter_chunks.mcmc.xml &
# build tree for visualization
finestructure -m T -x 10000 EUSToneSNPr8maf01p.fineRADpainter_chunks.out EUSToneSNPr8maf01p.fineRADpainter_chunks.mcmc.xml EUSToneSNPr8maf01p.fineRADpainter_chunks.mcmcTree.xml &

cp * /home/nrh44/EUST_RADseq_2018/fineRADstructure/EUSTallSNPr8maf01p
# plot results in R


### BayeScan 

# then run code below, prior odds 100 (at least)
/programs/BayeScan/BayeScan EUSTallSNPr8maf01p.bayescan -o EUSTallSNPr8maf01pbayescanLongRunOD10 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &

/programs/BayeScan/BayeScan EUSTallSNPr8.bayescan -o EUSTallSNPr8pbayescanLongRunOD10 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &


# make env.txt file with average of variable for each population
# note that order of populations is not alphabetical but by VCF 
# extracted WorldClim variables in R (see spatialstarling.R script)


/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bioPCA1 -env bioPCA1.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bioPCA1.chain2 -env bioPCA1.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &


# all variables
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio1mean -env bio1mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio2mean -env bio2mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio3mean -env bio3mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio4mean -env bio4mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio5mean -env bio5mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio6mean -env bio6mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio7mean -env bio7mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio8mean -env bio8mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio9mean -env bio9mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio10mean -env bio10mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio11mean -env bio11mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio12mean -env bio12mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio13mean -env bio13mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio14mean -env bio14mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio15mean -env bio15mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio16mean -env bio16mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio17mean -env bio17mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio18mean -env bio18mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio19mean -env bio19mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &

# running second chain to compare convergence
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio1mean.chain2 -env bio1mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio2mean.chain2 -env bio2mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio3mean.chain2 -env bio3mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio4mean.chain2 -env bio4mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio5mean.chain2 -env bio5mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio6mean.chain2 -env bio6mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio7mean.chain2 -env bio7mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio8mean.chain2 -env bio8mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio9mean.chain2 -env bio9mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio10mean.chain2 -env bio10mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio11mean.chain2 -env bio11mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio12mean.chain2 -env bio12mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio13mean.chain2 -env bio13mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio14mean.chain2 -env bio14mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio15mean.chain2 -env bio15mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio16mean.chain2 -env bio16mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio17mean.chain2 -env bio17mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio18mean.chain2 -env bio18mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &
/programs/bayescenv-1.1/bin/bayescenv EUSTallSNPr8.bayescan -o EUSTallSNPr8bayescenv.bio19mean.chain2 -env bio19mean.txt -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 &

### GPhoCS
bin/G-PhoCS-1-2-3 control-file-EUSTrad.ctl &> EUSTrad-GPhoCS1.log &

bin/readTrace mcmc.log -d1000 -b500 | awk '{print $5, $6, $7, $8, $9, $10, $11 }' &> readTrace.log

vcftools --vcf EUSTallSNPr8.vcf --012 --out EUSTallSNPr8
vcftools --vcf EUSTallSNPr8maf01p.vcf --012 --out EUSTallSNPr8maf01p


# Simon Martin scripts
# https://github.com/simonhmartin/genomics_general#diversity-and-divergence-analyses-in-sliding-windows, downloaded 9.11.18
python parseVCF.py -i EUSTallSNPr8.vcf --skipIndels --minQual 30 --gtf flag=DP min=5 | gzip > EUSTallSNPr8.geno.gz
# not actually phased, but in that format
# no outgroup provided. minor allele frequency used.
python sfs.py -g EUSTallSNPr8.geno.gz -f phased --dadiFormat --pref EUSTallSNPr8

# edit this to update below
python popgenWindows.py -w 50000 -m 5000 -g input.geno.gz -o output.csv.gz -f phased -T 5 -p popA A1,A2,A3,A4 -p popB B1,B2,B3,B4,B6,B6 -p popC -p popD --popsFile pops.txt 


# running stairway plot
cp -r /programs/stairway_plot/ /workdir/nrh44
cp -r /home/nrh44/EUST_RADseq_2018/stairway/* /workdir/nrh44/stairway_plot
cd stairway_plot/stairway_plot_es

unzip stairway_plot_v0.2.zip

export JAVA_HOME=/usr/local/jdk1.7.0_51
export PATH=$JAVA_HOME/bin:$PATH

PATH=/usr/local/jdk1.7.0/bin:$PATH
export PATH

# on personal computer
export CLASSPATH=~/Applications/stairway_plot_v2/stairway_plot_es
echo $CLASSPATH

java -cp stairway_plot_es Stairbuilder eustrad_fold.blueprint
bash eustrad_fold.blueprint.sh # in stairway_plot_v2beta2 directory!
bash eustrad_fold.blueprint.plot.sh


# see R for GEAs 
# identify genes near GEA candidates 
# bed files include variants under selection +/- 5000 bp up and downstream

# sed to replace genbank scaffold names with actual scaffolds
bedtools getfasta -fi EUSTref.fasta -bed LFMMcandidates.bed -fo LFMMcandidates.fa &
bedtools getfasta -fi EUSTref.fasta -bed RDAcandidates.bed -fo RDAcandidates.fa &

bedtools getfasta -fi EUSTref.fasta -bed LFMMcandidates.5000.bed -fo LFMMcandidates.5000.fa &
bedtools getfasta -fi EUSTref.fasta -bed RDAcandidates.5000.bed -fo RDAcandidates.5000.fa &

# blast fasta files with E = 0.001 and S.vulgaris records
# downloaded hit table to keep stats
# uploaded list of all genbank accessions to bioDBnet, output Gene Info and remove duplicates
# gives list of gene IDs plus description 
# just need to replace [ in txt file with tabs and then copy into hit table .csv




