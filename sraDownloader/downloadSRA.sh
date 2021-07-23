### IF SRA Toolkit is installed skip installation part start from export PATH
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2019/RNAseq/Supplementary_Materials/S1_Getting_raw_reads_from_SRA.html
# download the gzip file
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz
# unzip the file
tar -xzvf sratoolkit.2.10.9-ubuntu64.tar.gz
# add the 'bin' directory to the PATH - note the you will need to do this
# everytime you start a new terminal and wish to use the toolkit
export PATH=/data/apps/sratoolkit.2.10.9-ubuntu64/bin/:${PATH}
# create a directory to which to download the sra files
mkdir sra
# use the vdb-config tool to set the download directory
vdb-config -i 
# this pops up an interactive window instructions below
#Use the vdb-config window to set the Default Import Path to the new sra directory we just created. Use tab to navigate to Change under the Set  Default Import Path (the highlighting indicates the active field the arrow keys and tab to navigate to the correct directory. When you have changed the directory Save and Exit

#### Download SRA files
#https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP056956&o=acc_s%3Aa
prefetch SRR12458584 
prefetch SRR12458585 
prefetch SRR12458586 
prefetch SRR12458587 
prefetch SRR12458588 
prefetch SRR12458589 
prefetch SRR12458590 
prefetch SRR12458591 
prefetch SRR12458592 
prefetch SRR12458593 

#Extracting sra files
mkdir fastq
for sraFile in SRR*/*.sra; do
  echo "Extracting fastq from "${sraFile}
  fastq-dump \
     --origfmt \
     --gzip \
     --outdir fastq \
     ${sraFile}
 done
