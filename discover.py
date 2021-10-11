#Disccover v 1.1
import argparse
import subprocess
from pathlib import Path
import os
import glob

def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data

def getStxSubType():
    subprocess.call("pwd", shell=True)
    subprocess.call("blastn -query stx.fasta -db ../../discover/data/stx -task blastn -evalue 0.001 -out shigatoxin -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 1 -perc_identity 95.0", shell=True)
    file = open("shigatoxin")
    read = file.readlines()
    stx = open("shigatoxin_fc", "w")
    for c in read:
        line = c.split("\t")
        if float(line[2]) > 95 and float(line[3]) > 1200:
            stx.write(line[1] + "\t")
            stx.write(line[2] + "\t")
            stx.write(line[3] + "\t")
            stx.write(line[15] + "\t")
            stx.write("\n")
    stx.close()
    file.close()
    shigatoxin_typing = openFileAsTable("shigatoxin_fc")
    if len(shigatoxin_typing) == 0:
        str_shigatoxin_subtype = "No subtype match found"
    else:
        # get corresponding subtypes
        str_shigatoxin_subtype = ""
        shigatoxin_subtypes = []
        shigatoxin_subtypes_raw = []
        shigatoxin_types = openFileAsTable("../../discover/data/stx_subtypes")
        for subtype in shigatoxin_typing:
            blast_pident_100 = float(subtype[1]) == 100
            if (blast_pident_100):
                for item in shigatoxin_types:
                    if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                        shigatoxin_subtypes.append(item[1])
                        shigatoxin_subtypes_raw.append(item[1])
       # partial matches
        for subtype in shigatoxin_typing:
            for item in shigatoxin_types:
                if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                    if item[1][0:4] == "stx1":
                        shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                        shigatoxin_subtypes_raw.append(item[1])
                    if item[1][0:4] == "stx2":
                        shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                        shigatoxin_subtypes_raw.append(item[1])
        shigatoxin_subtypes.sort()
        str_shigatoxin_subtype = " ".join(shigatoxin_subtypes)
    return str_shigatoxin_subtype

def getSeroGroup(sero_file, antigen):
    subprocess.call("awk -F '\t' '$4>800 { print $2 FS $3 FS $4 FS $16 }' " + sero_file + " | sort -nrk 2 -nrk 3 > serogroup_fc", shell=True)
    subprocess.call("awk -F , '!seen[$0]++' serogroup_fc > serogroup_fcd", shell=True)
    sero_typing = openFileAsTable("serogroup_fcd")
    if len(sero_typing) == 0:
        serogroup = antigen + "?"
    else:
        serogroup = sero_typing[0][0][sero_typing[0][0].rfind(antigen):]
    return serogroup

parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='directory', help='relative path to directory with FASTQ data')
parser.add_argument('-i', dest='input', help='type of FASTQ data: ion (Ion Torrent), illumina (Illumina)')
args = parser.parse_args()

os.chdir('./'+args.directory)
os.system("ln -s " + os.popen("which trimmomatic").read().strip() + " trimmomatic.jar")
subprocess.call("mkdir Contigs/", shell=True)

if args.input =='illumina':
    for file in glob.glob("*_R1*"):
        sample=file.split('_')[0]
        p1_file=file
        p2_file=file.replace('_R1', '_R2')
        find_p2=os.listdir('./')
        if p2_file in find_p2:
            # TRIMMING
            subprocess.call(
            "trimmomatic PE -phred33 " + p1_file + " " + p2_file + " trimmed_" + sample + "_1.fq " + sample+"_1unpaired trimmed_" + sample + "_2.fq " + sample+"_2unpaired"+" SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36", shell=True)
            # ASSEMBLY
            subprocess.call(
            "perl " + "../discover/scripts/spades.pl "+sample+"_spades_contigs "+sample+"_spades_contig_stats "+sample+"_spades_scaffolds "+sample+"_spades_scaffold_stats "+sample+"_spades_log NODE spades.py --disable-gzip-output --isolate --pe1-ff --pe1-1"+ " trimmed_" + sample + "_1.fq "+ "--pe1-2 trimmed_" + sample + "_2.fq",
            shell=True)
            subprocess.call(
            "perl " + "../discover/scripts/filter_spades_repeats.pl -i "+sample+"_spades_contigs -t "+sample+"_spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o "+sample+"_output_with_repeats -u "+sample+"_contigs -n "+sample+"_repeat_sequences_only -e 5000 -f "+sample+"_discarded_sequences -s "+sample+"_summary",
            shell=True)
            if len(sample + "_contigs") != 0:
                Path(sample + "_disc").mkdir(exist_ok=True)
                subprocess.call("mv " + "trimmed_" + sample + "_1.fq " + sample + "_disc/", shell=True)
                subprocess.call("mv " + "trimmed_" + sample + "_2.fq " + sample + "_disc/", shell=True)
                subprocess.call("cp " + sample + "_contigs Contigs/", shell=True)
                subprocess.call("mv " + sample + "_contigs " + sample + "_disc/", shell=True)
                subprocess.call("rm -rf output_dir", shell=True)
        else:
            log=open(sample+'_log.txt', 'w')
            log.write(sample+'_R2: NOT FOUND')

    subprocess.call("rm *spades*", shell=True)
    subprocess.call("rm *unpaired", shell=True)
    subprocess.call("rm *summary", shell=True)
    subprocess.call("rm *sequences", shell=True)
    subprocess.call("rm *sequences_only", shell=True)
    subprocess.call("rm *repeats", shell=True)

    for file in glob.glob("*_disc/"):
        sample = file.split('_')[0]
        os.chdir(file)

    # VIRULOTYPER
        subprocess.call("abricate --db virulotyper "+sample + "_contigs > virulotyper_results.txt", shell=True)

    # SHIGATOXIN TYPER
        subprocess.call(
            "../../discover/scripts/stx_subtype_pe.sh ../.. " + "trimmed_" + sample + "_1.fq " + "trimmed_" + sample + "_2.fq " + sample + "_contigs",
            shell=True)
        results_stx = []
        results_stx.append(getStxSubType())
        stx = open("shigatoxin_results.txt", "w")
        stx.write("SHIGATOXIN TYPER: " + "\n")
        [stx.write(c + "\n") for c in results_stx]
        stx.close()

    # SEROTYPER O&H
        subprocess.call("../../discover/scripts/serotype.sh ../.. y "+"trimmed_"+sample+"_1.fq "+ "trimmed_"+sample+"_2.fq "+sample+"_contigs",
        shell=True)
        results_o = []
        results_h = []
        results_o.append(getSeroGroup('serogroup_O', "O"))
        results_h.append(getSeroGroup('serogroup_H', "H"))
        serotyper = open("serotyper_results.txt", "w")
        serotyper.write("SEROTYPER O&H " + "\n")
        serotyper.write("SEROTYPER O: ")
        [serotyper.write(c + ", ") for c in results_o]
        serotyper.write("\n")
        serotyper.write("SEROTYPER H: ")
        [serotyper.write(c + ", ") for c in results_h]
        serotyper.close()

    # MLST
        subprocess.call("mlst --scheme ecoli "+sample+"_contigs > mlst_results.txt", shell=True)

    # AMRGENES, abrichate
        subprocess.call("abricate -db ncbi "+sample+"_contigs > amr_abrichate_results.txt", shell=True)

    # CHEwBBACCA
        chewbbca=os.popen("which chewBBACA.py").read()
        subprocess.call("mkdir chewbbca", shell=True)
        subprocess.call("cp "+sample+"_contigs ./chewbbca/"+sample+"_contigs.fasta",shell=True)
        subprocess.call("python "+chewbbca+"chewBBACA.py AlleleCall -i ./chewbbca/ -g  ../../discover/chewBBACA_db/schema_chewBBACA_cgMLST_V4 -o chewbacca_results/ --cpu 6 --bsr 0.6 --ptf ../../discover/chewBBACA_db/trained_eColi.trn --fr", shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_statistics.tsv ./chewbbca_results.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_alleles.tsv ./results_alleles.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/RepeatedLoci.txt ./RepeatedLoci.txt", shell=True)
        os.chdir('../')

if args.input=='ion':
    list_ngs_data=glob.glob("*fastq")+glob.glob("*fastqsanger")
    for file in list_ngs_data:
        sample=file.split('_')[0]
        # TRIMMING
        subprocess.call("trimmomatic SE -phred33 " + file + " trimmed_" + sample + ".fq SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55",
        shell=True)
        # ASSEMBLY
        subprocess.call(
        "perl " + "../discover/scripts/spades.pl "+sample+"_spades_contigs "+sample+"_spades_contig_stats "+sample+"_spades_scaffolds "+sample+"_spades_scaffold_stats "+sample+"_spades_log NODE spades.py --disable-gzip-output --isolate --iontorrent -s trimmed_" + sample + ".fq",
        shell=True)
        subprocess.call(
        "perl " + "../discover/scripts/filter_spades_repeats.pl -i "+sample+"_spades_contigs -t "+sample+"_spades_contig_stats -c 0.33 -r 1.75 -l 1000 -o "+sample+"_output_with_repeats -u "+sample+"_contigs -n "+sample+"_repeat_sequences_only -e 5000 -f "+sample+"_discarded_sequences -s "+sample+"_summary",
        shell=True)
        if len(sample + "_contigs") != 0:
            Path(sample + "_disc").mkdir(exist_ok=True)
            subprocess.call("mv " + " trimmed_" + sample + ".fq " + sample + "_disc/", shell=True)
            subprocess.call("cp " + sample + "_contigs Contigs/", shell=True)
            subprocess.call("mv " + sample + "_contigs " + sample + "_disc/", shell=True)
            subprocess.call("rm -rf output_dir", shell=True)
    
    subprocess.call("rm *spades*", shell=True)
    subprocess.call("rm *summary", shell=True)
    subprocess.call("rm *sequences", shell=True)
    subprocess.call("rm *sequences_only", shell=True)
    subprocess.call("rm *repeats", shell=True)

    for file in glob.glob("*_disc/"):
        sample = file.split('_')[0]
        os.chdir(file)

    # VIRULOTYPER
        subprocess.call("abricate --db virulotyper "+sample + "_contigs > virulotyper_results.txt",
        shell=True)

    # SHIGATOXIN TYPER
        subprocess.call("../../discover/scripts/stx_subtype_se.sh ../.. "+"trimmed_"+sample+".fq "+sample+"_contigs",
        shell=True)
        results_stx = []
        results_stx.append(getStxSubType())
        stx=open("shigatoxin_results.txt", "w")
        stx.write("SHIGATOXIN TYPER: "+"\n")
        [stx.write(c+"\n") for c in results_stx]
        stx.close()

    # SEROTYPER O&H
        subprocess.call("../../discover/scripts/serotype.sh ../.. n trimmed_"+sample+".fq xxx "+sample+"_contigs",
        shell=True)
        results_o = []
        results_h = []
        results_o.append(getSeroGroup('serogroup_O', "O"))
        results_h.append(getSeroGroup('serogroup_H', "H"))
        serotyper = open("serotyper_results.txt", "w")
        serotyper.write("SEROTYPER O&H " + "\n")
        serotyper.write("SEROTYPER O: ")
        [serotyper.write(c + ", ") for c in results_o]
        serotyper.write("\n")
        serotyper.write("SEROTYPER H: ")
        [serotyper.write(c + ", ") for c in results_h]
        serotyper.close()

    # MLST
        subprocess.call("mlst --scheme ecoli "+sample+"_contigs > mlst_results.txt", shell=True)

    # AMRGENES, abrichate
        subprocess.call("abricate -db ncbi "+sample+"_contigs > amr_abrichate_results.txt", shell=True)

    #CHEwBBACCA
        chewbbca=os.popen("which chewBBACA.py").read()
        subprocess.call("mkdir chewbbca", shell=True)
        subprocess.call("cp "+sample+"_contigs ./chewbbca/"+sample+"_contigs.fasta",shell=True)
        subprocess.call("python "+chewbbca+"chewBBACA.py AlleleCall -i ./chewbbca/ -g  ../../discover/chewBBACA_db/schema_chewBBACA_cgMLST_V4 -o chewbacca_results/ --cpu 6 --bsr 0.6 --ptf ../../discover/chewBBACA_db/trained_eColi.trn --fr", shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_statistics.tsv ./chewbbca_results.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/results_alleles.tsv ./results_alleles.tsv", shell=True)
        subprocess.call("cp chewbacca_results/results_*/RepeatedLoci.txt ./RepeatedLoci.txt", shell=True)
        os.chdir('../')

subprocess.call("python ../discover_results.py", shell=True)
