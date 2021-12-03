#Disccover v 1.1
import os
import glob
import csv
import operator
from numpy import mean

#Tab_results.txt and  Tab_AMR.txt
list_gene_fine=['asta', 'ipad', 'saa', 'agg3a', 'cif', 'espj', 'toxb', 'aggd', 'lpfa', 'cfac', 'agg3d', 'katp', 'espa', 'ltca', 'cnf1', 'espc', 'espb', 'aata', 'agg4c', 'tccp', 'stx1a', 'ipah9', 'lnga', 'tir', 'cdtb', 'iss', 'espf', 'tsh', 'espp', 'virf', 'agg5a', 'sfas', 'k88ab', 'bfpa', 'efa1', 'capu', 'agga', 'eila', 'pic', 'cma', 'vat', 'agg4d', 'agg4a', 'cci', 'siga', 'aar', 'stx2a', 'ehxa', 'mchb', 'senb', 'aafd', 'iha', 'aggc', 'stx2b', 'fim41a', 'etpd', 'fasa', 'mcma', 'nleb', 'aafa', 'nlea', 'sat', 'feda', 'nfae', 'iron', 'suba', 'pet', 'pera', 'aafc', 'air', 'agg3c', 'aaic', 'eae', 'epea', 'gad', 'sepa', 'rpea', 'hlye', 'fedf', 'aap', 'orf4', 'aggr', 'f17a', 'nlec', 'orf3', 'cba', 'fana', 'stb', 'aafb', 'mchc', 'eata', 'sta1', 'agg3b', 'mchf', 'aggb', 'stx1b', 'f17g', 'celb', 'cofa', 'espi', 'agg4b', 'irea']

entries = os.listdir('./')
if 'Tab_results.txt' not in entries:
    tab1=open('Tab_results.txt', 'w')
    prima_riga = 'Sample' + '\t' + 'Avg Scaffold coverage' + '\t' + 'Burst size' + '\t' + 'MLST' + '\t' + 'STX subtype' + '\t' + 'Serotype O' + '\t'+ 'Serotype H' + '\t'
    for c in list_gene_fine:
        prima_riga += c
        prima_riga += '\t'
    prima_riga+='\n'
    tab1.write(prima_riga)
else:
    tab1 = open('Tab_results.txt', 'a')

if 'Tab_AMR.txt' not in entries:
    tab2 = open('Tab_AMR.txt', 'w')
    prima_riga = 'Sample' + '\t' + 'AMR genes'
    tab2.write(prima_riga)
else:
    tab2 = open('Tab_AMR.txt', 'a')

if 'Tab_cgMLST.txt' not in entries:
    tab3 = open('Tab_cgMLST.txt', 'w')
    for file in glob.glob("*_disc/"):
        prima_riga = ''
        os.chdir(file)
        prima_riga = 'Sample' + '\t'

        # extract info for MLST gene results
        csv_file = open("results_alleles.tsv")
        read_csv = list(csv.reader(csv_file, delimiter="\t"))

        for gene in read_csv[0][1:]:
            prima_riga += gene + '\t'
        prima_riga += "\n"

        tab3.write(prima_riga)
        os.chdir("../")
        break
else:
    tab3 = open('Tab_cgMLST.txt', 'a')

if 'Contamination_sheet.txt' not in entries:
    tab4 = open('Contamination_sheet.txt', 'w')
    prima_riga = 'Sample' + '\t' + 'Species'+ '\t' + 'Gene' + '\t' + 'PC' + '\t' + 'NDC' + '\n'
    tab4.write(prima_riga)
else:
    tab4 = open('Contamination_sheet.txt', 'a')


for file in glob.glob("*_disc/"):
    namesample = file.split('_')[0]
    os.chdir(file)
    row_sample=namesample+'\t'

    if os.path.isfile("./virulotyper_results.txt")==True and len("./virulotyper_results.txt")>1:
        #VIRULOTYPER RESULTS
        #read the file virulotyper_results.txt
        file1 = open('virulotyper_results.txt')
        read_csv = list(csv.reader(file1, delimiter="\t"))

        #list of genes with coverage >80.0
        geni_ED=[]
        if read_csv[0][4]!="STRAND":
        	for line in read_csv[1:]:
        		if float(line[8])>=80.0:
        			in_line=[]
        			in_line.append(line[12])
        			in_line.append(float(line[8]))
        			in_line.append(float(line[9]))
        			in_line.append(float(line[1].split('_')[5]))
        			geni_ED.append(in_line)
        else:
        	for line in read_csv[1:]:
        		if float(line[9])>=80.0:
        			in_line=[]
        			in_line.append(line[13])
        			in_line.append(float(line[9]))
        			in_line.append(float(line[10]))
        			in_line.append(float(line[1].split('_')[5]))
        			geni_ED.append(in_line)

        #order genes by coverage, identity and read mean coverage
        if geni_ED:
            ED_best_gene=[]
            tot_mean=[]
            for geni in list_gene_fine:
                list_gene=[]
                for line in geni_ED:
                    if line[0].split('_')[0]==geni:
                        list_gene.append(line)
                        tot_mean.append(line[3])
                s=sorted(list_gene, key=operator.itemgetter(1, 2, 3), reverse=True)
                if len(s)!=0:
                    ED_best_gene.append(s[0])

            #coverage mean
            mean_gene=str(round(mean(tot_mean)))

            #empty list from list_gene_fine
            find_gene=['0']*len(list_gene_fine)
            for x in ED_best_gene:
                gen = x[0].split('_')[0]
                for c in list_gene_fine:
                    if gen==c:
                        find_gene[list_gene_fine.index(c)]=x[0]
            file1.close()
        else:
            mean_gene = "ND"
            find_gene = ['ND'] * len(list_gene_fine)

    else:
        mean_gene = "ND"
        find_gene = ['ND'] * len(list_gene_fine)

    if os.path.isfile("./mlst_results.txt") == True:
        #MLST
        mlst=''
        species=''
        file2=open('mlst_results.txt').readlines()
        for line in file2:
            mlst+='ST'+line.split()[2]+'; '
            species+=line.split()[1]+'; '
    else:
        mlst="ND"

    if os.path.isfile("./chewbbca_results.tsv") == True and len("./chewbbca_results.tsv")>1:
        #extract EXC+INF from chewBBACA
        file3=open('chewbbca_results.tsv')
        file3_line=file3.readlines()
        exc=str(int(file3_line[1].split()[1])+int(file3_line[1].split()[2]))
        file3.close()
    else:
        exc="ND"

    if os.path.isfile("./shigatoxin_results.txt") == True and len("./shigatoxin_results.txt")>1:
        # SHIGATOXIN TYPER
        file4 = open("shigatoxin_results.txt")
        file4_lines=file4.readlines()
        list_of_stx = ''
        if file4_lines[1].find("No subtype match found")==-1:
            for line in file4_lines[1:]:
                for c in line.strip('\n').split(' '):
                    list_of_stx += c[0:4] + c[4].lower() + '; '
        else:
            list_of_stx+="ND  "
        file4.close()
    else:
        list_of_stx = "ND  "

    if os.path.isfile("./serotyper_results.txt") == True and len("./serotyper_results.txt")>1:
        # SEROTYPER- STX O&H
        file5 = open("serotyper_results.txt")
        file5_lines=file5.readlines()
        sero_o = file5_lines[1].strip('\n').split(':')[1:][0].replace(',', ';')
        sero_h = file5_lines[2].strip('\n').split(':')[1:][0].replace(',', ';')
        file5.close()
    else:
        sero_o ="ND  "
        sero_h ="ND  "

    #add info to sample row
    row_sample=namesample+'\t'
    row_sample+=mean_gene+'\t'
    row_sample+=exc+'\t'
    row_sample+=mlst[:-2]+'\t'
    row_sample += list_of_stx[:-2] + '\t'
    row_sample += sero_o[:-2] + '\t'
    row_sample +=sero_h[:-2] + '\t'

    for gene in find_gene:
        row_sample += str(gene) + '\t'

    row_sample+='\n'
    tab1.write(row_sample)

    if os.path.isfile("./amr_abrichate_results.txt") == True:
        #extract info for AMR results
        csv_file = open("amr_abrichate_results.txt")
        read_csv = list(csv.reader(csv_file, delimiter="\t"))
        riga2 = namesample + ' :' + '\t'

        #list of genes with coverage >80.0
        geni_amr = []
        list_geni_amr=[]
        if read_csv[0][4]!="STRAND":
        	for line in read_csv[1:]:
        		if float(line[8]) >= 80.0:
        			in_line = []
        			in_line.append(line[4])
        			in_line.append(float(line[8]))
        			in_line.append(float(line[9]))
        			in_line.append(float(line[1].split('_')[5]))
        			geni_amr.append(in_line)
        			list_geni_amr.append(line[4])
        else:
        	for line in read_csv[1:]:
        		if float(line[9]) >= 80.0:
        			in_line = []
        			in_line.append(line[5])
        			in_line.append(float(line[9]))
        			in_line.append(float(line[10]))
        			in_line.append(float(line[1].split('_')[5]))
        			geni_amr.append(in_line)
        			list_geni_amr.append(line[5])

        list_geni_amr=list(set(list_geni_amr))

        # order genes by coverage, identity and read mean covrage
        ED_best_amr = []
        for geni in list_geni_amr:
            list_gene = []
            for line in geni_amr:
                if line[0] == geni:
                    list_gene.append(line)

            s = sorted(list_gene, key=operator.itemgetter(1, 2, 3), reverse=True)
            if len(s) != 0:
                ED_best_amr.append(s[0])

        for x in ED_best_amr:
            riga2 += x[0] + "; "

    else:
        riga2=namesample + ' :' + '\tND'+"\n"

    csv_file.close()
    tab2.write('\n' + riga2)

    if os.path.isfile("./results_alleles.tsv") == True:
        # extract info for MLST gene results
        row_sample = namesample + '\t'
        csv_file = open("results_alleles.tsv")
        read_csv = list(csv.reader(csv_file, delimiter="\t"))

        for value in read_csv[1][1:]:
            row_sample += value + '\t'
        row_sample += "\n"
        csv_file.close()
        tab3.write(row_sample)
    else:
        row_sample = namesample + '\tND'+"\n"
        tab3.write(row_sample)

    if os.path.isfile("./RepeatedLoci.txt") == True:
        #MLST
        file4=open('RepeatedLoci.txt').readlines()
        row_sample = namesample + '\t' + species[:-2] + '\t'
        if len(file4)==1:
            row_sample+= '\t' + '\t' + '\t' + '\n'
            tab4.write(row_sample)
        else:
            list_loci = []
            for line in file4[1:]:
                repetition = line.split('\t')
                list_loci.append(repetition[0] + '\t' + repetition[1] + '\t' + repetition[2].replace('\n','')+ '\n')
            for rep in list_loci:
                locus=row_sample+rep
                tab4.write(locus)
    else:
	row_sample+= '\t' + '\t' + '\t' + '\n'
	tab4.write(row_sample)
    os.chdir('../')

tab1.close()
tab2.close()
tab3.close()
tab4.close()



