#change the "" names below here and then save - close and doible clikc
#dont forget has to be in the same folder as the two input files.

#INPUT fasta file of promoter regions of interest
#e.g.
#>Minc3s00001g00001
#GACAATAAAAAGTCGGGTTGCTTTAATTTTCATATTTTCTGATTCAGCTTGAAGAGACGCATCGAATGATAACAAAATTAAAAATATCTCACTTTCCTTTGTTGAATTAGAATTTTAAATAAGTTAGGGTTAATCACACTTGAAAACGGTTGGAAAAACAAAAAATAAAAAATATTTTTTTAAATTCACAATTATACATTTTTAAACTTTTTAGTTTATAAACTCATGTTAATTAAATATTTATTTATTTTTGAATGCCCATAAATAAAATTTCAAACGTGACAATTTTTTATTTTTAAATTTTGCCTCTTTTCTTCTTTTAATTTTATTCGTATCCTCGCGTTTTTTATTCATTAATTAAATTTTAAAAAATTAAAAAAATTTTTTTTATTTCAGATCCAAAGTCGTTACCCACAACCTCTCCTTACACGCGCATAATTTTCATCTTTTTGTCAACAAATCTGAAGCCCATTAGAAGGCAAAAAGTATTAGTGTTGACGACAGTGAGTAGGAGTAGGAAGACATTGAAAATCAGCAGCAGCCGAGTGGTAGAGTGAAGAAAAAAGTGTTTGGAAATTGAAACGGGTAAAAAAGGATTAAAACGAGTAGAAACTACTTAAAACGAGTAGTTGACAACAAATTCAAGCAGTAATTGAAAAAAAATCATCTCTAGTAGAAGAACAACTAGGAGAATAAACTCGAGTAGCTTTTTAGAACAACTAGAAAAAACATTGTTTAACCAAGCTACAGGCAACAACATTTTTTAATTAAAAATATTTTTTGTGTTTATTGTAAAAAATATTTTTCTACCAAAAAAAAGTGCTAAACAACTGACTATATTTTGAAAAAAAATTTTTTTGGACTTTGAAGGACACTGAAAAAAAAATTTTTTTAATTGCAAAAACTCGTTTTTTCCTTCTATCCACGTGAAGAAGTAGGAAAAAGGACCGAAAGCAATAACAACTTTTTAAAGAAGAAAAAAAATAGAAGAAGAAAAAAA
#etc...
name_of_in_fasta_file = "MINC.all.plus_effectors_1000_upstream.fasta.awk"

#INPUT list of names for only the promoters that correspond to secreted protein prediction.
#e.g.
#Minc3s00001g00003
#Minc3s00001g00024
#etc...

opened_sec_string = open("Minc3-SP-no-TM.acc")
read_sec_string = opened_sec_string.read()

#INOUT here motif of interest (TGCACTT in our case)
input_list = ['T','G','C','A','C','T','T']


#output
#will output three tables: % of sequences with given motif that are secreted, number secreted with given motif, total number with given motif
#each table contains 6 columns - the first contains the motif sequence, the following 5 show the results for 1-5 copies.
#
#e.g.
#TGCACTT	15	26	28	11	0
#GGCACTT .. .etc


options = ['A','T','G','C']
base = ""
string = ""
list_of_options = ""
for y in range(len(options)):
    for x in range(len(input_list)):
        for z in range(len(input_list)):
            if x is z:
                base = options[y]
            else:
                base = input_list[z]
            string = string + base
        #print (string)
        list_of_options = list_of_options + string + "\n"
        string = ""
        
        
list_of_variants = list_of_options.split("\n")

print(list_of_variants)
myset = set(list_of_variants)
print (myset)

opened_sec_string = open("Minc3-SP-no-TM.acc")
read_sec_string = opened_sec_string.read()

opened = open(name_of_in_fasta_file)
readed= opened.read()
splited = readed.split(">")
out = ""
sec_pct = ""
sec_number = ""
tot_num = ""
for variant in myset:
    if len(variant)>0:
        #print (variant)
        out = ""
        for line in splited:
            #print (line)
            if len (line)>0:
                linesplit = line.split("\n")
                name = linesplit[0]
                entry = linesplit[1]
                
                count_plus = entry.count(variant)

                complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                reverse_complement = "".join(complement.get(base, base) for base in reversed(variant))
                #print (reverse_complement)
                count_minus = entry.count(reverse_complement)
                count_total = count_plus+ count_minus
                out = out + name + "\t" + str(count_total) + "\n"
        outfile = open(name_of_in_fasta_file + "_" + variant + ".counts","w")
        outfile.write(out)
        outfile.close()
        #print (variant + " done")
        count_total = 0
        test_parse = out
        out = 0
        listx = ""
        listy = ""
        listz = ""
        for x in range(5):

            count_all_genes = 0
            count_all_genes_secreted= 0
            counter_genes = 0
            counter_secreted = 0
            split_test = test_parse.split("\n")
            for line in split_test:
                if len(line)>0:
                    #print (line)
                    linesplit = line.split("\t")
                    #print (linesplit[1])
                    count_all_genes = count_all_genes + 1
                    if name in read_sec_string:
                        count_all_genes_secreted = count_all_genes_secreted + 1
                    name = linesplit[0]
                    number_of_boxes = linesplit[1]
                    if int(number_of_boxes)>x:
                        #print (str(x) + "_" + line)
                        counter_genes = counter_genes + 1
                        
                        if name in read_sec_string:
                            counter_secreted = counter_secreted + 1
                        #print (str(x) + "_" + line + "secreted!")
            listx = listx + str(counter_genes) + "\t"
            listy = listy + str(counter_secreted) + "\t"
            do_nothing = 0
            if counter_genes >0 and counter_secreted>0:
                listz = listz + str(int((float(counter_secreted)/float(counter_genes))*100))+ "\t"
            else:
                listz = listz + "0\t"    
            #print ("number of boxes = " + str(x) + " #genes = " + str(count_all_genes) + " #secreted = " +str(count_all_genes_secreted)+  " #genes with "+str(x)+" or more boxes = " + str(counter_genes) + " number secreted = " + str(counter_secreted))
        #print (" #genes with motif =\t" +listx)
        #print (" # of those secreted =\t" +listy)
        #print (variant+"\t" +listz)#sec pct
        #print (variant+"\t" +listy)#secreted number
        print (variant+"\t" +listz+"\n")#sec pct
        #print (variant+"\t" +listx)#num genes
        sec_pct = sec_pct + variant+"\t" +listz+"\n" 
        sec_number = sec_number +variant+"\t" +listy + "\n"
        tot_num = tot_num + variant+"\t" +listx + "\n"

print(sec_pct)
print(sec_number)
print(tot_num)

opened.close()
        
