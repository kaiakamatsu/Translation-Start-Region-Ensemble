from Bio import SeqIO
import pandas as pd

#reads a text file of gene ids to check and returns a list of those ids (remember that each item has a "\n" at the end)
def id_tocheck():
    with open('check') as file:
        lines = file.read().split("\n")
    lines.pop(-1)
    return lines

#read dna file and create a list of nucleotide sequences for every gene
def read_dna() :
    filename = 'Mus_musculus.GRCm39.cdna.all.fa'
    dnasequences = list(SeqIO.parse(filename, 'fasta'))
    return dnasequences

#read protein file and create a list of amino acid sequences for every gene
def read_protein() :
    filename = 'Mus_musculus.GRCm39.pep.all.fa'
    proteinsequences = list(SeqIO.parse(filename, 'fasta'))
    return proteinsequences

#takes the element in the list and returns the gene id
def gene_id(element) :
    start = str(element.description).index("gene:")
    end = str(element.description).index("gene_biotype")
    geneid = str(element.description)[start : end]
    return geneid

#takes the element in the list and returns the protein id
def protein_id(element) :
    start = str(element.description).index("gene:")
    end = str(element.description).index("transcript")
    proteinid = str(element.description)[start: end]
    return proteinid

#takes the element from the pep file and returns just the transcript id
def protein_transcript(element):
    start = str(element.description).index("transcript:")
    end = str(element.description).index("gene_biotype:")
    proteintranscript = str(element.description)[start + 11 :end]
    return proteintranscript

#takes the element from the cdna file and returns just the transcript id
def cdna_transcript(element):
    #start = str(element.description).index("transcript:")
    end = str(element.description).index("cdna")
    cdnatranscript = str(element.description)[0:end]
    return cdnatranscript

#takes an element in the dna list and the entire protein list and finds the index of a protein that matches the genes id
def find_id (dnaelement, proteinlist) :
    gene = gene_id(dnaelement)
    index = 0
    for i in proteinlist :
        protein = str(protein_id(i))
        if gene == protein :
            return index
        index = index+1
    return -1

#takes the entire list of proteins and finds the one with the same transcript id as the given dna
def find_transcript (dnaelement, proteinlist) :
    transcript = cdna_transcript(dnaelement)
    index = 0
    for i in proteinlist:
        protein = str(protein_transcript(i))
        if transcript == protein:
            return index
        index = index + 1
    return -1

#takes a gene id and the entire cdna file and returns a list of possible transcript ids of that gene, along with their corresponding sequences
def list_transcript (gene, cdna):
    transcripts = []
    sequences = []
    for i in cdna:
        if (gene == gene_id(i)[5: gene_id(i).index(".")]):
            transcripts.append(cdna_transcript(i))
            sequences.append(i)
    return transcripts, sequences

#find the first occurrence of ATG in the cdna and reads codons starting from that start position, matches this amino acid sequence to the given peptide chain, if there is no match, retry with a later occurrence of ATG
#no ATG, return -1, otherwise return the index of the start codon
def find_ATG(dnasequence, peptidechain) :
    index = 0
    while index < len(dnasequence):
        dna = str(dnasequence)
        index = dna.find("ATG", index, len(dnasequence))
        if index != -1:
            dnasequence2 = dnasequence[index:len(dnasequence)]
            dnaaminoacid = dnasequence2.translate()
            check = match(dnaaminoacid, peptidechain)
            if check > -1 :
                return index
            index = index + 3
        else:
            index = len(dnasequence)
    return -1

#takes an amino acid letter and returns a list of codons that could be encoding for that amino acid
def aa_to_codon(letter) :
    if letter == "F":
        return ["TTT", "TTC"]
    elif letter == "L":
        return ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
    elif letter == "S":
        return ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
    elif letter == "Y":
        return ["TAT", "TAC"]
    elif letter == "X":
        return ["TAA", "TAG", "TGA"]
    elif letter == "C":
        return ["TGT", "TGC"]
    elif letter == "W":
        return ["TGG"]
    elif letter == "P":
        return ["CCT", "CCC", "CCA", "CCG"]
    elif letter == "H":
        return ["CAT", "CAC"]
    elif letter == "Q":
        return ["CAA", "CAG"]
    elif letter == "R":
        return ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]
    elif letter == "I":
        return ["ATT", "ATC", "ATA"]
    elif letter == "M":
        return ["ATG"]
    elif letter == "T":
        return ["ACT", "ACC", "ACA", "ACG"]
    elif letter == "N":
        return ["AAT", "AAC"]
    elif letter == "K":
        return ["AAA", "AAG"]
    elif letter == "V":
        return ["GTT", "GTC", "GTA", "GTG"]
    elif letter == "A":
        return ["GCT", "GCC", "GCA", "GCG"]
    elif letter == "D":
        return ["GAT", "GAC"]
    elif letter == "E":
        return ["GAA", "GAG"]
    elif letter == "G":
        return ["GGG", "GGT", "GGC", "GGA"]
    else:
        return ""

#takes the first amino acid letter - for each possible codon, find the first occurrence of it in the cdna and reads codons starting from that position, matches the peptide chain to amino acid sequence, if there is no match try again with a later occurence
#return -1 if there is no match at all, otherwise return the index of the start codon
def find_other_start(aa, dnasequence, peptidechain) :
    codons = aa_to_codon(aa)
    for i in codons:
        index = 0
        while index < len(dnasequence):
            dna = str(dnasequence)
            index = dna.find(i, index, len(dnasequence))
            if index != -1:
                dnasequence2 = dnasequence[index:len(dnasequence)]
                dnaaminoacid = dnasequence2.translate()
                check = match(dnaaminoacid, peptidechain)
                if check > -1:
                    return index
                index = index + 3
            else:
                index = len(dnasequence)
    return -1

#match amino acid sequences
def match(dnaaminoacids, peptidechain) :
    pep = str(peptidechain)
    dna = str(dnaaminoacids)
    return dna.find(pep)

#take the cdna sequence, the id, transcript and the start codon index, determines if it is weak or strong
def strong_or_weak (dnasequence, id, transcript, index) :
    infolist = []
    infolist.append(id[5:len(id)])
    infolist.append(transcript)
    print(transcript)
    if index == -1:
        infolist.append("")
        infolist.append("False Match")
    else:
        if len(dnasequence[0:index]) < 6:
            infolist.append(str(dnasequence[0:index+4]))
            infolist.append("Cannot Determine")
        else:
            kozak = "GCCACC"
            kozak2 = "GCCGCC"
            kozak3 = "ACCACC"
            kozak4 = "ACCGCC"
            sixbefore = str(dnasequence[index - 6 :index])
            oneafter=""
            start = str(dnasequence[index: index + 3])
            if len(dnasequence[index : len(dnasequence)]) > 3 :
                oneafter = str(dnasequence[index+3])
                infolist.append(str(dnasequence[index - 6: index + 4]))
            else:
                infolist.append(str(dnasequence[index - 6: index + 3]))
            if sixbefore == kozak and oneafter == "G" and start == "ATG":
                infolist.append("Strong")
            elif sixbefore == kozak2 and oneafter == "G" and start == "ATG":
                infolist.append("Strong")
            elif sixbefore == kozak3 and oneafter == "G" and start == "ATG":
                infolist.append("Strong")
            elif sixbefore == kozak4 and oneafter == "G" and start == "ATG":
                infolist.append("Strong")
            else:
                infolist.append("Weak")
    return infolist


dna = read_dna()  #reads the cdna file and stores a list of genes with their nucleotide sequences
pep = read_protein()  #reads the protein file and stores a list of peptide sequences
genesofinterest = id_tocheck()  #reads the file that stores the gene ids to investigate (make sure each gene id is from ensemble and is separated by a new line)

#below are four empty lists that will be used later to create a data frame for pandas
geneidlist = []   #list of gene ids
transcriptidlist = []  #list of transcript ids
startregionlist = []  #list of nucleotide sequences - start region of translation for different transcripts
resultlist = []  #list of results after the anaylsis of the start region


for gene in genesofinterest:  #iterating through each gene from the gene ids we want to investigate
    transcripts, sequences = list_transcript(gene, dna) #look for the gene id in the cdna file, return a list of possible transcripts for that gene along with their dna sequences (this is necessary because there are multiple transcripts versions for the same gene)
    if len(transcripts) == 0: #if there are no possible transcripts for that gene in the cdna file (geneid is not in the file), store that information
        geneidlist.append(gene)
        transcriptidlist.append(" ")
        startregionlist.append(" ")
        resultlist.append("gene not in file")
    else:
        for i in transcripts:  #for each element in the list of possible transcripts
            for c in sequences:  #for each element in the list of dna sequences that we stored before
                if i == cdna_transcript(c):  #make sure that the dna sequence has the matching transcript ids with the transcript in question
                    matchingprotein = find_transcript(c, pep)  #from the list of peptide sequences, find the protein with the same transcript id as the transcript in question
                    if matchingprotein == -1:  #if there is no protein with that transcript id, store that information
                        geneidlist.append(gene_id(c)[5:len(gene_id(c))])
                        transcriptidlist.append(i)
                        startregionlist.append(" ")
                        resultlist.append("no matching protein")
                    else:
                        firstaa = str(pep[matchingprotein].seq[0]) #retrieve the first amino acid of the peptide sequence for the protein
                        if firstaa == "M":  #if the protein starts with a methionine
                            startindex = find_ATG(c.seq, pep[matchingprotein].seq) #finds the start region by translating the cdna sequence starting from every occurence of ATG and matching that translation to the peptide sequence of the protein
                        else:
                            startindex = find_other_start(firstaa, c.seq, pep[matchingprotein].seq) #finds the start region by translating the cdna sequence from every possible codon that could encode for the first amino acid and matching that translation to the peptide sequence of the protein
                        #at this point, we know where the translation starts on the cdna sequence, or we know that the cdna sequence and peptide sequence didn't match
                        info = strong_or_weak(c.seq, gene_id(c), i, startindex)  #analyzes the 6 nucleotides before the start codon, the start codon itself, and one nucleotide after the start codon
                        #the "info" contains an ordered list of information: [geneid, transcriptid, 10 nucleotide start region, result of the analysis]
                        geneidlist.append(info[0])
                        transcriptidlist.append(info[1])
                        startregionlist.append(info[2])
                        resultlist.append(info[3])


#creates a dataframe using the information from above and stores it in a friendly spreadsheet using pandas
df = pd.DataFrame({'GeneID': geneidlist,
                       'Transcript': transcriptidlist,
                       'Start Region': startregionlist,
                       'Results': resultlist})
writer = pd.ExcelWriter('output.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='output1', index=False)
writer.save()