#for encoding peptides and proteins
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
from sklearn.model_selection import train_test_split
import pandas as pd
import os, re, math, platform
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import linalg as la
import argparse
from itertools import *
SPLITSEED = 810
from sklearn.decomposition import KernelPCA as KPCA

_AALetter = ['A', 'C', 'D', 'E', 'F', 'G', 'H',
             'I', 'K', 'L', 'M', 'N', 'P', 'Q',
             'R', 'S', 'T', 'V', 'W', 'Y']
_ftypes = ["AAC", "DiAAC", "PHYCs"]

word2int1DDict = {'G': 1, 'A': 2, 'V': 3, 'L': 4, 'I': 5, 'P': 6,
                               'F': 7, 'Y': 8, 'W': 9, 'S': 10, 'T': 11, 'C': 12,
                               'M': 13, 'N': 14, 'Q': 15, 'D': 16, 'E': 17,
                               'K': 18, 'R': 19, 'H': 20, 'X': 21, 'B': 22,
                               'J': 23, 'O': 24, 'U': 25, 'Z': 26}

int1D2wordDict = {1: 'G',2: 'A',3: 'V',4: 'L', 5: 'I',6: 'P',7: 'F', 8: 'Y',
                  9: 'W', 10: 'S', 11: 'T', 12: 'C', 13: 'M', 14: 'N', 15: 'Q',
                  16: 'D', 17: 'E', 18: 'K', 19: 'R', 20: 'H', 21: 'X', 22: 'B',
                  23: 'J', 24: 'O', 25: 'U', 26: 'Z'}

#读取文件
def read_fasta(fname):
    with open(fname, "rU") as f:
        seq_dict = [(record.id, record.seq._data) for record in SeqIO.parse(f, "fasta")]
    seq_df = pd.DataFrame(data=seq_dict, columns=["Id", "Sequence"])
    return seq_df

#AAC
def insert_AAC(seq_df):
    # Compute AAC for peptide in specific A.A
    def get_aac(seq, aa):
        return seq.count(aa) / len(seq) * 100

    # processing data_frame
    data_size = seq_df.size
    for ll in _AALetter:
        seq_df['AAC_{}'.format(ll)] = list(map(get_aac, seq_df['Sequence'], [ll] * data_size))
    return seq_df

#AAI
def AAI_1(fastas):
    encodings = []
    fileAAindex1 = open(R'Features/pre/AAindex_1.txt')
    fileAAindex2 = open(R'Features/pre/AAindex_2.txt')
    records1 = fileAAindex1.readlines()[1:]
    records2 = fileAAindex2.readlines()[1:]
    AAindex1 = []
    AAindex2 = []
    for i in records1:
        AAindex1.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
    for i in records2:
        AAindex2.append(i.rstrip().split()[1:] if i.rstrip() != '' else None)
    index = {}
    for i in range(len(_AALetter)):
        index[_AALetter[i]] = i
    fastas_len = len(fastas)
    for i in range(len(AAindex1)):
        total = 0
        for j in range(fastas_len):
            temp = AAindex1[i][index[fastas[j]]]
            total = total + float(temp)
        encodings.append(total / fastas_len)
    for i in range(len(AAindex2)):
        total = 0
        for j in range(fastas_len):
            temp = AAindex2[i][index[fastas[j]]]
            total = total + float(temp)
        encodings.append(total)
    return encodings

def AAI(seqs):
    encodings = []
    header = []
    for i in range(36):
        header.append('AAI_' + str(i))
    encodings.append(header)
    for fastas in seqs:
        fastas_NT5 = "%s" % fastas[:5]
        fastas_CT5 = "%s" % fastas[-5:]

        encodings_full = AAI_1(fastas)
        encodings_CT5 = AAI_1(fastas_CT5)
        encodings_NT5 = AAI_1(fastas_NT5)
        encodings.append(encodings_full + encodings_NT5 + encodings_CT5)
    return encodings

def insert_AAI(seq_df):
    enconding = AAI(seq_df['Sequence'])
    enconding = pd.DataFrame(enconding[1:], columns=enconding[0])
    seq_df = pd.concat([seq_df, enconding.iloc[:, :]], axis=1)
    return seq_df


#Kmer_KPCA
#定义Kmer函数
def TransDict_from_list(groups):
    transDict = dict()
    tar_list = ['0', '1', '2', '3', '4', '5', '6']
    result = {}
    index = 0
    for group in groups:
        g_members = sorted(group)  # Alphabetically sorted list
        for c in g_members:
            # print('c' + str(c))
            # print('g_members[0]' + str(g_members[0]))
            result[c] = str(tar_list[index])  # K:V map, use group's first letter as represent.
        index = index + 1
    return result

def translate_sequence(seq, TranslationDict):
    '''
    Given (seq) - a string/sequence to translate,
    Translates into a reduced alphabet, using a translation dict provided
    by the TransDict_from_list() method.
    Returns the string/sequence in the new, reduced alphabet.
    Remember - in Python string are immutable..

    '''
    import string
    from_list = []
    to_list = []
    for k, v in TranslationDict.items():#遍历TD
        from_list.append(k)
        to_list.append(v)
    # TRANS_seq = seq.translate(str.maketrans(zip(from_list,to_list)))
    TRANS_seq = seq.translate(str.maketrans(str(from_list), str(to_list)))
    # TRANS_seq = maketrans( TranslationDict, seq)
    return TRANS_seq

def get_3_protein_trids():
    nucle_com = []
    chars = ['0', '1', '2', '3', '4', '5', '6']
    base = len(chars)
    end = len(chars) ** 3
    for i in range(0, end):
        n = i
        ch0 = chars[n % base]#求模运算，相当于mod，也就是计算除法的余数
        n = n / base
        ch1 = chars[int(n % base)]
        n = n / base
        ch2 = chars[int(n % base)]
        nucle_com.append(ch0 + ch1 + ch2)
    return nucle_com

def get_4_nucleotide_composition_KPCA(tris, seq, pythoncount=True):
    seq_len = len(seq)
    tri_feature = [0] * len(tris)
    k = len(tris[0])
    note_feature = [[0 for cols in range(len(seq) - k + 1)] for rows in range(len(tris))]
    if pythoncount:
        for val in tris:
            num = seq.count(val)
            tri_feature.append(float(num) / seq_len)
    else:
        # tmp_fea = [0] * len(tris)
        for x in range(len(seq) + 1 - k):
            kmer = seq[x:x + k]
            if kmer in tris:
                ind = tris.index(kmer)
                # tmp_fea[ind] = tmp_fea[ind] + 1
                note_feature[ind][x] = note_feature[ind][x] + 1
        estimator = KPCA(n_components=1, kernel='rbf', gamma=15)
        tri_feature = estimator.fit_transform(note_feature).T
        tri_feature = tri_feature.flatten()
    # print tri_feature
        # pdb.set_trace()
    return tri_feature

def Kmer_KPCA(seq_df):
    encoding = []
    heard = []
    for i in range(343):
        heard.append("Kmer_KPCA_" + str(i))
    encoding.append(heard)
    protein_tris = get_3_protein_trids()
    groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    group_dict = TransDict_from_list(groups)
    for seq in seq_df:
        protein_seq = translate_sequence(seq, group_dict)
        protein_tri_fea = get_4_nucleotide_composition_KPCA(protein_tris, protein_seq, pythoncount =False)
        protein_tri_fea = list(protein_tri_fea)
        encoding.append(protein_tri_fea)
    return encoding

def insert_Kmer_KPCA(seq_df):
    encoding = Kmer_KPCA(seq_df["Sequence"])
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding],axis=1)
    return seq_df


#  CKSAAGP

def generateGroupPairs(groupKey):
    gPair = {}
    for key1 in groupKey:
        for key2 in groupKey:
            gPair[key1+'.'+key2] = 0
    return gPair

def minSequenceLength(fastas):#查看最小长度
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen

def cksaagp(fastas, gap = 3, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    if minSequenceLength(fastas) < gap+2:
        print('Error: all the sequence length should be greater than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
        return 0

    group = {'alphaticr': 'GAVLMI',
             'aromatic': 'FYW',
             'postivecharger': 'KRH',
             'negativecharger': 'DE',
             'uncharger': 'STCPNQ'}

    AA = 'ARNDCQEGHILKMFPSTWYV'

    groupKey = group.keys()

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    gPairIndex = []
    for key1 in groupKey:
        for key2 in groupKey:
            gPairIndex.append(key1+'.'+key2)

    encodings = []
    header = ['#']
    for g in range(gap + 1):
        for p in gPairIndex:
            header.append(p+'.gap'+str(g))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for g in range(gap + 1):
            gPair = generateGroupPairs(groupKey)
            sum = 0
            for p1 in range(len(sequence)):
                p2 = p1 + g + 1
                if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
                    gPair[index[sequence[p1]]+'.'+index[sequence[p2]]] = gPair[index[sequence[p1]]+'.'+index[sequence[p2]]] + 1
                    sum = sum + 1

            if sum == 0:
                for gp in gPairIndex:
                    code.append(0)
            else:
                for gp in gPairIndex:
                    code.append(gPair[gp] / sum)

        encodings.append(code)

    return encodings

def insert_CKSAAGP(seq_df, gap=2):
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = cksaagp(fastas, gap=gap)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]], axis=1)
    return seq_df

#PAAC
def minSequenceLengthWithNormalAA(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen
def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def paac(fastas, lambdaValue=30, w=0.05, **kw):
    if minSequenceLengthWithNormalAA(fastas) < lambdaValue + 1:
        print('Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
        return 0

    dataFile = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'\data\PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/PAAC.txt'
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
        AAProperty1.append([(j-meanI)/fenmu for j in i])

    encodings = []
    header = ['#']
    for aa in AA:
        header.append('Xc1.' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('Xc2.lambda' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(
                sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
                len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return encodings

def insert_PAAC(seq_df, lamb=3, w=0.4):
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = paac(fastas, lambdaValue=lamb, w=w)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]], axis=1)
    return seq_df

#CTD

# CTD_function
group1 = {'hydrophobicity_PRAM900101': 'RKEDQN', 'normwaalsvolume': 'GASTPDC', 'polarity': 'LIFWCMVY',
          'polarizability': 'GASDT', 'charge': 'KR', 'secondarystruct': 'EALMQKRH', 'solventaccess': 'ALFCGIVW'}
group2 = {'hydrophobicity_PRAM900101': 'GASTPHY', 'normwaalsvolume': 'NVEQIL', 'polarity': 'PATGS',
          'polarizability': 'CPNVEQIL', 'charge': 'ANCQGHILMFPSTWYV', 'secondarystruct': 'VIYCWFT',
          'solventaccess': 'RKQEND'}
group3 = {'hydrophobicity_PRAM900101': 'CLVIMFW', 'normwaalsvolume': 'MHKFRYW', 'polarity': 'HQRKNED',
          'polarizability': 'KMHFRYW', 'charge': 'DE', 'secondarystruct': 'GNPSD', 'solventaccess': 'MSPTHY'}
groups = [group1, group2, group3]
propertys = ('hydrophobicity_PRAM900101', 'normwaalsvolume', 'polarity', 'polarizability', 'charge', 'secondarystruct',
             'solventaccess')

def Count_C(sequence1, sequence2):
    sum = 0
    for aa in sequence1:
        sum = sum + sequence2.count(aa)
    return sum

def Count_D(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]
    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence))
                    break
        if myCount == 0:
            code.append(0)
    return code
def CTD(seqs):
    encodings = []
    header = []
    for i in range(147):
        header.append('CTD_' + str(i))
    encodings.append(header)
    for seq in seqs:
        code = []
        code2 = []
        CTDD1 = []
        CTDD2 = []
        CTDD3 = []
        aaPair = [seq[j:j + 2] for j in range(len(seq) - 1)]
        for p in propertys:
            c1 = Count_C(group1[p], seq) / len(seq)
            c2 = Count_C(group2[p], seq) / len(seq)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]

            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
            code2 = code2 + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
            CTDD1 = CTDD1 + [value / float(len(seq)) for value in Count_D(group1[p], seq)]
            CTDD2 = CTDD2 + [value / float(len(seq)) for value in Count_D(group2[p], seq)]
            CTDD3 = CTDD3 + [value / float(len(seq)) for value in Count_D(group3[p], seq)]
        encodings.append(code + code2 + CTDD1 + CTDD2 + CTDD3)
    return encodings

def insert_CTD(seq_df):
    enconding = CTD(seq_df['Sequence'])
    enconding = pd.DataFrame(enconding[1:], columns = enconding[0])
    seq_df = pd.concat([seq_df, enconding.iloc[:, :]],axis = 1)
    return seq_df

#AAE
def AAE_1(fastas):
    length = float(len(fastas))
    amino_acids = dict.fromkeys(_AALetter, 0)
    encodings = []
    for AA in amino_acids:
        hits = [a.start() for a in list(re.finditer(AA, fastas))]
        p_prev = 0
        p_next = 1
        sum = 0
        while p_next < len(hits):
            distance = (hits[p_next] - hits[p_prev]) / length
            sum += distance * math.log(distance, 2)
            p_prev = p_next
            p_next += 1
        amino_acids[AA] = -sum
        encodings.append(amino_acids[AA])
    return encodings

def AAE(seq):
    encodings = []
    header = []
    for i in range(60):
        header.append('AAE_' + str(i))
    encodings.append(header)
    for fastas in seq:
        fastas_NT5 = "%s" % fastas[:5]
        fastas_CT5 = "%s" % fastas[-5:]
        encodings_full = AAE_1(fastas)
        encodings_CT5 = AAE_1(fastas_CT5)
        encodings_NT5 = AAE_1(fastas_NT5)
        encodings.append(encodings_full + encodings_NT5 + encodings_CT5)
    return encodings

def insert_AAE(seq_df):
    encoding = AAE(seq_df['Sequence'])
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, :]], axis=1)
    return seq_df

#PHYCS
def insert_phycs(seq_df):
    #  Function for compute Isoelectric Point or net_charge of peptide
    def get_ieq_nc(seq, is_iep=True):
        protparam = PA(seq)
        return protparam.isoelectric_point() if is_iep else protparam.charge_at_pH(7.0)

    # Calculating IsoElectricPoints and NeutralCharge
    data_size = seq_df.size
    seq_df['IEP'] = list(map(get_ieq_nc, seq_df['Sequence'], [True] * data_size))  # IsoElectricPoints
    seq_df['Net Charge'] = list(map(get_ieq_nc, seq_df['Sequence'], [False] * data_size))  # Charge(Neutral)

    # Calculating hydrophobic moment (My assume all peptides are alpha-helix)
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'eisenberg')
    descrpt.calculate_moment(window=1000, angle=100, modality='max')
    seq_df['Hydrophobic Moment'] = descrpt.descriptor.reshape(-1)

    # Calculating "Hopp-Woods" hydrophobicity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'hopp-woods')
    descrpt.calculate_global()
    seq_df['Hydrophobicity'] = descrpt.descriptor.reshape(-1)

    # Calculating Energy of Transmembrane Propensity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'tm_tend')
    descrpt.calculate_global()
    seq_df['Transmembrane Propensity'] = descrpt.descriptor.reshape(-1)

    # Calculating Levitt_alpha_helical Propensity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'levitt_alpha')
    descrpt.calculate_global()
    seq_df['Alpha Helical Propensity'] = descrpt.descriptor.reshape(-1)

    # Calculating Aliphatic Index
    descrpt = GlobalDescriptor(seq_df['Sequence'].values)
    descrpt.aliphatic_index()
    seq_df['Aliphatic Index'] = descrpt.descriptor.reshape(-1)

    # Calculating Boman Index
    descrpt = GlobalDescriptor(seq_df['Sequence'].values)
    descrpt.boman_index()
    seq_df['Boman Index'] = descrpt.descriptor.reshape(-1)

    return seq_df

#GTPC
def GTPC(fastas):
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}

	groupKey = group.keys()
	baseNum = len(groupKey)
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key

	encodings = []
	header = ['#'] + triple
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])

		code = [name]
		myDict = {}
		for t in triple:
			myDict[t] = 0

		sum = 0
		for j in range(len(sequence) - 3 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
			sum = sum +1

		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
		encodings.append(code)

	return encodings

def insert_GTPC(seq_df):
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = GTPC(fastas)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]],axis=1)
    return seq_df

#DPC
diPeptides = [aa1 + aa2 for aa1 in _AALetter for aa2 in _AALetter]
def DPC(seqs):
    encodings = []
    hearder = []
    for i in range(400):
        hearder.append('DPC_' + str(i))
    encodings.append(hearder)
    for seq in seqs:
        AADict = {}
        for aa in range(len(_AALetter)):
            AADict[_AALetter[aa]] = aa

        tmpCode = [0] * 400
        for j in range(len(seq) - 2 + 1):
            tmpCode[AADict[seq[j]] * 20 + AADict[seq[j + 1]]] = tmpCode[AADict[seq[j]] * 20 + AADict[seq[j + 1]]] + 1
        if sum(tmpCode) != 0:
            tmpDPC = [i / sum(tmpCode) for i in tmpCode]
        encodings.append(tmpDPC)
    return encodings

def insert_DPC(seq_df):
    encodings = DPC(seq_df['Sequence'])
    encodings = pd.DataFrame(encodings[1:], columns=encodings[0])
    seq_df = pd.concat([seq_df, encodings.iloc[:,:]], axis=1)
    return seq_df

# QSO

def minSequenceLengthWithNormalAA(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen

def QSOrder(fastas, nlag=3, w=0.1, **kw):
    if minSequenceLengthWithNormalAA(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0

# 	dataFile = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'\data\Schneider-Wrede.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Schneider-Wrede.txt'
# 	dataFile1 = re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + r'\data\Grantham.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Grantham.txt'
    dataFile = 'data/Schneider-Wrede.txt'
    dataFile1 = 'data/Grantham.txt'
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    AA1 = 'ARNDCQEGHILKMFPSTWYV'

    DictAA = {}
    for i in range(len(AA)):
        DictAA[AA[i]] = i

    DictAA1 = {}
    for i in range(len(AA1)):
        DictAA1[AA1[i]] = i

    with open(dataFile) as f:
        records = f.readlines()[1:]
    AADistance = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance.append(array)
    AADistance = np.array([float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))

    with open(dataFile1) as f:
        records = f.readlines()[1:]
    AADistance1 = []
    for i in records:
        array = i.rstrip().split()[1:] if i.rstrip() != '' else None
        AADistance1.append(array)
    AADistance1 = np.array(
        [float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape((20, 20))

    encodings = []
    header = ['#']
    for aa in AA1:
        header.append('Schneider.Xr.' + aa)
    for aa in AA1:
        header.append('Grantham.Xr.' + aa)
    for n in range(1, nlag + 1):
        header.append('Schneider.Xd.' + str(n))
    for n in range(1, nlag + 1):
        header.append('Grantham.Xd.' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        arraySW = []
        arrayGM = []
        for n in range(1, nlag + 1):
            arraySW.append(sum([AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
            arrayGM.append(sum([AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
        myDict = {}
        for aa in AA1:
            myDict[aa] = sequence.count(aa)
        for aa in AA1:
            code.append(myDict[aa] / (1 + w * sum(arraySW)))
        for aa in AA1:
            code.append(myDict[aa] / (1 + w * sum(arrayGM)))
        for num in arraySW:
            code.append((w * num) / (1 + w * sum(arraySW)))
        for num in arrayGM:
            code.append((w * num) / (1 + w * sum(arrayGM)))
        encodings.append(code)
    return encodings

def insert_QSO(seq_df):
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = QSOrder(fastas)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]],axis=1)
    return seq_df

#NMBroto
def NMBroto(fastas,
            props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102', 'CHOC760101', 'BIGC670101', 'CHAM810101',
                   'DAYM780201'], nlag=4, **kw):
    if minSequenceLength(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0

    AA = 'ARNDCQEGHILKMFPSTWYV'
    fileAAidx = 'data/AAidx.txt'
    with open(fileAAidx) as f:
        records = f.readlines()[1:]
    myDict = {}
    for i in records:
        array = i.rstrip().split('\t')
        myDict[array[0]] = array[1:]

    AAidx = []
    AAidxName = []
    for i in props:
        if i in myDict:
            AAidx.append(myDict[i])
            AAidxName.append(i)
        else:
            print('"' + i + '" properties not exist.')
            return None

    AAidx1 = np.array([float(j) for i in AAidx for j in i])
    AAidx = AAidx1.reshape((len(AAidx), 20))
    pstd = np.std(AAidx, axis=1)
    pmean = np.average(AAidx, axis=1)

    for i in range(len(AAidx)):
        for j in range(len(AAidx[i])):
            AAidx[i][j] = (AAidx[i][j] - pmean[i]) / pstd[i]

    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i

    encodings = []
    header = ['#']
    for p in props:
        for n in range(1, nlag + 1):
            header.append(p + '.lag' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        N = len(sequence)
        for prop in range(len(props)):
            for n in range(1, nlag + 1):
                if len(sequence) > nlag:
                    # if key is '-', then the value is 0
                    rn = sum(
                        [AAidx[prop][index.get(sequence[j], 0)] * AAidx[prop][index.get(sequence[j + n], 0)] for j in
                         range(len(sequence) - n)]) / (N - n)
                else:
                    rn = 'NA'
                code.append(rn)
        encodings.append(code)
    return encodings

def insert_NMBroto(seq_df):
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = NMBroto(fastas)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]],axis=1)
    return seq_df
#将序列进行1,2,3.。编码
def word2int_1(fastas):
    encodings=[]
    j = len(fastas)
    for i in range(j):
        encodings.append(word2int1DDict[fastas[i]])
    return encodings

def word2int(seq):
    encodings = []
    header = []
    for i in range(100):
        header.append("word2int_"+ str(i))
    encodings.append(header)
    for fastas in seq:
        encoding = word2int_1(fastas)
        encodings.append(encoding)
    return encodings

def insert_word2int(seq_df):
    encodings = word2int(seq_df['Sequence'])
    encodings = pd.DataFrame(encodings[1:], columns=encodings[0])
    seq_df = pd.concat([seq_df, encodings.iloc[:, :]], axis=1)
    return seq_df


#ASDC
"""ASDC"""
Amino_acids = ['A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
'R', 'S','T','V', 'W','Y']
Amino_acids_ = list(product(Amino_acids,Amino_acids))
Amino_acids_ = [i[0]+i[1] for i in Amino_acids_]

def ASDC(seqs):
    header = []
    for i in range(400):
        header.append('ASDC_'+ str(i))
    seqs_ = []
    seqs_.append(header)
    for seq in seqs:
        ASDC_feature = []
        skip = 0
        for i in range(len(seq)): 
            ASDC_feature.extend(Skip(seq,skip)) 
            skip+=1
        seqs_.append([ASDC_feature.count(i)/len(ASDC_feature) for i in Amino_acids_])
    return seqs_

def Skip(seq,skip):
	element = []
	for i in range(len(seq)-skip-1):
		element.append(seq[i]+seq[i+skip+1])
	return element

def insert_ASDC(seq_df):
    encoding = ASDC(seq_df['Sequence'])
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, :]], axis=1)
    return seq_df

"""PSAAC"""
def PSAAC(seqs):
    header = []
    for i in range(40):
        header.append('ASDC_'+ str(i))
    seqs_ = []
    seqs_.append(header)
    PSAAC_profile_forward = []
    PSAAC_profile_backward = []
    forward_seq = []
    backward_seq = []
    i = 1
    for seq in seqs:
        forward_seq.append(list(seq[:5]))
        backward_seq.append(list(seq[-5:]))
    for position in range(5):
        PSAAC_profile_forward.append([list(np.array(forward_seq)[:,position]).count(amino)/len(seqs) for amino in Amino_acids])

    for position in range(5):
        PSAAC_profile_backward.append([list(np.array(backward_seq)[:,position]).count(amino)/len(seqs) for amino in Amino_acids])

    for seq in forward_seq:
        num = 0
        new_seq = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for amino in seq:
            index_ = Amino_acids.index(amino)
            new_seq[index_] = np.array(PSAAC_profile_forward)[num,index_]
            num+=1

        seqs_.append(new_seq)
    
    for seq in backward_seq:
        num = 0
        new_seq = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for amino in seq:
            index_ = Amino_acids.index(amino)
            new_seq[index_] = np.array(PSAAC_profile_backward)[num,index_]
            num+=1

        seqs_[i].extend(new_seq)
        i+=1
    return seqs_

def insert_PSAAC(seq_df):
    encodings = PSAAC(seq_df['Sequence'])
    encoding = pd.DataFrame(encodings[1:], columns=encodings[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, :]], axis=1)
    return seq_df