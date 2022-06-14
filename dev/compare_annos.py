import string
from nltk.corpus import stopwords
import Levenshtein as lev
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', None)


def remove_punctuations(txt, punct = string.punctuation):
    '''
    This function will remove punctuations from the input text
    '''
    return ''.join([c for c in txt if c not in punct])
  
def remove_stopwords(txt, sw = list(stopwords.words('english'))):
    '''
    This function will remove the stopwords from the input txt
    '''
    return ' '.join([w for w in txt.split() if w.lower() not in sw])

def clean_text(txt):
    '''
    This function will clean the text being passed by removing specific line feed characters
    like '\n', '\r', and '\'
    '''
    
    txt = txt.replace('\n', ' ').replace('\r', ' ').replace('\'', '')
    txt = remove_punctuations(txt)
    txt = remove_stopwords(txt)
    return txt.lower()

def calc_dist(anno1, anno2, th=0.50):

    dist = []
    for a1, a2 in zip(anno1, anno2):
        if ((a1 != '') & (a2 != '')):
            ld = lev.distance(a1, a2)
            ml = max(len(a1), len(a2)) 
            if ml != 0:
                d = 1 - ld/ml
            else:
                d = 0
        else:
            d = np.nan
        dist.append(d)

    return dist


anno_file = 'WGS000001.ORF_annotation_table.txt'
anno_df = pd.read_csv(anno_file, sep='\t', header=0)
anno_df.fillna('', inplace=True)

# Clean columns
# cog-2020
anno_df['cog-2020'] = [clean_text(x.split(' [')[0])
                       if '[' in x else clean_text(x)
                       for x in anno_df['cog-2020']
                       ]
# metacyc-26
anno_df['metacyc-26'] = [clean_text(x.split(' (')[0])
                         if '(' in x else clean_text(x)
                         for x in anno_df['metacyc-26']
                         ]
# refseq-protein-212
anno_df['refseq-protein-212'] = [clean_text(x.split(' [')[0].replace('MULTISPECIES: ', ''))
                                 if '[' in x else clean_text(x)
                                 for x in anno_df['refseq-protein-212']
                                 ]
# uniprot_sprot
anno_df['uniprot_sprot'] = [clean_text(x.split(' OS ')[0])
                            if ' OS ' in x else clean_text(x)
                            for x in anno_df['uniprot_sprot']
                            ]
# uniref100
anno_df['uniref100'] = [clean_text(x.split(' n ')[0])
                        if ' n ' in x else clean_text(x)
                        for x in anno_df['uniref100']
                        ]
# uniref50
anno_df['uniref50'] = [clean_text(x.split(' n ')[0])
                       if ' n ' in x else clean_text(x)
                       for x in anno_df['uniref50']
                       ]
# uniref90
anno_df['uniref90'] = [clean_text(x.split(' n ')[0])
                       if ' n ' in x else clean_text(x)
                       for x in anno_df['uniref90']
                       ]

anno_df['rs_ur50'] = calc_dist(anno_df['refseq-protein-212'], anno_df['uniref50'])
anno_df['rs_ur90'] = calc_dist(anno_df['refseq-protein-212'], anno_df['uniref90'])
anno_df['rs_ur100'] = calc_dist(anno_df['refseq-protein-212'], anno_df['uniref100'])
anno_df['ur50_ur90'] = calc_dist(anno_df['uniref50'], anno_df['uniref90'])
anno_df['ur50_ur100'] = calc_dist(anno_df['uniref50'], anno_df['uniref100'])
anno_df['ur90_ur100'] = calc_dist(anno_df['uniref90'], anno_df['uniref100'])

nan = np.nan
print('RefSeq vs UniRef50:')
print('\tAvg. Similarity:', anno_df['rs_ur50'].mean())
print('\tCount:', len(anno_df.query("`refseq-protein-212` != ''")['refseq-protein-212']), ':',
                  len(anno_df.query("uniref50 != ''")['uniref50'])
                  )

print('RefSeq vs UniRef90:')
print('\tAvg. Similarity:', anno_df['rs_ur90'].mean())
print('\tCount:', len(anno_df.query("`refseq-protein-212` != ''")['refseq-protein-212']), ':',
                  len(anno_df.query("uniref90 != ''")['uniref90'])
                  )

print('RefSeq vs UniRef100:')
print('\tAvg. Similarity:', anno_df['rs_ur100'].mean())
print('\tCount:', len(anno_df.query("`refseq-protein-212` != ''")['refseq-protein-212']), ':',
                  len(anno_df.query("uniref100 != ''")['uniref100'])
                  )

print('UniRef50 vs UniRef90:')
print('\tAvg. Similarity:', anno_df['ur50_ur90'].mean())
print('\tCount:', len(anno_df.query("uniref50 != ''")['uniref50']), ':',
                  len(anno_df.query("uniref90 != ''")['uniref90'])
                  )

print('UniRef50 vs UniRef100:')
print('\tAvg. Similarity:', anno_df['ur50_ur100'].mean())
print('\tCount:', len(anno_df.query("uniref50 != ''")['uniref50']), ':',
                  len(anno_df.query("uniref100 != ''")['uniref100'])
                  )

print('UniRef90 vs UniRef100:')
print('\tAvg. Similarity:', anno_df['ur90_ur100'].mean())
print('\tCount:', len(anno_df.query("uniref90 != ''")['uniref90']), ':',
                  len(anno_df.query("uniref100 != ''")['uniref100'])
                  )











