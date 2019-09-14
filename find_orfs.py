from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import pickle
import sys

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return (end1 >= start2) and (end2 >= start1)

def binarySearch (arr, l, r, x): 
  
    # Check base case 
    if r >= l: 
  
        mid = int(l + (r - l)/2)
  
        # If element is present at the middle itself 
        if overlap(arr[mid][0], arr[mid][1], x[0], x[1]):
            return mid 
          
        # If element is smaller than mid, then it can only 
        # be present in left subarray 
        elif arr[mid][0] > x[0]: 
            return binarySearch(arr, l, mid-1, x) 
  
        # Else the element can only be present in right subarray 
        else: 
            return binarySearch(arr, mid+1, r, x) 
  
    else: 
        # Element is not present in the array 
        return -1

# aoi is a pickled list of accessions of interest
def find_orf(email, df_in, soi):
    
    ct = 0
    rd = dict()

    print("reading infile")
    df_in = pd.read_pickle(df_in)

    df_sub = df_in[df_in['query_name'].isin(soi)] # limit dataframe to just sRNAs of interest
    df_sub.sort_values(by=['target_name', 'start'], inplace=True) # sort by coords

    with open(soi, 'rb') as infile:
        soi = pickle.load(infile)

    for acc in df_sub['target_name'].unique():

        print("fetching genbank", acc)

        Entrez.email = email
        foi = Entrez.efetch(db='nucleotide', id=acc, \
            rettype='gbwithparts', retmode='text')
        
        # read all the genbank features into a sorted list
        print("parsing genbank...")
            
        feat_array = []
        tup_array = []
        
        gb = SeqIO.read(foi, 'genbank')
        
        for feat in gb.features:
            if feat.type == 'CDS':
                feat_array.append(feat)
                tup_array.append((feat.location.start, feat.location.end))
        
        a = df_sub[df_sub['target_name'] == acc].sort_values(by='start')

        qtup = list(zip(a['start'].astype('int32'), a['end'].astype('int32')))

        nsp = 0
        for tup in qtup:
            
            result = binarySearch(tup_array, nsp, len(tup_array)-1, tup)
            if result != -1:
                qn = a['query_name'].where(a['start'] == tup[0]).dropna().values[0]
                qs = a['strand'].where(a['start'] == tup[0]).dropna().values[0]
                rd[ct] = {'query_name': qn, \
                    'qstart': tup[0], \
                    'qend': tup[1], \
                        'qstrand': qs, \
                    'target_genome': acc, \
                        'genome_desc': gb.description, \
                        'CDS_start': feat_array[result].location.start, \
                        'CDS_end': feat_array[result].location.end, \
                        'CDS_strand': feat_array[result].location.strand, \
                        'CDS_product': feat_array[result].qualifiers['product'][0]}
                nsp = result
                print(qn)
                print(acc)
                print(feat_array[result])
                ct += 1
                
    with open('mother_gene_isrK_OxyS_results.p', 'wb') as outfile:
        pickle.dump(rd, outfile)

    resdf = pd.DataFrame.from_dict(rd, orient='index')
    resdf.to_csv('mother_gene_isrK_OxyS_results.csv')

if __name__ == '__main__':
    if len(sys.argv) == 4:
         find_orf(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: find_orfs.py email@domain.com path_to_df.p sRNAs_of_interest.p")
         sys.exit(0)
