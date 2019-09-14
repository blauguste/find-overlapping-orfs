from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import pickle

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
def find_orf(df_in, aoi, path_to_gbs):
    
    ct = 0
    rd = dict()

    aoi = pickle.load(open(aoi, "rb"))

    print("reading infile")
    df_in = pd.read_pickle(df_in)
    df_sub = df_in[df_in['target_name'].isin(aoi)] # limit dataframe to just files with genbanks locally available

    # sort the dataframe by sRNA start/end position

    for acc in aoi:

        print("indexing", acc)
        tree_path =  path_to_gbs + acc + '.gb'
        
        with open(tree_path, "r") as gb_in:
        
            # read all the genbank features into a sorted list
            
            feat_array = []
            tup_array = []
            
            gb = SeqIO.read(gb_in, 'genbank')
            
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
                            'CDS_start': feat_array[result].location.start, \
                            'CDS_end': feat_array[result].location.end, \
                            'CDS_strand': feat_array[result].location.strand, \
                            'CDS_product': feat_array[result].qualifiers['product'][0]}
                    nsp = result
                    print(qn)
                    print(acc)
                    print(feat_array[result])
                    ct += 1
                    
    with open('mother_gene_90_results.p', 'wb') as outfile:
        pickle.dump(rd, outfile)

    resdf = pd.DataFrame.from_dict(rd, orient='index')
    resdf.to_csv('mother_gene_90_results.csv')






