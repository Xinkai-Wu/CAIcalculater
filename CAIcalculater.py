import pandas as pd
import numpy as np
import re

import sys

# load CDS 
if len(sys.argv) < 2:
    print("Please provide the coding sequence!")
else:
    cds_seq = sys.argv[1]

# load weigthed RSCU Score
ws_df = pd.read_csv("codon_adaptive_ws.txt",sep="\t",index_col=0)

def _cai_seq(seq, ws):
    """Calculate the CAI (float) for the provided cds sequence (string).
    This method uses the Index and returns the CAI for the given sequence.
    """
    seq = seq.upper()

    codon_id, codon_count = np.unique(re.findall('...', seq),return_counts=True)
    codon_count = pd.Series(codon_count, index=codon_id, name='codons')
    
    # remove start or stop codons
    df = pd.concat([ws, codon_count], axis=1, sort=True).loc[ws.index]
    
    # calculate cai
    codon_prop = df['codons'] / df['codons'].sum()
    return np.exp(np.sum(df.iloc[:,0] * codon_prop))

# cds_seq = 'ATGGTG'

print("CAI of coding sequence:")
print(_cai_seq(cds_seq,ws_df['cui_ln_ws']))
