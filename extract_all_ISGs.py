#-----------------------#
# A short script to extract all relevant ISGs from Shaw et al
# Threshold can be altered
#-----------------#

import pandas as pd

df = pd.read_csv('shaw_ISG_data.txt', sep='\t')

# Filtr stage - Can alter threshold here
filtered = df[df['Log2FC'] > 2]

filtered['Gene Name'].to_csv('filtered_ISG.txt', sep='\t', index=False)
