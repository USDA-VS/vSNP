#!/usr/bin/env python

import glob
import pandas as pd
import time
from datetime import datetime

'''
working directory to contain Excel files to simply merge to a single file
'''

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

files = glob.glob("*.xlsx")
df_list = []
for fff in files:
    df_list.append(pd.read_excel(fff))

result = pd.concat(df_list, sort=False)
result.to_excel(f'combined_excelworksheets-{st}.xlsx')
