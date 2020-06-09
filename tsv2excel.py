import sys
import pandas as pd
import glob
import os

targetDir = sys.argv[1]

for result in glob.glob(os.path.join(targetDir,"*.tsv")):
    filename = os.path.basename(".".join(result.split(".")[:-1]))
    print(filename)
    try:
        df = pd.read_csv(result, sep="\t")
        df.to_excel(os.path.join(targetDir,"{}.xlsx".format(filename)))
    except Exception as exp:
        print(exp)


    # print(result)
