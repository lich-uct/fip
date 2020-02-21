import os

from fiprofiling import relationmatrices as rm

with open(os.path.join(os.path.dirname(__file__), 'drugbank_ecfp6_1024.csv')) as csvfile:
    corm = rm.CORM.from_csv(csvfile)
coprm = rm.COPRM.from_CORM(corm)
pmirm = rm.PMIRM.from_COPRM(coprm)
zpmirm = rm.ZPMIRM.from_PMIRM(pmirm)
