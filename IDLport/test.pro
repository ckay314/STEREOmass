PRO testIt
fileA = '/Users/kaycd1/wombat/fits/20241028_002330_d4c2A.fts'
fileB = '/Users/kaycd1/wombat/fits/20241028_125330_d4c2A.fts'
files = [fileA, fileB] 
secchi_prep, files, headers, images  


dim = images[*,*,1]-images[*,*,0]
print, ''
print, ''
print, ''
mass=SCC_CALC_CME_MASS(dim, headers[1], /all)
print, mass[512,512] 


END

PRO testSP
fileA = '/Users/kaycd1/wombat/fits/20241028_002330_d4c2A.fts'
fileB = '/Users/kaycd1/wombat/fits/20241028_125330_d4c2A.fts'
secchi_prep, [fileA, fileB], hdr, im
END