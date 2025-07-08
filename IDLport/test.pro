PRO testIt

; COR2A
;fileA = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_002330_d4c2A.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_125330_d4c2A.fts'

; COR2B
fileA = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_162400_d4c2B.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_185400_d4c2B.fts'


;COR1A
fileA = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_001000_n4c1A.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_011000_n4c1A.fts'

;COR1B
fileA = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_001000_n4c1B.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_011000_n4c1B.fts'


files = [fileA, fileB] 
secchi_prep, files, headers, images  

dim = images[*,*,1]-images[*,*,0]
mass=SCC_CALC_CME_MASS(dim, headers[1], /all)
print, mass[420,96]

END




PRO testSP
fileA = '/Users/kaycd1/wombat/fits/20241028_002330_d4c2A.fts'
fileB = '/Users/kaycd1/wombat/fits/20241028_125330_d4c2A.fts'
secchi_prep, [fileA, fileB], hdr, im
END