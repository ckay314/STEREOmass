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
;fileA = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_001000_n4c1B.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/COR1B_20111114_011000_n4c1B.fts'

;HI1A
;fileA = '/Users/kaycd1/wombat/fits/testing/HI1A_20100501_000901_s4h1A.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/HI1A_20100501_164901_s4h1A.fts'

;HI1B
;fileA = '/Users/kaycd1/wombat/fits/testing/HI1B_20070801_012900_s4h1B.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/HI1B_20070801_164900_s4h1B.fts'

;HI2A
;fileA = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_000921_s4h2A.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_180921_s4h2A.fts'

;HI2B
fileA = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_000921_s4h2B.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_100921_s4h2B.fts'

files = [fileA, fileB] 
secchi_prep, files, headers, images  

dim = images[*,*,1]-images[*,*,0]
mass=SCC_CALC_CME_MASS(dim, headers[1], /all)
print, mass[140,196]

END




PRO testSP
fileA = '/Users/kaycd1/wombat/fits/20241028_002330_d4c2A.fts'
fileB = '/Users/kaycd1/wombat/fits/20241028_125330_d4c2A.fts'
secchi_prep, [fileA, fileB], hdr, im
END