PRO testIt

; COR2A
;fileA = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_002330_d4c2A.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/COR2A_20241028_125330_d4c2A.fts'

; COR2B
fileA = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_162400_d4c2B.fts'
fileB = '/Users/kaycd1/wombat/fits/testing/COR2B_20120712_185400_d4c2B.fts'


;COR1A
;fileA = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_001000_n4c1A.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/COR1A_20121231_011000_n4c1A.fts'

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
;fileA = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_000921_s4h2B.fts'
;fileB = '/Users/kaycd1/wombat/fits/testing/HI2B_20130201_100921_s4h2B.fts'

files = [fileA, fileB] 
secchi_prep, files, headers, images, outsize=1024  
dim = images[*,*,1]-images[*,*,0]
mass=SCC_CALC_CME_MASS(dim, headers[1], /all)
print, mass[444,111]
		;loadct, 1
		;window,1,xsize=1024,ysize=1024
		;tv, signum(mass) * alog10(abs(mass))
END




PRO testPSP
; PSP WISPR
fileA = '/Users/kaycd1/wombat/fits/testing/psp_L2_wispr_20250610T000025_V0_1221.fits'
fileB = '/Users/kaycd1/wombat/fits/testing/psp_L2_wispr_20250610T203026_V0_1221.fits'


;fileA = '/Users/kaycd1/wombat/fits/testing/psp_L3_wispr_20210121T000421_V1_1211.fits'
;fileB = '/Users/kaycd1/wombat/fits/testing/psp_L3_wispr_20210121T005221_V1_1211.fits'


files = [fileA, fileB] 
wispr_prep, files, headers, images  
;help, headers[0]
print, sd
; think this just views a processed image?
;im = wispr_mk_frame(fileA,bytscl=0)

dim = images[*,*,1];-images[*,*,0]
mass=SCC_CALC_CME_MASS(dim, headers[1], /all)


dn2msb1 = 3.93e-14 ; Inner conversion factor from dn to msb
dn2msb2 = 5.78e-14 ; Outer cf

print, max(mass) / dn2msb1
;print, max(dim)
;help, mass
if 0 then begin
	; force a decent range to plot
	y=histogram(mass, nbins=100, locations=x, /nan)
	mx=max(y, loc)
	;stop
	while mx/total(y) gt .5 do begin
		if loc eq 0 then hr=[x[0], x[loc+1]] else hr=[x[loc-1],x[loc+1]]
		y=histogram(mass, nbins=100, locations=x, min=hr[0], max=hr[1], /nan)
		mx=max(y, loc)
	endwhile

	cdf=total(y, /cumulative)/total(y)
	r1=where(cdf gt .45)
	r2=where(cdf gt .55)
	if r2[0] eq r1[0] then r2[0]=r2[1]
	range_im=[x[r1[0]], x[r2[0]]] 
	imbs=bytscl(mass, range_im[0], range_im[1])

	loadct, 0
	window,1,xsize=1024,ysize=960
	tv, imbs
endif

		
end



PRO testCoords
;fileA = '/Users/kaycd1/wombat/fits/testing/HI2A_20090301_000921_s4h2A.fts'
fileA = '/Users/kaycd1/wombat/fits/20120712_172400_d4c2A.fts'
secchi_prep, fileA, hdr, im

myWCS = fitshead2wcs(hdr)
pt = [-10873.3709868, -260.11013321]
print, wcs_get_pixel(myWCS, pt )



END







