PRO makeMass, im0, im, sat
	; Make the front part of the file path
	prefix = '/Volumes/SRP/vourla1/secchi/lz/L0/'
	IF sat eq 'A' THEN prefix = prefix + 'a/img/cor2/'
	IF sat eq 'B' THEN prefix = prefix + 'b/img/cor2/'

	; Clean up the image name and add to path
	im0 = REPSTR(im0, '.CK', '.fts')
	im = REPSTR(im, '.CK', '.fts')

	files = [prefix+im0, prefix+im]
	secchi_prep, files, headers, images, /smask_on, outsize=1024, /silent

	; Get the base difference
	dim = images[*,*,1]-images[*,*,0]

	; Convert dim to mass
	mass=SCC_CALC_CME_MASS(dim, headers[1], /all)

	; Save the file
	startIdx = STRPOS(im,'/', /REVERSE_SEARCH)
	endIdx = STRPOS(im,'.', /REVERSE_SEARCH)
	newstr = STRMID(im, startIdx+1, endIdx-startIdx-1) + '_mass'
	print, 'Saving ', newstr
		
	; Un secchi the header
	newhead = struct2fitshead(headers[1])
	
	; Path to the output folder
	outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'
	WRITEFITS, outpath+newstr+'.fts', mass, newhead
END



PRO runAll 
	; Open the info file
	fullFile = 'haveMassImage.txt'
	OPENR, lun, fullFile, /GET_LUN

	; Read one line at a time, saving the result into array
	infoF = ''
	line = ''
	counter = 0

	WHILE NOT EOF(lun) DO BEGIN 
		; Read the line and add to the array
		READF, lun, line 
	  	infoF = [infoF, line]
		counter++
	ENDWHILE

	FOR i = 248, counter DO BEGIN
		myLine = infoF[i]
		;print, myLine
		print, 'On case ', i, ' out of ', counter
		splitLine = STRSPLIT(myLine, /EXTRACT)
		if splitLine[1] ne 'Majumdar' THEN BEGIN
		 	;corID = splitLine[1]
			;corID = REPSTR(corID, 'c', '')
			Aim0 = splitLine[6]
			Aim  = splitLine[7]
			Bim0 = splitLine[11]
			Bim  = splitLine[12]
			if Aim ne 'None' THEN makeMass, Aim0, Aim, 'A'
			if Bim ne 'None' THEN makeMass, Bim0, Bim, 'B'
			
		ENDIF 
	ENDFOR
END