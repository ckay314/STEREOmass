PRO massIt, date, corID
; Syntax is massIt, YYYYMMDD, 'corsetID_0*'
; Command line calling needs corID to be a string to include _

; Make the path to the info file based on date/corID
basePath = '/Volumes/SRP/vourla1/laura/web_database/catalog/'

; Path to the output folder
outPath = '/Users/kaycd1/STEREO_Mass/MassFits/'

; Break date string into parts
date = STRTRIM(date,1)
yr = STRMID(date, 0,4)
date2 = STRMID(date, 2,6)
mn = STRMID(date, 4,2)

; Pull satellite from corID
corID = STRTRIM(corID,1)
sat = 'A'
print, corID
if STRMATCH(corID, '*B*') then sat = 'B'

; Build the full path
infoFile = basePath + sat + '/' + yr + '/' + mn + '/' + date2 + '/' + corID + '/' + corID +'.info'

; Check that file exists (some don't)
isFile = FILE_TEST(infoFile)
IF (isFile eq 1) THEN BEGIN

	print, 'Reading in CORSET file:'
	print, infoFile


	; Open the info file
	OPENR, lun, infoFile, /GET_LUN

	; Read one line at a time, saving the result into array
	infoF = ''
	line = ''
	; Look for the lines indicating the files we want while doing this
	baseIdx  = -1
	fileIdx0 = -1
	counter = 0
	WHILE NOT EOF(lun) DO BEGIN 
		counter += 1
	
		; Read the line and add to the array
		READF, lun, line 
	  	infoF = [infoF, line]
  
	  	; Look for the base image line
	  	if strmatch(line, 'Base image*') then baseIdx = counter + 1
	
		; Look for the start of the files
		if strmatch(line, 'Files used*') then fileIdx0 = counter + 1
	ENDWHILE

	; Close the file and free the file unit
	FREE_LUN, lun

	; Check that we got all indices
	IF baseIdx eq -1 or fileIdx0 eq -1 then begin
		print, 'Error in determining either base image or sequence images'
		print, 'baseIdx', baseIdx
		print, 'fileIdx0', fileIdx0
		STOP
	ENDIF

	fileIdxs = [baseIdx, indgen(counter-fileIdx0+1) + fileIdx0]
	files = []

	FOR i = 0, n_elements(fileIdxs)-1 DO BEGIN
		files = [files, infoF[fileIdxs[i]]]
	ENDFOR

	; Run everyone through secchi_prep
	secchi_prep, files, headers, images, /smask_on, outsize=1024, /silent

	FOR i = 1, n_elements(fileIdxs)-1 DO BEGIN
		; Get the base difference
		dim = images[*,*,i]-images[*,*,0]
	
		; Rotate it (or not?)
		;scc_roll_image, headers[i], dim, missing=0, cubic=1

		; Convert dim to mass
		mass=SCC_CALC_CME_MASS(dim, headers[i], /all)

		; Save the file
		myFile = infoF[fileIdxs[i]]
		startIdx = STRPOS(myFile,'/', /REVERSE_SEARCH)
		endIdx = STRPOS(myFile,'.', /REVERSE_SEARCH)
		newstr = corID + '_'+STRMID(myFile, startIdx+1, endIdx-startIdx-1) + '_mass'
		print, 'Saving ', newstr
		
		; Un secchi the header
		newhead = struct2fitshead(headers[i])
		WRITEFITS, outpath+newstr+'.fts', mass, newhead
	ENDFOR

	; Plot it, or not
	IF 0 THEN BEGIN
		; the best mass option so far
		maxMass = max(abs(mass)) / 40
		rngIt = (mass / maxMass)
		rngIt[where(rngIt LT -0.99)] = -0.99
		rngIt[where(rngIt GT 0.99)] = 0.99
		sclIt = rngIt * 128 + 128

		loadct, 0
		window,1,xsize=1024,ysize=1024
		tv, sclIt
		; CK has given up on writing time stamps because IDL is stupid and difficult
		; and the font looks like trash anyway
		;WRITE_PNG, 'temp.png', TVRD(/TRUE)
	ENDIF
ENDIF ELSE BEGIN
	print, "Cannot find file "+ infoFile
ENDELSE
END


PRO runFullYear, thisYr
; Open the info file
strYr = string(thisYr)
strYr = strYr.trim()
yrFile = 'allGood'+ strYr + '.txt'
OPENR, lun, yrFile, /GET_LUN

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

FOR i = 1, counter DO BEGIN
	print, ''
	print, '|--------------------------------|'
	print, '   On event'+ string(i) + ' of' + string(counter)
	print, '|--------------------------------|'
	myLine = infoF[i]
	splitLine = STRSPLIT(myLine, /EXTRACT)
	corID = splitLine[1]
	corID = REPSTR(corID, 'c', '')
	date  = splitLine[2]
	
	IF thisYr eq 2014 THEN BEGIN
		date = REPSTR(date, '/', '')
		date = STRMID(date, 4,4)+STRMID(date, 0,2)+STRMID(date, 2,2)
		print, date
	ENDIF ELSE BEGIN
		date = REPSTR(date, '-', '')
	ENDELSE
	print, date, corID
	massIt, date, corID
ENDFOR




END