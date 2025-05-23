pro testHDR
;im = readfits('MassFits/1895.1417_0B_20120102_175400_d4c2B_mass.fts', hdr)
;im = readfits('/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/20120102/20120102_172400_d4c2B.fts', hdr)
im = readfits('MassFits/test.fts', hdr)
print, im[1]
end

pro suckItSP
file = '/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/20120102/20120102_172400_d4c2B.fts'

secchi_prep, file, hdr, im, /smask_on, outsize=1024, /silent
print, hdr
ahdr = struct2fitshead(hdr)
fits_info, file
writefits, 'MassFits/test.fts', im, ahdr

end

pro test3
;file = '/Volumes/SRP/vourla1/secchi/lz/L0/b/img/cor2/20120102/20120102_172400_d4c2B.fts'
file = 'MassFits/test.fts'
hdr = headfits(file)
print, hdr
end