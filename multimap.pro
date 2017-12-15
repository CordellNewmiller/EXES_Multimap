; wishlist
; 
; sky subtraction with options
; prebinning along slit
; hypersampling
; medium mode
; 
; 

function blankexescube

  emptycube = {EXEScube}

  emptycube.und=ptr_new(/ALLOCATE_HEAP)
  emptycube.err=ptr_new(/ALLOCATE_HEAP)
  emptycube.wvm=ptr_new(/ALLOCATE_HEAP)
  emptycube.slit=ptr_new(/ALLOCATE_HEAP)
  emptycube.sky=ptr_new(/ALLOCATE_HEAP)



  emptycube.slitvec=ptr_new(/ALLOCATE_HEAP)
  emptycube.stepvec=ptr_new(/ALLOCATE_HEAP)

  emptycube.orderB=ptr_new(/ALLOCATE_HEAP)
  emptycube.orderE=ptr_new(/ALLOCATE_HEAP)
  emptycube.orderS=ptr_new(/ALLOCATE_HEAP)
  emptycube.orderT=ptr_new(/ALLOCATE_HEAP)

  emptycube.wavevec=ptr_new(/ALLOCATE_HEAP) ; vector

  emptycube.xpos=ptr_new(/ALLOCATE_HEAP)
  emptycube.ypos=ptr_new(/ALLOCATE_HEAP)

  emptycube.voxelind=ptr_new(/ALLOCATE_HEAP)



  return, emptycube

end

function exesmapinit, dir, undfilename, wvmfilename
  
  ;!null = {EXEScube}
  m = blankexescube()
    
  m.dir = dir
  m.undfilename = undfilename
  m.wvmfilename = wvmfilename
  
  skies = 3
  undraw = readfits(m.dir + m.undfilename, undheader) 
  *m.und = undraw[*,*,0:(((undraw.dim)[2]/2)-skies-1)]
  *m.err = undraw[*,*,(((undraw.dim)[2]/2)):(-skies-1)]
  
  *m.sky = undraw[*,*,(((undraw.dim)[2]/2)-skies):(((undraw.dim)[2]/2)-1)]
  
  m.dim = undraw.dim
  m.dim[2] = (((undraw.dim)[2]/2)-skies)
  
  wvmraw = readfits(m.dir + m.wvmfilename, wvmheader)
  *m.wvm = reform(wvmraw[*,*,0])

  *m.slit = reform(wvmraw[*,*,1])
  
  
  m.wn0 = double(sxpar(undheader,'WNO0'))
  m.object = strtrim(sxpar(undheader,'OBJECT'),2)
  m.instcfg = strtrim(sxpar(undheader,'INSTCFG'),2)
  m.slitwid = double(sxpar(undheader,'SLITWID'))
  
  m.ra = double(sxpar(undheader,'TELRA'))
  m.dec = double(sxpar(undheader,'TELDEC')) 
  m.angle = double(sxpar(undheader,'TELVPA'))
  
  
  *m.orderB = fix(strsplit(sxpar(undheader,'ORDR_B'),/extract,","))
  *m.orderE = fix(strsplit(sxpar(undheader,'ORDR_E'),/extract,","))
  *m.orderS = fix(strsplit(sxpar(undheader,'ORDR_S'),/extract,","))
  *m.orderT = fix(strsplit(sxpar(undheader,'ORDR_T'),/extract,","))
  
  m.step = double(sxpar(undheader,'MAPINTX'))
  m.pltscale = double(sxpar(undheader,'PLTSCALE'))
  
  
  return, m

  
end

function multimapinit, dir, undfilelist, wvmfile
  ; makes a hash of EXEScube objects
  
  n = undfilelist.length
  
  ;!null = blankexesmap()
  m = replicate({EXEScube}, n)
  for i = 0,(n-1) do begin
    m[i] = exesmapinit(dir,undfilelist[i],wvmfile)

  endfor
  
  return, m
  
end

pro exescubetrim, m, wavebound

  ; to do: find order, rotate if not medium mode
  order = 0

  n = m.length


  if strcmp(m[0].instcfg,'MEDIUM',/fold_case) then begin
    for i = 0,(n-1) do begin

      bottom=(*m[i].orderB)[order]
      ;ending=(*m[i].orderB)[order]
      ;start=(*m[i].orderB)[order]
      top=(*m[i].orderT)[order]



      *m[i].und = (*m[i].und)[(wavebound[0]):(wavebound[1]),bottom:top,*]
      *m[i].err = (*m[i].err)[(wavebound[0]):(wavebound[1]),bottom:top,*]
      *m[i].wvm = (*m[i].wvm)[(wavebound[0]):(wavebound[1]),bottom:top]
      *m[i].slit= (*m[i].slit)[(wavebound[0]):(wavebound[1]),bottom:top]

      m[i].dim = (*m[i].und).dim

      *m[i].slitvec = reform((*m[i].slit)[0,*])
      *m[i].stepvec = findgen((m[i].dim)[2]) * (m[i].step)
      *m[i].wavevec = reform((*m[i].wvm)[*,0])

    endfor
  endif else if strcmp(m[0].instcfg,'HIGH_MED',/fold_case) then begin
    for i = 0,(n-1) do begin

      bottom=(*m[i].orderB)[order]
      ending=(*m[i].orderB)[order]
      start=(*m[i].orderB)[order]
      top=(*m[i].orderT)[order]



;      *m[i].und = (*m[i].und)[(wavebound[0]):(wavebound[1]),bottom:top,*]
;      *m[i].err = (*m[i].err)[(wavebound[0]):(wavebound[1]),bottom:top,*]
;      *m[i].wvm = (*m[i].wvm)[(wavebound[0]):(wavebound[1]),bottom:top]
;      *m[i].slit= (*m[i].slit)[(wavebound[0]):(wavebound[1]),bottom:top]
;
;      m[i].dim = (*m[i].und).dim
;
;      *m[i].slitvec = reform((*m[i].slit)[0,*])
;      *m[i].stepvec = findgen((m[i].dim)[2]) * (m[i].step)
;      *m[i].wavevec = reform((*m[i].wvm)[*,0])
;
    endfor
  endif


end

function exesshrinkbin, data, newdim, error=error, duperror=duperror
  ;
  ; leave no keyword if input is data
  ; use /error if input is a variance plane
  ;       (gives equal weight to each measurement)
  ; use /duperror if input was resampled and is a variance plane
  ;       treats exactly equal measurements as the same measurement, weighted appropriately

  if newdim.length eq 2 then begin
    index = congrid(indgen(newdim),(data.dim)[0],(data.dim)[1])
  endif else if newdim.length eq 1 then begin
    index = congrid(indgen(newdim),(data.dim)[0])
  endif
  newarr=dblarr(newdim)

  if keyword_set(duperror) then begin
    for i=0,(newarr.length-1) do begin
      idat=data(where(i eq index))
      idatu=uniqc(idat,invec)
      newarr[i] = sqrt(total((invec*idatu)^2))/idat.length
      if 2*total(finite(data(where(i eq index)))) lt total(i eq index) then newarr[i] = !values.f_nan
    endfor
  endif else if keyword_set(error) then begin
    for i=0,(newarr.length-1) do begin
      newarr[i] = sqrt(total((data(where(i eq index)))^2,/nan))/total(i eq index)
      if 2*total(finite(data(where(i eq index)))) lt total(i eq index) then newarr[i] = !values.f_nan
    endfor
  endif else begin
    for i=0,(newarr.length-1) do begin
      newarr[i] = mean(data(where(i eq index)),/nan)
      if 2*total(finite(data(where(i eq index)))) lt total(i eq index) then newarr[i] = !values.f_nan
    endfor
  endelse


  return, newarr
end

pro exescubeprebin, m, binsize

  newdim=[und0crop.dim]*[1/binsize0,1,1]
  newund = dblarr(newdim)
  newerr = dblarr(newdim)

  for ei = 0,(m.length-1) do begin
    print,'Prebinning map ',i,' of ',(m.length-1)
    for i=0,((m[ei].dim)[2]-1) do begin
      (*m[ei].und)[*,*,i] = exesshrinkbin(reform((*m[ei].und)[*,*,i]),newdim[0:1])
      (*m[ei].err)[*,*,i] = exesshrinkbin(reform((*m[ei].err)[*,*,i]),newdim[0:1],/error)
    endfor
  endfor


end

pro exescubeplace, m

  n = m.length

  for i = 0,(n-1) do begin
    
    ; to add: account for multiple orders
    
    xcen=.5 * (*m[i].slitvec)[-1]
    ycen=.5 * (*m[i].stepvec)[-1]

    xvec = (*m[i].slitvec) -xcen
    yvec = (*m[i].stepvec) -ycen

    xunrot = xvec # replicate(1, n_elements(yvec))
    yunrot = replicate(1, n_elements(xvec)) # yvec
    
    intang=90.
    cosang = cos(-(m[i].angle +intang) * !pi/180.)
    sinang = sin(-(m[i].angle +intang) * !pi/180.)
    
    *m[i].xpos = xunrot * cosang + yunrot * sinang
    *m[i].ypos = yunrot * cosang - xunrot * sinang
        
    
  endfor

end

function pointercollapse, plist
  ;collapses arrays from pointers into a single dimension vector
  
  
  n = plist.length
  
  outvec = (*plist[0])[*]

  if n gt 1 then begin
    for i = 1,(n-1) do begin
      outvec=[outvec,(*plist[i])[*]]
    endfor
  endif
  
  return, outvec

end

;function exescubecoadd1, m
;  
;  n = m.length
;  
;  res = 2. ;m[0].slitwid
;  
;  xcollapse = pointercollapse(m.xpos)
;  ycollapse = pointercollapse(m.ypos)
;  
;  xmin= min(xcollapse,max=xmax)
;  ymin= min(ycollapse,max=ymax)
;  nx = ceil((xmax-xmin)/res)
;  ny = ceil((ymax-ymin)/res)
;  
;  
;  nw = (m[0].dim)[0]
;  
;  flnumerator=dblarr(nw,nx,ny)
;  varnumerator=dblarr(nw,nx,ny)
;  denominator=dblarr(nw,nx,ny)
;  
;  c = m[0]
;  *c.slitvec = indgen(nx)*res + xmin
;  *c.stepvec = indgen(ny)*res + ymin
;
;  ; throw out last map
;  for ei = 0,(n-1 -1) do begin
;    print, "Crunching exposure " + strtrim(string(ei+1),2)+ " of " + strtrim(string(n),2)
;    tic
;    for xi =0,(nx-1) do begin      
;      mx = ((xi*res+xmin) lt (*m[ei].xpos) ) and ((*m[ei].xpos) le ((xi+1)*res+xmin))
;      for yi = 0,(ny-1) do begin
;        my = ((yi*res+ymin) lt (*m[ei].ypos) ) and ((*m[ei].ypos) le ((yi+1)*res+ymin))
;        mask = where(mx and my)
;        for wi = 0,(nw-1) do begin
;          flnumerator[wi,xi,yi] = flnumerator[wi,xi,yi] + total(((reform((*m[ei].und)[wi,*,*]))[mask]) * ((reform((*m[ei].err)[wi,*,*]))[mask]))
;          varnumerator[wi,xi,yi] = varnumerator[wi,xi,yi] + total( ((reform((*m[ei].err)[wi,*,*]))[mask])^4  )
;          denominator[wi,xi,yi] = denominator[wi,xi,yi]+ total( ((reform((*m[ei].err)[wi,*,*]))[mask])  )
;        endfor
;      endfor
;      
;      time=toc()
;      timeleft = string( time*nx / (xi+1.) - time,FORMAT='(%"%3d")')
;      print,timeleft+'s',FORMAT='(%"\b\b\b\b%s",$)'
;
;    endfor
;    print,""
;
;  endfor
;  
;  
;  
;  
;  *c.und = (flnumerator / denominator) * (denominator ne 0)
;  *c.err = ((varnumerator^(0.5)) / denominator) * (denominator ne 0)
;  
;  
;
;  return, c
;end

function exescubecoadd, m, res

  n = m.length


  xcollapse = pointercollapse(m.xpos)
  ycollapse = pointercollapse(m.ypos)



  xmin= min(xcollapse,max=xmax)
  ymin= min(ycollapse,max=ymax)
  nx = ceil((xmax-xmin)/res)
  ny = ceil((ymax-ymin)/res)
  
  wmin = (*m[0].wavevec)[0]
  ;wres = wresp*((*m[0].wavevec)[1]-wmin)
  nw = (m[0].dim)[0]
  
  c = m[0]
  *c.slitvec = indgen(nx)*res + xmin
  *c.stepvec = indgen(ny)*res + ymin

  c.dim = [nw,nx,ny]
  
  newind = lindgen(c.dim)
  
  message, "Indexing...", /informational
  for ei = 0,(n-1) do begin ;loop over each exposure
    
    *m[ei].voxelind = lonarr(m[ei].dim)
    
    
    for syi = 0, ((m[ei].dim)[2]-1) do begin
      for sxi = 0, ((m[ei].dim)[1]-1) do begin
        for wi = 0,(nw-1) do begin
          (*m[ei].voxelind)[wi,sxi,syi] = newind(wi,floor(((*m[ei].xpos)[sxi,syi]-xmin)/res),floor(((*m[ei].ypos)[sxi,syi]-ymin)/res))
        endfor
      endfor
    endfor
    
;    
;    
;    for i = 0, (*m[ei].voxelind).length do begin
;      
;      (*m[ei].voxelind)[i] = newind(,floor(((*m[ei].xpos)[i]-xmin)/res),floor(((*m[ei].ypos)[i]-ymin)/res))
;      
;    endfor
;    
;    
;    for xi =0,(nx-1) do begin
;      mx = ((xi*res+xmin) lt (*m[ei].xpos) ) and ((*m[ei].xpos) le ((xi+1)*res+xmin))
;      for yi = 0,(ny-1) do begin
;        my = ((yi*res+ymin) lt (*m[ei].ypos) ) and ((*m[ei].ypos) le ((yi+1)*res+ymin))
;        mask = where(mx and my)
;        for wi = 0,(nw-1) do begin
;          flnumerator[wi,xi,yi] = flnumerator[wi,xi,yi] + total(((reform((*m[ei].und)[wi,*,*]))[mask]) * ((reform((*m[ei].err)[wi,*,*]))[mask]))
;          varnumerator[wi,xi,yi] = varnumerator[wi,xi,yi] + total( ((reform((*m[ei].err)[wi,*,*]))[mask])^4  )
;          denominator[wi,xi,yi] = denominator[wi,xi,yi]+ total( ((reform((*m[ei].err)[wi,*,*]))[mask])  )
;        endfor
;      endfor
;
;
;    endfor
;
  endfor
  
  
  flnumerator=dblarr(nw,nx,ny)
  varnumerator=dblarr(nw,nx,ny)
  denominator=dblarr(nw,nx,ny)

  for ei = 0,(n-1) do begin
    message, "Coadding map " + strtrim(string(ei+1),2)+ " of " + strtrim(string(n),2), /informational
    for j = 0, ((*m[ei].voxelind).length-1) do begin
          flnumerator[(*m[ei].voxelind)[j]] = flnumerator[(*m[ei].voxelind)[j]] + (*m[ei].und)[j] * ((*m[ei].err)[j])^(-2)
          ;varnumerator[(*m[ei].voxelind)[j]] = varnumerator[(*m[ei].voxelind)[j]] + ((*m[ei].err)[j])^4  
          denominator[(*m[ei].voxelind)[j]] = denominator[(*m[ei].voxelind)[j]] + ((*m[ei].err)[j])^(-2)
    endfor
  endfor



;  *c.und = (flnumerator / denominator) * (denominator ne 0)
;  *c.err = ((varnumerator^(0.5)) / denominator) * (denominator ne 0)
  *c.und = (flnumerator / denominator) * (denominator ne 0)
  *c.err = (denominator)^(-0.5) * (denominator ne 0)



  return, c
end

function exescubecoadd3, m, res

  n = m.length


  xcollapse = pointercollapse(m.xpos)
  ycollapse = pointercollapse(m.ypos)



  xmin= min(xcollapse,max=xmax)
  ymin= min(ycollapse,max=ymax)
  nx = ceil((xmax-xmin)/res)
  ny = ceil((ymax-ymin)/res)

  wmin = (*m[0].wavevec)[0]
  ;wres = wresp*((*m[0].wavevec)[1]-wmin)
  nw = (m[0].dim)[0]

  c = m[0]
  *c.slitvec = indgen(nx)*res + xmin
  *c.stepvec = indgen(ny)*res + ymin

  c.dim = [nw,nx,ny]

  newind = lindgen(c.dim)

  message, "Indexing...", /informational
  for ei = 0,(n-1) do begin ;loop over each exposure

    *m[ei].voxelind = lonarr(m[ei].dim)


    for syi = 0, ((m[ei].dim)[2]-1) do begin
      for sxi = 0, ((m[ei].dim)[1]-1) do begin
        for wi = 0,(nw-1) do begin
          (*m[ei].voxelind)[wi,sxi,syi] = newind(wi,floor(((*m[ei].xpos)[sxi,syi]-xmin)/res),floor(((*m[ei].ypos)[sxi,syi]-ymin)/res))
        endfor
      endfor
    endfor

  endfor


  flnumerator=dblarr(nw,nx,ny)
  varnumerator=dblarr(nw,nx,ny)
  denominator=dblarr(nw,nx,ny)

  for ei = 0,(n-1) do begin
    message, "Coadding map " + strtrim(string(ei+1),2)+ " of " + strtrim(string(n),2), /informational
    for j = 0, ((*m[ei].voxelind).length-1) do begin
      flnumerator[(*m[ei].voxelind)[j]] = flnumerator[(*m[ei].voxelind)[j]] + (*m[ei].und)[j] * ((*m[ei].err)[j])^(-2)
      denominator[(*m[ei].voxelind)[j]] = denominator[(*m[ei].voxelind)[j]] + ((*m[ei].err)[j])^(-2)
    endfor
  endfor



  *c.und = (flnumerator / denominator) * (denominator ne 0)
  *c.err = (denominator)^(-0.5) * (denominator ne 0)



  return, c
end

function exescubeinit, m, wavebound, res, skyoption
  
  exesskysubtract, m, skyoption
  
  
  
  exescubetrim, m, wavebound

  exescubeplace, m

  c=exescubecoadd(m,res)

  return, c
end

pro exesskysubtract, m, option
  n = m.length
  case option of
    0:        ; No subtraction
    1: begin  ; Median of 3 sky frames
      for i = 0,(n-1) do begin
        sky = median((*m[i].sky), dimension=3)
        for j = 0,((m[i].dim)[2]-1) do begin
          (*m[i].und)[*,*,j] = reform((*m[i].und)[*,*,j]) -sky
        endfor
      endfor
    end
    2: begin  ; Proportional average of start and end frames
      for i = 0,(n-1) do begin
        
        fn = (m[i].dim)[2]
        
        skyA = reform((*m[i].und)[*,*,0])
        skyB = reform((*m[i].und)[*,*,fn-1])
        
        for j = 0,(fn-1) do begin
          (*m[i].und)[*,*,j] = reform((*m[i].und)[*,*,j]) - (skyA*(fn-1-j)+skyB*(j))/(fn-1)
        endfor
      endfor
    end
    
  endcase
  
  
end

function m82init
  ; first map 64
  ; maps 67 to 76 have the same parameters
  ;

  dir = "~/Research/EXES/M82/"

  ;  if file_test(dir+'coaddcache.sav') then begin
  ;    restore, dir+'coaddcache.sav'
  ;  endif else begin
  undfilenames = ["F0363_EX_SPE_7500171_NONEEXEECHL_UND_0067.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0068.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0069.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0070.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0071.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0072.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0073.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0074.fits", $
    "F0363_EX_SPE_7500171_NONEEXEECHL_UND_0075.fits"]
  ;"F0363_EX_SPE_7500171_NONEEXEECHL_UND_0076.fits"]
  wvmfilename = "F0363_EX_SPE_7500171_NONEEXEECHL_WVM_0079-0082.fits"

  m82 = multimapinit(dir,undfilenames,wvmfilename)
  

  m82[*].object = 'Messier 82'
  m82[*].feature = 'Argon II'

  wavebound = [333,460]
  res = 2;1.5 ;arcsec  ;m[0].slitwid
  skyoption = 2

  coM82=exescubeinit(m82, wavebound, res,skyoption)


  return, coM82

  
end

function jupinit
  ; first map 64
  ; maps 67 to 76 have the same parameters
  ;

  dir = "/home/cordell/Research/jupiter/data/Repiped/"

  ;  if file_test(dir+'coaddcache.sav') then begin
  ;    restore, dir+'coaddcache.sav'
  ;  endif else begin
  undfilenames = ["F0158_EX_SPE_86000210_EXEELONEXEECHL_UND_0051.fits"]
  wvmfilename = "F0158_EX_SPE_86000210_EXEELONEXEECHL_WVM_0050.fits"

  JupS0 = multimapinit(dir,undfilenames,wvmfilename)

  JupS0[*].object = 'Jupiter'
  JupS0[*].feature = 'Hydrogen S0'
  
;
;  wavebound = [333,460]
;
;  exescubetrim, m82, wavebound
;
;  exescubeplace, m82
;
;  ;    c=exescubecoadd(m82)
;  ;
;  ;    save, c , filename = dir+'coaddcache.sav'
;  ;  endelse
;
;
;

  return, JupS0


end

function coaddinit, m
    if file_test(m[0].dir+'coaddcache.sav') then begin
      restore, m[0].dir+'coaddcache.sav'
    endif else begin
      c=exescubecoadd(m)
  
      save, c , filename = m[0].dir+'coaddcache.sav'
    endelse
    return, c

end

;pro multimapA
;    
;  
;  m = m82init()
;  
;
;  
;  n = m.length
;  ds=1
;  if !cpu.HW_NCPU ge 8 then ds=2
;  !p.charsize = 1.5*ds
;
;
;  window, 0, xs=ds*550,ys=ds*1000
;  for i = 0,(n-1) do begin
;    !p.multi=[i+1,2,(n+1)/2]
;    f = total((*m[i].und),1)
;    xs=(*m[i].slitvec)
;    ys=abs(*m[i].stepvec)
;    aspect = (m[i].pltscale)/abs(m[i].step)
;    display, f, xs, ys, min=-4, max=7,aspect=aspect,/noerase,title=string(i)
;    ;wait, 1
;  endfor
;  !p.multi=0
;
;  write_png, m[0].dir+'m82totals.png',tvrd(/true)
;  ;print, m
;  
;  i=4
;  
;  s = reform((*m[i].und)[*,*,18])
;  
;  resolve_routine, 'jupmaps',/compile_full_file
;  
;  sr= shrinkbin(s,[16,30])
;  wv = shrinkbin((*m[i].wavevec),[16])
;  
;  pl=rainbowplot(wv,sr[*,10:20],/cbar,thick=3*ds,xtit='Wavenumber',ytit='Intensity',title='Single Slit, Binned')
;
;  
;  
;  end
;pro multimapC
;
;
;  m = m82init()
;
;  
;  farr = total((*m[0].und),1)
;  xarr = (*m[0].xpos)
;  yarr = (*m[0].ypos)
;  
;  for i = 1,8 do begin
;    farr = [farr,total((*m[i].und),1)]
;    xarr = [xarr,(*m[i].xpos)]
;    yarr = [yarr,(*m[i].ypos)]
;    
;  endfor
;  
;  
;  ds=1
;  if !cpu.HW_NCPU ge 8 then ds=2
;  !p.charsize = 1.5*ds
;
;  
;  
;  
;  window, 0, xs=ds*550,ys=ds*1000
;  !p.multi=[0,1,2]
;  
;  
;  r = min_curve_surf(farr,xarr,yarr)
;  
;  contour,r
;  
;  
;  ;display, f, min=-1, max=1, xs, ys, title="Coadd"
;  ;!p.multi=[1,1,2]
;
;
;
;
;
;  !p.multi=0
;
;  ;  write_png, m[0].dir+'m82totals.png',tvrd(/true)
;  ;  ;print, m
;  ;
;  ;  i=4
;  ;
;  ;  s = reform((*m[i].und)[*,*,18])
;  ;
;  ;  resolve_routine, 'jupmaps',/compile_full_file
;  ;
;  ;  sr= shrinkbin(s,[16,30])
;  ;  wv = shrinkbin((*m[i].wavevec),[16])
;  ;
;  ;  pl=rainbowplot(wv,sr[*,10:20],/cbar,thick=3*ds,xtit='Wavenumber',ytit='Intensity',title='Single Slit, Binned')
;
;  ;mwrfits, (*c.und), c.dir+'coadded.fits',/create
;
;
;end

function exesaddstack, c, outstack

  wavemax=max((*c.wavevec),min=wavemin)

;  f = total((*c.und)/(*c.err)^2,1) / (total(1/(*c.err)^2,1)*(wavemax-wavemin))
  f = total((*c.und),1) / ((wavemax-wavemin))

  if n_params() gt 1 then begin
    e = (total((*c.err)^2,1))^(0.5)
    outstack = dblarr([f.dim,2])
    outstack[*,*,0] = f
    outstack[*,*,1] = e
  endif
  
  return, f
end

function multimapauto, _EXTRA=ex

  dir = "./" ; Working directory
  
  undfilenames = file_search(dir,'*_UND_*.fits')

  wvmfilename = (file_search(dir,'*_WVM_*.fits'))[0]

  cubes = multimapinit(dir,undfilenames,wvmfilename)


;  cubes[*].object = 'Messier 82'
;  cubes[*].feature = 'Argon II'

  wavebound = [333,460]
  res = 2;1.5 ;arcsec  ;m[0].slitwid
  skyoption = 2

  coaddedcube=exescubeinit(cubes, wavebound, res,skyoption)


  return, coaddedcube
end

pro multimapM82


  c = m82init()

  ds=1
  if !cpu.HW_NCPU ge 8 then ds=2
  !p.charsize = 1.5*ds
  !p.charthick = 1.5*ds
  
  
  fmin=-1
  fmax=2

  f = exesaddstack(c,fstack)

  window, 1, xs=ds*650,ys=ds*1000
  !p.COLOR = '000000'xL
  !p.BACKGROUND = 'FFFFFF'xL
  !p.thick = ds
  !p.multi=[0,1,2]
    
  ;f = rotate(f,1)
  xs=(*c.slitvec)
  ys=(*c.stepvec)
  display, f, min=fmin, max=fmax, xs, ys, aspect=1, title="Coadd",xtit='Arcsec',ytit='Arcsec'

  intensity = 'Intensity (erg/s/cm^2/sr)'
  cgcolorbar,title=intensity,minrange=fmin,maxrange=fmax,/vertical,/fit,/right
  !p.multi=[1,1,2]



  contour, f, xs, ys, levels=0.01*indgen(12)^2,    xrange=[-20,31],yrange=[-16,19], xstyle=1,ystyle=1, title='M82 ArII',/isotropic,thick=ds,xthick=ds,ythick=ds,xtit='Arcsec',ytit='Arcsec'




  !p.multi=0

  ;  write_png, m[0].dir+'m82totals.png',tvrd(/true)
  ;  ;print, m
  ;
  ;  i=4
  ;
  ;  s = reform((*m[i].und)[*,*,18])
  ;
  ;  resolve_routine, 'jupmaps',/compile_full_file
  ;
  ;  sr= shrinkbin(s,[16,30])
  ;  wv = shrinkbin((*m[i].wavevec),[16])
  ;
  ;  pl=rainbowplot(wv,sr[*,10:20],/cbar,thick=3*ds,xtit='Wavenumber',ytit='Intensity',title='Single Slit, Binned')

  arcsectodeg=2.778e-4 ;deg/arcsec
  !null=min(abs(xs),xmloc)
  !null=min(abs(ys),ymloc)
  res = xs[1]-xs[0]
    
  mkhdr, outhdr, fstack
  make_astr, astr, crpix=[xmloc,ymloc]+1,crval=[c.ra,c.dec],delt=[res,res]*arcsectodeg
  putast, outhdr, astr

  
  mwrfits, (*c.und), c.dir+'m82_coadded_cube.fits',/create
  mwrfits, fstack, c.dir+'m82_integrated.fits',outhdr,/create
  
  write_png, c.dir+'m82contour.png',tvrd(/true)

  


end

pro multimapdisplay, c

  
  ds=1
  if !cpu.HW_NCPU ge 8 then ds=2
  
  !p.charsize = 1.5*ds
  !p.charthick = 1.5*ds


  fmin=-1
  fmax=2

  f = exesaddstack(c,fstack)

  window, 1, xs=ds*650,ys=ds*1000
  !p.COLOR = '000000'xL
  !p.BACKGROUND = 'FFFFFF'xL
  !p.thick = ds
  !p.multi=[0,1,2]

  ;f = rotate(f,1)
  xs=(*c.slitvec)
  ys=(*c.stepvec)
  display, f, min=fmin, max=fmax, xs, ys, aspect=1, title="Coadd",xtit='Arcsec',ytit='Arcsec'

  intensity = 'Intensity (erg/s/cm^2/sr)'
  cgcolorbar,title=intensity,minrange=fmin,maxrange=fmax,/vertical,/fit,/right
  !p.multi=[1,1,2]



  contour, f, xs, ys, levels=0.01*indgen(12)^2,    xrange=[-20,31],yrange=[-16,19], xstyle=1,ystyle=1, title='M82 ArII',/isotropic,thick=ds,xthick=ds,ythick=ds,xtit='Arcsec',ytit='Arcsec'
  
  !p.multi=0

end

pro multimap,_EXTRA=ex

  ;
  ; Optional Keywords
  ; 
  ; DIRECTORY = directory
  ;   Default: present working directory
  ;   Type: string, e.g. '/home/joe/data/andromeda/' or './andromeda/'
  ;   Input data directory containing UND map files to be coadded, as well as a corresponding WVM flat file, as processed from the EXES pipeline.
  ; 
  ; OUTDIR = outdir
  ;   Default: "Multimap_Output" directory inside DIRECTORY
  ;   Type: string
  ;   Location to save the coadded results as FITS files. Outdir is created if it doesn't already exist.
  ;   
  ; FILENAME = filename
  ;   Default: Object name + '_Coadded_Cube'  
  ;   Type: String
  ;   
  ; WAVEBOUND = wavebound
  ;   Default: full input range
  ;   
  ;   
  ; 
  ; 
  ;
  
  
  
  
  c = multimapauto(_EXTRA=ex)
  
  ds=1
;  if !cpu.HW_NCPU ge 8 then ds=2
  
  !p.charsize = 1.5*ds
  !p.charthick = 1.5*ds


  fmin=-1
  fmax=2

  f = exesaddstack(c,fstack)

  window, 1, xs=ds*650,ys=ds*1000
  !p.COLOR = '000000'xL
  !p.BACKGROUND = 'FFFFFF'xL
  !p.thick = ds
  !p.multi=[0,1,2]

  ;f = rotate(f,1)
  xs=(*c.slitvec)
  ys=(*c.stepvec)
  display, f, min=fmin, max=fmax, xs, ys, aspect=1, title="Coadd",xtit='Arcsec',ytit='Arcsec'

  intensity = 'Intensity (erg/s/cm^2/sr)'
  cgcolorbar,title=intensity,minrange=fmin,maxrange=fmax,/vertical,/fit,/right
  !p.multi=[1,1,2]



  contour, f, xs, ys, levels=0.01*indgen(12)^2,    xrange=[-20,31],yrange=[-16,19], xstyle=1,ystyle=1, title='M82 ArII',/isotropic,thick=ds,xthick=ds,ythick=ds,xtit='Arcsec',ytit='Arcsec'
  
  !p.multi=0
  
  if not file_test(c.dir+'Multimap_Output',/directory) then begin
    file_mkdir, c.dir+'Multimap_Output'
  endif
  
  outdir = c.dir+'Multimap_Output/'
  outname = strcompress(c.object,/remove_all)
  


  ;  write_png, m[0].dir+'m82totals.png',tvrd(/true)
  ;  ;print, m
  ;
  ;  i=4
  ;
  ;  s = reform((*m[i].und)[*,*,18])
  ;
  ;  resolve_routine, 'jupmaps',/compile_full_file
  ;
  ;  sr= shrinkbin(s,[16,30])
  ;  wv = shrinkbin((*m[i].wavevec),[16])
  ;
  ;  pl=rainbowplot(wv,sr[*,10:20],/cbar,thick=3*ds,xtit='Wavenumber',ytit='Intensity',title='Single Slit, Binned')

  arcsectodeg=2.778e-4 ;deg/arcsec
  !null=min(abs(xs),xmloc)
  !null=min(abs(ys),ymloc)
  res = xs[1]-xs[0]

  mkhdr, outhdr, fstack
  make_astr, astr, crpix=[xmloc,ymloc]+1,crval=[c.ra,c.dec],delt=[res,res]*arcsectodeg
  putast, outhdr, astr


  mwrfits, (*c.und), outdir+outname+'_Coadded_Cube.fits',/create
  mwrfits, fstack, outdir+outname+'_Integrated.fits',outhdr,/create

;  write_png, +'m82contour.png',tvrd(/true)




end




