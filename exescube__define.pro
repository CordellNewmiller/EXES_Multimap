
;function exescube::init
;    
;    
;    self.und = ptr_new(/ALLOCATE_HEAP)
;    self.err = ptr_new(/ALLOCATE_HEAP)
;    self.wvm = ptr_new(/ALLOCATE_HEAP)
;
;
;
;    self.xvec = ptr_new(/ALLOCATE_HEAP)
;    self.yvec = ptr_new(/ALLOCATE_HEAP)
;
;    self.wavevec = ptr_new(/ALLOCATE_HEAP)
;  
;  print, 'hello world'
;  
;end


function exescube::hello
  print, 'hello world'
end

pro exescube__define

  void = {EXEScube, $
    object: '', $ ; e.g. Jupiter
    feature: '', $ e.g. S(0)
    dir: '', $ e.g. ~/Research/jupiter/data/
    undfilename: '', $
    wvmfilename:'', $
    instcfg:'', $

    wn0:0D, $
    ra:0D, $
    dec:0D, $
    angle:0D, $
    step:0D, $
    pltscale:0D,$
    slitwid:0D,$

    und:ptr_new(/ALLOCATE_HEAP), $  science frames [spectral, slit, step]
    err:ptr_new(/ALLOCATE_HEAP), $  uncertainty frames [spectral, slit, step]
    wvm:ptr_new(/ALLOCATE_HEAP), $  wavemap [spectral, slit]
    slit:ptr_new(/ALLOCATE_HEAP), $ 
    sky:ptr_new(/ALLOCATE_HEAP), $  sky frames [spectral, slit, frame]

    dim:[0,0,0], $


    slitvec:ptr_new(/ALLOCATE_HEAP), $
    stepvec:ptr_new(/ALLOCATE_HEAP), $


    orderB:ptr_new(/ALLOCATE_HEAP), $
    orderE:ptr_new(/ALLOCATE_HEAP), $
    orderS:ptr_new(/ALLOCATE_HEAP), $
    orderT:ptr_new(/ALLOCATE_HEAP), $


    wavevec:ptr_new(/ALLOCATE_HEAP), $ ; vector
      
    xpos:ptr_new(/ALLOCATE_HEAP), $
    ypos:ptr_new(/ALLOCATE_HEAP), $
    
    voxelind:ptr_new(/ALLOCATE_HEAP) $


      

  }
end