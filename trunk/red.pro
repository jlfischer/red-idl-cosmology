;+
; NAME:
;   RED
;
; PURPOSE:
; compute cosmographic quantities...
; 
; Available cosmological routines:
;  
;    DHUBBLE()              - Hubble distance
;    THUBBLE()              - Hubble time
;    DANGULAR(z)            - angular diameter distance
;    DANGULARDIFF(z1,z2)    - difference in angular diameter
;    DCOMOVINGLOS(z)        - line-of-sight co-moving distance
;    DCOMOVINGTRANSVERSE(z) - transverse co-moving distance
;    DLUMINOSITY(z)         - luminosity distance
;    DMODULUS(z)            - distance modulus
;    DVCOMOVING(z)          - differential co-moving volume
;    VCOMOVING(z)           - integrated co-moving volume
;    GETAGE(z)              - look-back time
;    GETREDSHIFT(age)       - redshift for a given age
;  
;    REDOMEGAL()            - return current omegalambda
;    REDOMEGA0()            - return current omega0
;    REDh100()              - return current h100
;  
; Available units:
;    cm    - centimeter          s    - second
;    meter - meter               yr   - year
;    pc    - parsec [default]    Myr  - megayear
;    kpc   - kiloparsec          Gyr  - gigayear [default]
;    Mpc   - megaparsec
; ----------------------------------------------------------------
; 
;
; CALLING SEQUENCE:
;  red, z=z, omega0=omega0, omegalambda=omegalambda, h100=h100, $
;       help=help, default=default, verbose=verbose, _extra=extra
;
; INPUTS:
; 
;
; OPTIONAL INPUTS:
;
;	
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; COMMENTS:
;
;
; EXAMPLE:
;   red,/default,z=1.0,/verbose
;
; PROCEDURES USED:
;	COSMOLOGY_ROUTINES
;
; MODIFICATION HISTORY:
;	Moustakas & Moustakas 2002
;-

pro red, z=z, omega0=omega0, omegalambda=omegalambda, h100=h100, $
         help=help, default=default, verbose=verbose, _extra=extra

    common cosmology, olh

    forward_function asinh, cube, cuberoot, epeebles, dangular, dangulardiff, $
     dlosfunc, dcomovinglos, dcomovingtransverse, dhubble, thubble, $
     dluminosity,dmodulus, dvcomoving, vcomoving, $
     agefunc, getage, getredshift

; if no cosmology has been defined then define it
    
    if n_elements(olh) eq 0L then begin

       resolve_routine, 'cosmology_routines', /compile_full_file

       IF n_elements(default) NE 0 THEN BEGIN 
          omega0 = 0.3D
          omegalambda = 0.7D
          h100 = 0.70D
       ENDIF 

       if n_elements(omega0) EQ 0 then omega0 = 0.3
       if n_elements(omegalambda) EQ 0 then omegalambda = 0.7
       if n_elements(h100) EQ 0 then h100 = 0.70
       
       olh = {name: 'Cosmology Parameters', $ ; cosmology structure
              omega0: omega0, $
              omegalambda: omegalambda, $
              h100: h100}

    endif else begin

       IF n_elements(default) NE 0 THEN BEGIN 
          omega0 = 0.3
          omegalambda = 0.7
          h100 = 0.70
       ENDIF 

;       if n_elements(omega0) EQ 0 then omega0 = 0.3D
;       if n_elements(omegalambda) EQ 0 then omegalambda = 0.7D
;       if n_elements(h100) EQ 0 then h100 = 0.70D
       
       IF n_elements(omega0)      NE 0 THEN olh.omega0 = omega0
       IF n_elements(omegalambda) NE 0 THEN olh.omegalambda = omegalambda
       IF n_elements(h100)        NE 0 THEN olh.h100 = h100       

    endelse

; default cosmology
    
    IF n_elements(default) NE 0 THEN BEGIN 
       omega0 = 0.3
       omegalambda = 0.7
       h100 = 0.70
    ENDIF 

    if keyword_set(verbose) then begin
       print, 'Omega Matter = '+$
        strcompress(string(olh.omega0,format='(F7.4)'),/remove)
       print, 'Omega Lambda = '+$
        strcompress(string(olh.omegalambda,format='(F7.4)'),/remove)
       print, 'H_0 / 100    = '+$
        strcompress(string(olh.h100,format='(F7.4)'),/remove)
    endif
    
; print to the screen all the available cosmology routines
    
    if keyword_set(help) then begin

       print,''
       print,'red, z=z, omega0=omega0, omegalambda=omegalambda, h100=h100, $'
       print,'help=help, default=default, verbose=verbose, _extra=extra'
       print,''
       print, '----------------------------------------------------------------'
       print, 'Available cosmological routines:'
       print, ' '
       print, '   DHUBBLE()              - Hubble distance'
       print, '   THUBBLE()              - Hubble time'
       print, '   DANGULAR(z)            - angular diameter distance'
       print, '   DANGULARDIFF(z1,z2)    - difference in angular diameter'
       print, '   DCOMOVINGLOS(z)        - line-of-sight co-moving distance'
       print, '   DCOMOVINGTRANSVERSE(z) - transverse co-moving distance'
       print, '   DLUMINOSITY(z)         - luminosity distance'
       print, '   DMODULUS(z)            - distance modulus'
       print, '   DVCOMOVING(z)          - differential co-moving volume'
       print, '   VCOMOVING(z)           - integrated co-moving volume'
       print, '   GETAGE(z)              - look-back time'
       print, '   GETREDSHIFT(age)       - redshift for a given age'
       print, ' '
       print, '   REDOMEGAL()            - return current omegalambda'
       print, '   REDOMEGA0()            - return current omega0'
       print, '   REDh100()              - return current h100'
       print, ' '
       print, 'Available units:'
       print, '   cm    - centimeter          s    - second'
       print, '   meter - meter               yr   - year'
       print, '   pc    - parsec [default]    Myr  - megayear'
       print, '   kpc   - kiloparsec          Gyr  - gigayear [default]'
       print, '   Mpc   - megaparsec'
       print, '----------------------------------------------------------------'

    endif

    nz = n_elements(z)
    if nz ne 0L then for i = 0L, nz-1L do begin

       agez0 = getage(0.0)
;      print, 'Cosmology (OLH) = '
       print, 'Age (z=0.000) = '+strcompress(agez0,/remove)+' Gyr'

       if z[i] gt 0.0 then begin

          agez = getage(z[i])
          dmod = dmodulus(z[i])
          dang = dangular(z[i],/mpc)
          dlum = dluminosity(z[i],/mpc)
          xscl = dangular(z[i],/kpc)/206265.0D

          print, 'Age (z='+strcompress(string(z,format='(F5.3)'),/remove)+') = '+$
            strcompress(agez,/remove)+' Gyr'
          print, 'Lookback time = '+strcompress(agez0-agez,/remove)+' Gyr'
          print, 'DModulus      = '+strcompress(dmod,/remove)+' mag'
          print, 'DAngular      = '+strcompress(dang,/remove)+' Mpc'
          print, 'DLuminosity   = '+strcompress(dlum,/remove)+' Mpc'
          print, 'Scale         = '+strcompress(xscl,/remove)+' kpc/arcsec'

       endif
       
    endfor

return 
end 
