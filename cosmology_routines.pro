pro cosmology_routines, _extra=extra
return
end
; ----------------------------------------------------------------------
; units conversion programs
; ----------------------------------------------------------------------
function convert_distance, cm=cm, meter=meter, kpc=kpc, Mpc=Mpc, Gpc=Gpc

    conversion = 1.0D ; default [pc]

; distance units
    
    if keyword_set(cm) then conversion = 3.086D18     ; [cm/pc]    
    if keyword_set(meter) then conversion = 3.086D16  ; [m/pc]    
    if keyword_set(kpc) then conversion = 1.0D-3      ; [kpc/pc]
    if keyword_set(Mpc) then conversion = 1.0D-6      ; [Mpc/pc]   
    if keyword_set(Gpc) then conversion = 1.0D-9      ; [Gpc/pc]   

return, conversion
end

function convert_time, s=s, yr=yr, Myr=Myr

    conversion = 1.0D ; default [Gyr]

; time units
    
    secyr = 365.25*24.0*3600D                       ; [s/yr]
    if keyword_set(s) then conversion = 1.0D9*secyr ; [s/Gyr]
    if keyword_set(yr) then conversion = 1.0D9      ; [yr/Gyr]
    if keyword_set(Myr) then conversion = 1.0D3     ; [Myr/Gyr]
    
return, conversion
end
; ----------------------------------------------------------------------
; numerical integration routine from C. Markwardt
; ----------------------------------------------------------------------

; ----------------------------------------------------------------------
; cosmological parameter access programs
; ----------------------------------------------------------------------
FUNCTION redomegal
   COMMON cosmology
   return, olh.omegalambda
END 
FUNCTION redomega0
   COMMON cosmology
   return, olh.omega0
END 
FUNCTION redh100
   COMMON cosmology
   return, olh.h100
END 
; ----------------------------------------------------------------------
; functionality programs
; ----------------------------------------------------------------------
function asinh, x
    x = double(x)
    return, alog(x+sqrt(1.d0+x^2))
end 
function cube, x 
    return, x*x*x
end 
function cuberoot, x
    return, double(x)^(1.0/3.0d0)
end 
function sqr, x
    return, x*x
end    
function epeebles, z
; E(z) is the time derivative of the logarithm of the scale factor
; a(t), which is proportional to (1+z); ie (1+z)*e(z) is proportional
; to the derivative of z wrt the lookback time (Peebles, p. 312, Hogg
; Eq. 13)
    common cosmology

    omegaR = 1.0d0 - olh.omega0 - olh.omegalambda ; curvature density parameter
    ez = sqrt(olh.omega0*cube(1.0+z)+omegaR*sqr(1.0+z)+olh.omegalambda)
    return, ez

end 

; ----------------------------------------------------------------------
; cosmological measurements
; ----------------------------------------------------------------------
function dhubble, _extra=extra
; return the Hubble distance in pc (Hogg Eq. 4)
    common cosmology

    speedoflight = 2.99792458d5 ; (km/s)
    dh = speedoflight / olh.h100 / 100.d0 ; Mpc
    conv = convert_distance(_extra=extra)
    return, conv * dh * 1.0d6
    
end 
function thubble, _extra=extra
; return the Hubble time in Gyr; recall, 1./(1 km/s/mpc) = 9.77813 Gyr
; (Hogg Eq. 3)
    common cosmology
    conv = convert_time(_extra=extra)
    return, conv * 9.77813d0 / olh.h100
end

function dangular, z, _extra=extra
; angular diameter distance (Hogg Eq. 17) (pc)
    common cosmology
    conv = convert_distance(_extra=extra)
    return, conv * dcomovingtransverse(z) / (1. + z)
end 
function dangulardiff, z1, z2, _extra=extra
; angular diameter distance difference between two objects at
; redshifts z1 and z2 (Hogg Eq. 18)
    common cosmology

    omegaR = 1.0d0 - olh.omega0 - olh.omegalambda ; curvature density parameter

    conv = convert_distance(_extra=extra)
    if omegaR ge 0.0 then begin 
       dm1 = dcomovingtransverse(z1)
       dm2 = dcomovingtransverse(z2)
       dh = dhubble()
       da12 = (1./(1.+z2)) * $
         ( dm2 * sqrt(1.+omegaR*sqr(dm1)/sqr(dh)) - $
           dm1 * sqrt(1.+omegaR*sqr(dm2)/sqr(dh)) )
       return, conv * da12
    endif else begin 
       message,'No formula available for a negative curvature density parameter.'
       return, -1.0
    endelse 

end 
function dlosfunc, z
; construct the appropriate function based on peeble's E(z) to
; determine the total line-of-sight comoving distance
    common cosmology
    return,1.0d0/(epeebles(z))
end 
function dcomovinglos, z, _extra=extra
; total line-of-sight comoving distance in pc (Hogg Eq. 14)
    common cosmology
    conv = convert_distance(_extra=extra)
    return, conv * dhubble() * qromb('dlosfunc',0.d0,z,/double)
;   return, conv * dhubble() * qpint1d('dlosfunc',0.d0,z)
end 
function dcomovingtransverse, z, _extra=extra
; transverse comoving distance in pc (Hogg Eq. 15) (pc)
    common cosmology

    conv = convert_distance(_extra=extra)
    omegaR = 1.0d0 - olh.omega0 - olh.omegalambda ; curvature density parameter

    if omegaR gt 0.0 then $
      dct = dhubble() / sqrt(omegaR) * $
      sinh( sqrt(omegaR) * dcomovinglos(z) / $
            dhubble() ) else $
      if omegaR eq 0.0 then $
      dct = dcomovinglos(z) else $
      if omegaR lt 0.0 then $
      dct = dhubble() / sqrt(abs(omegaR)) * $
      sin( sqrt(abs(omegaR)) * dcomovinglos(z) / $
           dhubble() )
    return, conv * dct

end 
function dluminosity, z, _extra=extra
 ; luminosity distance in pc (Hogg Eq. 20)
    common cosmology
    conv = convert_distance(_extra=extra)
    return, conv * dcomovingtransverse(z) * (1.0 + z)
end 
function dmodulus, z, _extra=extra
; distance modulus (Hogg Eq. 24)
    common cosmology
    conv = convert_distance(_extra=extra)
    return, conv * 5.0d0 * alog10(dluminosity(z,/pc) / 10.0d0)
end 
function dvcomoving, z, _extra=extra
; differential comoving volume.  returns the comoving volume element
; per steradian, per redshift interval dz (Hogg Eq. 27) (pc^3)
    common cosmology
    conv = convert_distance(_extra=extra)
    return, conv^3.0 * dhubble() * sqr(1.+z) * sqr(dangular(z)) / epeebles(z)
end 
function vcomoving, z, _extra=extra
; total comoving volume, all sky, to redshift z (obtained by
; integrating the comoving volume element from the present to redshift
; z).  (Hogg Eq. 28)  c.f. luminosity densities, number counts, space
; densities, etc.
    common cosmology
    
    omegaR = 1.0d0 - olh.omega0 - olh.omegalambda ; curvature density parameter
    dh = dhubble()
    dm = dcomovingtransverse(z)

    if omegaR gt 0.0 then $
      vc = (4.*!dpi*cube(dh))/(2.*omegaR) * $
      ( dm/dh*sqrt(1.+omegaR*sqr(dm)/sqr(dh)) - $
        1./sqrt(abs(omegaR))*asinh(sqrt(abs(omegaR))*dm/dh) ) else $
      if omegaR eq 0.0 then $
      vc = (4.*!dpi/3.)*cube(dm) else $
      if omegaR lt 0.0 then $
      vc = (4.*!dpi*cube(dh))/(2.*omegaR) * $
      ( dm/dh*sqrt(1.+omegaR*sqr(dm)/sqr(dh)) - $
        1./sqrt(abs(omegaR))*asin(sqrt(abs(omegaR))*dm/dh) )
    
    conv = convert_distance(_extra=extra)

return, conv^3.0 * vc
end 
function agefunc, z
; construct the appropriate function based on peeble's E(z) for
; the integration of the age of the universe at some redshift
    common cosmology
    return, 1.0d0/(epeebles(z)*(1.+z))
end 
function getage, z, _extra=extra
; lookback time to an object (the difference between the age of the
; universe now and the age of the universe at the time the photons
; were emitted according to that object) (Hogg Eq. 29)
    common cosmology
    conv = convert_time(_extra=extra)
;   return, conv * thubble() * qromb('agefunc',z,1000.d0,/double)
    nz = n_elements(z) & age = fltarr(nz)
    for i = 0L, nz-1L do age[i] = conv*thubble()*qromb('agefunc',z[i],1000D0,/double)
;   for i = 0L, nz-1L do age[i] = conv*thubble()*qpint1d('agefunc',z[i],!values.f_infinity)
    if (nz eq 1L) then age = age[0]
    return, age
end    
function getredshift, age, _extra=extra
; jm06apr25uofa - given an age, return the redshift; note: linear
;                 interpolation works better than spline  

    common cosmology

    conv = convert_time(_extra=extra)

    zmax = 10.0 & zmin = 1E-2 & dlogz = 0.1
    biglogz = findgen((alog10(zmax)-alog10(zmin))/dlogz+1)*dlogz+alog10(zmin)
    bigz = 10.0^biglogz
    bigage = getage(bigz)

return, conv*(interpol(bigz,bigage,age)>0)
end
;function getredshift, age
;; adapted with permission from ignacio ferreras' c code
;; /* secant method used to compute the root of age-get_age */
;; /* see numerical recipes cv2 pg.354 for details */
;    common cosmology
;
;    tu = getage(0.0)
;
;    if age gt tu then return, -1.0
;    if (tu - age) lt 1.0e-4 then return, 0.0
;    z = 0.1
;    t = getage(z)
;    while age lt t do begin 
;       t = getage(z)
;       z = z+0.1
;    endwhile 
;    x2 = z
;    x1 = z - 0.1
;    fl = age - getage(x1)
;    f = age - getage(x2)
;    if (abs(fl) lt abs(f)) then begin 
;       rts = x1
;       xl = x2
;       swap = fl
;       fl = f
;       f = swap
;    endif else begin 
;       xl = x1
;       rts = x2
;    endelse 
;    for j=1,50 do begin 
;       dx = (xl-rts)*f/(f-fl)
;       xl = rts
;       fl = f
;       rts = rts + dx
;       f = age - getage(rts)
;       if (abs(dx) lt 0.001 or abs(f) lt 1.0e-7) then return, rts 
;    endfor 
;    print, 'Convergence error finding z for age = '+string(age)+'.'
;    return, -1.0
;end 
