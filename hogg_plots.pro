;+
; NAME:
;	HOGG_PLOTS
;
; PURPOSE:
;	Generate Figures 1-7 from Hogg 2000, astro-ph/9905116.
;
; CALLING SEQUENCE:
;                   hogg_plots,figure=figure,postscript=postscript
;
; INPUTS:
;        figure number to reproduce, default=fig 1
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
;                cosmology
;
; COMMENTS:
;
;
; EXAMPLE:
;         
;
; PROCEDURES USED:
;
;
; MODIFICATION HISTORY:
;
;-

pro hogg_plots,figure=figure,postscript=postscript

    IF n_elements(figure) EQ 0 THEN figure = 1
    
    IF n_elements(postscript) NE 0 THEN BEGIN 
       outps = 'red_hogg_fig'+strcompress(string(figure),/remove_all)+'.ps'
       set_plot,'ps'
       device,file=outps,xsize=7,ysize=7,/inches,/portrait,/times
       !p.font = 0              ; postscript fonts
       oldcharsize = !p.charsize
       !p.charsize = 1.2
    ENDIF ELSE window,0,xsize=500,ysize=500
    
    zmin = 0.0
    zmax = 5.0
    dz = 0.1

    zarray = dz*findgen((zmax-zmin+dz)/dz)
    nz = n_elements(zarray)

    omega_0 = [1.0,0.05,0.2]
    omega_lambda = [0.0,0.0,0.8]

    ncosmology = n_elements(omega_0)
    
    FOR i = 0L, ncosmology-1L DO BEGIN 
;       red, omega_0=omega_0[i], omega_lambda=omega_lambda[i];,/verb
       red, omega0=omega_0[i], omegalambda=omega_lambda[i],/verb
       dh = dhubble()           ; hubble distance
       th = thubble()           ; hubble time

                                ; the proper motion distance, aka the
                                ; comoving transverse distance
       IF figure EQ 1 THEN BEGIN 
          cdist = dcomovingtransverse(zarray)/dh
          if i eq 0L then plot, zarray, cdist, line=i, $
           xr=[0,5],yr=[0,3],/xst,/yst,$
           xtit='redshift z',$
           ytit='proper motion distance D!dM!n/D!dH' else $ 
           oplot, zarray, cdist, line=i

                                ; angular diameter distance
       ENDIF ELSE IF figure EQ 2 THEN BEGIN 
          cdist = dangular(zarray)/dh
          if i eq 0L then plot, zarray, cdist, line=i, $
           xr=[0,5],yr=[0,0.5],/xst,/yst,$
           xtit='redshift z',$
           ytit='angular diameter distance D!dA!n/D!dH' else $ 
           oplot, zarray, cdist, line=i
          
                                ; luminosity distance
       ENDIF ELSE IF figure EQ 3 THEN BEGIN 
          cdist = dluminosity(zarray)/dh
          if i eq 0L then plot, zarray, cdist, line=i, $
           xr=[0,5],yr=[0,16],/xst,/yst,$
           xtit='redshift z',$
           ytit='luminosity distance D!dL!n/D!dH' else $ 
           oplot, zarray, cdist, line=i
          
                                ; distance modulus
       ENDIF ELSE IF figure EQ 4 THEN BEGIN 
          cdist = dmodulus(zarray)+5.d0*alog10(redh100())
          if i eq 0L then plot, zarray, cdist, line=i, $
           xr=[0,5],yr=[40,50],/xst,/yst,$
           xtit='redshift z',$
           ytit='distance modulus DM + 5 log h (mag)' else $ 
           oplot, zarray, cdist, line=i
          
                                ; comoving volume element
       ENDIF ELSE IF figure EQ 5 THEN BEGIN 
          cvol = dvcomoving(zarray)/(dhubble()^3)
          if i eq 0L then plot, zarray, cvol, line=i, $
           xr=[0,5],yr=[0,1.1],/xst,/yst,$
           xtit='redshift z',$
           ytit='comoving volume element [1/D!dH!n]!u3!n dV!d!nc/dz/dO' else $ 
           oplot, zarray, cvol, line=i

                                ; lookback time tl/th and age t/th
       ENDIF ELSE IF figure EQ 6 THEN BEGIN 
          cage = getage(zarray)/th
          clook = getage(0.d0)/th-cage
          if i eq 0L then plot, zarray, cage, line=i, $
           xr=[0,5],yr=[0,1.2],/xst,/yst,$
           xtit='redshift z',$
           ytit='lookback time t!dL!n/t!dH!n and age t/t!dH!n' else $ 
           oplot, zarray, cage, line=i
          oplot,zarray,clook,line=i

                                ; dimensionless intersection probability
       ENDIF ELSE IF figure EQ 7 THEN BEGIN 
          dpdz = (1.d0+zarray)^2/epeebles(zarray)
          if i eq 0L then plot, zarray, dpdz, line=i, $
           xr=[0,5],yr=[0,6],/xst,/yst,$
           xtit='redshift z',$
           ytit='dimensionless intersection probability dP/dz' else $ 
           oplot, zarray, dpdz, line=i
          
       ENDIF 
       
    ENDFOR 

    IF n_elements(postscript) NE 0 THEN BEGIN
       device,/close
       set_plot,'x'
       !p.charsize = oldcharsize
    ENDIF 

return 
END 
