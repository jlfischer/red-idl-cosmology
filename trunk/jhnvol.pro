function redvol,z
  return,dvcomoving(z,/Mpc)
end

function jhnvol,z1,z2
; comoving volume (integrated) in given redshift interval per square arcsec
     vol=dblarr(n_elements(z1))
     for i=0l,n_elements(z1)-1 do begin
         if n_elements(z1) gt 1 and i mod 100 eq 0 then message,strcompress(i),/inf
         vol[i] = 3.05d-4/3600.0d0*qpint1d('redvol',z1[i],z2[i])
     endfor
  return,vol
end
 
