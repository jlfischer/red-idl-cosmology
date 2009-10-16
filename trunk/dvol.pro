;;;;;;;;;;;;;;;;;
function redvol,z
  return,dvcomoving(z,/Mpc)
end

function nvol,z1,z2
; comoving volume (integrated) in given redshift interval per square arcsec
     deg2sr    = (180.d0/!dpi)^2 ; deg^2/sr
     return,1./deg2sr/3600.*qromb('redvol',z1,z2,/double)
end

pro dvol, z1, z2, nz, _extra=_extra
     np = 200l
     zarr = (findgen(np)/(np-1)*(1-z1/z2)+z1/z2)*z2
     varr = fltarr(np)
     for i=0,np-1 do varr[i]=nvol(zarr[i],z2)
     totvol = nvol(z1,z2)

     nzarr = reverse(lindgen(nz)+1)
     divz = fltarr(nz)
     for i=0,nz-1 do divz[i]=interpol(zarr,varr,totvol/nz*nzarr[i])

     intdivz = [divz,z2]
     for i=1,n_elements(intdivz)-1 do print,intdivz[i],nvol(intdivz[i-1],intdivz[i])

 end 
