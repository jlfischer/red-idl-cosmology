%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% idl cosmology package 'red', v1.0, 02jun13 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  this package consists of two programs, the driver routine 'red.pro'
  and the file 'cosmology_routines.pro' which contains all the
  cosmological functions. 

  comments, corrections, and suggestions (all welcome) to: 
  leonidas moustakas (oxford)   l.moustakas1[at]physics.ox.ac.uk
  john moustakas (steward, UA)  ioannis[at]as.arizona.edu

  to set up, put directory ./red/ from this distribution in a
  convenient place (eg under ~/idl/), and add directory location to
  the IDL_PATH in your idl-startup script or ~/.login file:
  setenv IDL_PATH $IDL_PATH{:}+/home/user/idl/red

  to start, initialize once by running red (could place this in a 
  startup file)

IDL> red

  which sets up the default cosmology (check it with /verbose)

IDL> red,/default,/verbose
Omega Matter = 0.3000
Omega Lambda = 0.7000
H_0 / 100    = 0.7500

  the parameters can be set explicitely, e.g. 

IDL> red,omega0=0.2,omegalambda=0.8,h100=0.66,/verb
Omega Matter = 0.2000
Omega Lambda = 0.8000
H_0 / 100    = 0.6600

  these parameters can be accessed on command or in code, e.g. by 

IDL> print,redh100(),redomegal(),redomega0()
      0.66000003      0.80000001      0.20000000

  once any particular cosmology is set, e.g.

IDL> red,/def

  it is remembered.  look at hogg_plots for examples...


IDL> red,/help
----------------------------------------------------------------
Available cosmological routines:
 
   DHUBBLE()              - Hubble distance
   THUBBLE()              - Hubble time
   DANGULAR(z)            - angular diameter distance
   DANGULARDIFF(z1,z2)    - difference in angular diameter
   DCOMOVINGLOS(z)        - line-of-sight co-moving distance
   DCOMOVINGTRANSVERSE(z) - transverse co-moving distance
   DLUMINOSITY(z)         - luminosity distance
   DMODULUS(z)            - distance modulus
   DVCOMOVING(z)          - differential co-moving volume
   VCOMOVING(z)           - integrated co-moving volume
   GETAGE(z)              - look-back time
   GETREDSHIFT(age)       - redshift for a given age
 
   REDOMEGAL()            - return current omegalambda
   REDOMEGA0()            - return current omega0
   REDh100()              - return current h100

Available units:
   cm    - centimeter          s    - second
   meter - meter               yr   - year
   pc    - parsec [default]    Myr  - megayear
   kpc   - kiloparsec          Gyr  - gigayear [default]
   Mpc   - megaparsec
----------------------------------------------------------------

  get quick summary of useful quantities at a specified redshift

IDL> red,z=1.0
Age (z=0.000) = 12.568951 Gyr
Age (z=1.000) = 5.3678168 Gyr
Lookback time = 7.2011342 Gyr
DModulus      = 43.950422 mag
DAngular      = 1541.7868 Mpc
DLuminosity   = 6167.1471 Mpc
Scale         = 7.4747862 kpc/arcsec

  change output units by switch

IDL> print,dluminosity(1.0)
   6.1671471e+09
IDL> print,dluminosity(1.0,/Mpc)
       6167.1471
IDL> print,dluminosity(1.0,/cm)
   1.9031816e+28

  input-redshift arrays are okay.

IDL> zz=findgen(6)
IDL> print,zz,dvcomoving(zz)/dhubble()^3
 0.00000    1.00000     2.00000     3.00000    4.00000     5.00000
 0.0000000  0.33799393  0.49311672  0.49368840 0.45352966  0.40726417

  for more examples, look at other README.* files in the distribution, 
  and at the program hogg_plots.pro

  enjoy!

