RED is a set of routines written in IDL for performing cosmological calculations based on David Hogg's extraordinarily useful [astro-ph article](http://arxiv.org/abs/astro-ph/9905116).  There are two main programs: the driver routine 'RED.PRO' and the file 'COSMOLOGY\_ROUTINES.PRO' which contains all the cosmological functions.

Downloading and installation are very quick and straight-forward.  First, grab the latest version of the code using [subversion](http://subversion.tigris.org):

```
svn co http://red-idl-cosmology.googlecode.com/svn/trunk/ red
```

or follow the instructions on the [Source](http://code.google.com/p/red-idl-cosmology/source/checkout) tab, above.  Next, move the 'red' directory to a local IDL sub-directory (e.g., ${HOME}/idl); if this directory is not part of your IDL path then add the line "setenv IDL\_PATH $IDL\_PATH{:}+${HOME}/idl/red" to your '.idlenv' file.  And then compute away!  To get started, see the README.red file for specific examples.

Comments, corrections, and suggestions (all welcome) should be directed to Leonidas Moustakas (JPL) or John Moustakas (University of California, San Diego).