Notes to fix `renv` on an M1 make that used Macports:

 * do the following:

        mkdir -p ~/.R
        nano ~/.R/Makevars

   Then add the lines:

        CPPFLAGS=-I/opt/local/include
        LDFLAGS=-L/opt/local/lib

 * install gfortran-14.2-universal.pkg from [here](https://mac.r-project.org/tools/)
