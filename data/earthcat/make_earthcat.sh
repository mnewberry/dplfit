#!/bin/bash
# Download Southern California Seismographic Network data catalog

FILELIST="1984.catalog 1987.catalog 1990.catalog 1993.catalog 1996.catalog 1999.catalog 1985.catalog 1988.catalog 1991.catalog 1994.catalog 1997.catalog 2000.catalog 1986.catalog 1989.catalog 1992.catalog 1995.catalog 1998.catalog"
# To include 1932-1959 and 1960-1974:
#FILELIST="193259.catalog 196074.catalog $FILELIST"

for FI in $FILELIST 
do wget http://web.archive.org/web/20020607105011/http://www.scecdc.scec.org/ftp/catalogs/SCSN/`basename $FI`
done

cat *.catalog | cut -c 47-49 > mags.all
