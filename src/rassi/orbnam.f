************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      CHARACTER*8 FUNCTION ORBNAM(ISORB,LORB)
      IMPLICIT NONE
#include "WrkSpc.fh"
      CHARACTER*8 STRING8
      CHARACTER*2 ORBTYP
      INTEGER ISORB,IPART,ISMLAB,ISOIND,LORB,NPART
      INTEGER LOINFO,KOINFO
      NPART= IWORK(LORB+5)
      KOINFO=19
      LOINFO=LORB-1+KOINFO
      ISMLAB= IWORK(LOINFO+ 1+(ISORB-1)*8)
      IPART= IWORK(LOINFO+ 4+(ISORB-1)*8)
      ISOIND= IWORK(LOINFO+ 2+(ISORB-1)*8)
      ORBTYP='De'
      IF(IPART.EQ.NPART-1) ORBTYP='Fr'
      IF(IPART.EQ.NPART-2) ORBTYP='Se'
      IF(IPART.EQ.NPART-3) ORBTYP='In'
      IF(IPART.LE.NPART-4) ORBTYP='Ac'
      WRITE(STRING8,'(A2,I1,A1,I3.3,1X)')ORBTYP,ISMLAB,':',ISOIND
      ORBNAM=STRING8
      RETURN
      END
