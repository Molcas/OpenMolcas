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
      CHARACTER(LEN=8) FUNCTION ORBNAM(ISORB,ORBTAB)
      IMPLICIT NONE
      Integer ISORB, ORBTAB(*)
      CHARACTER(LEN=8) STRING8
      CHARACTER(LEN=2) ORBTYP
      INTEGER IPART,ISMLAB,ISOIND,NPART
      INTEGER KOINFO
      NPART= ORBTAB(6)
      KOINFO=19
      ISMLAB= ORBTAB(KOINFO+ 1+(ISORB-1)*8)
      IPART=  ORBTAB(KOINFO+ 4+(ISORB-1)*8)
      ISOIND= ORBTAB(KOINFO+ 2+(ISORB-1)*8)
      ORBTYP='De'
      IF(IPART.EQ.NPART-1) ORBTYP='Fr'
      IF(IPART.EQ.NPART-2) ORBTYP='Se'
      IF(IPART.EQ.NPART-3) ORBTYP='In'
      IF(IPART.LE.NPART-4) ORBTYP='Ac'
      WRITE(STRING8,'(A2,I1,A1,I3.3,1X)')ORBTYP,ISMLAB,':',ISOIND
      ORBNAM=STRING8

      END FUNCTION ORBNAM
