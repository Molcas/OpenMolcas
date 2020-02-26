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
      SUBROUTINE MKEPS(FIFA,DREF)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      INTEGER I, ID, IEPS, IEPSA, IEPSE, IEPSI, ISTLT, ISYM
      INTEGER ITOT, NA, NI, NO
      REAL*8 E
      REAL*8 FIFA(*),DREF(*)

      CALL QENTER('MKEPS')

c   Orbital energies, EPS, EPSI,EPSA,EPSE:
      IEPS=0
      IEPSI=0
      IEPSA=0
      IEPSE=0
      ISTLT=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NO=NORB(ISYM)
        DO I=1,NI
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSI=IEPSI+1
          EPSI(IEPSI)=E
        END DO
        DO I=NI+1,NI+NA
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSA=IEPSA+1
          EPSA(IEPSA)=E
        END DO
        DO I=NI+NA+1,NO
          E=FIFA(ISTLT+(I*(I+1))/2)
          IEPS=IEPS+1
          EPS(IEPS)=E
          IEPSE=IEPSE+1
          EPSE(IEPSE)=E
        END DO
        ISTLT=ISTLT+(NO*(NO+1))/2
      END DO

C EASUM=CONTRACT EPSA WITH DIAGONAL OF ACTIVE DENS
      EASUM=0.0D0
      DO ISYM=1,NSYM
        NA=NASH(ISYM)
        DO I=1,NA
          ITOT=NAES(ISYM)+I
          ID=(ITOT*(ITOT+1))/2
          EASUM=EASUM+EPSA(ITOT)*DREF(ID)
        END DO
      END DO

#ifdef _DEBUG_
      WRITE(6,*)
      WRITE(6,*)'      ORBITAL ENERGIES, EPS:'
      WRITE(6,'(1X,5F12.6)')(EPS(I),I=1,NORBT)
      WRITE(6,*)'      INACTIVE ORBITAL ENERGIES, EPSI:'
      WRITE(6,'(1X,5F12.6)')(EPSI(I),I=1,NISHT)
      WRITE(6,*)'        ACTIVE ORBITAL ENERGIES, EPSA:'
      WRITE(6,'(1X,5F12.6)')(EPSA(I),I=1,NASHT)
      WRITE(6,*)'      EXTERNAL ORBITAL ENERGIES, EPSE:'
      WRITE(6,'(1X,5F12.6)')(EPSE(I),I=1,NSSHT)
      WRITE(6,*)'      Active energies contr. w. DREF:'
      WRITE(6,'(1X,5F12.6)') EASUM
#endif

      E=0.0D0
      DO I=1,NASHT
       E=E+EPSA(I)*DREF((I*(I+1))/2)
      END DO

#ifdef _DEBUG_
      WRITE(6,*)' Compare to E=',E
      WRITE(6,*)' DREF:'
      WRITE(6,'(1x,5F16.8)')(DREF(I),I=1,(NASHT*(NASHT+1))/2)
#endif

      CALL QEXIT('MKEPS')

      RETURN
      END
