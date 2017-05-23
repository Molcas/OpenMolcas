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
      SUBROUTINE PRFSBTAB(IFSBTAB)
      IMPLICIT NONE
      INTEGER IFSBTAB(*)
C     INTEGER NFSB,NASPRT,LPOS,ISPART,IFSB
      INTEGER NFSB,NASPRT,     ISPART,IFSB
      INTEGER NDET,ISTA,NHEAD,KPOS
      IF(IFSBTAB(2).NE.73) THEN
        WRITE(6,*)' PRFSBTAB error: Not a Fock Sector Block Table.'
        WRITE(6,*)' Table type code  =',IFSBTAB(2)
        CALL ABEND()
      END IF
      WRITE(6,*)
      WRITE(6,*)'============================================='
      WRITE(6,*)' Fock Sector Table printout'
      WRITE(6,'(a,i9)')'               Table size:',IFSBTAB(1)
      WRITE(6,'(a,i9)')'          Table type code:',IFSBTAB(2)
      WRITE(6,'(a,i9)')' Nr of Fock Sector Blocks:',IFSBTAB(3)
      WRITE(6,'(a,i9)')'      Nr of Subpartitions:',IFSBTAB(4)
      WRITE(6,'(a,i9)')' Total nr of Determinants:',IFSBTAB(5)
      WRITE(6,'(a,i9)')'        Hash Map Capacity:',IFSBTAB(6)
      WRITE(6,'(a,i9)')'        Hash Map 1st word:',IFSBTAB(7)
      WRITE(6,*)
      WRITE(6,*)'FS Block   BlkSiz    Start indx     Substring Types'
      NFSB=IFSBTAB(3)
      NASPRT=IFSBTAB(4)
      NHEAD=7
      KPOS=NHEAD+1
      DO IFSB=1,NFSB
        NDET=IFSBTAB(KPOS+NASPRT)
        ISTA=IFSBTAB(KPOS+NASPRT+1)
        WRITE(6,'(1X,I6,3X,I7,3x,I10,5X,10I4)')
     &   IFSB,NDET,ISTA,(IFSBTAB(KPOS-1+ISPART),ISPART=1,NASPRT)
        KPOS=KPOS+NASPRT+2
      END DO
      WRITE(6,*)'============================================='
      RETURN
      END
