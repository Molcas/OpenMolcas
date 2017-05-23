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
      SUBROUTINE PRPRTTAB(IPRTTAB)
      IMPLICIT NONE
      INTEGER IPRTTAB(*)
      INTEGER IPART,NPART
      INTEGER NSYM,ISYM
C     INTEGER MXMEM

C Executable statements
      WRITE(6,*)
      WRITE(6,*)' Partition table printout'
      WRITE(6,'(A,I5)')'Table size        NSIZE=',IPRTTAB(1)
      WRITE(6,'(A,I5)')'Table type ID     ITYPE=',IPRTTAB(2)
      WRITE(6,'(A,I5)')'Nr of partitions  NPART=',IPRTTAB(3)
      WRITE(6,'(A,I5)')'Nr of symm labels NSYM =',IPRTTAB(4)
      NPART=IPRTTAB(3)
      NSYM =IPRTTAB(4)
      WRITE(6,'(8X,I5,5X,8I5)')(IPRTTAB(5+ISYM),ISYM=0,NSYM)
      DO IPART=1,NPART
       WRITE(6,'(I3,5X,I5,5X,8I5)')IPART,
     &      (IPRTTAB(5+ISYM+(NSYM+1)*IPART),ISYM=0,NSYM)
      END DO
      RETURN
      END
