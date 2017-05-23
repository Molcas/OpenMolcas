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
      REAL*8 FUNCTION OVERLAP_RASSI(IFSBTAB1,IFSBTAB2,PSI1,PSI2)
      IMPLICIT NONE
      REAL*8 PSI1(*),PSI2(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER ISSTARR(50)
      INTEGER NFSB1,NASPRT1,NDETS1,KSTARR1
      INTEGER NFSB2,NASPRT2,NDETS2,NHSH2,KHSH2,KSTARR2
      INTEGER IFSB1,KPOS1,ISPART
      INTEGER IFSB2,KPOS2
      INTEGER IBLKPOS1,NBLKSIZ1
      INTEGER IBLKPOS2,NBLKSIZ2
      real*8 ddot_
      external ddot_

C Purpose: Compute the overlap of two wave functions
C The FS blocks of the two wave functions:
      NFSB1  =IFSBTAB1(3)
      NASPRT1=IFSBTAB1(4)
      NDETS1 =IFSBTAB1(5)
      KSTARR1=8
      NFSB2  =IFSBTAB2(3)
      NASPRT2=IFSBTAB2(4)
      NDETS2 =IFSBTAB2(5)
      NHSH2  =IFSBTAB2(6)
      KHSH2  =IFSBTAB2(7)
      KSTARR2=8
      OVERLAP_RASSI=0.0D0
      IF (NFSB1.EQ.0) RETURN
      IF (NFSB2.EQ.0) RETURN
      IF (NASPRT1.NE.NASPRT2) GOTO 901
      IF (NDETS1.EQ.0) RETURN
      IF (NDETS2.EQ.0) RETURN

C Loop over FS blocks of the PSI1 wave function
      DO IFSB1=1,NFSB1
        KPOS1=KSTARR1+(NASPRT1+2)*(IFSB1-1)
        DO ISPART=1,NASPRT1
          ISSTARR(ISPART)=IFSBTAB1(KPOS1-1+ISPART)
        END DO
        NBLKSIZ1 =IFSBTAB1(KPOS1+NASPRT1  )
        IBLKPOS1 =IFSBTAB1(KPOS1+NASPRT1+1)
C Find this block in the PSI2 structure.
        CALL HSHGET(ISSTARR,NASPRT2,NASPRT2+2,IFSBTAB2(KSTARR2),
     &                NHSH2,IFSBTAB2(KHSH2),IFSB2)
        IF(IFSB2.EQ.0) GOTO 100
        KPOS2=KSTARR2+(NASPRT2+2)*(IFSB2-1)
        NBLKSIZ2 =IFSBTAB2(KPOS2+NASPRT2  )
        IF(NBLKSIZ1.NE.NBLKSIZ2) THEN
          WRITE(6,*)' OVERLAP Error: The same FS block has not'
          WRITE(6,*)' the same size in PSI1 and PSI2.'
          CALL ABEND()
        END IF
        IBLKPOS2 =IFSBTAB2(KPOS2+NASPRT2+1)
        OVERLAP_RASSI=OVERLAP_RASSI+
     &         DDOT_(NBLKSIZ1,PSI1(IBLKPOS1),1,PSI2(IBLKPOS2),1)
 100    CONTINUE
      END DO
      RETURN
 901  CONTINUE
      WRITE(6,*)' OVERLAP Error: The two wave function structures'
      WRITE(6,*)' have different nr of subpartitions!'
      CALL ABEND()
      END
