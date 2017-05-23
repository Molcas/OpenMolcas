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
      SUBROUTINE PRCMAT(NSS,XMATR,XMATI)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)
C Write out matrix elements over states as a complex matrix
C in square format
      DO JSTA=1,NSS,2
       JEND=MIN(NSS,JSTA+1)
       WRITE(6,*)
       WRITE(6,'(1X,A8,12X,I3,35X,I3)')' STATE  ',(JSS,JSS=JSTA,JEND)
       DO ISS=1,NSS
       WRITE(6,'(1X,I4,2x,2(A1,F10.6,A1,F10.6,A1,3x))')
     &           ISS,('(',XMATR(ISS,JSS),',',XMATI(ISS,JSS),
     &           ')',JSS=JSTA,JEND)
       END DO
      END DO
      RETURN
      END
      SUBROUTINE MULMAT(NSS,XMATR,XMATI,ee,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)
      COMPLEX*16 Z(NSS,NSS)
      ee=0.d0
      DO ISS=1,NSS
      DO JSS=1,NSS
      Z(ISS,JSS)=(0.0d0,0.0d0)
      enddo
      enddo
      DO ISS=1,NSS
      DO JSS=1,NSS
      ee=ee+XMATR(ISS,JSS)*XMATR(ISS,JSS)+
     & XMATI(ISS,JSS)*XMATI(ISS,JSS)
      Z(ISS,JSS)=Z(ISS,JSS)+
     &DCMPLX(XMATR(ISS,JSS),XMATI(ISS,JSS))
      enddo
      enddo
      RETURN
      END
C THE ORIGI : WRITE(6,'(1X,I4,2x,2(A1,G16.9,A1,G16.9,A1,3x))')
