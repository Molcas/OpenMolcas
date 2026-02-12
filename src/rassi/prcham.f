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
      SUBROUTINE PRCHAM(NSS,CHAMR,CHAMI)
      use definitions, only: iwp, wp, u6

      IMPLICIT NONE

      integer(kind=iwp), intent(in):: NSS
      REAL(kind=wp), intent(in):: CHAMR(NSS,NSS),CHAMI(NSS,NSS)

      integer(kind=iwp) JSTA,JEND,ISS,JSS
C Write out a complex Hamiltonian (or other hermitian matrix)
C in a triangular format.
      DO JSTA=1,NSS,2
       JEND=MIN(NSS,JSTA+1)
       WRITE(u6,*)
       WRITE(u6,'(1X,A8,11X,I4,33X,I4)')'SO-STATE',(JSS,JSS=JSTA,JEND)
       DO ISS=JSTA,NSS
       WRITE(u6,'(1X,I4,2x,2(A1,F15.11,A1,F15.11,A1,3x))')
     &           ISS,('(',CHAMR(ISS,JSS),',',CHAMI(ISS,JSS),
     &           ')',JSS=JSTA,MIN(ISS,JEND))
       END DO
      END DO

      END SUBROUTINE PRCHAM
