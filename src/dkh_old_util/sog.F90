!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SOG(N,SS,SINV,P,G,A1)
!     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
!     ORTHOGONALIZED BASIS.
!     N              DIMENSION OF MATRICES. ISIZE=N*(N+1)/2
!     SS(ISIZE)      ORIGINAL OVERLAP MATRIX (LOWER TRIANGULAR)
!                    WILL NOT BE DESTROYED
!     P (ISIZE)      OUTPUT TRANSFORMATION MATRIX
!     G (ISIZE)      SCRATCH
!     A1(N)          SCRATCH

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: SS(N*(N+1)/2)
real(kind=wp), intent(out) :: SINV(N,N), P(N*(N+1)/2), G(N*(N+1)/2), A1(N)
integer(kind=iwp) :: I, I_F, IH, IJ, IL, IQ, J, J1, JL, JQ, K, L, LG
real(kind=wp) :: ETOT, S1KK, SUM_

JL = 0
IQ = 0
do J=1,N
  IL = JL
  JQ = IQ
  S1KK = SS(IQ+J)
  G(IL+J) = One
  if (J > 1) then
    J1 = J-1
    JL = 0
    do K=1,J1
      LG = JQ
      ETOT = Zero
      do L=1,K
        LG = LG+1
        JL = JL+1
        ETOT = ETOT+SS(LG)*G(JL)
      end do
      S1KK = S1KK-ETOT*ETOT
      A1(K) = ETOT
    end do
    I_F = 1
    JL = IL
    do K=1,J1
      SUM_ = Zero
      JL = JL+1
      I_F = I_F+K-1
      IH = I_F
      do L=K,J1
        IH = IH+L-1
        SUM_ = SUM_+A1(L)*G(IH)
      end do
      G(JL) = -SUM_
    end do
  end if
  S1KK = One/sqrt(S1KK)
  JL = IL
  do K=1,J
    JL = JL+1
    IQ = IQ+1
    G(JL) = G(JL)*S1KK
    P(IQ) = G(JL)
  end do
end do

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    SINV(I,J) = Zero
    SINV(J,I) = P(IJ)
  end do
end do

return

end subroutine SOG
