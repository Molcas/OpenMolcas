!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine Sogr(idbg,N,SS,SINV,P,G,A1)
!     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
!     ORTHOGONALIZED BASIS.
!     N              DIMENSION OF MATRICES. ISIZE=N*(N+1)/2
!     SS(ISIZE)      ORIGINAL OVERLAP MATRIX (LOWER TRIANGULAR)
!                    WILL NOT BE DESTROYED
!     P (ISIZE)      OUTPUT TRANSFORMATION MATRIX
!     G (ISIZE)      SCRATCH
!     A1(N)          SCRATCH

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: idbg, N
real(kind=wp), intent(in) :: SS(*)
real(kind=wp), intent(out) :: SINV(N,N)
real(kind=wp), intent(_OUT_) :: P(*), G(*), A1(*)
integer(kind=iwp) :: I, ierr, IH, IJ, IL, IQ, J, J1, JF, JL, JQ, K, L, LG
real(kind=wp) :: ETOT, S1KK, SUM_

if (iDbg > 0) call PrMat(idbg,SS,n,0,'SS')
ierr = 0
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
    JF = 1
    JL = IL
    do K=1,J1
      SUM_ = Zero
      JL = JL+1
      JF = JF+K-1
      IH = JF
      do L=K,J1
        IH = IH+L-1
        SUM_ = SUM_+A1(L)*G(IH)
      end do
      G(JL) = -SUM_
    end do
  end if
  if (s1kk <= 1.0e-16_wp) then
    write(u6,*) '    Sogr| j=',j,' s1kk=',s1kk
    ierr = ierr+1
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
if (ierr > 0) then
  write(u6,'(a)') 'function has negative norm'
  call Abend()
end if
if (iDbg > 0) call PrMat(idbg,P,n,0,'P')

return

end subroutine Sogr
