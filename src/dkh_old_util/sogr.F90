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

implicit real*8(A-H,O-Z)
dimension SS(*), P(*), G(*), A1(*), SINV(N,N)
integer ierr

if (iDbg > 0) call PrMat(idbg,SS,n,0,'SS')
ierr = 0
JL = 0
IQ = 0
do J=1,N
  IL = JL
  JQ = IQ
  S1KK = SS(IQ+J)
  G(IL+J) = 1.d0
  if (J == 1) GO TO 341
  J1 = J-1
  JL = 0
  do K=1,J1
    LG = JQ
    ETOT = 0.d0
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
    SUM = 0.d0
    JL = JL+1
    JF = JF+K-1
    IH = JF
    do L=K,J1
      IH = IH+L-1
      SUM = SUM+A1(L)*G(IH)
    end do
    G(JL) = -SUM
  end do
341 continue
  if (s1kk <= 1.D-16) then
    write(6,*) '    Sogr| j=',j,' s1kk=',s1kk
    ierr = ierr+1
  end if
  S1KK = 1.d0/sqrt(S1KK)
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
    SINV(I,J) = 0.d0
    SINV(J,I) = P(IJ)
  end do
end do
if (ierr > 0) call errex_rel('function has negative norm')
if (iDbg > 0) call PrMat(idbg,P,n,0,'P')

return

end subroutine Sogr
