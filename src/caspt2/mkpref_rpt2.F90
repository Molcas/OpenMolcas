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
! Copyright (C) 2006, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2006  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKPREF_RPT2(N,G2,PREF,NPREF)
! Compute PREF(PQRS) = <0| 0.5*Epqrs |0>
! from G2(P,Q,R,S) = <0| Epqrs |0>
! Storage differs: PREF is triangular
! in the Fortran-like indices PQ, RS.

use Index_Functions, only: iTri
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: N, NPREF
real(kind=wp), intent(in) :: G2(N,N,N,N)
real(kind=wp), intent(out) :: PREF(NPREF)
integer(kind=iwp) :: I, IJ, IJKL, IJKLT, IJLK, IJT, J, JI, JIKL, JILK, K, KL, KLT, L, LK
real(kind=wp) :: P1, P2

IJT = 0
IJKLT = 0
do I=1,N
  do J=1,I
    IJT = IJT+1
    IJ = I+N*(J-1)
    JI = J+N*(I-1)
    KLT = 0
    outer: do K=1,N
      do L=1,K
        KLT = KLT+1
        if (KLT > IJT) exit outer
        IJKLT = IJKLT+1
        KL = K+N*(L-1)
        LK = L+N*(K-1)

        P1 = Half*G2(I,J,K,L)
        P2 = Half*G2(I,J,L,K)
        IJKL = iTri(IJ,KL)
        IJLK = iTri(IJ,LK)
        JIKL = iTri(JI,KL)
        JILK = iTri(JI,LK)
        PREF(IJKL) = P1
        PREF(IJLK) = P2
        PREF(JIKL) = P2
        PREF(JILK) = P1
      end do
    end do outer
  end do
end do

end subroutine MKPREF_RPT2
