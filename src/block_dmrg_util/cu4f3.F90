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
! Compute F3 using cumulant reconstruction except for G3-dependent terms

function CU4F3H(NAC,E,ES,G1,G2,F1,F2,iP,iQ,jP,jQ,kP,kQ)

use Constants, only: Zero, Two, Three, Six, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CU4F3H
integer(kind=iwp), intent(in) :: NAC, iP, iQ, jP, jQ, kP, kQ
real(kind=wp), intent(in) :: E(NAC), ES, G1(NAC,NAC), G2(NAC,NAC,NAC,NAC), F1(NAC,NAC), F2(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: lT

CU4F3H = Zero

! Compute these 3PDM contributions in other place
!
!      // G(I,i)G(JKL,jkl)[4]
!CU4F3H = CU4F3H &
!         + G1(lT,lT)*G3(iP,iQ,jP,jQ,kP,kQ)
!      // G(I,j)G(JKL,ikl)[12]
!CU4F3H = CU4F3H &
!         -Half*G1(iP,lT)*G3(lT,iQ,jP,jQ,kP,kQ) &
!         -Half*G1(jP,lT)*G3(lT,jQ,iP,iQ,kP,kQ) &
!         -Half*G1(kP,lT)*G3(lT,kQ,iP,iQ,jP,jQ) &
!         -Half*G1(lT,iQ)*G3(iP,lT,jP,jQ,kP,kQ) &
!         -Half*G1(lT,jQ)*G3(jP,lT,iP,iQ,kP,kQ) &
!         -Half*G1(lT,kQ)*G3(kP,lT,iP,iQ,jP,jQ)

!     // G(I,i)G(JKL,jkl)[4]
CU4F3H = CU4F3H &
         +G1(iP,iQ)*F2(jP,jQ,kP,kQ) &
         +G1(jP,jQ)*F2(iP,iQ,kP,kQ) &
         +G1(kP,kQ)*F2(iP,iQ,jP,jQ)
!     // G(I,j)G(JKL,ikl)[12]
CU4F3H = CU4F3H &
         -Half*G1(iP,jQ)*F2(jP,iQ,kP,kQ) &
         -Half*G1(iP,kQ)*F2(kP,iQ,jP,jQ) &
         -Half*G1(jP,iQ)*F2(iP,jQ,kP,kQ) &
         -Half*G1(jP,kQ)*F2(kP,jQ,iP,iQ) &
         -Half*G1(kP,iQ)*F2(iP,kQ,jP,jQ) &
         -Half*G1(kP,jQ)*F2(jP,kQ,iP,iQ)
!     // G(IJ,ij)G(KL,kl)[3]
CU4F3H = CU4F3H &
         +G2(iP,iQ,jP,jQ)*F1(kP,kQ) &
         +G2(iP,iQ,kP,kQ)*F1(jP,jQ) &
         +F1(iP,iQ)*G2(jP,jQ,kP,kQ)
!     // G(IJ,ik)G(KL,jl)[12]
CU4F3H = CU4F3H &
         -Half*G2(iP,iQ,jP,kQ)*F1(kP,jQ) &
         -Half*G2(jP,jQ,iP,kQ)*F1(kP,iQ) &
         -Half*G2(kP,kQ,iP,jQ)*F1(jP,iQ) &
         -Half*F1(iP,jQ)*G2(jP,iQ,kP,kQ) &
         -Half*F1(iP,kQ)*G2(kP,iQ,jP,jQ) &
         -Half*F1(jP,kQ)*G2(kP,jQ,iP,iQ)
!     // G(I,i)G(J,j)G(KL,kl)[6]
CU4F3H = CU4F3H &
         -Two*G1(iP,iQ)*G1(jP,jQ)*F1(kP,kQ) &
         -Two*G1(iP,iQ)*G1(kP,kQ)*F1(jP,jQ) &
         -Two*G1(jP,jQ)*G1(kP,kQ)*F1(iP,iQ) &
         -Two*G1(iP,iQ)*ES*G2(jP,jQ,kP,kQ) &
         -Two*G1(jP,jQ)*ES*G2(iP,iQ,kP,kQ) &
         -Two*G1(kP,kQ)*ES*G2(iP,iQ,jP,jQ)
!     // G(I,j)G(J,i)G(KL,kl)[6]
CU4F3H = CU4F3H &
         +G1(iP,jQ)*G1(jP,iQ)*F1(kP,kQ) &
         +G1(iP,kQ)*G1(kP,iQ)*F1(jP,jQ) &
         +G1(jP,kQ)*G1(kP,jQ)*F1(iP,iQ)
!     // G(I,i)G(J,k)G(KL,jl)[24]
CU4F3H = CU4F3H &
         +G1(iP,iQ)*G1(jP,kQ)*F1(kP,jQ) &
         +G1(iP,iQ)*G1(kP,jQ)*F1(jP,kQ) &
         +G1(jP,jQ)*G1(iP,kQ)*F1(kP,iQ) &
         +G1(jP,jQ)*G1(kP,iQ)*F1(iP,kQ) &
         +G1(kP,kQ)*G1(iP,jQ)*F1(jP,iQ) &
         +G1(kP,kQ)*G1(jP,iQ)*F1(iP,jQ) &
         +ES*G1(iP,jQ)*G2(jP,iQ,kP,kQ) &
         +ES*G1(jP,iQ)*G2(iP,jQ,kP,kQ) &
         +ES*G1(iP,kQ)*G2(kP,iQ,jP,jQ) &
         +ES*G1(kP,iQ)*G2(iP,kQ,jP,jQ) &
         +ES*G1(jP,kQ)*G2(kP,jQ,iP,iQ) &
         +ES*G1(kP,jQ)*G2(jP,kQ,iP,iQ)
!     // G(I,j)G(J,k)G(KL,il)[24]
CU4F3H = CU4F3H &
         -Half*G1(iP,jQ)*G1(jP,kQ)*F1(kP,iQ) &
         -Half*G1(iP,kQ)*G1(kP,jQ)*F1(jP,iQ) &
         -Half*G1(jP,iQ)*G1(iP,kQ)*F1(kP,jQ) &
         -Half*G1(jP,kQ)*G1(kP,iQ)*F1(iP,jQ) &
         -Half*G1(kP,iQ)*G1(iP,jQ)*F1(jP,kQ) &
         -Half*G1(kP,jQ)*G1(jP,iQ)*F1(iP,kQ)
!     // G(I,i)G(J,j)G(K,k)G(L,l)[1]
CU4F3H = CU4F3H &
         +Six*G1(iP,iQ)*G1(jP,jQ)*G1(kP,kQ)*ES
!     // G(I,j)G(J,i)G(K,k)G(L,l)[6]
CU4F3H = CU4F3H &
         -Three*G1(iP,jQ)*G1(jP,iQ)*G1(kP,kQ)*ES &
         -Three*G1(iP,kQ)*G1(kP,iQ)*G1(jP,jQ)*ES &
         -Three*G1(jP,kQ)*G1(kP,jQ)*G1(iP,iQ)*ES
!     // G(I,j)G(J,k)G(K,i)G(L,l)[8]
CU4F3H = CU4F3H &
         +OneHalf*G1(iP,jQ)*G1(jP,kQ)*G1(kP,iQ)*ES &
         +OneHalf*G1(iP,kQ)*G1(kP,jQ)*G1(jP,iQ)*ES

! Loop over lT
do lT=1,NAC
  !     // G(IJ,ik)G(KL,jl)[12]
  CU4F3H = CU4F3H &
           -Half*G2(iP,iQ,jP,lT)*G2(lT,jQ,kP,kQ)*E(lT) &
           -Half*G2(iP,iQ,kP,lT)*G2(lT,kQ,jP,jQ)*E(lT) &
           -Half*G2(jP,jQ,iP,lT)*G2(lT,iQ,kP,kQ)*E(lT) &
           -Half*G2(jP,jQ,kP,lT)*G2(lT,kQ,iP,iQ)*E(lT) &
           -Half*G2(kP,kQ,iP,lT)*G2(lT,iQ,jP,jQ)*E(lT) &
           -Half*G2(kP,kQ,jP,lT)*G2(lT,jQ,iP,iQ)*E(lT)
  !     // G(IJ,kl)G(KL,ij)[3], G(IJ,lk)G(KL,ji)[3]
  !     // G(IJ,kl)G(KL,ji)[3], G(IJ,lk)G(KL,ij)[3]
  CU4F3H = CU4F3H &
           +( &
             +Two*G2(iP,kQ,jP,lT)*G2(kP,iQ,lT,jQ)*E(lT) &
             +Two*G2(iP,jQ,kP,lT)*G2(jP,iQ,lT,kQ)*E(lT) &
             +Two*G2(iP,jQ,lT,kQ)*G2(jP,iQ,kP,lT)*E(lT) &
             +Two*G2(iP,lT,jP,kQ)*G2(kP,jQ,lT,iQ)*E(lT) &
             +Two*G2(iP,lT,kP,jQ)*G2(jP,kQ,lT,iQ)*E(lT) &
             +Two*G2(iP,kQ,lT,jQ)*G2(jP,lT,kP,iQ)*E(lT) &
             +G2(iP,kQ,jP,lT)*G2(kP,jQ,lT,iQ)*E(lT) &
             +G2(iP,jQ,kP,lT)*G2(jP,kQ,lT,iQ)*E(lT) &
             +G2(iP,jQ,lT,kQ)*G2(jP,lT,kP,iQ)*E(lT) &
             +G2(iP,lT,jP,kQ)*G2(kP,iQ,lT,jQ)*E(lT) &
             +G2(iP,lT,kP,jQ)*G2(jP,iQ,lT,kQ)*E(lT) &
             +G2(iP,kQ,lT,jQ)*G2(jP,iQ,kP,lT)*E(lT) &
           )/Six
  !     // G(I,j)G(J,i)G(KL,kl)[6]
  CU4F3H = CU4F3H &
           +G1(iP,lT)*G1(lT,iQ)*G2(jP,jQ,kP,kQ)*E(lT) &
           +G1(jP,lT)*G1(lT,jQ)*G2(iP,iQ,kP,kQ)*E(lT) &
           +G1(kP,lT)*G1(lT,kQ)*G2(iP,iQ,jP,jQ)*E(lT)
  !     // G(I,i)G(J,k)G(KL,jl)[24]
  CU4F3H = CU4F3H &
           +G1(iP,iQ)*G1(jP,lT)*G2(lT,jQ,kP,kQ)*E(lT) &
           +G1(iP,iQ)*G1(kP,lT)*G2(lT,kQ,jP,jQ)*E(lT) &
           +G1(jP,jQ)*G1(iP,lT)*G2(lT,iQ,kP,kQ)*E(lT) &
           +G1(jP,jQ)*G1(kP,lT)*G2(lT,kQ,iP,iQ)*E(lT) &
           +G1(kP,kQ)*G1(iP,lT)*G2(lT,iQ,jP,jQ)*E(lT) &
           +G1(kP,kQ)*G1(jP,lT)*G2(lT,jQ,iP,iQ)*E(lT) &
           +G1(iP,iQ)*G1(lT,jQ)*G2(jP,lT,kP,kQ)*E(lT) &
           +G1(iP,iQ)*G1(lT,kQ)*G2(kP,lT,jP,jQ)*E(lT) &
           +G1(jP,jQ)*G1(lT,iQ)*G2(iP,lT,kP,kQ)*E(lT) &
           +G1(jP,jQ)*G1(lT,kQ)*G2(kP,lT,iP,iQ)*E(lT) &
           +G1(kP,kQ)*G1(lT,iQ)*G2(iP,lT,jP,jQ)*E(lT) &
           +G1(kP,kQ)*G1(lT,jQ)*G2(jP,lT,iP,iQ)*E(lT)
  !     // G(I,j)G(J,k)G(KL,il)[24]
  CU4F3H = CU4F3H &
           -Half*G1(iP,jQ)*G1(jP,lT)*G2(lT,iQ,kP,kQ)*E(lT) &
           -Half*G1(iP,kQ)*G1(kP,lT)*G2(lT,iQ,jP,jQ)*E(lT) &
           -Half*G1(jP,iQ)*G1(iP,lT)*G2(lT,jQ,kP,kQ)*E(lT) &
           -Half*G1(jP,kQ)*G1(kP,lT)*G2(lT,jQ,iP,iQ)*E(lT) &
           -Half*G1(kP,iQ)*G1(iP,lT)*G2(lT,kQ,jP,jQ)*E(lT) &
           -Half*G1(kP,jQ)*G1(jP,lT)*G2(lT,kQ,iP,iQ)*E(lT) &
           -Half*G1(iP,lT)*G1(lT,jQ)*G2(jP,iQ,kP,kQ)*E(lT) &
           -Half*G1(iP,lT)*G1(lT,kQ)*G2(kP,iQ,jP,jQ)*E(lT) &
           -Half*G1(jP,lT)*G1(lT,iQ)*G2(iP,jQ,kP,kQ)*E(lT) &
           -Half*G1(jP,lT)*G1(lT,kQ)*G2(kP,jQ,iP,iQ)*E(lT) &
           -Half*G1(kP,lT)*G1(lT,iQ)*G2(iP,kQ,jP,jQ)*E(lT) &
           -Half*G1(kP,lT)*G1(lT,jQ)*G2(jP,kQ,iP,iQ)*E(lT) &
           -Half*G1(lT,iQ)*G1(iP,jQ)*G2(jP,lT,kP,kQ)*E(lT) &
           -Half*G1(lT,iQ)*G1(iP,kQ)*G2(kP,lT,jP,jQ)*E(lT) &
           -Half*G1(lT,jQ)*G1(jP,iQ)*G2(iP,lT,kP,kQ)*E(lT) &
           -Half*G1(lT,jQ)*G1(jP,kQ)*G2(kP,lT,iP,iQ)*E(lT) &
           -Half*G1(lT,kQ)*G1(kP,iQ)*G2(iP,lT,jP,jQ)*E(lT) &
           -Half*G1(lT,kQ)*G1(kP,jQ)*G2(jP,lT,iP,iQ)*E(lT)
  !     // G(I,k)G(J,l)G(KL,ij)[12]
  CU4F3H = CU4F3H &
           -Half*G1(iP,kQ)*G1(jP,lT)*G2(kP,iQ,lT,jQ)*E(lT) &
           -Half*G1(iP,jQ)*G1(kP,lT)*G2(jP,iQ,lT,kQ)*E(lT) &
           -Half*G1(jP,iQ)*G1(kP,lT)*G2(iP,jQ,lT,kQ)*E(lT) &
           -Half*G1(iP,lT)*G1(jP,kQ)*G2(lT,iQ,kP,jQ)*E(lT) &
           -Half*G1(iP,lT)*G1(kP,jQ)*G2(lT,iQ,jP,kQ)*E(lT) &
           -Half*G1(jP,lT)*G1(kP,iQ)*G2(lT,jQ,iP,kQ)*E(lT) &
           -Half*G1(iP,jQ)*G1(lT,kQ)*G2(jP,iQ,kP,lT)*E(lT) &
           -Half*G1(iP,kQ)*G1(lT,jQ)*G2(kP,iQ,jP,lT)*E(lT) &
           -Half*G1(jP,iQ)*G1(lT,kQ)*G2(iP,jQ,kP,lT)*E(lT) &
           -Half*G1(jP,kQ)*G1(lT,iQ)*G2(kP,jQ,iP,lT)*E(lT) &
           -Half*G1(kP,iQ)*G1(lT,jQ)*G2(iP,kQ,jP,lT)*E(lT) &
           -Half*G1(kP,jQ)*G1(lT,iQ)*G2(jP,kQ,iP,lT)*E(lT)
  !     // G(I,j)G(J,i)G(K,k)G(L,l)[6]
  CU4F3H = CU4F3H &
           -Three*G1(iP,lT)*G1(lT,iQ)*G1(jP,jQ)*G1(kP,kQ)*E(lT) &
           -Three*G1(jP,lT)*G1(lT,jQ)*G1(iP,iQ)*G1(kP,kQ)*E(lT) &
           -Three*G1(kP,lT)*G1(lT,kQ)*G1(iP,iQ)*G1(jP,jQ)*E(lT)
  !     // G(I,j)G(J,k)G(K,i)G(L,l)[8]
  CU4F3H = CU4F3H &
           +OneHalf*G1(iP,jQ)*G1(jP,lT)*G1(lT,iQ)*G1(kP,kQ)*E(lT) &
           +OneHalf*G1(iP,kQ)*G1(kP,lT)*G1(lT,iQ)*G1(jP,jQ)*E(lT) &
           +OneHalf*G1(jP,kQ)*G1(kP,lT)*G1(lT,jQ)*G1(iP,iQ)*E(lT) &
           +OneHalf*G1(iP,lT)*G1(lT,jQ)*G1(jP,iQ)*G1(kP,kQ)*E(lT) &
           +OneHalf*G1(iP,lT)*G1(lT,kQ)*G1(kP,iQ)*G1(jP,jQ)*E(lT) &
           +OneHalf*G1(jP,lT)*G1(lT,kQ)*G1(kP,jQ)*G1(iP,iQ)*E(lT)
  !     // G(I,k)G(J,l)G(K,i)G(L,j)[3]
  CU4F3H = CU4F3H &
           +OneHalf*G1(iP,kQ)*G1(jP,lT)*G1(kP,iQ)*G1(lT,jQ)*E(lT) &
           +OneHalf*G1(iP,jQ)*G1(kP,lT)*G1(jP,iQ)*G1(lT,kQ)*E(lT) &
           +OneHalf*G1(iP,lT)*G1(jP,kQ)*G1(lT,iQ)*G1(kP,jQ)*E(lT)
  !     // G(I,j)G(J,k)G(K,l)G(L,i)[6]
  CU4F3H = CU4F3H &
           -0.75_wp*G1(iP,jQ)*G1(jP,kQ)*G1(kP,lT)*G1(lT,iQ)*E(lT) &
           -0.75_wp*G1(iP,kQ)*G1(kP,jQ)*G1(jP,lT)*G1(lT,iQ)*E(lT) &
           -0.75_wp*G1(iP,jQ)*G1(jP,lT)*G1(lT,kQ)*G1(kP,iQ)*E(lT) &
           -0.75_wp*G1(iP,kQ)*G1(kP,lT)*G1(lT,jQ)*G1(jP,iQ)*E(lT) &
           -0.75_wp*G1(iP,lT)*G1(lT,jQ)*G1(jP,kQ)*G1(kP,iQ)*E(lT) &
           -0.75_wp*G1(iP,lT)*G1(lT,kQ)*G1(kP,jQ)*G1(jP,iQ)*E(lT)
end do

return

end function CU4F3H
