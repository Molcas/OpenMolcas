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

subroutine DoReadBPT2(iAOV,iAOst,iCmpa,kBasn,lBasn)
! Read back-transformed density elements of the Kth and Lth shells
! All elements of the Jth shell (auxiliary functions) are read
! Only for C1 symmetry

use Index_Functions, only: iTri
use SOAO_Info, only: iAOtSO
use pso_stuff, only: B_PT2, LuGamma2, nBasA, ReadBPT2
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iAOV(4), iAOSt(4), iCmpa(4), kBasn, lBasn
integer(kind=iwp) :: i3, i4, kAOk, kSO, kSO0, kSOk, lAOl, loc, lSO, lSO0, lSOl

B_PT2(:,:,:) = Zero

lSO0 = iAOtSO(iAOV(4)+1,0)+iAOst(4)-1
kSO0 = iAOtSO(iAOV(3)+1,0)+iAOst(3)-1
do i4=1,iCmpa(4)
  lSO = iAOtSO(iAOV(4)+i4,0)+iAOst(4)
  do i3=1,iCmpa(3)
    kSO = iAOtSO(iAOV(3)+i3,0)+iAOst(3)
    do lAOl=0,lBasn-1
      lSOl = lSO+lAOl
      do kAOk=0,kBasn-1
        kSOk = kSO+kAOk
        loc = iTri(kSOk,lSOl)
        read(unit=LuGAMMA2,rec=loc) B_PT2(1:nBasA,kSOk-kSO0,lSOl-lSO0)
      end do
    end do
  end do
end do
ReadBPT2 = .false.

end subroutine DoReadBPT2
