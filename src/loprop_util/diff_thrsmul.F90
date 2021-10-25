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

subroutine Diff_ThrsMul(MP,ThrsMul,ThrsMul_Clever,nAt,nij)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt, nij
real(kind=wp), intent(in) :: MP(nij,*), ThrsMul
real(kind=wp), intent(out) :: ThrsMul_Clever
integer(kind=iwp) :: iAtom, jAtom, k, kaunt, kauntA, kComp, l
real(kind=wp) :: dM, dMMax

dMMax = Zero
kauntA = 1
do iAtom=1,nAt
  do jAtom=1,iAtom
    kaunt = 1
    do l=0,1
      kComp = (l+1)*(l+2)/2
      do k=1,kComp
        dM = MP(kauntA,kaunt)
        if (abs(dM) > dMMax) dMMax = abs(dM)
        kaunt = kaunt+1
      end do
    end do
    kauntA = kauntA+1
  end do
end do

ThrsMul_Clever = dMMax*ThrsMul

return

end subroutine Diff_ThrsMul
