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

subroutine mk_G(G,GInv,nDimBC)

use Slapaf_Info, only: Curvilinear, Degen, dMass, Smmtrc, User_Def
use Constants, only: Zero, One, uToau
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDimBC
real(kind=wp), intent(out) :: G(nDimBC,nDimBC), GInv(nDimBC,nDimBC)
integer(kind=iwp) :: i, ii, ix, nsAtom
logical(kind=iwp) :: Auto

Auto = .not. User_Def
nsAtom = size(Smmtrc,2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the mass tensor

G(:,:) = Zero
GInv(:,:) = Zero
ii = 0
do i=1,nsAtom
  do ix=1,3
    if (Smmtrc(ix,i)) then
      ii = ii+1
      if (Auto .and. (.not. Curvilinear)) then
        G(ii,ii) = Degen(ix,i)/dMass(i)
      else
        G(ii,ii) = One/(Degen(ix,i)*dMass(i))
      end if
      GInv(ii,ii) = One/(G(ii,ii)*uToau)
    end if
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('G (cartesian)',' ',G,nDimBC,nDimBC)
call RecPrt('G-1 (cartesian)',' ',GInv,nDimBC,nDimBC)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine mk_G
