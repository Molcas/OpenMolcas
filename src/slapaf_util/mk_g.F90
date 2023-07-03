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

use Slapaf_Info, only: dMass, Degen, Smmtrc
use Slapaf_Parameters, only: Curvilinear, User_Def

implicit real*8(a-h,o-z)
#include "real.fh"
#include "constants2.fh"
real*8 G(nDimBC,nDimBC), GInv(nDimBC,nDimBC)
logical Auto

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
      GInv(ii,ii) = One/(G(ii,ii)*UTOAU)
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
