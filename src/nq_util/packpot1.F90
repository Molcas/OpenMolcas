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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine PackPot1(Packed,Full,nPack,Factor)

use nq_Info, only: mIrrep, mOrb, OffOrb2, OffOrbTri, nPot1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPack
real(kind=wp), intent(out) :: Packed(nPack)
real(kind=wp), intent(in) :: Full(nPot1), Factor
integer(kind=iwp) :: iIrrep, iOff1, IOff2, nOrbs, p, q

do iIrrep=0,mIrrep-1
  nOrbs = mOrb(iIrrep)
  IOff1 = OffOrbTri(iIrrep)
  IOff2 = OffOrb2(iIrrep)
  do P=1,nOrbs
    do Q=1,P
      Packed(IOff1+(P-1)*P/2+Q) = Factor*(Full(IOff2+(P-1)*nOrbs+Q)+Full(IOff2+(Q-1)*nOrbs+P))
    end do
  end do
end do

return

end subroutine PackPot1
