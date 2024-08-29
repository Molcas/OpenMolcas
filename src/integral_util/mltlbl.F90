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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

function MltLbl(Lbl1,Lbl2)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MltLbl
integer(kind=iwp), intent(in) :: Lbl1, Lbl2
integer(kind=iwp) :: iIrrep, ijSym, jIrrep

MltLbl = 0
do iIrrep=0,nIrrep-1
  if (iand(Lbl1,2**iIrrep) == 0) cycle
  do jIrrep=0,nIrrep-1
    if (iand(Lbl2,2**jIrrep) == 0) cycle
    ijSym = ieor(iIrrep,jIrrep)
    if (iand(MltLbl,2**ijSym) == 0) MltLbl = MltLbl+2**ijSym
  end do
end do

return

end function MltLbl
