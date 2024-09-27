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

use Symmetry_Info, only: Mul, nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MltLbl
integer(kind=iwp), intent(in) :: Lbl1, Lbl2
integer(kind=iwp) :: iIrrep, ijSym, jIrrep

MltLbl = 0
do iIrrep=1,nIrrep
  if (.not. btest(Lbl1,iIrrep-1)) cycle
  do jIrrep=1,nIrrep
    if (.not. btest(Lbl2,jIrrep-1)) cycle
    ijSym = Mul(iIrrep,jIrrep)
    if (.not. btest(MltLbl,ijSym-1)) MltLbl = MltLbl+2**(ijSym-1)
  end do
end do

return

end function MltLbl
