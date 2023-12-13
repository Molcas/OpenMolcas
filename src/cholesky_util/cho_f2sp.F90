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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

function Cho_F2SP(iSP)
!
! Thomas Bondo Pedersen, March 2006.
!
! Purpose: return reduced shell pair index for full shell pair index
!          iSP. If not found, 0 is returned.
!          Note: nnShl_SP is used to avoid problems in parallel runs
!          when swapping nnShl and nnShl_G. If properly set,
!          nnShl_SP = nnShl_G.

use Cholesky, only: iSP2F, nnShl_SP
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Cho_F2SP
integer(kind=iwp), intent(in) :: iSP
integer(kind=iwp) :: jSP

Cho_F2SP = 0
do jSP=1,nnShl_SP
  if (iSP2F(jSP) == iSP) then
    Cho_F2SP = jSP
    exit
  end if
end do

end function Cho_F2SP
