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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_P_Distrib_SP(iOpt,mySP,N_mySP)
!
! Thomas Bondo Pedersen, June 2007.
!
! Determine distribution of ShellPairs.
! iOpt=1: each node has same number of ShellPairs.
! iOpt=2: each node has same dimension.

use Cholesky, only: nnShl
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp), intent(_OUT_) :: mySP(*)
integer(kind=iwp), intent(out) :: N_mySP

N_mySP = 0
if (iOpt == 1) then
  call Cho_P_Distrib_Vec(1,nnShl,mySP,N_mySP)
else
  call Cho_P_Distrib_SP_byDim(mySP,N_mySP)
end if

end subroutine Cho_P_Distrib_SP
