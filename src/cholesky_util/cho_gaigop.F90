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

subroutine Cho_GAiGOp(X,n,Op)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: n, X(n)
character(len=*) :: Op
#include "cho_para_info.fh"
integer(kind=iwp) :: iv, kv

if (Cho_Real_Par) then
  iv = 0
  do while (iv < n)
    kv = min(n-iv,32000000)
    call GAiGOp(X(iv+1),kv,Op)
    iv = iv+kv
  end do
end if

end subroutine Cho_GAiGOp
