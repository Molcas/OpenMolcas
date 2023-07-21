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

subroutine Cho_GAdGOp(X,n,Op)

implicit none
integer n, iv, kv
real*8 X(n)
character*(*) Op
#include "cho_para_info.fh"

if (Cho_Real_Par) then
  iv = 0
  do while (iv < n)
    kv = min(n-iv,32000000)
    call GAdGOp(X(iv+1),kv,Op)
    iv = iv+kv
  end do
end if

end subroutine Cho_GAdGOp
