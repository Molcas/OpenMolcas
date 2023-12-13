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

subroutine WelMem( &
#                 define _CALLING_
#                 include "mem_interface.fh"
                 )

use Definitions, only: iwp

implicit none
#include "mem_interface.fh"
integer(kind=iwp) :: i, jsum, k

#include "macros.fh"
unused_var(lr)

k = la+lb
jsum = 1
do i=1,k
  jsum = jsum+3**i
end do
nHer = 1
Mem = jsum+max((k+1)*(k/2+1)*(k/4+1)+1,9+3**k,5)

return

end subroutine WelMem
