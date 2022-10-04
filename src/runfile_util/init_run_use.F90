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

subroutine init_run_use()

use Definitions, only: iwp

implicit none
#include "run_use_common.fh"
integer(kind=iwp) :: i

do i=1,nTocCA
  i_run_CA_used(i) = 0
end do
do i=1,nTocDA
  i_run_DA_used(i) = 0
end do
do i=1,nTocDS
  i_run_DS_used(i) = 0
end do
do i=1,nTocIA
  i_run_IA_used(i) = 0
end do
do i=1,nTocIS
  i_run_IS_used(i) = 0
end do

return

end subroutine init_run_use
