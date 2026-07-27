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
<<<<<<< HEAD
=======
<<<<<<<< HEAD:src/casvb_util/gethess_cvb.F90
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine gethess_cvb(hess)
========
>>>>>>> upstream-openmolcas/master
! Copyright (C) 2023, Roland Lindh                                     *
!***********************************************************************

module Breit
<<<<<<< HEAD

use Constants, only: Zero
use Definitions, only: iwp, wp

implicit none
=======
>>>>>>>> upstream-openmolcas/master:src/rys_util/breit.F90

use casvb_global, only: nfr
use Definitions, only: wp, iwp

implicit none
<<<<<<<< HEAD:src/casvb_util/gethess_cvb.F90
real(kind=wp), intent(out) :: hess(nfr,nfr)
integer(kind=iwp) :: ivar

call unitmat(hess,nfr)
do ivar=1,nfr
  call hess_cvb(hess(:,ivar))
end do
========
>>>>>>> upstream-openmolcas/master
private

integer(kind=iwp) :: nComp = 1, nOrdOp = 0
real(kind=wp) :: D_tensor(3,3) = Zero
logical(kind=iwp) :: Do_BP_integrals = .false.
<<<<<<< HEAD

public :: D_tensor, Do_BP_integrals, nComp, nOrdOp

end module Breit
=======
>>>>>>>> upstream-openmolcas/master:src/rys_util/breit.F90

public :: D_tensor, Do_BP_integrals, nComp, nOrdOp

<<<<<<<< HEAD:src/casvb_util/gethess_cvb.F90
end subroutine gethess_cvb
========
end module Breit
>>>>>>>> upstream-openmolcas/master:src/rys_util/breit.F90
>>>>>>> upstream-openmolcas/master
