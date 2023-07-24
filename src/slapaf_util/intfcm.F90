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

subroutine IntFcm(lOld_Implicit)
!***********************************************************************
!                                                                      *
! Object: to initialize the Hessian matrix for the first iteration.    *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May 1991                                                 *
!***********************************************************************

use Slapaf_Info, only: lOld
use Slapaf_procedures, only: OldFCM
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(inout) :: lOld_Implicit
integer(kind=iwp) :: nHess, nQQ
real(kind=wp) :: rDum(1)
logical(kind=iwp) :: Found_IRC, Hess_Found
real(kind=wp), allocatable :: Hess(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Read force constant matrix from old interphase

if (lOld) then

  ! Explicit request to use an old force constant matrix stored
  ! on an old runfile.

  call OldFcm(Hess,nQQ,'RUNOLD')

else

  ! If this is not an IRC calculation explore if the runfile
  ! contains a Hessian. If so, pull it off the runfile.

  call qpg_iScalar('IRC',Found_IRC)

  if (.not. Found_IRC) then
    call qpg_dArray('Hess',Hess_Found,nHess)

    if (Hess_Found .and. (nHess > 0)) then
      lOld_Implicit = .true.
      call OldFcm(Hess,nQQ,'RUNFILE')
    end if

  end if

end if

if ((.not. lOld) .and. lOld_Implicit) lOld = .true.

if (lOld) then
# ifdef _DEBUGPRINT_
  call RecPrt('IntFcm: Final Hessian',' ',Hess,nQQ,nQQ)
# endif
  call Put_dArray('Hss_Q',Hess,nQQ**2)
  call Put_dArray('Hss_upd',rDum,0)
  call mma_deallocate(Hess)
end if

return

end subroutine IntFcm
