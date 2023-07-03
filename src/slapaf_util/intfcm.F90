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

use Slapaf_Parameters, only: lOld

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
real*8 rDum(1)
logical lOld_Implicit, Hess_Found, Found_IRC
real*8, allocatable :: Hess(:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine OldFcm(Hess,nQQ,Lbl)
    real*8, allocatable :: Hess(:)
    integer nQQ
    character*(*) Lbl
  end subroutine OldFcm
end interface

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
