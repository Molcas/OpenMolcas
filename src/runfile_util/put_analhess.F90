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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Put_AnalHess
!
!> @brief
!>   Write the symmetry blocked nuclear Hessian in Cartesian coordinates on the runfile
!> @author R. Lindh
!>
!> @details
!> The utility will write the symmetry blocked nuclear Hessian in Cartesian coordinates on the runfile.
!>
!> @param[in] AnalHess  Array with the symmetry blocked nuclear Hessian in Cartesian coordinates
!> @param[in] nAnalHess Size of the array \p AnalHess
!***********************************************************************

subroutine Put_AnalHess(AnalHess,nAnalHess)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAnalHess
real(kind=wp), intent(in) :: AnalHess(nAnalHess)
integer(kind=iwp) :: inLoop, irderr, iSI1(7), iTmp
logical(kind=iwp) :: Found
character(len=88) :: Label

call Put_dArray('Analytic Hessian',AnalHess,nAnalHess)

! Add the iteration number corresponding to this Hessian
iSI1(2) = 0
call Qpg_iArray('Slapaf Info 1',Found,iTmp)
if (Found) call Get_iArray('Slapaf Info 1',iSI1,7)
call Getenvf('MOLCAS_ITER',Label)
read(Label,'(I80)') iTmp
call Getenvf('EMIL_InLoop',Label)
read(Label,*,IOStat=irderr) inLoop
if (irderr /= 0) inLoop = 0
if (inLoop < 1) iTmp = 0
if (iTmp == 0) then
  call Put_iScalar('HessIter',iTmp)
else
  call Put_iScalar('HessIter',iSI1(2)+1)
end if

return

end subroutine Put_AnalHess
