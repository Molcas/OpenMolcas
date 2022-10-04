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
!  Get_AnalHess
!
!> @brief
!>   Read the the symmetry-blocked nuclear Hessian from the runfile and return a
!>   pointer to the array's location in \c Work and the length of the array
!> @author R. Lindh
!>
!> @details
!> The utility will read the symmetry-blocked nuclear Hessian from the runfile and
!> return a pointer and the length of the array.
!>
!> @param[out] Hess  Array with the symmetry-blocked nuclear Hessian
!>                   in Cartesian coordinates
!> @param[out] nHess Size of the array of the symmetry-blocked nuclear Hessian
!***********************************************************************

subroutine Get_AnalHess(Hess,nHess)

implicit real*8(A-H,O-Z)
real*8 Hess(nHess)
character*24 Label
logical Found

Label = 'Analytic Hessian'
call qpg_dArray(Label,Found,nAnalHess)
if ((.not. Found) .or. (nAnalHess == 0)) then
  write(6,*) 'Get_AnalHess: Hessian not found!'
  call Abend()
end if
if (nAnalHess /= nHess) then
  write(6,*) 'Get_AnalHess: nAnalHess/=nHess'
  write(6,*) 'nAnalHess=',nAnalHess
  write(6,*) 'nHess=',nHess
  call Abend()
end if
call Get_dArray(Label,Hess,nHess)

return

end subroutine Get_AnalHess
