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

subroutine DecideOnDirect(CanDoDirect,FoundTwoEls,DoDirect,DoCholesky)
!-----------------------------------------------------------------------
!                                                                      -
!  All modules using two-electron integrals should call this routine   -
!  to decide whether to do the calculation integral-direct.            -
!                                                                      -
!  On input:     CanDoDirect : Direct capability of module (T/F).      -
!                FoundTwoEls : Two-electron integral file available.   -
!                                                                      -
!  On exit:      DoDirect:     Do calculation integral-direct.         -
!                DoCholesky:   Do Cholesky calculation                 -
!                                                                      -
!-----------------------------------------------------------------------

use Definitions, only: iwp, u6

implicit none
logical(kind=iwp), intent(in) :: CanDoDirect, FoundTwoEls
logical(kind=iwp), intent(out) :: DoDirect, DoCholesky
integer(kind=iwp) :: iOptSeward
logical(kind=iwp) :: AlwaysDirect, isDirect, Expert, NeverDirect

! Read option variable set in Seward
call Get_iScalar('System BitSwitch',iOptSeward)

call DecideOnCholesky(DoCholesky)
if (DoCholesky) then
  DoDirect = btest(iOptSeward,13)
  return
end if

isDirect = btest(iOptSeward,0)
Expert = btest(iOptSeward,1)
AlwaysDirect = isDirect .and. (.not. Expert)
NeverDirect = (.not. isDirect) .and. (.not. Expert)

if (AlwaysDirect) then
  if (.not. CanDoDirect) then
    write(u6,'(A)') ' Error, cannot do integral-direct calculation!'
    write(u6,'(A)') ' Turn off DIRECT option in SEWARD input.'
    call Abend()
  else
    DoDirect = .true.
  end if
else if ((.not. FoundTwoEls) .and. (NeverDirect .or. (.not. CanDoDirect))) then
  ! No integrals, no direct calculation (allowed) - we have to crash!
  write(u6,'(2A)') ' Two-electron integral file was not found!'
  if (CanDoDirect) write(u6,'(A)') ' Try keyword DIRECT in SEWARD.'
  call Abend()
else if (.not. FoundTwoEls) then
  DoDirect = .true.
else
  DoDirect = .false.
end if

return

end subroutine DecideOnDirect
