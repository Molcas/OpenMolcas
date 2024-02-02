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

subroutine mean_field(EXCH,N,H,X,Y,Z,zJ,T,W,thrs,dM,sM,ST,dbg)

implicit none
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: EXCH, N
real(kind=8), intent(in) :: H, X, Y, Z, zJ, T
real(kind=8), intent(in) :: thrs
real(kind=8), intent(in) :: W(EXCH)
complex(kind=8), intent(in) :: DM(3,EXCH,EXCH)
complex(kind=8), intent(in) :: SM(3,EXCH,EXCH)
logical, intent(in) :: dbg
! output
real(kind=8), intent(out) :: ST(3)
integer :: iopt

iopt = 2

if (iopt == 1) then

  if (dbg) write(6,'(A)') 'mean_field: enter mean_field_exch'
  call mean_field_exch(N,H,X,Y,Z,zJ,T,thrs,W(1:N),DM(1:3,1:N,1:N),SM(1:3,1:N,1:N),ST)
  if (dbg) write(6,'(A)') 'mean_field: exit mean_field_exch'

else if (iopt == 2) then

  if (dbg) write(6,'(A)') 'mean_field: enter mean_field_all'
  call mean_field_all(EXCH,N,H,X,Y,Z,zJ,T,thrs,W,DM,SM,ST)
  if (dbg) write(6,'(A)') 'mean_field: exit mean_field_all'

else

  write(6,'(A)') 'MEAN_FIELD:  iopt is not defined:',iopt

end if

return

end subroutine mean_field
