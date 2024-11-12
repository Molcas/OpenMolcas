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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N
real(kind=wp), intent(in) :: H, X, Y, Z, zJ, T, W(EXCH), thrs
complex(kind=wp), intent(in) :: DM(3,EXCH,EXCH), SM(3,EXCH,EXCH)
real(kind=wp), intent(out) :: ST(3)
logical(kind=iwp), intent(in) :: dbg
integer(kind=iwp) :: iopt
complex(kind=wp), allocatable :: DM_TMP(:,:,:), SM_TMP(:,:,:)

iopt = 2

if (iopt == 1) then

  if (dbg) write(u6,'(A)') 'mean_field: enter mean_field_exch'
  call mma_allocate(DM_TMP,3,N,N,label='DM_TMP')
  call mma_allocate(SM_TMP,3,N,N,label='SM_TMP')
  DM_TMP(:,:,:) = DM(:,1:N,1:N)
  SM_TMP(:,:,:) = SM(:,1:N,1:N)
  call mean_field_exch(N,H,X,Y,Z,zJ,T,W(1:N),thrs,DM_TMP,SM_TMP,ST)
  call mma_deallocate(DM_TMP)
  call mma_deallocate(SM_TMP)
  if (dbg) write(u6,'(A)') 'mean_field: exit mean_field_exch'

else if (iopt == 2) then

  if (dbg) write(u6,'(A)') 'mean_field: enter mean_field_all'
  call mean_field_all(EXCH,N,H,X,Y,Z,zJ,T,W,thrs,DM,SM,ST)
  if (dbg) write(u6,'(A)') 'mean_field: exit mean_field_all'

else

  write(u6,'(A)') 'MEAN_FIELD:  iopt is not defined:',iopt

end if

return

end subroutine mean_field
