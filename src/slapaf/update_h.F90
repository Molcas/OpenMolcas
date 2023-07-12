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

subroutine Update_H(nWndw,H,nInter,nIter,iOptC,dq,g,iOptH,jPrint,GNrm,nsAtom,Store,AllowFindTS)

use Slapaf_Info, only: IRC, MF, Mode
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nWndw, nInter, nIter, iOptH, jPrint, nsAtom
real(kind=wp), intent(inout) :: H(nInter,nInter)
integer(kind=iwp), intent(inout) :: iOptC
real(kind=wp), intent(in) :: dq(nInter,nIter), g(nInter,nIter), GNrm
logical(kind=iwp), intent(in) :: Store, AllowFindTS
integer(kind=iwp) :: ierr, IterHess
logical(kind=iwp) :: Old_MF
real(kind=wp), allocatable :: Tmp(:,:)
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
! Update the Hessian. The procedure here loops over gradients
! of the nWndw last iteration. The reason for this is that the
! initial Hessian is reset to new values at each iteration.
! The anharmonic constants used here is the most recently
! updated version.

call DrvUpH(nWndw,nIter,H,nInter,dq,g,iOptH,IterHess)

call Chk4NAN(nInter*nInter,H,ierr)

if (ierr /= 0) call SysAbendMsg('Update_H','NaNs in Hessian','')
if (Store) call Put_dArray('Hss_upd',H,nInter**2)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Check if we have an old reaction mode

Old_MF = DDot_(3*nsAtom,MF,1,MF,1) /= Zero
Old_MF = Old_MF .and. (Mode /= 0) .and. (IRC == 0)
call mma_allocate(Tmp,3,nsAtom,Label='Tmp')
if (Old_MF) then
  if (jPrint >= 6) write(u6,*) ' Reading old reaction mode from disk'
  Tmp(:,:) = MF(:,:)
  Mode = 1  ! any number between 1 and nInter
  iOptC = ibset(iOptC,13)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Massage the new Hessian

call FixHess(H,nInter,iOptC,Tmp,GNrm,nsAtom,(IterHess == nIter),AllowFindTS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Store the new reaction mode if any!

if ((Mode > 0) .and. (Mode <= nInter)) then
  if (jPrint >= 6) write(u6,*) ' Storing new reaction mode on disk'
  MF(:,:) = Tmp(:,:)
end if
call mma_deallocate(Tmp)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (jPrint >= 99) call RecPrt('Update_H: Updated Hessian',' ',H,nInter,nInter)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Update_H
