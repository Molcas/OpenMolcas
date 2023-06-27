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

use Slapaf_Info, only: MF
use Slapaf_Parameters, only: IRC, Mode

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
#include "real.fh"
real*8 H(nInter,nInter), dq(nInter,nIter), g(nInter,nIter)
logical Old_MF, Store, AllowFindTS
real*8, allocatable :: Tmp(:,:)

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
  if (jPrint >= 6) write(6,*) ' Reading old reaction mode from disk'
  Tmp(:,:) = MF(:,:)
  Mode = 1  ! any number between 1 and nInter
  iOptC = ior(iOptC,8192)
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
  if (jPrint >= 6) write(6,*) ' Storing new reaction mode on disk'
  MF(:,:) = Tmp(:,:)
end if
call mma_deallocate(Tmp)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (jPrint >= 99) then
  call RecPrt('Update_H: Updated Hessian',' ',H,nInter,nInter)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Update_H
