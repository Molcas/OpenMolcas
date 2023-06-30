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

subroutine DrvUpH(nWndw,nIter,H,nInter,dq,g,iOptH,IterHess)

use Slapaf_Info, only: mRowH
use NewH_mod, only: DiagMM, UpdMask
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nWndw, nIter, nInter, iOptH
real(kind=wp), intent(inout) :: H(nInter,nInter)
real(kind=wp), intent(in) :: dq(nInter,nIter), g(nInter,nIter+1)
integer(kind=iwp), intent(inout) :: IterHess
integer(kind=iwp) :: i, iSt, lIter
logical(kind=iwp) :: DoMask, Found
!#define _DEBUGPRINT_

!                                                                      *
!***********************************************************************
!                                                                      *
iSt = max(2,nIter-(nWndw-1))
call Qpg_iScalar('HessIter',Found)
if (Found) then
  call Get_iScalar('HessIter',IterHess)
  iSt = max(iSt,IterHess+1)
else
  IterHess = 0
end if
if (allocated(mRowH)) iSt = max(iSt,size(mRowH)+2)
#ifdef _DEBUGPRINT_
Lu = u6
write(Lu,*) 'DrvUpH: iSt,kIter=',iSt,nIter
call RecPrt('DrvUpH: Initial Hessian',' ',H,nInter,nInter)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. btest(iOptH,3)) then
  write(Lu,*)
  if (nIter < iSt) then
    write(Lu,*) 'No update of Hessian on the first iteration'
  else
    write(Lu,'(A,30I3)') 'Hessian update from points:',(lIter,lIter=iSt-1,nIter)
  end if
  write(Lu,*)
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
DoMask = .false.
if (allocated(UpdMask)) then
  if (size(UpdMask) == nInter) DoMask = .true.
end if
if (DoMask) then
  do i=1,nInter
    if (UpdMask(i) /= 0) then
      H(i,:) = Zero
      H(:,i) = Zero
      H(i,i) = DiagMM
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Update the Hessian over the window

#ifdef _DEBUGPRINT_
call RecPrt('DrvUpH: Initial Hessian',' ',H,nInter,nInter)
#endif
do lIter=iSt,nIter
# ifdef _DEBUGPRINT_
  write(Lu,*) 'DrvUpH: Call NewH, lIter=',lIter
# endif
  call NewH(nInter,lIter,dq,g,H,iOptH,nIter)
end do
#ifdef _DEBUGPRINT_
call RecPrt('DrvUpH: Updated Hessian',' ',H,nInter,nInter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine DrvUpH
