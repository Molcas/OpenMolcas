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

subroutine ValANM(nAtom,nInter,nIter,Bmx,Degen,rInt,Cx,Label,nWndw)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAtom, nInter, nIter, nWndw
real(kind=wp) :: BMx(3*nAtom,3*nAtom), Degen(3*nAtom), rInt(nInter,nIter), Cx(3*nAtom,nIter)
character(len=*) :: Label
integer(kind=iwp) :: iEnd, iIter, ij, iSt, j, M, N, NRHS
real(kind=wp), allocatable :: ScrC(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Values:    q=BuX
! Gradients: g=Bu(dE/dX)
!                                                                      *
!***********************************************************************
!                                                                      *
iSt = nIter
iEnd = iSt-min(nIter,nWndw+1)+1
if (Label == 'Values') then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Values:    q=BX

  call mma_allocate(ScrC,3*nAtom*(iSt-iEnd+1),Label='ScrC')

  do iIter=iSt,iEnd,-1
    do j=1,3*nAtom
      ij = (iIter-iEnd)*3*nAtom+j
      ScrC(ij) = Cx(j,iIter)*Degen(j)
    end do
  end do

  call DGEMM_('T','N',nInter,iSt-iEnd+1,3*nAtom,One,BMx,3*nAtom,ScrC,3*nAtom,Zero,rInt(1,iEnd),nInter)

  call mma_deallocate(ScrC)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  M = 3*nAtom
  N = nInter
  !NRHS = nIter
  NRHS = iSt-iEnd+1
  call Eq_Solver('N',M,N,NRHS,BMx,.false.,Degen,Cx(1,iEnd),rInt(1,iEnd))
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,'(A)') ' In ValANM: New '
call RecPrt(Label,' ',rInt(1,iEnd),nInter,iSt-iEnd+1)
#endif

return

end subroutine ValANM
