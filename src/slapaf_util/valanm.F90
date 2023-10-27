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
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nAtom, nInter, nIter, nWndw
real(kind=wp), intent(in) :: BMx(3*nAtom,3*nAtom), Degen(3*nAtom), Cx(3*nAtom,nIter)
real(kind=wp), intent(inout) :: rInt(nInter,nIter)
character(len=*), intent(in) :: Label
integer(kind=iwp) :: iEnd, iIter, ij, iSt, j, NRHS
real(kind=wp), allocatable :: ScrC(:)

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
!NRHS = nIter
NRHS = iSt-iEnd+1
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

  call DGEMM_('T','N',nInter,NRHS,3*nAtom,One,BMx,3*nAtom,ScrC,3*nAtom,Zero,rInt(:,iEnd),nInter)

  call mma_deallocate(ScrC)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Eq_Solver('N',3*nAtom,nInter,NRHS,BMx,.false.,Degen,Cx(:,iEnd),rInt(:,iEnd))
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
