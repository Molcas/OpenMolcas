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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
real*8 BMx(3*nAtom,3*nAtom), rInt(nInter,nIter), Degen(3*nAtom), Cx(3*nAtom,nIter)
character Label*(*)
real*8, allocatable :: ScrC(:)

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

  call DGEMM_('T','N',nInter,iSt-iEnd+1,3*nAtom,1.0d0,BMx,3*nAtom,ScrC,3*nAtom,0.0d0,rInt(1,iEnd),nInter)

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
write(6,'(A)') ' In ValANM: New '
call RecPrt(Label,' ',rInt(1,iEnd),nInter,iSt-iEnd+1)
#endif

return

end subroutine ValANM
