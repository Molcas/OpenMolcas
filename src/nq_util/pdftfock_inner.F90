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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************
subroutine PDFTFock_Inner(Fock,Kern,MO1,MO2,mGrid)

use nq_Info

! Input
integer mGrid
real*8, dimension(mGrid*nOrbt) :: MO1, MO2
real*8, dimension(mGrid) :: Kern
! Output
real*8, dimension(nPot1) :: Fock
! Intermediate
real*8, dimension(mGrid*nOrbt) :: KernMO
integer iGrid, iIrrep, iOff1, iOff2

call dcopy_(mGrid*nOrbt,MO1,1,KernMO,1)

do iGrid=1,mGrid
  call DScal_(nOrbt,Kern(iGrid),KernMO(iGrid),mGrid)
end do

do iIrrep=0,mIrrep-1
  IOff1 = OffOrb(iIrrep)*mGrid+1
  IOff2 = OffOrb2(iIrrep)+1
  call DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,1.0d0,KernMO(IOff1),mGrid,MO2(IOff1),mGrid,1.0d0,Fock(iOff2),mOrb(iIrrep))
end do

return

end subroutine PDFTFock_Inner
