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

use nq_Info, only: mIrrep, mOrb, nOrbt, nPot1, OffOrb, OffOrb2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Fock(nPot1)
integer(kind=iwp), intent(in) :: mGrid
real(kind=wp), intent(in) :: Kern(mGrid), MO1(mGrid,nOrbt), MO2(mGrid,nOrbt)
integer(kind=iwp) :: iGrid, iIrrep, iOff1, iOff2
real(kind=wp), allocatable :: KernMO(:,:)

call mma_allocate(KernMO,mGrid,nOrbt,label='KernMO')

do iGrid=1,mGrid
  KernMO(iGrid,:) = MO1(iGrid,:)*Kern(iGrid)
end do

do iIrrep=0,mIrrep-1
  IOff1 = OffOrb(iIrrep)+1
  IOff2 = OffOrb2(iIrrep)+1
  call DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,One,KernMO(:,IOff1:),mGrid,MO2(:,IOff1:),mGrid,One,Fock(iOff2),mOrb(iIrrep))
end do

call mma_deallocate(KernMO)

return

end subroutine PDFTFock_Inner
