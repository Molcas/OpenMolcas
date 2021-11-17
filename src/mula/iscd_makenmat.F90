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

subroutine ISCD_MakenMat(n_max,nOsc,lNMAT,lnTabDim,Graph2,nTabDim,nMat0)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n_max, nOsc, lNMAT, lnTabDim, Graph2(n_max+1,n_max+1,nOsc)
integer(kind=iwp), intent(out) :: nTabDim(0:lnTabDim), nMat0(nOsc)
integer(kind=iwp) :: i, iDet, iDNR, iIndex, iQ, iQuanta, nd, nD_0, nvTabDim
integer(kind=iwp), allocatable :: iVec(:)
integer(kind=iwp), external :: iDetnr

! Initialize.
nMat0(:) = 0
iIndex = 0
nTabDim(0) = iIndex
call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
call mma_allocate(iVec,nOsc,label='iVec')

! Macrocycle on iQuanta
nD_0 = 0
do iQuanta=1,n_max
  iVec(:) = 0
  iQ = -1
  iVec(1) = -1
  call TabDim(iQuanta,nOsc,nd)
  call TabDim(iQuanta-1,nOsc,nvTabDim)
  nd = nd-nvTabDim

  ! Microcycle on iDet
  do iDet=1,nD
    iVec(1) = iVec(1)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) exit
        iQ = iQ-iVec(i)+1
        iVec(i) = 0
        iVec(i+1) = iVec(i+1)+1
      end do
    end if
    iVec(nOsc) = iQuanta-iq
    iDNR = iDetnr(iVec,Graph2,nosc,n_max)
    iDNR = iDNR-nD_0
    nMat0(:) = iVec
    nTabDim(iDNR+nD_0) = iIndex
    call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
  end do
  nD_0 = nD_0+nD
end do
call mma_deallocate(iVec)

return

end subroutine ISCD_MakenMat
