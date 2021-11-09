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

implicit real*8(a-h,o-z)
implicit integer(i-n)
integer Graph2(n_max+1,n_max+1,nOsc)
integer nTabDim(0:lnTabDim), nMat0(nOsc)
integer nvTabDim
#include "WrkSpc.fh"
#include "io_mula.fh"

! Initialize.
do i=1,nOsc
  nMat0(i) = 0
end do
iIndex = 0
nTabDim(0) = iIndex
call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
call GetMem('iVec','Allo','INTE',ipiVec,nOsc)

! Macrocycle on iQuanta
nD_0 = 0
do iQuanta=1,n_max
  do jv=1,nOsc
    iWork(ipiVec+jv-1) = 0
  end do
  iQ = -1
  iWork(ipiVec) = -1
  call TabDim2_drv(iQuanta,nOsc,nd)
  call TabDim2_drv(iQuanta-1,nOsc,nvTabDim)
  nd = nd-nvTabDim

  ! Microcycle on iDet
  do iDet=1,nD
    iWork(ipiVec) = iWork(ipiVec)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) exit
        iQ = iQ-iWork(ipiVec+i-1)+1
        iWork(ipiVec+i-1) = 0
        iWork(ipiVec+i) = iWork(ipiVec+i)+1
      end do
    end if
    iWork(ipiVec+nOsc-1) = iQuanta-iq
    iDNR = iDetnr(iWork(ipiVec),Graph2,nosc,n_max)
    iDNR = iDNR-nD_0
    do iv=1,nOsc
      nMat0(iv) = iWork(ipiVec+iv-1)
    end do
    nTabDim(iDNR+nD_0) = iIndex
    call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
  end do
  nD_0 = nD_0+nD
end do
call GetMem('iVec','Free','INTE',ipiVec,nOsc)

return

end subroutine ISCD_MakenMat
