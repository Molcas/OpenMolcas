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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine Tw_corr_drv(EOrb,nEO,CMO,nCMO,Ecorr)

use InfSCF, only: nnOc, nSym, nOcc, nDel, nOrb, nFro, nBas
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer nEO, nCMO
real*8 EOrb(nEO), CMO(nCMO), Ecorr
integer i, iOff, ipEOkk, ipEVir, iRC, iSym, jOff, jOkk, jOrb, jVir, kOff, nExt, nOkk
real*8, allocatable :: Eov(:)

call mma_Allocate(Eov,nEO,Label='Eov')

ipEOkk = 1
ipEVir = ipEOkk+nnOc
iOff = 0
jOff = 0
kOff = 0
do iSym=1,nSym
  nOkk = nFro(iSym)+nOcc(iSym,1)
  nExt = nBas(iSym)-nDel(iSym)-nOkk
  jOrb = 1+jOff
  jOkk = ipEOkk+iOff
  do i=0,nOkk-1
    Eov(jOkk+i) = EOrb(jOrb+i)
  end do
  jOrb = jOrb+nOkk
  jVir = ipEVir+kOff
  do i=0,nExt-1
    Eov(jVir+i) = EOrb(jOrb+i)
  end do
  iOff = iOff+nOkk
  jOff = jOff+nOrb(iSym)
  kOff = kOff+nExt
end do

call Tw_corr(irc,Ecorr,CMO,Eov(:ipEVir-1),Eov(ipEVir:))

call mma_deallocate(Eov)

return

end subroutine Tw_corr_drv
