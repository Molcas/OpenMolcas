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

subroutine Get_PUVXLen(NPUVX)
! Rewritten from mcpdft/alloc.f

use input_mclr, only: nSym, nOrb, nAsh

implicit none
integer NPUVX
integer iSp, iSq, iSr, iSs, nAq, iSpq, iSpqr, nAr, nAs, nOp, nRS

NPUVX = 0
do iSp=1,nSym
  nOp = NORB(iSp)
  do iSq=1,nSym
    nAq = NASH(iSq)
    iSpq = ieor(iSp-1,iSq-1)
    do iSr=1,nSym
      iSpqr = ieor(iSpq,iSr-1)+1
      nAr = NASH(iSr)
      do iSs=1,iSr
        if (iSpqr /= iSs) GO TO 11
        nAs = NASH(iSs)
        nRS = nAr*nAs
        if (iSs == iSr) nRS = (nAr+nAr**2)/2
        NPUVX = NPUVX+nOp*nAq*nRS
11      continue
      end do
    end do
  end do
end do

end subroutine Get_PUVXLen
