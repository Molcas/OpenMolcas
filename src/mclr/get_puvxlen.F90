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

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use input_mclr, only: nAsh, nOrb, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NPUVX
integer(kind=iwp) :: iSp, iSpq, iSpqr, iSq, iSr, iSs, nAq, nAr, nAs, nOp, nRS

NPUVX = 0
do iSp=1,nSym
  nOp = NORB(iSp)
  do iSq=1,nSym
    nAq = NASH(iSq)
    iSpq = Mul(iSp,iSq)
    do iSr=1,nSym
      iSpqr = Mul(iSpq,iSr)
      nAr = NASH(iSr)
      do iSs=1,iSr
        if (iSpqr /= iSs) exit
        nAs = NASH(iSs)
        nRS = nAr*nAs
        if (iSs == iSr) nRS = nTri_Elem(nAr)
        NPUVX = NPUVX+nOp*nAq*nRS
      end do
    end do
  end do
end do

end subroutine Get_PUVXLen
