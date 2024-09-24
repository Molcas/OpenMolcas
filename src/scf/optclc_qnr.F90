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

subroutine OptClc_QNR(CInter,nCI,nD,Grd1,Xnp1,mOV,Ind,MxOptm,kOptim,kOV)

use LnkLst, only: LLGrad, LLx

implicit none
integer nCI, nD, mOV, MxOptm, kOptim, kOV(2)
real*8 CInter(nCI,nD), Grd1(mOV), Xnp1(mOV)
integer Ind(MxOptm)
interface
  subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,kOptim,kOV,LL,DD)
    implicit none
    integer nCI, nD, mOV, MxOptm, kOptim, kOV(2), LL
    real*8 CInter(nCI,nD), Array(mOV)
    integer Ind(MxOptm)
    real*8, optional :: DD
  end subroutine OptClc_X
end interface

call OptClc_X(CInter,nCI,nD,Grd1,mOV,Ind,MxOptm,kOptim,kOV,LLGrad)
call OptClc_X(CInter,nCI,nD,Xnp1,mOV,Ind,MxOptm,kOptim,kOV,LLx)

return

end subroutine OptClc_QNR
