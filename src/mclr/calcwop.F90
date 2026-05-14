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
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CalcWop(Wop,D,PUVX,NPUVX,IndTUVX,Coeff,Off_Ash)

use MCLR_Data, only: ipMat, nDens, nNA
use input_mclr, only: nAsh, nBas, nIsh, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Wop(nDens)
integer(kind=iwp), intent(in) :: NPUVX, IndTUVX(nnA,nnA,nnA,nnA), Off_Ash(nSym)
real(kind=wp), intent(in) :: D(nnA,nnA), PUVX(NPUVX), Coeff
integer(kind=iwp) :: iLoc1, it, iu, jAsh, jSym, pt, qu, t, u, v, x
real(kind=wp) :: tempd1

do jSym=1,nSym
  jAsh = nAsh(jSym)
  if (jAsh == 0) cycle
  do iu=1,jAsh
    u = iu+off_Ash(jSym)
    qu = iu+nIsh(jSym)
    iLoc1 = (qu-1)*nBas(jSym)+ipMat(jSym,jSym)-1
    do it=1,jAsh
      t = it+off_Ash(jSym)
      pt = it+nIsh(jSym)
      tempd1 = Zero
      do v=1,nnA
        do x=1,nnA
          if (IndTUVX(t,u,v,x) /= 0) tempd1 = tempd1+D(x,v)*PUVX(IndTUVX(t,u,v,x))
        end do
      end do
      Wop(iLoc1+pt) = tempd1
    end do
  end do
end do

Wop(:) = Coeff*Wop(:)

end subroutine CalcWop
