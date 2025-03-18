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

use MCLR_Data, only: nNA, nDens2, ipMat
use input_mclr, only: nSym, nAsh, nBas, nIsh

implicit none
! Input
integer NPUVX
real*8 Coeff
real*8, dimension(nnA**2) :: D
real*8, dimension(NPUVX) :: PUVX
integer, dimension(nnA,nnA,nnA,nnA) :: IndTUVX
integer, dimension(nSym) :: Off_Ash
! Output
real*8, dimension(nDens2) :: Wop
! Auxiliaries
integer jSym, it, iu, t, u, v, x, pt, qu, iLoc1, iLoc2, jAsh
real*8 tempd1

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
      tempd1 = 0.0d0
      do v=1,nnA
        iLoc2 = (v-1)*nnA
        do x=1,nnA
          if (IndTUVX(t,u,v,x) /= 0) tempd1 = tempd1+D(iLoc2+x)*PUVX(IndTUVX(t,u,v,x))
        end do
      end do
      Wop(iLoc1+pt) = tempd1
    end do
  end do
end do

call DScal_(nDens2,Coeff,Wop,1)

end subroutine CalcWop
