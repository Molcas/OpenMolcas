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

subroutine CalcW(W,GDMat,PUVX,NPUVX,IndTUVX)

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: nNA
use input_mclr, only: nRoots, ntAsh
use Constants, only: Zero

implicit none
! Output
real*8, dimension(nTri_Elem(nRoots),nTri_Elem(nRoots)) :: W
! Input
integer NPUVX
real*8, dimension(nTri_Elem(nRoots),nnA,nnA) :: GDMat
real*8, dimension(NPUVX) :: PUVX
integer, dimension(ntAsh,ntAsh,ntAsh,ntAsh) :: IndTUVX
! Auxiliary Quantities
integer K, L, M, N, IKL, IMN, it, iu, iv, ix

do K=1,nRoots
  do L=1,K
    IKL = iTri(K,L)
    do M=1,nRoots
      do N=1,M
        IMN = iTri(M,N)
        W(IKL,IMN) = Zero
        do it=1,nnA
          do iu=1,nnA
            do iv=1,nnA
              do ix=1,nnA
                if (IndTUVX(it,iu,iv,ix) /= 0) W(IKL,IMN) = W(IKL,IMN)+GDMat(IKL,it,iu)*GDMat(IMN,iv,ix)*PUVX(IndTUVX(it,iu,iv,ix))
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

end subroutine CalcW
