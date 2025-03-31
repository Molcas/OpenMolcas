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

subroutine CMSRHSGDMat(GDMat)

use ipPage, only: W
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nNA, n2Dens, ipCI, n1Dens
use MCLR_Data, only: XISPSM
use input_mclr, only: State_Sym, nRoots, nCSF

implicit none
! Input
! Output
real*8, dimension(nRoots*(nRoots+1)/2,nnA,nnA) :: GDMat
! Auxiliary quantities
real*8, dimension(:), allocatable :: GDArray
real*8, dimension(n2dens) :: rdum
integer I, J, IOrb, JOrb, NIJ
! I:state index, "excited state" of state J when I /= J
! IOrb:row index,    orbital index for state I
! JOrb:column index, orbital index for state J
integer nConfL, nConfR, iL, iR
real*8, allocatable :: CIL(:), CIR(:)

! (D^IJ)_pq = <I|E_pq|J>, setting I>=J
!  <I|E_pq|J>=<J|E_qp|I>
call mma_allocate(GDArray,n1dens)
iL = state_sym
iR = state_sym
nConfR = max(ncsf(iR),nint(xispsm(iR,1)))
nConfL = max(ncsf(iL),nint(xispsm(iL,1)))
call mma_allocate(CIR,nConfR)
call mma_allocate(CIL,nConfL)
do I=1,nRoots
  call CSF2SD(W(ipCI)%A(1+(I-1)*ncsf(iL)),CIL,iL)
  do J=1,I
    call CSF2SD(W(ipCI)%A(1+(J-1)*ncsf(iR)),CIR,iR)
    call Densi2_mclr(1,GDArray,rdum,CIL,CIR,0,0,0,n1dens,n2dens)
    NIJ = I*(I-1)/2+J
    do IOrb=1,nnA
      do JOrb=1,nnA
        GDMat(NIJ,IOrb,JOrb) = GDArray((JOrb-1)*nnA+IOrb)
      end do
    end do
  end do
end do
call mma_deallocate(GDArray)
call mma_deallocate(CIL)
call mma_deallocate(CIR)

end subroutine CMSRHSGDMat
