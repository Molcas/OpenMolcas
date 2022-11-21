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

subroutine InitDB(nMult,natom,nAtQM,nGrdPt,Cord,Grid,T,TT,TTT,Ext,DB,IsMM)
! Compute DB = d(ExtPot[TtT^-1]Tt)/dq

use espf_global, only: MxExtPotComp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nMult, natom, nAtQM, nGrdPt, IsMM(natom)
real(kind=wp), intent(in) :: Cord(3,natom), Grid(3,nGrdPt), T(nMult,nGrdPt), TT(nMult,nMult), TTT(nGrdPt,nMult), &
                             Ext(MxExtPotComp,natom)
real(kind=wp), intent(out) :: DB(nGrdPt,3,nAtQM)
real(kind=wp), allocatable :: DT(:,:,:,:), DTT(:,:,:,:), DTTT(:,:,:,:), DTTTT(:,:,:,:)

! Various derivatives

call mma_allocate(DT,nMult,nGrdPt,3,nAtQM,label='DT')
call mma_allocate(DTT,nMult,nMult,3,nAtQM,label='DTT')
call mma_allocate(DTTT,nMult,nGrdPt,3,nAtQM,label='DTTT')
call mma_allocate(DTTTT,nMult,nMult,3,nAtQM,label='DTTTT')

call CalcDT(nMult,nGrdPt,natom,nAtQM,IsMM,Cord,Grid,T,TT,DT,DTT,DTTTT,DTTT)

call mma_deallocate(DT)
call mma_deallocate(DTT)
call mma_deallocate(DTTTT)

! Finally dB = dTTT * V_ext + TTT * dV_ext

call CalcDB(nMult,nGrdPt,natom,nAtQM,IsMM,TTT,DTTT,Ext,DB)
call mma_deallocate(DTTT)

return

end subroutine InitDB
