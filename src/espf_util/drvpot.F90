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
! Copyright (C) 1991, Roland Lindh                                     *
!               2001, Hans-Joachim Werner                              *
!***********************************************************************

subroutine Drvpot(CCoor,opnuc,ncmp,ptchrg,ngrid,iaddpot)
!***********************************************************************
!                                                                      *
! Object: driver for computation of one-electron property matrices     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '91                                              *
!                                                                      *
!     Modified for Properties only by HJW Aug 2001                     *
!     Restricted to POT: Ignacio Fdez. Galvan, March 2019              *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nBas, nCnttp
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Integral_interfaces, only: int_kernel, int_mem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ncmp, ngrid, iaddpot
real(kind=wp), intent(in) :: CCoor(3,ngrid)
real(kind=wp), intent(_OUT_) :: opnuc(*)
real(kind=wp), intent(inout) :: ptchrg(*)
integer(kind=iwp) :: i, iIrrep, iopadr(1), jCnt, jCnttp, jxyz, mCnt, nc, nComp, ndc, nOrdOp, nSym, ntdg
real(kind=wp) :: dummy(1), rHrmt
logical(kind=iwp) :: Do_ESPF
character(len=8) :: Label
integer(kind=iwp), allocatable :: ip(:), kOper(:), lOper(:)
real(kind=wp), allocatable :: Centr(:,:), Dens(:), Nuc(:)
procedure(int_kernel) :: PotInt
procedure(int_mem) :: NAMem

!                                                                      *
!***********************************************************************
!                                                                      *
call IniSewM('mltpl',0)

call Set_Basis_Mode('Valence')
call Setup_iSD()
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
ntdg = 0
do iIrrep=0,nIrrep-1
  ntdg = ntdg+nTri_Elem(nBas(iIrrep))
end do
call DecideOnESPF(Do_ESPF)

call mma_allocate(Centr,3,S%mCentr)
ndc = 0
nc = 1
do jCnttp=1,nCnttp
  mCnt = dbsc(jCnttp)%nCntr
  if (dbsc(jCnttp)%Aux) mCnt = 0
  do jCnt=1,mCnt
    ndc = ndc+1
    do i=0,nIrrep/dc(ndc)%nStab-1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Centr(1:3,nc))
      nc = nc+1
    end do
    jxyz = jxyz+3
  end do
end do
nc = nc-1

nComp = 1
nOrdOp = 0
call mma_allocate(ip,nComp,label='ip')
call mma_allocate(lOper,nComp,label='lOper')
call mma_allocate(kOper,nComp,label='kOper')
Label = 'Pot '
if ((iaddpot <= 0) .and. (.not. Do_ESPF)) then
  call mma_allocate(Nuc,ngrid,label='Nuc')
  call Pot_nuc(CCoor,Nuc,ngrid)
else
  call mma_allocate(Nuc,ncmp,label='Nuc')
  Nuc(:) = Zero
end if
if (iaddpot < 0) then

  call mma_allocate(Dens,ntdg,Label='Dens')
  if (iaddpot == -1) then
    call Get_D1ao_Var(Dens,ntdg)
  else
    call Get_dArray_chk('D1ao',Dens,ntdg)
  end if
  call Drv1_Pot(Dens,CCoor,ptchrg,ngrid,1,0)
  call mma_deallocate(Dens)

  if (.not. Do_ESPF) then
    ptchrg(1:ngrid) = ptchrg(1:ngrid)+Nuc
    opnuc(1:ngrid) = Nuc
  end if
else
  lOper(1) = 2**nirrep-1
  kOper(1) = 0
  call OneEl(PotInt,NAMem,Label,ip,lOper,ncmp,CCoor,nOrdOp,Nuc,rHrmt,kOper,dummy,1,opnuc,iopadr,1,1,ptchrg,ngrid,iaddpot)
  if ((iaddpot == 0) .and. (.not. Do_ESPF)) opnuc(1) = Nuc(1)
end if
call mma_deallocate(ip)
call mma_deallocate(lOper)
call mma_deallocate(kOper)
call mma_deallocate(Nuc)

call mma_deallocate(Centr)
call Free_iSD()

return

end subroutine Drvpot
