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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_ComputeCoulombIntermediates(Timing,nD,ip_DBlocks,ip_V,CNorm)
! Thomas Bondo Pedersen, October 2010.
!
! Purpose: Compute Coulomb intermediates
!
!      V(J) = sum_uv C(uv,J)*D(uv)
!
!      using LDF fitting coefficients.
!      Frobenius norms of the fitting
!      coefficients are computed as a byproduct, stored as:
!      CNorm(4*(AB-1)+1)
!          =sum_uAvBJ sqrt[C(uAvB,J)**2]       {all J}
!      CNorm(4*(AB-1)+2)
!          =sum_uAvBJA sqrt[C(uAvB,JA)**2]     {J on A}
!      CNorm(4*(AB-1)+3)
!          =sum_uAvBJB sqrt[C(uAvB,JB)**2]     {J on B}
!      CNorm(4*(AB-1)+4)
!          =sum_uAvBJAB sqrt[C(uAvB,JAB)**2]   {J on AB (2CF)}
!
! FOR BLOCKED VERSION: Call LDF_ComputeCoulombIntermediates0

#ifdef _MOLCAS_MPP_
use Para_Info, only: nProcs, Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "ldf_atom_pair_info.fh"
logical(kind=iwp), intent(in) :: Timing
integer(kind=iwp), intent(in) :: nD, ip_DBlocks(nD), ip_V(nD)
real(kind=wp), intent(out) :: CNorm(4*NumberOfAtomPairs)
logical(kind=iwp) :: doNorm
real(kind=wp) :: tC1, tC2, tW1, tW2
integer(kind=iwp) :: TaskListID, iD, l_C, iAtomPair, iAtom, jAtom, nAtom, nuv, M, ipD, ipV, ipC
real(kind=wp), allocatable :: LDFCBlk(:)
logical(kind=iwp), external :: Rsv_Tsk
integer(kind=iwp), external :: LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD, LDF_nAtom
real(kind=wp), external :: ddot_
#include "WrkSpc.fh"

if (Timing) call CWTIme(tC1,tW1)

! Initialize V arrays
do iD=1,nD
  call LDF_ZeroAuxBasVector(ip_V(iD))
end do

! Allocate array for storing coefficients
l_C = 0
do iAtomPair=1,NumberOfAtomPairs
  iAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+1)
  jAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+2)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  M = LDF_nBasAux_Pair_wLD(iAtomPair)
  l_C = max(l_C,nuv*M)
end do
call mma_allocate(LDFCBlk,l_C,label='LDFCBlk')

! compute norm?
doNorm = .true.
#ifdef _MOLCAS_MPP_
! Init norm array
if ((nProcs > 1) .and. Is_Real_Par()) then
  if (doNorm) CNorm(:) = Zero
end if
#endif

! Compute V
nAtom = LDF_nAtom()
call Init_Tsk(TaskListID,NumberOfAtomPairs)
do while (Rsv_Tsk(TaskListID,iAtomPair))
  call LDF_CIO_ReadC_wLD(iAtomPair,LDFCBlk,l_C)
  iAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+1)
  jAtom = iWork(ip_AP_Atoms-1+2*(iAtomPair-1)+2)
  nuv = LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
  ipC = 1
  M = LDF_nBasAux_Atom(iAtom)
  if (doNorm) then
    CNorm(4*(iAtomPair-1)+1) = sqrt(dDot_(nuv*LDF_nBasAux_Pair_wLD(iAtomPair),LDFCBlk,1,LDFCBlk,1))
    CNorm(4*(iAtomPair-1)+2) = sqrt(dDot_(nuv*M,LDFCBlk(ipC),1,LDFCBlk(ipC),1))
  end if
  do iD=1,nD
    ipD = iWork(ip_DBlocks(iD)-1+iAtomPair)
    ipV = iWork(ip_V(iD)-1+iAtom)
    call dGeMV_('T',nuv,M,One,LDFCBlk(ipC),nuv,Work(ipD),1,One,Work(ipV),1)
  end do
  if (jAtom /= iAtom) then
    ipC = ipC+nuv*M
    M = LDF_nBasAux_Atom(jAtom)
    if (doNorm) then
      CNorm(4*(iAtomPair-1)+3) = sqrt(dDot_(nuv*M,LDFCBlk(ipC),1,LDFCBlk(ipC),1))
    end if
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+iAtomPair)
      ipV = iWork(ip_V(iD)-1+jAtom)
      call dGeMV_('T',nuv,M,One,LDFCBlk(ipC),nuv,Work(ipD),1,One,Work(ipV),1)
    end do
  else
    if (doNorm) then
      CNorm(4*(iAtomPair-1)+3) = CNorm(4*(iAtomPair-1)+2)
    end if
  end if
  if (iWork(ip_AP_2CFunctions-1+2*(iAtomPair-1)+1) > 0) then
    ipC = ipC+nuv*M
    M = iWork(ip_AP_2CFunctions-1+2*(iAtomPair-1)+1)
    if (doNorm) then
      CNorm(4*(iAtomPair-1)+4) = sqrt(dDot_(nuv*M,LDFCBlk(ipC),1,LDFCBlk(ipC),1))
    end if
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+iAtomPair)
      ipV = iWork(ip_V(iD)-1+nAtom+iAtomPair)
      call dGeMV_('T',nuv,M,One,LDFCBlk(ipC),nuv,Work(ipD),1,One,Work(ipV),1)
    end do
  else
    if (doNorm) then
      CNorm(4*(iAtomPair-1)+4) = Zero
    end if
  end if
end do
call Free_Tsk(TaskListID)
if (Timing) then
  call CWTIme(tC2,tW2)
  write(u6,'(A,2(1X,F12.2),A)') 'Time spent computing Coulomb (V) intermediates:   ',tC2-tC1,tW2-tW1,' seconds'
end if

! Deallocation
call mma_deallocate(LDFCBlk)

#ifdef _MOLCAS_MPP_
! Get complete V and norm on all nodes
if ((nProcs > 1) .and. Is_Real_Par()) then
  if (Timing) call CWTime(tC1,tW1)
  do iD=1,nD
    call LDF_P_AddAuxBasVector(ip_V(iD))
  end do
  if (doNorm) then
    call GAdGOp(CNorm,4*NumberOfAtomPairs,'+')
  end if
  if (Timing) then
    call CWTime(tC2,tW2)
    write(u6,'(A,2(1X,F12.2),A)') 'Parallel overhead for Coulomb (V) intermediates:  ',tC2-tC1,tW2-tW1,' seconds'
  end if
end if
#endif

end subroutine LDF_ComputeCoulombIntermediates
