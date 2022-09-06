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
! Copyright (C) 1990,1991,1993,1998,2005, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

subroutine Post_2Center_RI(A_Diag)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!             Modified to 2-center ERIs for RI June '05                *
!***********************************************************************

use Basis_Info, only: nBas_Aux
use Wrj12, only: Lu_A, Lu_Q, nChV
use Gateway_global, only: force_out_of_core
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: A_Diag(:)
integer nDmA(0:7), nDmB(0:7)
logical Out_of_Core
character Name_Q*6
real*8, allocatable :: Scr(:), X(:), Z(:)
integer, allocatable :: iDiag(:)
real*8, allocatable, target :: Am(:), Qm(:), A_k(:), Q_k(:)
real*8, pointer :: A_l(:) => null(), Q_l(:) => null()
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine SORT_mat(irc,nDim,nVec,iD_A,nSym,lu_A0,mode,lScr,Scr,Diag)
    integer irc
    integer nSym
    integer nDim(nSym)
    integer nVec(nSym)
    integer iD_A(*)
    integer lu_A0(nSym)
    character(LEN=7) mode
    integer lScr
    real*8 Scr(lScr)
    real*8, optional :: Diag(*)
  end subroutine SORT_mat
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *

nScr = 0
nBfn2 = 0
nBfnTot = 0
do iIrrep=0,nIrrep-1
  lJ = nBas_Aux(iIrrep)
  if (iIrrep == 0) lJ = lJ-1
  nDmA(iIrrep) = lJ
  nDmB(iIrrep) = 0
  nScr = max(nScr,3*lJ)
  nBfn2 = nBfn2+lJ**2
  nBfnTot = nBfnTot+lJ
end do
nA_Diag = nBfnTot

call mma_maxDBLE(MaxMem)
!                                                                      *
!***********************************************************************
!                                                                      *
! Fill in the lower part of the A matrix as it is stored on disk.

do iIrrep=0,nIrrep-1
  nB = nBas_Aux(iIrrep)
  if (iIrrep == 0) nB = nB-1 ! subtract dummy af
  call Square_A(Lu_A(iIrrep),nB,MaxMem,Force_Out_of_Core)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Pivoting of the A matrix

call mma_allocate(iDiag,nA_Diag,Label='iDiag')
call mma_maxDBLE(MaxMem2)

if (Force_Out_of_Core) MaxMem2 = 3*nBfnTot
!lScr = Min(MaxMem2,nScr)
lScr = max(MaxMem2-(nScr/3),nScr)
call mma_allocate(Scr,lScr,Label='Scr')

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_A,'GePivot',lScr,Scr,Diag=A_Diag)
ichk = 0
do iIrrep=0,nIrrep-1
  nChV(iIrrep) = nDmB(iIrrep)
  ichk = ichk+min(1,nDmA(iIrrep)-nDmB(iIrrep))
end do
if (ichk /= 0) then
  write(6,*)
  write(6,*) 'Post_2Center_RI'
  write(6,*) 'Detected lin. dependences in the auxiliary basis.'
  write(6,'(A,8I6)') ' # of AuxBas before l. d. removal: ',(nDmA(i),i=0,nIrrep-1)
  write(6,'(A,8I6)') ' # of AuxBas after  l. d. removal: ',(nDmB(i),i=0,nIrrep-1)
  write(6,*)
end if

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_A,'DoPivot',lScr,Scr)

! Note: after the 'DoPivot' call to Sort_mat, the A-matrix is
!       no longer stored as squared but as upper-triangular

call mma_deallocate(Scr)
call mma_deallocate(A_Diag)

!***********************************************************************
!     A-vectors are now on disk. Go ahead and compute the Q-vectors!
!***********************************************************************

ThrQ = 1.0d-14 ! Threshold for Inv_Cho_Factor

do iIrrep=0,nIrrep-1
  !nB = nBas_Aux(iIrrep)
  !if (iIrrep == 0) nB = nB-1
  nB = nDmB(iIrrep)
  if (nB == 0) Go To 777
  nQm = nB*(nB+1)/2

  nXZ = nB
  nQm_full = nB*(nB+1)/2

  if (Force_Out_of_Core) MaxMem = (8*(2*nQm_full+5*nXZ))/10
  Out_of_Core = 2*nQm_full+5*nXZ > MaxMem

  if (Out_Of_Core) then
    mQm = (nQm*MaxMem-5*nXZ)/(2*nQm_full)
    a = One
    b = -Two*dble(mQm)
    mB = int(-a/Two+sqrt((a/Two)**2-b))
    kQm = mB*(mB+1)/2
    if (kQm > mQm) then
      call WarningMessage(2,'Error in Post_2Center_RI')
      write(6,*) 'kQm > mQm!'
      write(6,*) 'MaxMem=',MaxMem
      write(6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
      write(6,*) 'nB,mB=',nB,mB
      call Abend()
    end if
  else
    mB = nB
    kQm = nQm
  end if

  lQm = kQm
  lAm = lQm

  if (lQm < 1) then
    call WarningMessage(2,'Error in Post_2Center_RI')
    write(6,*) 'lQm < 1'
    call Abend()
  end if

  ! Some of memory for scratch arrays for Inv_Cho_Factor
  ! Allocate memory for the A- and Q-vectors and initialize.

  lScr = nXZ
  call mma_allocate(Scr,lScr,Label='Scr')
  call mma_allocate(Z,nXZ,Label='Z')
  call mma_allocate(X,nXZ,Label='X')
  call mma_allocate(Am,lAm,Label='Am')
  call mma_allocate(Qm,lQm,Label='Qm')
  call mma_allocate(A_k,nXZ,Label='A_k')
  call mma_allocate(Q_k,nXZ,Label='Q_k')

  Am(:) = Zero
  Qm(:) = Zero
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Process the A_ks to generate Q_ks.

  iSeed = 55+iIrrep
  Lu_Q(iIrrep) = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'QMAT',iIrrep
  call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)

  iAddr = 0
  nMem = mB
  do kCol=1,nB

    iAddr_ = iAddr
    if (kCol <= nMem) then
      ! Point to A_k in Am
      iOff = (kCol-1)*kCol/2
      A_l(1:kCol) => Am(iOff+1:iOff+kCol)
      if (kCol == 1) then
        nAm = nMem*(nMem+1)/2
        call dDaFile(Lu_A(iIrrep),2,Am,nAm,iAddr_)
      end if
      ! Point to Q_k in Qm
      Q_l(1:kCol) => Qm(iOff+1:iOff+kCol)
    else if (kCol > nMem) then
      ! Use special scratch for A_k
      A_l(1:kCol) => A_k(1:kCol)
      call dDaFile(Lu_A(iIrrep),2,A_l,kCol,iAddr_)
      ! Use special scratch for Q_k
      Q_l(1:kCol) => Q_k(1:kCol)
    end if

    LinDep = 2
    call Inv_Cho_Factor(A_l,kCol,Am,Qm,nMem,Lu_A(iIrrep),Lu_Q(iIrrep),Scr,lScr,Z,X,ThrQ,Q_l,LinDep)

    if (LinDep /= 0) then
      call WarningMessage(2,'Error in Post_2Center_RI')
      write(6,*) 'Inv_Cho_Factor found linear dependence!'
      call Abend()
    end if

    ! Write the new A/Q-vector to file

    iAddr_ = iAddr
    if (kCol == nMem) then
      nQm = kCol*(kCol+1)/2
      call dDaFile(Lu_Q(iIrrep),1,Qm,nQm,iAddr)
      call dDaFile(Lu_A(iIrrep),1,Am,nQm,iAddr_)
    else if (kCol > nMem) then
      call dDaFile(Lu_Q(iIrrep),1,Q_l,kCol,iAddr)
      call dDaFile(Lu_A(iIrrep),1,A_l,kCol,iAddr_)
    end if

  end do

  Q_l => null()
  A_l => null()
  call mma_deallocate(Q_k)
  call mma_deallocate(A_k)
  call mma_deallocate(Qm)
  call mma_deallocate(Am)
  call mma_deallocate(X)
  call mma_deallocate(Z)
  call mma_deallocate(Scr)
  call DaClos(Lu_A(iIrrep))
777 continue
end do ! iIrrep
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Sort the Q-matrix back to the original order.

call mma_maxDBLE(MaxMem2)

if (Force_Out_of_Core) MaxMem2 = 2*nBfnTot
lScr = min(MaxMem2,max(nBfn2,2*nBfnTot))
call mma_allocate(Scr,lScr,Label='Scr')

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_Q,'Restore',lScr,Scr)

! Note: after the 'Restore' call to Sort_mat, the Q-matrix is
!       no longer stored as upper-triangular but as RECTANGULAR
!       (nDmA,nDmB). The column index is still pivoted.

call mma_deallocate(Scr)
call mma_deallocate(iDiag)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
return

end subroutine Post_2Center_RI
