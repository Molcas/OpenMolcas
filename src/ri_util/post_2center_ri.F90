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

use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas_Aux
use RI_glob, only: Lu_A, Lu_Q, nChV
use Gateway_global, only: force_out_of_core
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: A_Diag(*)
integer(kind=iwp) :: i, iAddr, iAddr_, ichk, iIrrep, iOff, irc, iSeed, kCol, kQm, lAm, LinDep, lJ, lQm, lScr, MaxMem, MaxMem2, mB, &
                     mQm, nA_Diag, nAm, nB, nBfn2, nBfnTot, nDmA(0:7), nDmB(0:7), nMem, nQm, nQm_full, nScr, nXZ
real(kind=wp) :: a, b, dum(1), ThrQ
logical(kind=iwp) :: Out_of_Core
character(len=6) :: Name_Q
integer(kind=iwp), allocatable :: iDiag(:)
real(kind=wp), allocatable :: Scr(:), X(:), Z(:)
real(kind=wp), allocatable, target :: A_k(:), Am(:), Q_k(:), Qm(:)
real(kind=wp), pointer :: A_l(:), Q_l(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
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

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_A,'GePivot',lScr,Scr,A_Diag)
ichk = 0
do iIrrep=0,nIrrep-1
  nChV(iIrrep) = nDmB(iIrrep)
  ichk = ichk+min(1,nDmA(iIrrep)-nDmB(iIrrep))
end do
if (ichk /= 0) then
  write(u6,*)
  write(u6,*) 'Post_2Center_RI'
  write(u6,*) 'Detected lin. dependences in the auxiliary basis.'
  write(u6,'(A,8I6)') ' # of AuxBas before l. d. removal: ',(nDmA(i),i=0,nIrrep-1)
  write(u6,'(A,8I6)') ' # of AuxBas after  l. d. removal: ',(nDmB(i),i=0,nIrrep-1)
  write(u6,*)
end if

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_A,'DoPivot',lScr,Scr,dum)

! Note: after the 'DoPivot' call to Sort_mat, the A-matrix is
!       no longer stored as squared but as upper-triangular

call mma_deallocate(Scr)

!***********************************************************************
!     A-vectors are now on disk. Go ahead and compute the Q-vectors!
!***********************************************************************

ThrQ = 1.0e-14_wp ! Threshold for Inv_Cho_Factor

do iIrrep=0,nIrrep-1
  !nB = nBas_Aux(iIrrep)
  !if (iIrrep == 0) nB = nB-1
  nB = nDmB(iIrrep)
  if (nB == 0) cycle
  nQm = nTri_Elem(nB)

  nXZ = nB
  nQm_full = nTri_Elem(nB)

  if (Force_Out_of_Core) MaxMem = (8*(2*nQm_full+5*nXZ))/10
  Out_of_Core = 2*nQm_full+5*nXZ > MaxMem

  if (Out_Of_Core) then
    mQm = (nQm*MaxMem-5*nXZ)/(2*nQm_full)
    a = One
    b = -Two*real(mQm,kind=wp)
    mB = int(-a*Half+sqrt((a*Half)**2-b))
    kQm = nTri_Elem(mB)
    if (kQm > mQm) then
      call WarningMessage(2,'Error in Post_2Center_RI')
      write(u6,*) 'kQm > mQm!'
      write(u6,*) 'MaxMem=',MaxMem
      write(u6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
      write(u6,*) 'nB,mB=',nB,mB
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
    write(u6,*) 'lQm < 1'
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
      iOff = nTri_Elem(kCol-1)
      A_l(1:kCol) => Am(iOff+1:iOff+kCol)
      if (kCol == 1) then
        nAm = nTri_Elem(nMem)
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
      write(u6,*) 'Inv_Cho_Factor found linear dependence!'
      call Abend()
    end if

    ! Write the new A/Q-vector to file

    iAddr_ = iAddr
    if (kCol == nMem) then
      nQm = nTri_Elem(kCol)
      call dDaFile(Lu_Q(iIrrep),1,Qm,nQm,iAddr)
      call dDaFile(Lu_A(iIrrep),1,Am,nQm,iAddr_)
    else if (kCol > nMem) then
      call dDaFile(Lu_Q(iIrrep),1,Q_l,kCol,iAddr)
      call dDaFile(Lu_A(iIrrep),1,A_l,kCol,iAddr_)
    end if

  end do

  nullify(Q_l,A_l)
  call mma_deallocate(Q_k)
  call mma_deallocate(A_k)
  call mma_deallocate(Qm)
  call mma_deallocate(Am)
  call mma_deallocate(X)
  call mma_deallocate(Z)
  call mma_deallocate(Scr)
  call DaClos(Lu_A(iIrrep))
end do ! iIrrep
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Sort the Q-matrix back to the original order.

call mma_maxDBLE(MaxMem2)

if (Force_Out_of_Core) MaxMem2 = 2*nBfnTot
lScr = min(MaxMem2,max(nBfn2,2*nBfnTot))
call mma_allocate(Scr,lScr,Label='Scr')

call SORT_mat(irc,nDmA,nDmB,iDiag,nIrrep,Lu_Q,'Restore',lScr,Scr,dum)

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
