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

subroutine CD_AInv_Inner(n,m,ADiag,Lu_A,Lu_Q,Thr_CD)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, Lu_A
integer(kind=iwp), intent(out) :: m
real(kind=wp), intent(inout) :: ADiag(n)
integer(kind=iwp), intent(inout) :: Lu_Q
real(kind=wp), intent(in) :: Thr_CD
integer(kind=iwp) :: iAddr, iAddr_, iOff, kCol, kQm, lAm, LinDep, lQm, lScr, Lu_Z, MaxMem, MaxMem2, mb, mQm, nAm, nB, nBfn2, &
                     nBfnTot, nDim, nMem, nQm, nQm_full, nScr, nVec, nXZ
real(kind=wp) :: a, b, Thr, ThrQ
logical(kind=iwp) :: Out_of_Core
integer(kind=iwp), allocatable :: iADiag(:)
real(kind=wp), allocatable :: Scr(:), Z(:), X(:)
real(kind=wp), allocatable, target :: Qm(:), Am(:), Q_k(:), A_k(:)
real(kind=wp), pointer :: Q_l(:), A_l(:)

nScr = 3*n
call mma_maxDBLE(MaxMem)
lScr = min(MaxMem,nScr)
call mma_allocate(Scr,lScr,Label='Scr')
call mma_allocate(iADiag,n,Label='iADiag')

nDim = n

Thr = Thr_CD*0.1_wp
Lu_Z = 7
call DaName_MF_WA(Lu_Z,'ZMAT09')
call Get_Pivot_idx(ADiag,nDim,nVec,Lu_A,Lu_Z,iADiag,Scr,lScr,Thr)
m = nVec
if (nDim /= nVec) then
  write(u6,*)
  write(u6,*) 'Detected lin. dep. in the auxiliary basis'
  write(u6,'(A,I6)') ' # of aux. bfns before lin. dep. removal: ',nDim
  write(u6,'(A,I6)') ' # of aux. bfns after  lin. dep. removal: ',nVec
end if

call Pivot_Mat(nDim,nVec,Lu_A,Lu_Z,iADiag,Scr,lScr)

call mma_deallocate(Scr)

!***********************************************************************
!     A-vectors are now on disk. Go ahead and compute the Q-vectors!
!***********************************************************************

ThrQ = Thr_CD*0.1_wp ! Threshold for Inv_Cho_Factor

nB = nVec
if (nB /= 0) then
  nQm = nTri_Elem(nB)

  nXZ = nB
  nQm_full = nTri_Elem(nB)

  Out_of_Core = 2*nQm_full+5*nXZ > MaxMem

  if (Out_Of_Core) then
    mQm = (nQm*MaxMem-5*nXZ)/(2*nQm_full)
    a = One
    b = -Two*real(mQm,kind=wp)
    mB = int(-a*Half+sqrt((a*Half)**2-b))
    kQm = nTri_Elem(mB)
    if (kQm > mQm) then
      call WarningMessage(2,'Error in CD_AInv_Inner')
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
    call WarningMessage(2,'Error in CD_AInv_Inner')
    write(u6,*) 'lQm < 1'
    call Abend()
  end if

  ! Some of memory for scratch arrays for Inv_Cho_Factor
  ! Allocate memory for the A- and Q-vectors and initialize.

  lScr = nXZ
  call mma_allocate(Scr,lScr,Label='Scr')
  call mma_allocate(Qm,lQm,Label='Qm')
  call mma_allocate(Am,lAm,Label='Am')
  call mma_allocate(A_k,nXZ,Label='A_k')
  call mma_allocate(Q_k,nXZ,Label='Q_k')
  call mma_allocate(X,nXZ,Label='X')
  call mma_allocate(Z,nXZ,Label='Z')

  Am(:) = Zero
  Qm(:) = Zero
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Process the A_ks to generate Q_ks.

  iAddr = 0
  nMem = mB
  do kCol=1,nB

    iAddr_ = iAddr
    if (kCol <= nMem) then
      iOff = nTri_Elem(kCol-1)
      ! Point to A_k in Am
      A_l(1:kCol) => Am(iOff+1:iOff+kCol)
      if (kCol == 1) then
        nAm = nTri_Elem(nMem)
        call dDaFile(Lu_Z,2,Am,nAm,iAddr_)
      end if
      ! Point to Q_k in Qm
      Q_l(1:kCol) => Qm(iOff+1:iOff+kCol)
    else if (kCol > nMem) then
      ! Use special scratch for A_k
      A_l(1:kCol) => A_k(1:kCol)
      call dDaFile(Lu_Z,2,A_l,kCol,iAddr_)
      ! Use special scratch for Q_k
      Q_l(1:kCol) => Q_k(1:kCol)
    end if

    LinDep = 2
    call Inv_Cho_Factor(A_l,kCol,Am,Qm,nMem,Lu_Z,Lu_Q,Scr,lScr,Z,X,ThrQ,Q_l,LinDep)

    if (LinDep /= 0) then
      call WarningMessage(2,'Error in CD_AInv_Inner')
      write(u6,*) 'Inv_Cho_Factor found linear dependence!'
      call Abend()
    end if

    ! Write the new A/Q-vector to file

    iAddr_ = iAddr
    if (kCol == nMem) then
      nQm = nTri_Elem(kCol)
      call dDaFile(Lu_Q,1,Qm,nQm,iAddr)
      call dDaFile(Lu_Z,1,Am,nQm,iAddr_)
    else if (kCol > nMem) then
      call dDaFile(Lu_Q,1,Q_l,kCol,iAddr)
      call dDaFile(Lu_Z,1,A_l,kCol,iAddr_)
    end if

  end do

  nullify(Q_l,A_l)
  call mma_deallocate(X)
  call mma_deallocate(Z)
  call mma_deallocate(Q_k)
  call mma_deallocate(A_k)
  call mma_deallocate(Am)
  call mma_deallocate(Qm)
  call mma_deallocate(Scr)
end if
call DaEras(Lu_Z)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Sort the Q-matrix back to the original order.

call mma_maxDBLE(MaxMem2)

nBfnTot = n
nBfn2 = n**2
lScr = min(MaxMem2,max(nBfn2,2*nBfnTot))
call mma_allocate(Scr,lScr,Label='Scr')

call Restore_Mat(nDim,nVec,Lu_Q,Lu_A,iADiag,Scr,lScr,.true.)
call DaEras(Lu_Q)
Lu_Q = Lu_A

! Note: after the 'Restore' call to Sort_mat, the Q-matrix is
!       no longer stored as upper-triangular but as squared
!       (zeros have been added because the corresponding argument
!        is set to .true.). The column index is still pivoted.

call mma_deallocate(Scr)
call mma_deallocate(iADiag)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
return

end subroutine CD_AInv_Inner
