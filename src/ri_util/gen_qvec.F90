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

subroutine Gen_QVec(nIrrep,nBas_Aux)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
integer nBas_Aux(0:nIrrep-1), Lu_Q(0:7), Lu_A(0:7)
logical Out_Of_Core
character*6 Name_Q
real*8, allocatable, target :: Mem(:)
real*8, pointer :: A_l(:) => null(), Q_l(:) => null()
real*8, pointer :: A_k(:) => null(), Q_k(:) => null()
real*8, pointer :: Am(:) => null(), Qm(:) => null()
integer, allocatable :: iDiag(:)
real*8, allocatable :: Scr(:), X(:), Z(:)
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
ThrQ = 1.0D-14 ! Threshold for Inv_Cho_Factor

mB = 0
nA_Diag = 0
do iIrrep=0,nIrrep-1
  nB = nBas_Aux(iIrrep)
  nA_Diag = nA_Diag+nB
  mB = max(mB,nB)

  iSeed = 55+iIrrep
  Lu_Q(iIrrep) = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'QMAT',iIrrep
  call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)

  iSeed = 63+iIrrep
  Lu_A(iIrrep) = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'AVEC',iIrrep
  call DaName_MF_WA(Lu_A(iIrrep),Name_Q)
end do
nBfn2 = mB**2

call mma_allocate(Z,mB,Label='Z')
call mma_allocate(X,mB,Label='X')
lScr = 3*mB
call mma_allocate(Scr,lScr,Label='Scr')
call mma_maxDBLE(Mem_Max)
call mma_allocate(Mem,Mem_Max,Label='Mem')

do iIrrep=0,nIrrep-1
  nB = nBas_Aux(iIrrep)
  nQm = nB*(nB+1)/2

  Out_Of_Core = 2*nQm > Mem_Max
  if (Out_Of_Core) then
    MaxMem = Mem_Max-2*nB
    mQm = MaxMem/2
    a = One
    b = -Two*dble(mQm)
    mB = int(-a/Two+sqrt((a/Two)**2-b))
    kQm = mB*(mB+1)/2
    if (kQm > mQm) then
      call WarningMessage(2,'Error in Gen_QVec')
      write(6,*) 'kQm > mQm!'
      write(6,*) 'MaxMem=',MaxMem
      write(6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
      write(6,*) 'nB,mB=',nB,mB
      call Abend()
    end if
    iE = 2*kQm
    iS = iE+1
    iE = iE+mB
    Q_k(1:mB) => Mem(iS:iE)
    iS = iE+1
    iE = iE+mB
    A_k(1:mB) => Mem(iS:iE)
  else
    mB = nB
    kQm = nQm
  end if

  iS = 1
  iE = kQm
  Qm(1:kQm) => Mem(iS:iE)
  iS = iE+1
  iE = iE+kQm
  Am(1:kQm) => Mem(iS:iE)

  iAddr = 0
  do kCol=1,nB

    if (kCol <= mB) then
      iOff = (kCol-1)*kCol/2
      A_l(1:) => Am(1+iOff:)
    else
      A_l(1:) => A_k(1:)
    end if

    iAddr_ = iAddr
    if ((kCol <= mB) .and. (kCol == 1)) then
      call dDaFile(Lu_A(iIrrep),2,Am,kQm,iAddr_)
    else if (kCol > mB) then
      call dDaFile(Lu_A(iIrrep),2,A_l,kCol,iAddr_)
    end if
#   ifdef _DEBUGPRINT_
    write(6,*) 'kCol=',kCol
    call TriPrt('Am',' ',Am,mB)
    call RecPrt('Al',' ',A_l,1,kCol)
#   endif

    if (kCol <= mB) then
      iOff = (kCol-1)*kCol/2
      Q_l(1:) => Qm(1+iOff:)
    else
      Q_l(1:) => Q_k(1:)
    end if

    LinDep = 2
    call Inv_Cho_Factor(A_l,kCol,Am,Qm,mB,Lu_A(iIrrep),Lu_Q(iIrrep),Scr,lScr,Z,X,ThrQ,Q_l,LinDep)

    if (LinDep /= 0) then
      call WarningMessage(2,'Error in Gen_QVec')
      write(6,*) 'Inv_Cho_Factor found linear dependence!'
      call Abend()
    end if
#   ifdef _DEBUGPRINT_
    call TriPrt('Qm',' ',Qm,min(mB,kCol))
    call RecPrt('Ql',' ',Q_l,1,kCol)
#   endif

    ! Write the new A/Q-vector to file

    iAddr_ = iAddr
    if (kCol == mB) then
      lQm = kCol*(kCol+1)/2
      call dDaFile(Lu_Q(iIrrep),1,Qm,lQm,iAddr)
      call dDaFile(Lu_A(iIrrep),1,Am,lQm,iAddr_)
    else if (kCol > mB) then
      nQ_k = kCol
      call dDaFile(Lu_Q(iIrrep),1,Q_l,nQ_k,iAddr)
      call dDaFile(Lu_A(iIrrep),1,A_l,nQ_k,iAddr_)
    end if

  end do    ! kCol
  call DaClos(Lu_A(iIrrep))
end do      ! iIrrep

A_l => null()
Q_l => null()
A_k => null()
Q_k => null()
Am => null()
Qm => null()
call mma_deallocate(Mem)
call mma_deallocate(Scr)
call mma_deallocate(X)
call mma_deallocate(Z)

! Sort the Q-matrix to square storage.

call mma_allocate(iDiag,nA_Diag,Label='iDiag')
ik = 0
do iIrrep=0,nIrrep-1
  do k=1,nBas_Aux(iIrrep)
    ik = ik+1
    iDiag(ik) = k  ! dummy assignement
  end do
end do
call mma_maxDBLE(MaxMem2)
lScr = min(MaxMem2,nBfn2)
call mma_allocate(Scr,lScr,Label='Scr')

call SORT_Mat(irc,nBas_Aux,nBas_Aux,iDiag,nIrrep,Lu_Q,'Restore',lScr,Scr)

! Note: after the 'Restore' call to Sort_mat, the Q-matrix is
!       no longer stored as upper-triangular but as squared
!       (zeros have been added).

call mma_deallocate(Scr)
call mma_deallocate(iDiag)

do iIrrep=0,nIrrep-1
  call DaClos(Lu_Q(iIrrep))
end do

return

end subroutine Gen_QVec
