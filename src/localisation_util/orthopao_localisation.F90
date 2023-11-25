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
! Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine OrthoPAO_Localisation(X,nBas,nFro,nOrb2Loc,nSym,nPass,Test)
! Thomas Bondo Pedersen, December 2005.
! - revised January 2006.
!
! Purpose: orthonormalization of Cholesky PAOs according to
!
!          V = X^T*S*X
!          X <- X*V^(-1/2)
!
!          where S is the AO overlap matrix.
!          The orthonormalization is carried out nPass times.
!          After this routine, X will satisfy X^T*S*X=1.
!
! NOTE: X is assumed to contain all orbitals!!

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: X(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nOrb2Loc(nSym), nPass
logical(kind=iwp), intent(in) :: Test
integer(kind=iwp) :: i, iPass, iSym, iTask, kOffX, l_Scr, nB, nBasMax, nErr, nO2L, nO2LMax
real(kind=wp) :: xNrm
type(DSBA_Type) :: S, X2
real(kind=wp), allocatable :: Scr(:)
real(kind=wp), allocatable, target :: Vt(:,:)
real(kind=wp), pointer :: V(:,:), VISqrt(:,:), VSqrt(:,:)
real(kind=wp), parameter :: Tol = 1.0e-10_wp
character(len=*), parameter :: SecNam = 'OrthoPAO_Localisation'
real(kind=wp), external :: ddot_

! Check for quick return.
! -----------------------

if (nPass < 1) return

! Read S from disk, stored as full square.
! ----------------------------------------

call Allocate_DT(X2,nBas,nBas,nSym,label='X',Ref=X)
call Allocate_DT(S,nBas,nBas,nSym,label='S')
call GetOvlp_Localisation(S%A0,'Sqr',nBas,nSym)

! Allocations.
! ------------

nBasMax = nBas(1)
nO2LMax = nOrb2Loc(1)
do iSym=2,nSym
  nBasMax = max(nBasMax,nBas(iSym))
  nO2LMax = max(nO2LMax,nOrb2Loc(iSym))
end do

l_Scr = 2*(nBasMax**2)+nBasMax*(nBasMax+1)/2  ! needed in SqrtMt
call mma_allocate(Vt,nO2LMax**2,3,label='V')
call mma_allocate(Scr,l_Scr,label='Scr')

! Orthonormalization passes.
! --------------------------

do iPass=1,nPass
  do iSym=1,nSym

    V(1:nOrb2Loc(iSym),1:nOrb2Loc(iSym)) => Vt(1:nOrb2Loc(iSym)**2,1)
    VSqrt(1:nOrb2Loc(iSym),1:nOrb2Loc(iSym)) => Vt(1:nOrb2Loc(iSym)**2,2)
    VISqrt(1:nOrb2Loc(iSym),1:nOrb2Loc(iSym)) => Vt(1:nOrb2Loc(iSym)**2,3)

    ! Set pointer to PAO part of X.
    ! ------------------------------

    kOffX = nFro(iSym)+1

    ! Compute V = X^T*S*X.
    ! --------------------

    call GetUmat_Localisation(V,X2%SB(iSym)%A2(:,kOffX:),S%SB(iSym)%A2,X2%SB(iSym)%A2(:,kOffX:),Scr,nBas(iSym),nOrb2Loc(iSym))

    ! Compute V^(-1/2).
    ! -----------------

    iTask = 2 ! compute sqrt as well as inverse sqrt
    call SqrtMt(V,nOrb2Loc(iSym),iTask,VSqrt,VISqrt,Scr)

    ! Compute orthonormal X <- X*V^(-1/2).
    ! ------------------------------------

    nB = max(nBas(iSym),1)
    nO2L = max(nOrb2Loc(iSym),1)
    call dCopy_(nBas(iSym)*nOrb2Loc(iSym),X2%SB(iSym)%A2(:,kOffX:),1,Scr,1)
    call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),One,Scr,nB,VISqrt,nO2L,Zero,X2%SB(iSym)%A2(:,kOffX:),nB)

  end do
end do

nullify(V)
nullify(VSqrt)
nullify(VISqrt)

! Test orthonormalization (i.e. V=1?).
! ------------------------------------

if (Test) then
  nErr = 0
  do iSym=1,nSym
    V(1:nOrb2Loc(iSym),1:nOrb2Loc(iSym)) => Vt(1:nOrb2Loc(iSym)**2,1)
    kOffX = nFro(iSym)+1
    call GetUmat_Localisation(V,X2%SB(iSym)%A2(:,kOffX:),S%SB(iSym)%A2,X2%SB(iSym)%A2(:,kOffX:),Scr,nBas(iSym),nOrb2Loc(iSym))
    do i=1,nOrb2Loc(iSym)
      V(i,i) = V(i,i)-One
    end do
    xNrm = sqrt(dDot_(nOrb2Loc(iSym)**2,V,1,V,1))
    if (xNrm > Tol) then
      write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm,' (sym.',iSym,')'
      nErr = nErr+1
    end if
  end do
  if (nErr /= 0) then
    write(u6,*) SecNam,': failure after ',nPass,' passes'
    call SysAbendMsg(SecNam,'Orthonormalization failure!',' ')
  end if
  nullify(V)
end if

! De-allocations.
! ---------------

call mma_deallocate(Vt)
call mma_deallocate(Scr)
call Deallocate_DT(S)
call Deallocate_DT(X2)

end subroutine OrthoPAO_Localisation
