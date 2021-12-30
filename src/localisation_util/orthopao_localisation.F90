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

implicit real*8(a-h,o-z)
real*8 X(*)
integer nBas(nSym), nFro(nSym), nOrb2Loc(nSym)
logical Test
#include "WrkSpc.fh"
character*21 SecNam
parameter(SecNam='OrthoPAO_Localisation')
parameter(Tol=1.0d-10)
external ddot_

! Check for quick return.
! -----------------------

if (nPass < 1) return

! Read S from disk, stored as full square.
! ----------------------------------------

l_S = nBas(1)**2
do iSym=2,nSym
  l_S = l_S+nBas(iSym)**2
end do
call GetMem('S','Allo','Real',ip_S,l_S)
call GetOvlp_Localisation(Work(ip_S),'Sqr',nBas,nSym)

! Allocations.
! ------------

nBasMax = nBas(1)
nO2LMax = nOrb2Loc(1)
do iSym=2,nSym
  nBasMax = max(nBasMax,nBas(iSym))
  nO2LMax = max(nO2LMax,nOrb2Loc(iSym))
end do

l_V = nO2LMax**2
l_VSqrt = l_V
l_VISqrt = l_V
l_Scr = 2*(nBasMax**2)+nBasMax*(nBasMax+1)/2  ! needed in SqrtMt
call GetMem('V','Allo','Real',ip_V,l_V)
call GetMem('VSqrt','Allo','Real',ip_VSqrt,l_VSqrt)
call GetMem('VISqrt','Allo','Real',ip_VISqrt,l_VISqrt)
call GetMem('Scr','Allo','Real',ip_Scr,l_Scr)

! Orthonormalization passes.
! --------------------------

do iPass=1,nPass
  kX = 1
  kS = ip_S
  do iSym=1,nSym

    ! Set pointer to PAO part of X.
    ! ------------------------------

    kOffX = kX+nBas(iSym)*nFro(iSym)

    ! Compute V = X^T*S*X.
    ! --------------------

    call GetUmat_Localisation(Work(ip_V),X(kOffX),Work(kS),X(kOffX),Work(ip_Scr),l_Scr,nBas(iSym),nOrb2Loc(iSym))

    ! Compute V^(-1/2).
    ! -----------------

    iTask = 2 ! compute sqrt as well as inverse sqrt
    call SqrtMt(Work(ip_V),nOrb2Loc(iSym),iTask,Work(ip_VSqrt),Work(ip_VISqrt),Work(ip_Scr))

    ! Compute orthonormal X <- X*V^(-1/2).
    ! ------------------------------------

    nB = max(nBas(iSym),1)
    nO2L = max(nOrb2Loc(iSym),1)
    call dCopy_(nBas(iSym)*nOrb2Loc(iSym),X(kOffX),1,Work(ip_Scr),1)
    call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),1.0d0,Work(ip_Scr),nB,Work(ip_VISqrt),nO2L,0.0d0,X(kOffX),nB)

    ! Update pointers.
    ! ----------------

    kX = kX+nBas(iSym)**2
    kS = kS+nBas(iSym)**2

  end do
end do

! Test orthonormalization (i.e. V=1?).
! ------------------------------------

if (Test) then
  kX = 1
  kS = ip_S
  nErr = 0
  do iSym=1,nSym
    kOffX = kX+nBas(iSym)*nFro(iSym)
    call GetUmat_Localisation(Work(ip_V),X(kOffX),Work(kS),X(kOffX),Work(ip_Scr),l_Scr,nBas(iSym),nOrb2Loc(iSym))
    kOff = ip_V-1
    do i=1,nOrb2Loc(iSym)
      Work(kOff+nOrb2Loc(iSym)*(i-1)+i) = Work(kOff+nOrb2Loc(iSym)*(i-1)+i)-1.0d0
    end do
    xNrm = sqrt(dDot_(nOrb2Loc(iSym)**2,Work(ip_V),1,Work(ip_V),1))
    if (xNrm > Tol) then
      write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm,' (sym.',iSym,')'
      nErr = nErr+1
    end if
    kX = kX+nBas(iSym)**2
    kS = kS+nBas(iSym)**2
  end do
  if (nErr /= 0) then
    write(6,*) SecNam,': failure after ',nPass,' passes'
    call SysAbendMsg(SecNam,'Orthonormalization failure!',' ')
  end if
end if

! De-allocations.
! ---------------

call GetMem('Scr','Free','Real',ip_Scr,l_Scr)
call GetMem('VISqrt','Free','Real',ip_VISqrt,l_VISqrt)
call GetMem('VSqrt','Free','Real',ip_VSqrt,l_VSqrt)
call GetMem('V','Free','Real',ip_V,l_V)
call GetMem('S','Free','Real',ip_S,l_S)

end subroutine OrthoPAO_Localisation
