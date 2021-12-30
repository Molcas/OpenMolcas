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

subroutine Ortho_Orb(Xmo,Smat,nBas,nOrb2Loc,nPass,Test)
! Purpose: orthonormalization of orbitals according to
!
!          V = X^T*S*X
!          X <- X*V^(-1/2)
!
!          where S is the AO overlap matrix.
!          The orthonormalization is carried out nPass times.
!          After this routine, X will satisfy X^T*S*X=1.

implicit real*8(a-h,o-z)
real*8 Xmo(*), Smat(*)
integer nBas, nOrb2Loc
logical Test
#include "WrkSpc.fh"
character*9 SecNam
parameter(SecNam='Ortho_Orb')
parameter(Tol=1.0d-10)
external ddot_

! Check for quick return.
! -----------------------

if (nPass < 1) return

! Allocations.
! ------------

l_V = nOrb2Loc**2
l_VSqrt = l_V
l_VISqrt = l_V
l_Scr = 2*(nBas**2)+nBas*(nBas+1)/2 ! needed in SqrtMt
call GetMem('V','Allo','Real',ip_V,l_V)
call GetMem('VSqrt','Allo','Real',ip_VSqrt,l_VSqrt)
call GetMem('VISqrt','Allo','Real',ip_VISqrt,l_VISqrt)
call GetMem('Scr','Allo','Real',ip_Scr,l_Scr)

! Orthonormalization passes.
! --------------------------

do iPass=1,nPass

  ! Compute V = X^T*S*X.
  ! --------------------

  call GetUmat_Localisation(Work(ip_V),Xmo(1),Smat(1),Xmo(1),Work(ip_Scr),l_Scr,nBas,nOrb2Loc)

  ! Compute V^(-1/2).
  ! -----------------

  iTask = 2 ! compute sqrt as well as inverse sqrt
  call SqrtMt(Work(ip_V),nOrb2Loc,iTask,Work(ip_VSqrt),Work(ip_VISqrt),Work(ip_Scr))

  ! Compute orthonormal X <- X*V^(-1/2).
  ! ------------------------------------

  nB = max(nBas,1)
  nO2L = max(nOrb2Loc,1)
  call dCopy_(nBas*nOrb2Loc,Xmo(1),1,Work(ip_Scr),1)
  call DGEMM_('N','N',nBas,nOrb2Loc,nOrb2Loc,1.0d0,Work(ip_Scr),nB,Work(ip_VISqrt),nO2L,0.0d0,Xmo(1),nB)

end do

! Test orthonormalization (i.e. V=1?).
! ------------------------------------

if (Test) then
  nErr = 0
  call GetUmat_Localisation(Work(ip_V),Xmo(1),Smat(1),Xmo(1),Work(ip_Scr),l_Scr,nBas,nOrb2Loc)
  kOff = ip_V-1
  do i=1,nOrb2Loc
    Work(kOff+nOrb2Loc*(i-1)+i) = Work(kOff+nOrb2Loc*(i-1)+i)-1.0d0
  end do
  xNrm = sqrt(dDot_(nOrb2Loc**2,Work(ip_V),1,Work(ip_V),1))
  if (xNrm > Tol) then
    write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm
    nErr = nErr+1
  end if
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

return

end subroutine Ortho_Orb
