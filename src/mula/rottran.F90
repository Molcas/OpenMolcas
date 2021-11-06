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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine RotTranRem(Sinv,S,Mass,AtCoord,NumOfAt,NumInt)
!  Purpose:
!    Project total translation and total rotation out of inverted
!    S matrix.
!
!  Input:
!    S        : Real*8 three dimensional array.
!    Mass     : Real*8 array - masses of the atoms.
!    AtCoord  : Real*8 two dimensional array - coordinates
!               of the atoms.
!
!  Output:
!    Sinv     : Real*8 three dimensional array - Inverted
!               S matrix with rotation and translation projected out.
!
!  Calls:
!    Daxpy  (ESSL)
!    Dcopy  (ESSL)
!    Ddot_  (ESSL)
!
!  Uses:
!    LinAlg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!use Linalg
!implicit none
#include "Constants_mula.fh"
integer NumOfAt, NumInt, nFree
real*8 AtCoord(3,NumofAt)
real*8 S(3,NumOfAt,NumInt)
real*8 Sinv(3,NumOfAt,NumInt)
real*8 Mass(NumOfAt)
real*8 Det, X
integer iAtom, n, k, i, j
#include "WrkSpc.fh"

! Initialize.
nFree = 6
n = 3*NumOfAt
call GetMem('Amat','Allo','Real',ipAmat,3*NumOfAt*nFree)
lAmat = n
!Amat = 0.0d0
call dcopy_(3*NumOfAt*nFree,[0.0d0],0,Work(ipAmat),1)

! Invert S matrix.
call GetMem('Temp2','Allo','Real',ipTemp2,NumInt*Numint)

call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,1.0d0,S,3*NumOfAt,S,3*NumOfAt,0.0d0,Work(ipTemp2),NumInt)
! Invert, by solving eq Temp2*X=Temp1. Solution computed in-place.
! Thus, Temp1=Unit matrix(in) and contains solution (out).
call GetMem('Temp1','Allo','Real',ipTemp1,NumInt*Numint)

call dcopy_(NumInt**2,[0.0d0],0,Work(ipTemp1),1)
call dcopy_(NumInt,[1.0d0],0,Work(ipTemp1),NumInt+1)
call Dool_MULA(Work(ipTemp2),NumInt,NumInt,Work(ipTemp1),NumInt,NumInt,det)
!PAM01 Replacement for Dool_MULA, if superstable solution is wanted:
!Eps = 1.0D-8
!call SymSolve(Temp2,Temp1,Temp3,Temp4,Eps)
!call dcopy_(NumInt**2,Temp3,1,Temp1,1)
!PAM01 End of replacement code.

call DGEMM_('N','N',3*NumOfAt,NumInt,NumInt,1.0d0,S,3*NumOfAt,Work(ipTemp1),NumInt,0.0d0,Sinv,3*NumOfAt)
call GetMem('Temp1','Free','Real',ipTemp1,NumInt*Numint)
call GetMem('Temp2','Free','Real',ipTemp2,NumInt*Numint)

! Pure translation.
k = 1
do iAtom=1,NumOfAt
  Work(ipAmat+k-1) = 1.0d0
  Work(ipAmat+k+lAmat) = 1.0d0
  Work(ipAmat+k+1+lAmat*2) = 1.0d0
  k = k+3
end do

! Pure rotation.
k = 1
do iAtom=1,NumOfAt
  Work(ipAmat+k+lAmat*3-1) = -AtCoord(2,iAtom)
  Work(ipAmat+k+1+lAmat*3-1) = AtCoord(1,iAtom)
  Work(ipAmat+k+2+lAmat*3-1) = 0.0d0
  Work(ipAmat+k+lAmat*4-1) = AtCoord(3,iAtom)
  Work(ipAmat+k+1+lAmat*4-1) = 0.0d0
  Work(ipAmat+k+2+lAmat*4-1) = -AtCoord(1,iAtom)
  Work(ipAmat+k+lAmat*5-1) = 0.0d0
  Work(ipAmat+k+1+lAmat*5-1) = -AtCoord(3,iAtom)
  Work(ipAmat+k+2+lAmat*5-1) = AtCoord(2,iAtom)
  k = k+3
end do

! Scale with the mass of the atom.
n = 3*NumOfAt
call GetMem('AmatMass','Allo','Real',ipAmatMass,3*NumOfAt*nFree)
!AmatMass = Amat
call dcopy_(lAmat*nFree,Work(ipAmat),1,Work(ipAmatMass),1)
!PAM04: Replace the following code section...
!do i=1,nFree
!  Acol => AmatMass(:,i)
!  k = 1
!  do j=1,NumOfAt
!    jMass = (k+2)/3
!    Acol(k) = uToAu*Mass(jMass)*Acol(k)
!    Acol(k+1) = uToAu*Mass(jMass)*Acol(k+1)
!    Acol(k+2) = uToAu*Mass(jMass)*Acol(k+2)
!    k = k+3
!  end do
!end do
!PAM04: ... with the following...
do i=1,nFree
  do j=1,NumOfAt
    X = uToAu*Mass(j)
    Work(ipAmatMass+1+3*(j-1)+lAmat*(i-1)-1) = X*Work(ipAmatMass+1+3*(j-1)+lAmat*(i-1)-1)
    Work(ipAmatMass+2+3*(j-1)+lAmat*(i-1)-1) = X*Work(ipAmatMass+2+3*(j-1)+lAmat*(i-1)-1)
    Work(ipAmatMass+3+3*(j-1)+lAmat*(i-1)-1) = X*Work(ipAmatMass+3+3*(j-1)+lAmat*(i-1)-1)
  end do
end do
!PAM04: ... until here.

! Project rotation and translation out of S matrix.
call GetMem('Temp1','Allo','Real',ipTemp1,nFree*nFree)
call DGEMM_('T','N',nFree,nFree,3*NumOfAt,1.0d0,Work(ipAmatMass),3*NumOfAt,Work(ipAmat),3*NumOfAt,0.0d0,Work(ipTemp1),nFree)
call GetMem('Ainv','Allo','Real',ipAinv,nFree*nFree)
call dcopy_(nFree**2,[0.0d0],0,Work(ipAinv),1)
call dcopy_(nFree,[1.0d0],0,Work(ipAinv),nFree+1)
call Dool_MULA(Work(ipTemp1),nFree,nFree,Work(ipAinv),nFree,nFree,det)
call GetMem('Temp1','Free','Real',ipTemp1,nFree*nFree)

call GetMem('Temp2','Allo','Real',ipTemp2,nFree*Numint)
call DGEMM_('T','N',nFree,NumInt,3*NumOfAt,1.0d0,Work(ipAmatMass),3*NumOfAt,Sinv,3*NumOfAt,0.0d0,Work(ipTemp2),nFree)
call GetMem('AmatMass','Free','Real',ipAmatMass,3*NumOfAt*nFree)

call GetMem('Temp3','Allo','Real',ipTemp3,nFree*Numint)

call DGEMM_('N','N',nFree,NumInt,nFree,1.0d0,Work(ipAinv),nFree,Work(ipTemp2),nFree,0.0d0,Work(ipTemp3),nFree)
call GetMem('Ainv','Free','Real',ipAinv,nFree*nFree)
call GetMem('Temp2','Free','Real',ipTemp2,nFree*Numint)

call GetMem('Stemp','Allo','Real',ipStemp,3*NumOfAt*Numint)
call DGEMM_('N','N',3*NumOfAt,NumInt,nFree,1.0d0,Work(ipAmat),3*NumOfAt,Work(ipTemp3),nFree,0.0d0,work(ipStemp),3*NumOfAt)
call GetMem('Amat','Free','Real',ipAmat,3*NumOfAt*nFree)
call GetMem('Temp3','Free','Real',ipTemp3,nFree*Numint)

!Sinv = Sinv-Stemp
call daxpy_(3*NumOfAt*Numint,-1.0d0,Work(ipStemp),1,Sinv,1)

call GetMem('Stemp','Free','Real',ipStemp,3*NumOfAt*Numint)

end subroutine RotTranRem
