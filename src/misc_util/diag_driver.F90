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

subroutine Diag_Driver(JobZ,Range,UpLo,nDim,Triangular,Aux,lDimAux,vLower,vUpper,iLower,iUpper,EigVal,EigVec,lDimVec,iUnit_Matrix, &
                       iSort,Method,nFound,Ierr)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character :: JobZ, Range, UpLo, Method
integer(kind=iwp) :: nDim, lDimAux, iLower, iUpper, lDimVec, iUnit_Matrix, iSort, nFound, iErr
real(kind=wp) :: Triangular(*), Aux(*), vLower, vUpper, EigVal(*), EigVec(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: iMethod, iSize, iSuppZ, liErr(2), liScr, liW(2), liWork, lScr, lWork, nDim2
real(kind=wp) :: Tolerance, Work_L(1)
real(kind=wp), external :: dLamCh_
logical(kind=iwp), external :: lSame

! Sizes for arrays

nDim2 = nDim*nDim
iSize = nDim*(nDim+1)/2

! Determine which algorithm to use

iMethod = 0
if (lSame(Method,'A')) then
  iMethod = 0
else if (lSame(Method,'Q')) then
  iMethod = 1
else if (lSame(Method,'J')) then
  iMethod = 2
else
  write(u6,*) '!!! Diag_Driver called with an unknown method: ',Method
  write(u6,*) '!!! Supported methods: Q, J, and A'
  write(u6,*) '    Method = ''',Method,''''
  call Abend()
end if

if (iMethod <= 1) then

  ! Use the QL algorithm (dSyevR)

  call Square(Triangular,Aux,lDimAux,1,nDim)
  call dCopy_(nDim2,[Zero],0,EigVec,1)
  call dCopy_(nDim,[One],0,EigVec,nDim+1)

  ! Determine safe tolerance for dSyevR

  Tolerance = dLamCh_('Safe minimum')

  ! Determine optimal sizes of scratch arrays

  call GetMem('ISUPPZ  ','ALLO','INTE',iSuppZ,2*nDim)
  lWork = -1
  liWork = -1
  !C AOM 03.08.2005 - Added LiW(2), otherwise on Opteron crashed
  !C AOM 04.08        Also added liErr for the same reason
  call dsyevr_(JobZ,Range,UpLo,nDim,Aux,lDimAux,vLower,vUpper,iLower,iUpper,Tolerance,nFound,EigVal,EigVec,lDimVec,iWork(iSuppZ), &
               Work_L,lWork,liW,liWork,liErr(1))
  lWork = int(Work_L(1))
  liWork = liW(1)

  ! Allocate scratch arrays

  call GetMem('SCRATCH ','ALLO','REAL',lScr,lWork)
  call GetMem('ISCRATCH','ALLO','INTE',liScr,liWork)

  ! Run actual QL routine

  call dsyevr_(JobZ,Range,UpLo,nDim,Aux,lDimAux,vLower,vUpper,iLower,iUpper,Tolerance,nFound,EigVal,EigVec,lDimVec,iWork(iSuppZ), &
               Work(lScr),lWork,iWork(liScr),liWork,liErr(1))
  iErr = liErr(1)

  ! Free scratch

  call GetMem('SCRATCH ','FREE','REAL',lScr,lWork)
  call GetMem('ISCRATCH','FREE','INTE',liScr,liWork)
  call GetMem('ISUPPZ  ','FREE','INTE',iSuppZ,2*nDim)

  ! Check for convergence

  if (iErr /= 0) then
    write(u6,*) '!!! No Convergence in the QL algorithm.'
    if (lSame(Method,'A')) then
      write(u6,*) '!!! Trying Jacobi instead.'
      write(u6,*) '!!! Warning: This might be very slow.'
      iMethod = 2
    else
      call Abend()
    end if
  else
    call Chk4NAN(nDim**2,EigVec,Ierr)
    if (iErr > 0) then
      write(u6,*) 'At least one of the eigenvectors found with'
      write(u6,*) 'DSYEVR contained a NAN.'
      if (lSame(Method,'A')) then
        write(u6,*) 'Trying Jacobi instead.'
        write(u6,*) 'Warning: This might be very slow.'
        iMethod = 2
      else
        call Abend()
      end if
    end if
  end if
end if
if (iMethod == 2) then

  ! Use the Jacobi algorithm

  call dCopy_(iSize,Triangular,1,Aux,1)
  if (iUnit_Matrix == 1) then
    call dCopy_(nDim2,[Zero],0,EigVec,1)
    call dCopy_(nDim,[One],0,EigVec,nDim+1)
  end if
  call Jacob(Aux,EigVec,nDim,lDimVec)
  call vEig(nDim,Aux,EigVal)
end if

! Sort the eigenvalues and eigenvectors?

if (iSort == 1) then
  call JacOrd2(EigVal,EigVec,nDim,lDimVec)
else if (iSort == -1) then
  call SortEig(EigVal,EigVec,nDim,lDimVec)
end if

return

end subroutine Diag_Driver
