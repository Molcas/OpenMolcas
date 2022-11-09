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

subroutine Diag_Driver(JobZ,Rng,UpLo,nDim,Triangular,Aux,lDimAux,vLower,vUpper,iLower,iUpper,EigVal,EigVec,lDimVec,iUnit_Matrix, &
                       iSort,Method,nFound,Ierr)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, BLASR8

#include "intent.fh"

implicit none
character, intent(in) :: JobZ, Rng, UpLo, Method
integer(kind=iwp), intent(in) :: nDim, lDimAux, iLower, iUpper, lDimVec, iUnit_Matrix, iSort
real(kind=wp), intent(in) :: Triangular(*), vLower, vUpper
real(kind=wp), intent(_OUT_) :: Aux(*), EigVal(*), EigVec(*)
integer(kind=iwp), intent(out) :: nFound, iErr
integer(kind=iwp) :: iMethod, iSize, liErr(2), liW(2), liWork, lWork
real(kind=wp) :: Tolerance, Work_L(1)
integer(kind=iwp), allocatable :: iScr(:), iSuppZ(:)
real(kind=wp), allocatable :: Scr(:)
real(kind=BLASR8), external :: dLamCh
logical(kind=iwp), external :: lSame

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
  call unitmat(EigVec,nDim)

  ! Determine safe tolerance for dSyevR

  Tolerance = dLamCh('Safe minimum')

  ! Determine optimal sizes of scratch arrays

  call mma_allocate(iSuppZ,2*nDim,label='ISUPPZ')
  lWork = -1
  liWork = -1
  !C AOM 03.08.2005 - Added LiW(2), otherwise on Opteron crashed
  !C AOM 04.08        Also added liErr for the same reason
  call dsyevr_(JobZ,Rng,UpLo,nDim,Aux,lDimAux,vLower,vUpper,iLower,iUpper,Tolerance,nFound,EigVal,EigVec,lDimVec,iSuppZ,Work_L, &
               lWork,liW,liWork,liErr(1))
  lWork = int(Work_L(1))
  liWork = liW(1)

  ! Allocate scratch arrays

  call mma_allocate(Scr,lWork,label='SCRATCH')
  call mma_allocate(iScr,liWork,label='ISCRATCH')

  ! Run actual QL routine

  call dsyevr_(JobZ,Rng,UpLo,nDim,Aux,lDimAux,vLower,vUpper,iLower,iUpper,Tolerance,nFound,EigVal,EigVec,lDimVec,iSuppZ,Scr,lWork, &
               iScr,liWork,liErr(1))
  iErr = liErr(1)

  ! Free scratch

  call mma_deallocate(iSuppZ)
  call mma_deallocate(Scr)
  call mma_deallocate(iScr)

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
else if (iMethod == 2) then

  ! Use the Jacobi algorithm

  iSize = nTri_Elem(nDim)
  Aux(1:iSize) = Triangular(1:iSize)
  if (iUnit_Matrix == 1) call unitmat(EigVec,nDim)
  call Jacob(Aux,EigVec,nDim,lDimVec)
  call vEig(nDim,Aux,EigVal)
end if

! Sort the eigenvalues and eigenvectors?

if (iSort /= 0) call SortEig(EigVal,EigVec,nDim,lDimVec,1,iSort < 0)

return

end subroutine Diag_Driver
