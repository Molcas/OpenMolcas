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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine CalcGprime(Gprime,Mass,xvec,InterVec,AtCoord,NumOfAt,h,NumInt)
!  Purpose:
!    Calculate first derivatives of G.
!
!  Input:
!    Mass     : Real array - the mass of the atoms.
!    xvec     : Real array - the geometry in internal coordinates.
!    InterVec : Integer array.
!    NumOfAt  : Integer - the number of atoms.
!
!  Output:
!    Gprime   : Real two dimensional array - first derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(out) :: Gprime(ngdim,ngdim,ngdim)
real(kind=wp), intent(in) :: Mass(NumOfAt), xvec(NumInt), h
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: icoord, ih, iterm
real(kind=wp) :: xtmp(NumInt)
real(kind=wp), allocatable :: Gtemp(:,:,:), Stemp(:,:,:)

! Initialize.
call mma_allocate(Stemp,3,NumOfAt,NumInt,label='Stemp')
call mma_allocate(Gtemp,NumInt,NumInt,4,label='Gtemp')

do icoord=1,NumInt
  xtmp(:) = xvec
  Gtemp(:,:,:) = Zero
  iterm = 1
  do ih=-3,3,2
    xtmp(icoord) = xvec(icoord)+real(ih,kind=wp)*h
    !call Int_To_Cart(InterVec,xtmp,AtCoord,NumOfAt,NumInt,Mass)
    call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)
    Stemp(:,:,:) = Zero

    call CalcS(AtCoord,InterVec,Stemp,NumInt,NumOfAt)
    call CalcG(Gtemp(:,:,iterm),Mass,Stemp,NumInt,NumOfAt)

    iterm = iterm+1
  end do
  Gprime(1:NumInt,1:NumInt,icoord) = (Gtemp(:,:,1)-27.0_wp*Gtemp(:,:,2)+27.0_wp*Gtemp(:,:,3)-Gtemp(:,:,4))/(48.0_wp*h)
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)

call mma_deallocate(Stemp)
call mma_deallocate(Gtemp)

end subroutine CalcGprime
!####
subroutine CalcGdbleprime(Gdbleprime,Mass,xvec,InterVec,AtCoord,NumOfAt,h,NumInt)
!  Purpose:
!    Calculate second derivatives of G.
!
!  Input:
!    Mass       : Real array - the mass of the atoms.
!    xvec       : Real array - the geometry in internal coordinates.
!    InterVec   : Integer array.
!    NumOfAt    : Integer - the number of atoms.
!
!  Output:
!    Gdbleprime : Real two dimensional array - second derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: ngdim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Three
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(out) :: Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real(kind=wp), intent(in) :: Mass(NumOfAt), xvec(NumInt), h
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: icoord, jcoord
real(kind=wp), allocatable :: Gprime1(:,:,:), Gprime2(:,:,:), Gprime3(:,:,:), Gprime4(:,:,:), xtmp(:)

! Initialize.
call mma_allocate(xtmp,NumInt,label='xtmp')
call mma_allocate(Gprime1,NumInt,NumInt,NumInt,label='Gprime1')
call mma_allocate(Gprime2,NumInt,NumInt,NumInt,label='Gprime2')
call mma_allocate(Gprime3,NumInt,NumInt,NumInt,label='Gprime3')
call mma_allocate(Gprime4,NumInt,NumInt,NumInt,label='Gprime4')

do jcoord=1,NumInt
  xtmp(:) = xvec
  xtmp(jcoord) = xvec(jcoord)-Three*h
  call CalcGprime(Gprime1,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)-h
  call CalcGprime(Gprime2,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)+h
  call CalcGprime(Gprime3,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  xtmp(jcoord) = xvec(jcoord)+Three*h
  call CalcGprime(Gprime4,Mass,xtmp,InterVec,AtCoord,NumOfAt,h,NumInt)
  do icoord=1,NumInt
    Gdbleprime(1:NumInt,1:NumInt,icoord,jcoord) = &
      (Gprime1(:,:,icoord)-27.0_wp*Gprime2(:,:,icoord)+27.0_wp*Gprime3(:,:,icoord)-Gprime4(:,:,icoord))/(48.0_wp*h)
  end do
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_To_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)

call mma_deallocate(xtmp)
call mma_deallocate(Gprime1)
call mma_deallocate(Gprime2)
call mma_deallocate(Gprime3)
call mma_deallocate(Gprime4)

end subroutine CalcGdbleprime
