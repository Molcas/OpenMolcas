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
!    Mass     : Real*8 array - the mass of the atoms.
!    xvec     : Real*8 array - the geometry in internal
!               coordinates.
!    InterVec : Integer array.
!    NumOfAt  : Integer - the number of atoms.
!
!  Output:
!    Gprime   : Real*8 two dimensional array - first
!               derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp

!implicit none
#include "Constants_mula.fh"
#include "dims.fh"
real*8 h
real*8 Gprime(ngdim,ngdim,ngdim)
integer NumInt, NumOfAt
real*8 Mass(NumOfAt)
real*8 xvec(NumInt)
real*8 xtmp(NumInt)
integer InterVec(*)
real*8 AtCoord(3,NumOfAt)
integer icoord, iterm
integer ih
real*8, allocatable :: Gtemp(:,:,:), Stemp(:,:,:)

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
!    Mass       : Real*8 array - the mass of the atoms.
!    xvec       : Real*8 array - the geometry in internal
!                 coordinates.
!    InterVec   : Integer array.
!    NumOfAt    : Integer - the number of atoms.
!
!  Output:
!    Gdbleprime : Real*8 two dimensional array - second
!                 derivative of the inverse mass tensor G.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Three
use Definitions, only: wp

!implicit none
#include "Constants_mula.fh"
#include "dims.fh"
real*8 h
integer icoord, jcoord
integer NumInt, NumOfAt
real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real*8 Mass(NumOfAt)
real*8 xvec(NumInt)
integer InterVec(*)
real*8 AtCoord(3,NumOfAt)
real*8, allocatable :: Gprime1(:,:,:), Gprime2(:,:,:), Gprime3(:,:,:), Gprime4(:,:,:), xtmp(:)

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
