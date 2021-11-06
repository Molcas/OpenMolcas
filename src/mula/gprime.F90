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
integer ih, k, j
#include "WrkSpc.fh"

! Initialize.
call GetMem('Stemp','Allo','Real',ipStemp,3*NumOfAt*NumInt)
call GetMem('Gtemp','Allo','Real',ipGtemp,NumInt*NumInt*4)

do icoord=1,NumInt
  do iv=1,NumInt
    xtmp(iv) = xvec(iv)
  end do
  !Gtemp = 0.0d0
  call dcopy_(NumInt*NumInt*4,[0.0d0],0,Work(ipGtemp),1)
  iterm = 1
  do ih=-3,3,2
    xtmp(icoord) = xvec(icoord)+dble(ih)*h
    !call Int_To_Cart(InterVec,xtmp,AtCoord,NumOfAt,NumInt,Mass)
    call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)
    !Stemp = 0.0d0
    call dcopy_(3*NumOfAt*NumInt,[0.0d0],0,Work(ipStemp),1)

    call CalcS(AtCoord,InterVec,Work(ipStemp),NumInt,NumOfAt)
    !Gtemp(1,1,iterm)= 1+NumInt*(0+NumInt*(iterm-1))-1
    !call CalcG(Gtemp(1,1,iterm),Mass,Work(ipStemp))
    call CalcG(Work(ipGtemp+NumInt*NumInt*(iterm-1)),Mass,Work(ipStemp),NumInt,NumOfAt)

    iterm = iterm+1
  end do
  do k=1,NumInt
    do j=1,NumInt
      !Gtemp(j,k,1) = j+NumInt*(k-1+NumInt*(1-1))-1
      !Gtemp(j,k,2) = j+NumInt*(k-1+NumInt*(2-1))-1
      !Gtemp(j,k,3) = j+NumInt*(k-1+NumInt*(3-1))-1

      !Gprime(j,k,icoord) = (Gtemp(j,k,1)-27.0d0*Gtemp(j,k,2)+27.0d0*Gtemp(j,k,3)-Gtemp(j,k,4))/(48.0d0*h)
      Gprime(j,k,icoord) = (Work(ipGtemp+j+NumInt*(k-1)-1)-27.0d0*Work(ipGtemp+j+NumInt*(k-1+NumInt)-1)+ &
                           27.0d0*Work(ipGtemp+j+NumInt*(k-1+NumInt*2)-1)-Work(ipGtemp+j+NumInt*(k-1+NumInt*3)-1))/(48.0d0*h)
    end do
  end do
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)

call GetMem('Stemp','Free','Real',ipStemp,3*NumOfAt*NumInt)
call GetMem('Gtemp','Free','Real',ipGtemp,NumInt*NumInt*4)

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

!implicit none
#include "Constants_mula.fh"
#include "dims.fh"
real*8 h
integer icoord, jcoord, k, j
integer NumInt, NumOfAt
real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
real*8 Mass(NumOfAt)
real*8 xvec(NumInt)
integer InterVec(*)
real*8 AtCoord(3,NumOfAt)
#include "WrkSpc.fh"

! Initialize.
call GetMem('xtmp','Allo','Real',ipxtmp,NumInt)
NumInt3 = NumInt*NumInt*NumInt
call GetMem('Gprime1','Allo','Real',ipGprime1,NumInt3)
call GetMem('Gprime2','Allo','Real',ipGprime2,NumInt3)
call GetMem('Gprime3','Allo','Real',ipGprime3,NumInt3)
call GetMem('Gprime4','Allo','Real',ipGprime4,NumInt3)

do jcoord=1,NumInt
  call dcopy_(NumInt,xvec,1,Work(ipxtmp),1)
  !xtmp = xvec
  Work(ipxtmp+jcoord-1) = xvec(jcoord)-3.0d0*h
  call CalcGprime(Work(ipGprime1),Mass,Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
  Work(ipxtmp+jcoord-1) = xvec(jcoord)-1.0d0*h
  call CalcGprime(Work(ipGprime2),Mass,Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
  Work(ipxtmp+jcoord-1) = xvec(jcoord)+1.0d0*h
  call CalcGprime(Work(ipGprime3),Mass,Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
  Work(ipxtmp+jcoord-1) = xvec(jcoord)+3.0d0*h
  call CalcGprime(Work(ipGprime4),Mass,Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
  do icoord=1,NumInt
    do k=1,NumInt
      do j=1,NumInt
        !Gprime(j,k,icoord) = j+NumInt*(k-1+NumInt*(icoord-1))-1
        ivv = j+NumInt*(k-1+NumInt*(icoord-1))-1
        Gdbleprime(j,k,icoord,jcoord) = (Work(ipGprime1+ivv)-27.0d0*Work(ipGprime2+ivv)+27.0d0*Work(ipGprime3+ivv)- &
                                        Work(ipGprime4+ivv))/(48.0d0*h)
      end do
    end do
  end do
end do
!call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
call Int_To_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)

call GetMem('xtmp','Free','Real',ipxtmp,NumInt)
call GetMem('Gprime1','Free','Real',ipGprime1,NumInt3)
call GetMem('Gprime2','Free','Real',ipGprime2,NumInt3)
call GetMem('Gprime3','Free','Real',ipGprime3,NumInt3)
call GetMem('Gprime4','Free','Real',ipGprime4,NumInt3)

end subroutine CalcGdbleprime
