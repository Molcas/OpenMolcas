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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************
! Version of Oct 21

subroutine SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
! Calculate square root of positive definite symmetric matrix A
! if (ITASK  /=  2 ) Inverted square root matrix is also calculated
! In case of singularities in A A -1/2 is defined to have the same
! singularity

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NDIM, ITASK
real(kind=wp), intent(in) :: A(NDIM,NDIM)
real(kind=wp), intent(out) :: ASQRT(NDIM,NDIM), AMSQRT(NDIM,NDIM)
real(kind=wp), intent(_OUT_) :: SCR(*)
integer(kind=iwp) :: I, KLASYM, KLAVAL, KLAVEC, KLFREE, NTEST

! Length of SCR should at least be 2 * NDIM ** 2 + NDIM*(NDIM+1)/2
KLFREE = 1

KLASYM = KLFREE
KLAVAL = KLASYM
KLFREE = KLASYM+NDIM*(NDIM+1)/2

KLAVEC = KLFREE
KLFREE = KLFREE+NDIM**2

NTEST = 0

! TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
call TRIPAK(A,SCR(KLASYM),1,NDIM,NDIM)
call unitmat(SCR(KLAVEC),NDIM)
call NIDiag(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM)
call JACORD(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM)
do I=2,NDIM
  SCR(KLAVAL-1+I) = SCR(KLASYM-1+I*(I+1)/2)
end do
if (NTEST >= 1) then
  write(u6,*) ' Eigenvalues of matrix : '
  call WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
end if
! Check for negative eigenvalues
do I=1,NDIM
  if (abs(SCR(KLAVAL-1+I)) < 1.0e-14_wp) SCR(KLAVAL) = Zero
  if (SCR(KLAVAL-1+I) < Zero) then
    !write(u6,*) ' SQRTMT : Negative eigenvalue ',SCR(KLAVAL-1+I)
    !write(u6,*) ' SQRTMT : I will STOP '
    !stop ' SQRTMT : Negative eigenvalue '
    call SYSABENDMSG('sqrtmt','Internal error','Negative eigenvalue')
  end if
end do

do I=1,NDIM
  SCR(KLAVAL-1+I) = sqrt(SCR(KLAVAL-1+I))
end do
! XDIAXT(XDX,X,DIA,NDIM,SCR)
call XDIAXT(ASQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))

if (ITASK == 2) then
  do I=1,NDIM
    if (SCR(KLAVAL-1+I) > 1.0e-13_wp) then
      SCR(KLAVAL-1+I) = One/SCR(KLAVAL-1+I)
    else
      SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
    end if
  end do
  call XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
end if

if (NTEST >= 1) then
  write(u6,*) ' Info from SQRTMT '
  write(u6,*) ' ================='
  write(u6,*) ' Input matrix to SQRTMT '
  call WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
  write(u6,*) ' Square root of matrix '
  call WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
  if (ITASK == 2) then
    write(u6,*) ' Inverse square root of matrix '
    call WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
  end if
end if

return

end subroutine SQRTMT
