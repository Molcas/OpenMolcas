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

implicit real*8(A-H,O-Z)
dimension A(NDIM,NDIM)
dimension ASQRT(NDIM,NDIM), AMSQRT(NDIM,NDIM)
dimension SCR(*)

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
call DCopy_(NDIM**2,[0.0d0],0,SCR(KLAVEC),1)
call DCopy_(NDIM,[1.0d0],0,SCR(KLAVEC),1+NDIM)
call NIDiag(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM)
call JACORD(SCR(KLASYM),SCR(KLAVEC),NDIM,NDIM)
call COPDIA(SCR(KLASYM),SCR(KLAVAL),NDIM,1)
if (NTEST >= 1) then
  write(6,*) ' Eigenvalues of matrix : '
  call WRTMAT(SCR(KLAVAL),NDIM,1,NDIM,1)
end if
! Check for negative eigenvalues
do I=1,NDIM
  if (abs(SCR(KLAVAL-1+I)) < 1.0D-14) SCR(KLAVAL) = 0.0d0
  if (SCR(KLAVAL-1+I) < 0.0d0) then
    ! WRITE(6,*) ' SQRTMT : Negative eigenvalue ', SCR(KLAVAL-1+I)
    ! WRITE(6,*) ' SQRTMT : I will STOP '
    ! STOP       ' SQRTMT : Negative eigenvalue '
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
    if (SCR(KLAVAL-1+I) > 1.0D-13) then
      SCR(KLAVAL-1+I) = 1.0d0/SCR(KLAVAL-1+I)
    else
      SCR(KLAVAL-1+I) = SCR(KLAVAL-1+I)
    end if
  end do
  call XDIAXT(AMSQRT,SCR(KLAVEC),SCR(KLAVAL),NDIM,SCR(KLFREE))
end if

if (NTEST >= 1) then
  write(6,*) ' Info from SQRTMT '
  write(6,*) ' ================='
  write(6,*) ' Input matrix to SQRTMT '
  call WRTMAT(A,NDIM,NDIM,NDIM,NDIM)
  write(6,*) ' Square root of matrix '
  call WRTMAT(ASQRT,NDIM,NDIM,NDIM,NDIM)
  if (ITASK == 2) then
    write(6,*) ' Inverse square root of matrix '
    call WRTMAT(AMSQRT,NDIM,NDIM,NDIM,NDIM)
  end if
end if

return

end subroutine SQRTMT
