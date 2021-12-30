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

subroutine XDIAXT(XDX,X,DIA,NDIM,SCR)
! Obtain XDX = X * DIA * X(Transposed)
! where DIA is an diagonal matrix stored as a vector

implicit real*8(A-H,O-Z)
dimension XDX(NDIM,NDIM)
dimension X(NDIM,NDIM), DIA(NDIM)
dimension SCR(NDIM,NDIM)

! DIA * X(transposed)
do I=1,NDIM
  call COPVEC(X(1,I),SCR(1,I),NDIM)
  call SCALVE(SCR(1,I),DIA(I),NDIM)
end do
! X * DIA * X(Transposed)
call MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)

return

end subroutine XDIAXT
