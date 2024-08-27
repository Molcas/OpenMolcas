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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine ConMax(A,nPrim,mPrim,B,nCont,C,mCont)
!***********************************************************************
!                                                                      *
! Object: to find the largest element in the contraction matrix  for   *
!         each primitive index.                                        *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             July '91                                                 *
!***********************************************************************

implicit none
integer, intent(In) :: nPrim, nCont, mPrim, mCont
real*8, intent(In) :: B(nPrim,nCont), C(mPrim,mCont)
real*8, intent(Out) :: A(nPrim,mPrim)
integer iPrim, jPrim
real*8, external :: DDot_
real*8 Temp

do iPrim=1,nPrim
  Temp = DDot_(nCont,B(iPrim,1),nPrim,B(iPrim,1),nPrim)
  do jPrim=1,mPrim
    A(iPrim,jPrim) = Temp
  end do
end do

do jPrim=1,mPrim
  Temp = DDot_(mCont,C(jPrim,1),mPrim,C(jPrim,1),mPrim)
  do iPrim=1,nPrim
    A(iPrim,jPrim) = sqrt(A(iPrim,jPrim)*Temp)
  end do
end do

return

end subroutine ConMax
