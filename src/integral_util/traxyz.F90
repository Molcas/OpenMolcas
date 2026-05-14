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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Traxyz(nZeta,la,WInt,Scr,A)
!***********************************************************************
!                                                                      *
! Object: to transform the well-integrals from the local coordinate    *
!         system to the global one. The transformation matrix is       *
!         stored in A.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!              October '92.                                            *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la
real(kind=wp), intent(inout) :: WInt(nZeta*3**(la-1),3)
real(kind=wp), intent(out) :: Scr(nZeta*3**la)
real(kind=wp), intent(in) :: A(nZeta,3,3)
integer(kind=iwp) :: i, iLen, iOff, iVec, iZeta, jVec, mLen, nLen

#ifdef _DEBUGPRINT_
call RecPrt(' Enter Traxyz: WInt',' ',Wint,nZeta,3**la)
call RecPrt(' The transformation matrix',' ',A,nZeta,9)
#endif

! Transform

nLen = 3
mLen = 3**(la-1)
!write(u6,*) ' nLen, mLen=',nLen, mLen
! Loop over all angular indices to be transformed.
do i=1,la
  iVec = 0
  do iLen=1,mLen
    iOff = (iLen-1)*nLen
    do iZeta=1,nZeta
      iVec = iVec+1

      jVec = iOff*nZeta+iZeta
      Scr(jVec) = A(iZeta,1,1)*WInt(iVec,1)+A(iZeta,1,2)*WInt(iVec,2)+A(iZeta,1,3)*WInt(iVec,3)
      jVec = (iOff+1)*nZeta+iZeta
      Scr(jVec) = A(iZeta,2,1)*WInt(iVec,1)+A(iZeta,2,2)*WInt(iVec,2)+A(iZeta,2,3)*WInt(iVec,3)
      jVec = (iOff+2)*nZeta+iZeta
      Scr(jVec) = A(iZeta,3,1)*WInt(iVec,1)+A(iZeta,3,2)*WInt(iVec,2)+A(iZeta,3,3)*WInt(iVec,3)

    end do
  end do
  !call RecPrt(' Partially transformed WInt',' ',Scr,nZeta,mLen*3)
  WInt(:,:) = reshape(Scr(:),[nZeta*mLen,3])
end do

#ifdef _DEBUGPRINT_
call RecPrt('Exit Traxyz :Global well integrals',' ',WInt,nZeta,mLen*3)
#endif

end subroutine Traxyz
