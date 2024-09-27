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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Setup1(ExpA,nPrim,ExpB,mPrim,A,B,rKappa,Pcoor,ZInv)
!***********************************************************************
!                                                                      *
!     Object : to compute some data which is needed for the one-       *
!              electron integrals.                                     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nPrim, mPrim
real(kind=wp), intent(in) :: ExpA(nPrim), ExpB(mPrim), A(3), B(3), ZInv(nPrim,mPrim)
real(kind=wp), intent(out) :: rKappa(nPrim,mPrim), Pcoor(nPrim,mPrim,3)
integer(kind=iwp) :: jPrim
real(kind=wp) :: ab

#ifdef _DEBUGPRINT_
call RecPrt(' *** ExpA ***',' ',ExpA,1,nPrim)
call RecPrt(' *** ExpB ***',' ',ExpB,1,mPrim)
call RecPrt(' *** ZInv ***',' ',ZInv,nPrim,mPrim)
write(u6,*) 'A(:)=',A(:)
write(u6,*) 'B(:)=',B(:)
#endif
ab = (A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2

if (ab /= Zero) then
  do jPrim=1,mPrim
    rKappa(:,jPrim) = exp(-ExpA(:)*ExpB(jPrim)*ab*ZInv(:,jPrim))
    Pcoor(:,jPrim,1) = (ExpA(:)*A(1)+ExpB(jPrim)*B(1))*ZInv(:,jPrim)
    Pcoor(:,jPrim,2) = (ExpA(:)*A(2)+ExpB(jPrim)*B(2))*ZInv(:,jPrim)
    Pcoor(:,jPrim,3) = (ExpA(:)*A(3)+ExpB(jPrim)*B(3))*ZInv(:,jPrim)
  end do
else
  rKappa(:,:) = One
  PCoor(:,:,1) = A(1)
  PCoor(:,:,2) = A(2)
  PCoor(:,:,3) = A(3)
end if
#ifdef _DEBUGPRINT_
call RecPrt(' *** Kappa ***',' ',rKappa,nPrim,mPrim)
call RecPrt(' ***   Px  ***',' ',Pcoor(1,1,1),nPrim,mPrim)
call RecPrt(' ***   Py  ***',' ',Pcoor(1,1,2),nPrim,mPrim)
call RecPrt(' ***   Pz  ***',' ',Pcoor(1,1,3),nPrim,mPrim)
#endif

end subroutine Setup1
