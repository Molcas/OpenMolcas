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

subroutine ReaInd(iGP,iDT,iDistIm,iCNum,indma,ncparm,Sum1,s90um)

use qmstat_global, only: CordIm, Cordst, DipIm, nCent, nCha, nPart, nPol, Qimp, Qsta
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iGP(3), iDT(3), iDistIm, iCNum, indma, ncparm
real(kind=wp) :: Sum1, s90um
#include "WrkSpc.fh"
integer(kind=iwp) :: i, Inc, Inc2, j, k, l
real(kind=wp) :: D1x, D1y, D1z, Q1, Q2, x, X1, y,z

Sum1 = 0
!irekn = 0
!xled = 0
!yled = 0
!zled = 0
do i=1+(nPol*iCnum),indma
  ! The energy of the induced dipoles (iDT) in the field from the real charges. Yes, this
  ! is the polarization energy in a system without permanent dipoles, see good old
  ! Bottcher, eq. (3.129). Observe also that we here include effects of the reaction
  ! field on the induced dipoles, see the polarization loop.
  do j=1,3  !The energy of the induced dipoles (iDT) in
    Sum1 = Sum1+Work(iGP(j)+i-1)*Work(iDT(j)+i-1)
  end do
  ! IF WE WISH TO MONITOR THE INDUCED DIPOLES, UNCOMMENT THIS, AND THE COMMENTED THING ABOVE.
  !irekn = irekn+1
  !xled = xled+Work(iDt(1)+i-1)
  !yled = yled+Work(iDt(2)+i-1)
  !zled = zled+Work(iDt(3)+i-1)
  !if (irekn == 3) then
  !  irekn = 0
  !  TOT = sqrt(xled**2+yled**2+zled**2)
  !  write(u6,*) 'HHH',TOT
  !  xled = 0
  !  yled = 0
  !  zled = 0
  !end if
end do
Sum1 = Sum1*Half
S90um = Zero
!Energy of charge distribution in the reaction field to the induced dipoles.
!Once more, see Bottcher eq. (4.69): we are computing the product between charges and the
!potential connected with the reaction field.
do i=iCnum+1,nPart
  do j=1,nPol
    Q1 = Qimp((i-1)*nPol+j)
    D1x = DipIm(1,(i-1)*nPol+j)
    D1y = DipIm(2,(i-1)*nPol+j)
    D1z = DipIm(3,(i-1)*nPol+j)
    x = CordIm(1,(i-1)*nCent+j)
    y = CordIm(2,(i-1)*nCent+j)
    z = CordIm(3,(i-1)*nCent+j)
    Inc = ncparm*nCent*(i-(iCnum+1))+(j-1)*ncparm
    do l=nCent-nCha+1,nCent
      Inc2 = Inc+l
      Q2 = Qsta(l-nCent+nCha)
      do k=iCnum+1,nPart
        X1 = (X-Cordst(1,l+(k-1)*nCent))*D1x
        X1 = (Y-Cordst(2,l+(k-1)*nCent))*D1y+X1
        X1 = (Z-Cordst(3,l+(k-1)*nCent))*D1z+X1
        !Change sign on Q2 since we are in the backwards land, while Q1 and X1 already are backward.
        S90um = S90um-(Q1+X1*Work(iDistIm-1+inc2+(k-(iCnum+1))*nCent)**2)*Q2*Work(iDistIm-1+inc2+(k-(iCnum+1))*nCent)
      end do
    end do
  end do
end do

return

end subroutine ReaInd
