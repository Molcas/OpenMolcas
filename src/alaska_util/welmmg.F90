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

subroutine WelMmG(nHer,MmWelG,la,lb,lr)

! Statement function
nElem(i) = (i+1)*(i+2)/2

k = la+lb+1
jsum = 1
do i=1,k
  jsum = jsum+3**i
end do
nHer = 1
MmWelG = 2*jsum+max((k+1)*(k/2+1)*(k/4+1)+1,9+3**k,5)

! Add memory for contributions to the derivative

MmWelG = MmWelG+nElem(la+1)*nElem(lb)
if (la >= 1) MmWelG = MmWelG+nElem(la-1)*nElem(lb)
MmWelG = MmWelG+nElem(la)*nElem(lb+1)
if (lb >= 1) MmWelG = MmWelG+nElem(la)*nElem(lb-1)
MmWelG = MmWelG+2

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lr)

end subroutine WelMmG
