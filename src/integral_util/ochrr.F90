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

subroutine OCHRR(target,nPrim,nTrgt,la,lb,ipRs)
!***********************************************************************
! Object: this is a One Center HRR routine.                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

implicit none
integer, intent(In) :: nPrim, nTrgt, la, lb
integer, intent(Out) :: ipRs
real*8, intent(InOut) :: target(nPrim,nTrgt)
integer i, ixyz, ix, iz, nElem, Ind
integer iout, ixb, iyb, izb, ixyzb, iybMax, ixa, iyaMax, ixab, iya, iza, izab, ixyza, iTo, iFrom, iab
! Statment functions
nElem(i) = (i+1)*(i+2)/2
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

if ((la == 0) .or. (lb == 0)) then
  ipRs = 1
  return
end if
iab = 0
iout = nElem(la+lb)
ipRs = iout*nPrim+1
do ixb=0,lb
  iybMax = lb-ixb
  do iyb=0,iybMax
    izb = iybMax-iyb
    ixyzb = Ind(lb,ixb,izb)

    do ixa=0,la
      iyaMax = la-ixa
      ixab = ixa+ixb
      do iya=0,iyaMax
        iza = iyaMax-iya
        izab = iza+izb
        ixyza = Ind(la,ixa,iza)
        iTo = iout+nElem(la)*(ixyzb-1)+ixyza
        iFrom = iab+Ind(la+lb,ixab,izab)

        target(:,iTo) = target(:,iFrom)

      end do
    end do
  end do
end do

return

end subroutine OCHRR
