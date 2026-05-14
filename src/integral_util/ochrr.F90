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

subroutine OCHRR(tgt,nPrim,nTrgt,la,lb,ipRs)
!***********************************************************************
! Object: this is a One Center HRR routine.                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             May '90                                                  *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrim, nTrgt, la, lb
real(kind=wp), intent(inout) :: tgt(nPrim,nTrgt)
integer(kind=iwp), intent(out) :: ipRs
integer(kind=iwp) :: iab, iFrom, iout, iTo, ixa, ixab, ixb, ixyza, ixyzb, iya, iyaMax, iyb, iybMax, iza, izab, izb

if ((la == 0) .or. (lb == 0)) then
  ipRs = 1
else
  iab = 0
  iout = nTri_Elem1(la+lb)
  ipRs = iout*nPrim+1
  do ixb=0,lb
    iybMax = lb-ixb
    do iyb=0,iybMax
      izb = iybMax-iyb
      ixyzb = C_Ind(lb,ixb,izb)

      do ixa=0,la
        iyaMax = la-ixa
        ixab = ixa+ixb
        do iya=0,iyaMax
          iza = iyaMax-iya
          izab = iza+izb
          ixyza = C_Ind(la,ixa,iza)
          iTo = iout+nTri_Elem1(la)*(ixyzb-1)+ixyza
          iFrom = iab+C_Ind(la+lb,ixab,izab)

          tgt(:,iTo) = tgt(:,iFrom)

        end do
      end do
    end do
  end do
end if

return

end subroutine OCHRR
