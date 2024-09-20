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

!#define _DEBUGPRINT_
subroutine HRR1(ab1,nab1,a1b,na1b,cffAB,ab,nab,na,nb,na1,nb1,nPrim,la,lb)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             June '91                                                 *
!***********************************************************************

use Index_Functions, only: C_Ind3, nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nab1, na1b, nab, na, nb, na1, nb1, nPrim, la, lb
real(kind=wp), intent(out) :: ab1(nPrim,nab1)
real(kind=wp), intent(in) :: a1b(nPrim,na1b), cffAB(3), ab(nPrim,nab)
integer(kind=iwp) :: ipA1B, ipAB, ipAB1, ipxyz, ixa, ixb, ixyza, ixyza1, ixyzb, ixyzb1, iya, iyb, iza, izb
#ifdef _DEBUGPRINT_
character(len=72) :: Label
#endif

#ifdef _DEBUGPRINT_
write(Label,'(A,i1,A,i1,A)') ' Source: (',na1,',',nb,'|'
call RecPrt(Label,' ',a1b,nPrim,na1b)
call RecPrt(' Coordinates (A-B)',' ',cffAB,1,3)
write(Label,'(A,i1,A,i1,A)') ' Source: (',na,',',nb,'|'
call RecPrt(Label,' ',ab,nPrim,nab)
#endif

! Loop over indices of the target batch

do ixb=nb1,0,-1
  do iyb=nb1-ixb,0,-1
    izb = nb1-ixb-iyb
    ixyzb1 = C_Ind3(ixb,iyb,izb)

    do ixa=na,0,-1
      do iya=na-ixa,0,-1
        iza = na-ixa-iya

        ixyza = C_Ind3(ixa,iya,iza)

        ! Find a angular index which can be decremented

        if (ixb /= 0) then
          ipxyz = 1
          ixyza1 = C_Ind3(ixa,iya,iza)
          ixyzb = C_Ind3(ixb,iyb,izb)
        else if (iyb /= 0) then
          ipxyz = 2
          ixyza1 = C_Ind3(ixa,iya+1,iza)
          ixyzb = C_Ind3(ixb,iyb-1,izb)
        else
          ipxyz = 3
          ixyza1 = C_Ind3(ixa,iya,iza+1)
          ixyzb = C_Ind3(ixb,iyb,izb-1)
        end if
        if (la >= lb) then
          ipab1 = ixyza+nTri_Elem1(na)*(ixyzb1-1)
          ipa1b = ixyza1+nTri_Elem1(na1)*(ixyzb-1)
          ipab = ixyza+nTri_Elem1(na)*(ixyzb-1)
        else
          ipab1 = ixyzb1+nTri_Elem1(nb1)*(ixyza-1)
          ipa1b = ixyzb+nTri_Elem1(nb)*(ixyza1-1)
          ipab = ixyzb+nTri_Elem1(nb)*(ixyza-1)
        end if
        if (cffAB(ipxyz) /= Zero) then
          ab1(:,ipab1) = a1b(:,ipa1b)+cffAB(ipxyz)*ab(:,ipab)
        else
          ab1(:,ipab1) = a1b(:,ipa1b)
        end if
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(Label,'(A,i1,A,i1,A)') ' Target: (',na,',',nb1,'|'
call RecPrt(Label,' ',ab1,nPrim,nab1)
#endif

return

end subroutine HRR1
