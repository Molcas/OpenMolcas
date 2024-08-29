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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nab1, na1b, nab, na, nb, na1, nb1, nPrim, la, lb
real(kind=wp), intent(out) :: ab1(nPrim,nab1)
real(kind=wp), intent(in) :: a1b(nPrim,na1b), cffAB(3), ab(nPrim,nab)
integer(kind=iwp) :: i, ipA1B, ipAB, ipAB1, ipxyz, ixa, ixb, ixyza, ixyza1, ixyzb, ixyzb1, iya, iyb, iza, izb
#ifdef _DEBUGPRINT_
character(len=72) :: Label
#endif
! Statement functions
integer(kind=iwp) :: iy, iz, ixyz, Ind, nElem
Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2+iz+1
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

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
    ixyzb1 = Ind(iyb,izb)

    do ixa=na,0,-1
      do iya=na-ixa,0,-1
        iza = na-ixa-iya

        ixyza = Ind(iya,iza)

        ! Find a angular index which can be decremented

        if (ixb /= 0) then
          ipxyz = 1
          ixyza1 = Ind(iya,iza)
          ixyzb = Ind(iyb,izb)
        else if (iyb /= 0) then
          ipxyz = 2
          ixyza1 = Ind(iya+1,iza)
          ixyzb = Ind(iyb-1,izb)
        else
          ipxyz = 3
          ixyza1 = Ind(iya,iza+1)
          ixyzb = Ind(iyb,izb-1)
        end if
        if (la >= lb) then
          ipab1 = ixyza+nElem(na)*(ixyzb1-1)
          ipa1b = ixyza1+nElem(na1)*(ixyzb-1)
          ipab = ixyza+nElem(na)*(ixyzb-1)
        else
          ipab1 = ixyzb1+nElem(nb1)*(ixyza-1)
          ipa1b = ixyzb+nElem(nb)*(ixyza1-1)
          ipab = ixyzb+nElem(nb)*(ixyza-1)
        end if
        if (cffAB(ipxyz) /= Zero) then
          call DZaXpY(nPrim,cffAB(ipxyz),ab(1,ipab),1,a1b(1,ipa1b),1,ab1(1,ipab1),1)
        else
          !call dcopy_(nPrim,a1b(1,ipa1b),1,ab1(1,ipab1),1)
          do i=1,nPrim
            ab1(i,ipab1) = a1b(i,ipa1b)
          end do
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
