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
!               1996, Hans-Joachim Werner                              *
!***********************************************************************

subroutine Assg1(Temp,PAO,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,IfGrad,Indx,mVec)
!***********************************************************************
!                                                                      *
! Object: to assemble the gradients of the ERI's.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91; modified by H.-J. Werner, Mai 1996          *
!***********************************************************************

use Index_Functions, only: C_Ind3_Rev, nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Temp(9)
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, Indx(3,4)
real(kind=wp), intent(in) :: PAO(nT,nTri_Elem1(la),nTri_Elem1(lb),nTri_Elem1(lc),nTri_Elem1(ld)), &
                             xyz2D0(nRys,nT,0:la+1,0:lb+1,0:lc+1,0:ld+1,3), xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp), intent(out) :: mVec
#include "itmax.fh"
integer(kind=iwp) :: i, iCent, icir(3), Ind1(3,3), Ind2(3,3), ipa, ipb, ipc, ipd, ixa, ixabcd, ixb, ixbcd, ixc, ixcd, ixd, iya, &
                     iyabcd, iyb, iybcd, iyc, iycd, iyd, iza, izb, izc, izd, nVec(3)

Temp(:) = Zero

mVec = 0
nVec(:) = 0
do i=1,3       ! Cartesian directions
  do iCent=1,4 ! Centers of integral
    if (IfGrad(i,iCent)) then
      mVec = mVec+1
      nVec(i) = nVec(i)+1
      Ind1(nVec(i),i) = 3*(Indx(i,iCent)-1)+i
      Ind2(nVec(i),i) = mVec
    end if
  end do
end do

do ipd=1,nTri_Elem1(ld)
  icir(:) = C_Ind3_Rev(ipd,ld)
  ixd = icir(1)
  iyd = icir(2)
  izd = icir(3)

  do ipc=1,nTri_Elem1(lc)
    icir(:) = C_Ind3_Rev(ipc,lc)
    ixc = icir(1)
    iyc = icir(2)
    izc = icir(3)

    ixcd = ixc+ixd
    iycd = iyc+iyd

    do ipb=1,nTri_Elem1(lb)
      icir(:) = C_Ind3_Rev(ipb,lb)
      ixb = icir(1)
      iyb = icir(2)
      izb = icir(3)

      ixbcd = ixcd+ixb
      iybcd = iycd+iyb

      do ipa=1,nTri_Elem1(la)
        icir(:) = C_Ind3_Rev(ipa,la)
        ixa = icir(1)
        iya = icir(2)
        iza = icir(3)

        ixabcd = ixbcd+ixa
        iyabcd = iybcd+iya

        ! Compute all desired gradients with respect to an x-component.

        if (iyabcd /= 0) then

          select case (nVec(1))
            case (3)
              call ass3(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)), &
                        xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(3,1)),PAO(1,ipa,ipb,ipc,ipd), &
                        Temp(Ind2(1,1)),Temp(Ind2(2,1)),Temp(Ind2(3,1)),nT,nRys)
            case (2)
              call ass2(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)), &
                        xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),Temp(Ind2(2,1)),nT,nRys)
            case (1)
              call ass1(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)), &
                        PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),nT,nRys)
          end select

        else

          select case (nVec(1))
            case (3)
              call ass3a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)), &
                         xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(3,1)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,1)),Temp(Ind2(2,1)),Temp(Ind2(3,1)),nT,nRys)
            case (2)
              call ass2a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)), &
                         xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(2,1)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,1)),Temp(Ind2(2,1)),nT,nRys)
            case (1)
              call ass1a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,ixa,ixb,ixc,ixd,Ind1(1,1)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,1)),nT,nRys)
          end select

        end if

        ! Compute all desired gradients with respect to a y-component.

        if (ixabcd /= 0) then

          select case (nVec(2))
            case (3)
              call ass3(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)), &
                        xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(3,2)),PAO(1,ipa,ipb,ipc,ipd), &
                        Temp(Ind2(1,2)),Temp(Ind2(2,2)),Temp(Ind2(3,2)),nT,nRys)
            case (2)
              call ass2(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)), &
                        xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),Temp(Ind2(2,2)),nT,nRys)
            case (1)
              call ass1(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)), &
                        PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),nT,nRys)
          end select

        else

          select case (nVec(2))
            case (3)
              call ass3a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)), &
                         xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(3,2)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,2)),Temp(Ind2(2,2)),Temp(Ind2(3,2)),nT,nRys)
            case (2)
              call ass2a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)), &
                         xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(2,2)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,2)),Temp(Ind2(2,2)),nT,nRys)
            case (1)
              call ass1a(xyz2D0(1,1,iza,izb,izc,izd,3),xyz2D1(1,1,iya,iyb,iyc,iyd,Ind1(1,2)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,2)),nT,nRys)
          end select

        end if

        ! Compute all desired gradients with respect to a z-component.

        if (ixabcd*iyabcd /= 0) then

          select case (nVec(3))
            case (3)
              call ass3(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                        xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),PAO(1,ipa,ipb,ipc,ipd), &
                        Temp(Ind2(1,3)),Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
            case (2)
              call ass2(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                        xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),Temp(Ind2(2,3)),nT,nRys)
            case (1)
              call ass1(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                        PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),nT,nRys)
          end select

        else if ((ixabcd == 0) .and. (iyabcd /= 0)) then

          select case (nVec(3))
            case (3)
              call ass3a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                         xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,3)),Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
            case (2)
              call ass2a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                         xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),Temp(Ind2(2,3)),nT,nRys)
            case (1)
              call ass1a(xyz2D0(1,1,iya,iyb,iyc,iyd,2),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,3)),nT,nRys)
          end select

        else if ((iyabcd == 0) .and. (ixabcd /= 0)) then

          select case (nVec(3))
            case (3)
              call ass3a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                         xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,3)),Temp(Ind2(2,3)),Temp(Ind2(3,3)),nT,nRys)
            case (2)
              call ass2a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)), &
                         xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),Temp(Ind2(2,3)),nT,nRys)
            case (1)
              call ass1a(xyz2D0(1,1,ixa,ixb,ixc,ixd,1),xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,3)),nT,nRys)
          end select

        else

          select case (nVec(3))
            case (3)
              call ass3b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)), &
                         xyz2D1(1,1,iza,izb,izc,izd,Ind1(3,3)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),Temp(Ind2(2,3)), &
                         Temp(Ind2(3,3)),nT,nRys)
            case (2)
              call ass2b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),xyz2D1(1,1,iza,izb,izc,izd,Ind1(2,3)),PAO(1,ipa,ipb,ipc,ipd), &
                         Temp(Ind2(1,3)),Temp(Ind2(2,3)),nT,nRys)
            case (1)
              call ass1b(xyz2D1(1,1,iza,izb,izc,izd,Ind1(1,3)),PAO(1,ipa,ipb,ipc,ipd),Temp(Ind2(1,3)),nT,nRys)
          end select

        end if

      end do

    end do

  end do

end do

return

end subroutine Assg1
