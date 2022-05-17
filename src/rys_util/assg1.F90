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

subroutine Assg1(Temp,PAO,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,IfGrad,Index,mVec)
!***********************************************************************
!                                                                      *
! Object: to assemble the gradients of the ERI's.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91; modified by H.-J. Werner, Mai 1996          *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "iavec.fh"
real*8 PAO(nT,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,(lc+1)*(lc+2)/2,(ld+1)*(ld+2)/2), xyz2D0(nRys,nT,0:la+1,0:lb+1,0:lc+1,0:ld+1,3), &
       xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9), Temp(9)
logical IfGrad(3,4)
integer Ind1(3,3), Ind2(3,3), index(3,4), nVec(3)
! Statement function
nElem(i) = (i+1)*(i+2)/2

call dcopy_(9,[Zero],0,Temp,1)

ii = la*(la+1)*(la+2)/6
jj = lb*(lb+1)*(lb+2)/6
kk = lc*(lc+1)*(lc+2)/6
ll = ld*(ld+1)*(ld+2)/6

mVec = 0
do i=1,3       ! Cartesian directions
  nVec(i) = 0
  do iCent=1,4 ! Centers of integral
    if (IfGrad(i,iCent)) then
      mVec = mVec+1
      nVec(i) = nVec(i)+1
      Ind1(nVec(i),i) = 3*(index(i,iCent)-1)+i
      Ind2(nVec(i),i) = mVec
    end if
  end do
end do

do ipd=1,nElem(ld)
  ixd = ixyz(1,ll+ipd)
  iyd = ixyz(2,ll+ipd)
  izd = ixyz(3,ll+ipd)

  do ipc=1,nElem(lc)
    ixc = ixyz(1,kk+ipc)
    iyc = ixyz(2,kk+ipc)
    izc = ixyz(3,kk+ipc)

    ixcd = ixc+ixd
    iycd = iyc+iyd

    do ipb=1,nElem(lb)
      ixb = ixyz(1,jj+ipb)
      iyb = ixyz(2,jj+ipb)
      izb = ixyz(3,jj+ipb)

      ixbcd = ixcd+ixb
      iybcd = iycd+iyb

      do ipa=1,nElem(la)
        ixa = ixyz(1,ii+ipa)
        iya = ixyz(2,ii+ipa)
        iza = ixyz(3,ii+ipa)

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
