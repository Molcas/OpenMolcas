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

subroutine Assg1_mck(g1,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,IfGrad,Index,mVec,Index2)
!***********************************************************************
!                                                                      *
! Object: to assemble the gradients of the ERI's.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

implicit real*8(A-H,O-Z)
!#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "iavec.fh"
real*8 g1(nT,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,(lc+1)*(lc+2)/2,(ld+1)*(ld+2)/2,9), xyz2D0(nRys,nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), &
       xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9)
logical IfGrad(3,4)
integer Ind1(3), Ind2(3), index(3,4), Index2(3,4)
! Statement function
nElem(i) = (i+1)*(i+2)/2

ka = (la+1)*(la+2)/2
kb = (lb+1)*(lb+2)/2
kc = (lc+1)*(lc+2)/2
kd = (ld+1)*(ld+2)/2
nG1 = nT*9*ka*kb*kc*kd
call dcopy_(nG1,[Zero],0,G1,1)
call ICOPY(12,[0],0,Index2,1)

ii = la*(la+1)*(la+2)/6
jj = lb*(lb+1)*(lb+2)/6
kk = lc*(lc+1)*(lc+2)/6
ll = ld*(ld+1)*(ld+2)/6

do ipa=1,nElem(la)
  ipaii = ipa+ii
  ixa = ixyz(1,ipaii)
  iya = ixyz(2,ipaii)
  iza = ixyz(3,ipaii)

  do ipb=1,nElem(lb)
    ipbjj = ipb+jj
    ixb = ixyz(1,ipbjj)
    iyb = ixyz(2,ipbjj)
    izb = ixyz(3,ipbjj)

    ixab = ixa+ixb
    iyab = iya+iyb

    do ipc=1,nElem(lc)
      ipckk = ipc+kk
      ixc = ixyz(1,ipckk)
      iyc = ixyz(2,ipckk)
      izc = ixyz(3,ipckk)

      ixabc = ixab+ixc
      iyabc = iyab+iyc

      do ipd=1,nElem(ld)
        ipdll = ipd+ll
        ixd = ixyz(1,ipdll)
        iyd = ixyz(2,ipdll)
        izd = ixyz(3,ipdll)

        ixabcd = ixabc+ixd
        iyabcd = iyabc+iyd

        ! Compute all desired gradients with respect to an x-component.

        mVec = 0
        nVec = 0
        do iCent=1,4
          if (IfGrad(1,iCent)) then
            mVec = mVec+1
            nVec = nVec+1
            Ind1(nVec) = 3*(index(1,iCent)-1)+1
            Ind2(nVec) = mVec
            Index2(1,iCent) = mVec
          end if
        end do

        if (iyabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        else

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,ixa,ixb,ixc,ixd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        end if

        ! Compute all desired gradients with respect to a y-component.

        nVec = 0
        do iCent=1,4
          if (IfGrad(2,iCent)) then
            mVec = mVec+1
            nVec = nVec+1
            Ind1(nVec) = 3*(index(2,iCent)-1)+2
            Ind2(nVec) = mVec
            Index2(2,iCent) = mVec
          end if
        end do

        if (ixabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        else

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iza,izb,izc,izd,3)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iya,iyb,iyc,iyd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        end if

        ! Compute all desired gradients with respect to a z-component.

        nVec = 0
        do iCent=1,4
          if (IfGrad(3,iCent)) then
            mVec = mVec+1
            nVec = nVec+1
            Ind1(nVec) = 3*(index(3,iCent)-1)+3
            Ind2(nVec) = mVec
            Index2(3,iCent) = mVec
          end if
        end do

        if (ixabcd*iyabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)*xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        else if ((ixabcd == 0) .and. (iyabcd /= 0)) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,iya,iyb,iyc,iyd,2)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        else if ((iyabcd == 0) .and. (ixabcd /= 0)) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp = xyz2D0(iRys,iT,ixa,ixb,ixc,ixd,1)
                  tmp1 = tmp1+tmp*xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        else

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                tmp3 = 0.0d0
                do iRys=1,nRys
                  tmp1 = tmp1+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                  tmp3 = tmp3+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(3))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
                g1(iT,ipa,ipb,ipc,ipd,Ind2(3)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(3))+tmp3
              end do
            case (2)
              do iT=1,nT
                tmp1 = 0.0d0
                tmp2 = 0.0d0
                do iRys=1,nRys
                  tmp1 = tmp1+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = 0.0d0
                do iRys=1,nRys
                  tmp1 = tmp1+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
              end do
          end select

        end if

      end do

    end do

  end do

end do

return

end subroutine Assg1_mck
