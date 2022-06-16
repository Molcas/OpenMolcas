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

subroutine Assg1_mck(g1,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,IfGrad,Indx,mVec,Indx2)
!***********************************************************************
!                                                                      *
! Object: to assemble the gradients of the ERI's.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Index_Functions, only: C_Ind3_Rev, nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, Indx(3,4)
real(kind=wp), intent(out) :: g1(nT,nTri_Elem1(la),nTri_Elem1(lb),nTri_Elem1(lc),nTri_Elem1(ld),9)
real(kind=wp), intent(in) :: xyz2D0(nRys,nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp), intent(out) :: mVec, Indx2(3,4)
#include "itmax.fh"
integer(kind=iwp) :: iCent, icir(3), Ind1(3), Ind2(3), ipa, ipb, ipc, ipd, iRys, iT, ixa, ixab, ixabc, ixabcd, ixb, ixc, ixd, iya, &
                     iyab, iyabc, iyabcd, iyb, iyc, iyd, iza, izb, izc, izd, nVec
real(kind=wp) :: tmp, tmp1, tmp2, tmp3

g1(:,:,:,:,:,:) = Zero
Indx2(:,:) = 0

do ipa=1,nTri_Elem1(la)
  icir(:) = C_Ind3_Rev(ipa,la)
  ixa = icir(1)
  iya = icir(2)
  iza = icir(3)

  do ipb=1,nTri_Elem1(lb)
    icir(:) = C_Ind3_Rev(ipb,lb)
    ixb = icir(1)
    iyb = icir(2)
    izb = icir(3)

    ixab = ixa+ixb
    iyab = iya+iyb

    do ipc=1,nTri_Elem1(lc)
      icir(:) = C_Ind3_Rev(ipc,lc)
      ixc = icir(1)
      iyc = icir(2)
      izc = icir(3)

      ixabc = ixab+ixc
      iyabc = iyab+iyc

      do ipd=1,nTri_Elem1(ld)
        icir(:) = C_Ind3_Rev(ipd,ld)
        ixd = icir(1)
        iyd = icir(2)
        izd = icir(3)

        ixabcd = ixabc+ixd
        iyabcd = iyabc+iyd

        ! Compute all desired gradients with respect to an x-component.

        mVec = 0
        nVec = 0
        do iCent=1,4
          if (IfGrad(1,iCent)) then
            mVec = mVec+1
            nVec = nVec+1
            Ind1(nVec) = 3*(Indx(1,iCent)-1)+1
            Ind2(nVec) = mVec
            Indx2(1,iCent) = mVec
          end if
        end do

        if (iyabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
            Ind1(nVec) = 3*(Indx(2,iCent)-1)+2
            Ind2(nVec) = mVec
            Indx2(2,iCent) = mVec
          end if
        end do

        if (ixabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
            Ind1(nVec) = 3*(Indx(3,iCent)-1)+3
            Ind2(nVec) = mVec
            Indx2(3,iCent) = mVec
          end if
        end do

        if (ixabcd*iyabcd /= 0) then

          select case (nVec)
            case (3)
              do iT=1,nT
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
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
                tmp1 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                tmp3 = Zero
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
                tmp1 = Zero
                tmp2 = Zero
                do iRys=1,nRys
                  tmp1 = tmp1+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(1))
                  tmp2 = tmp2+xyz2D1(iRys,iT,iza,izb,izc,izd,Ind1(2))
                end do
                g1(iT,ipa,ipb,ipc,ipd,Ind2(1)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(1))+tmp1
                g1(iT,ipa,ipb,ipc,ipd,Ind2(2)) = g1(iT,ipa,ipb,ipc,ipd,Ind2(2))+tmp2
              end do
            case (1)
              do iT=1,nT
                tmp1 = Zero
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
