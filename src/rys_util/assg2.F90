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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Assg2(g2,nT,nRys,la,lb,lc,ld,xyz2D0,xyz2D1,xyz2D2,IfHss,Index1,Index2,ng,nh,PAO)
!***********************************************************************
!                                                                      *
! Object: to assemble the gradients of the ERI's.                      *
!                                                                      *
!     Author: Anders Bernhardsson                                      *
!             Dept. of Theoretical Chemistry,                          *
!             University of Lund, SWEDEN                               *
!             March '95                                                *
!***********************************************************************

use Index_Functions, only: C_Ind3_Rev, nTri_Elem, nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: g2(78)
integer(kind=iwp), intent(in) :: nT, nRys, la, lb, lc, ld, Index1(3,3), Index2(2,6,3), ng(3), nh(3)
real(kind=wp), intent(in) :: xyz2D0(nRys,nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9), &
                             xyz2D2(nRys,nT,0:la,0:lb,0:lc,0:ld,18), &
                             PAO(nT,nTri_Elem1(la),nTri_Elem1(lb),nTri_Elem1(lc),nTri_Elem1(ld))
logical(kind=iwp), intent(in) :: IfHss(4,3,4,3)
#include "itmax.fh"
integer(kind=iwp) :: I, ia1, ia2, ia3, ib1, ib2, ib3, ic1, ic2, ic3, iCar, iCent, icir(3), id1, id2, id3, iDer, ipa, ipb, ipc, &
                     ipd, iRys, it, ix1, ix2, jCar, jCent, jDer, kCar
real(kind=wp) :: tmp

!define _DEBUGPRINT_
!iRout = 248
!iPrint = nPrint(iRout)
g2(:) = Zero
#ifdef _DEBUGPRINT_
call RecPrt('Assg2: g2(0)',' ',g2,1,78)
#endif
kCar = 0 ! dummy initialize

! First we construct the non diagonal derivatives

do iCar=1,3
  do jCar=1,3

    ! Determine the permutation of the cartesian indices

    if (jCar == iCar) cycle
    if (iCar*jCar == 2) Kcar = 3
    if (iCar*jCar == 6) kCar = 1
    if (iCar*jCar == 3) KCar = 2

    ! Loop over the atomic centre, in order
    ! that the integrals are calculated in.

    do iDer=1,ng(iCar)
      do jDer=1,ng(jCar)
        iCent = Index1(iDer,iCar)
        jCent = Index1(jDer,jCar)

        if (IfHss(iCent,iCar,jCent,jCar)) then
          I = nTri_Elem((iCent-1)*3+iCar-1)+(jCent-1)*3+jCar
          ix1 = (iDer-1)*3+iCar
          ix2 = (jDer-1)*3+jCar

          ! Loop over angular momentas

          do ipd=1,nTri_Elem1(ld)
            icir(:) = C_Ind3_Rev(ipd,ld)
            id1 = icir(iCar)
            id2 = icir(jCar)
            id3 = icir(kCar)
            do ipc=1,nTri_Elem1(lc)
              icir(:) = C_Ind3_Rev(ipc,lc)
              ic1 = icir(iCar)
              ic2 = icir(jCar)
              ic3 = icir(kCar)
              do ipb=1,nTri_Elem1(lb)
                icir(:) = C_Ind3_Rev(ipb,lb)
                ib1 = icir(iCar)
                ib2 = icir(jCar)
                ib3 = icir(kCar)
                do ipa=1,nTri_Elem1(la)
                  icir(:) = C_Ind3_Rev(ipa,la)
                  ia1 = icir(iCar)
                  ia2 = icir(jCar)
                  ia3 = icir(kCar)

                  ! Loop over Rys-polynomia and exponents of the basis set!

                  tmp = Zero
                  do it=1,nt
                    do iRys=1,nRys
                      tmp = tmp+PAO(iT,ipa,ipb,ipc,ipd)*xyz2D0(iRys,iT,ia3,ib3,ic3,id3,kCar)*xyz2D1(iRys,iT,ia1,ib1,ic1,id1,ix1)* &
                            xyz2D1(iRys,iT,ia2,ib2,ic2,id2,ix2)
                    end do
                  end do
                  g2(I) = g2(I)+tmp
                end do
              end do
            end do
          end do
        end if
      end do
    end do
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Assg2: g2 non-diagonal',' ',g2,1,78)
#endif

! Then construct the diagonal derivatives

do iCar=1,3
  jCar = mod(iCar,3)+1
  kCar = mod(jCar,3)+1
  do iDer=1,nh(iCar)
    ix1 = (iDer-1)*3+iCar
    iCent = Index2(1,iDer,iCar)
    jCent = Index2(2,iDer,iCar)
    if (IfHss(iCent,iCar,jCent,iCar)) then
      I = nTri_Elem((iCent-1)*3+iCar-1)+(jCent-1)*3+iCar

      do ipd=1,nTri_Elem1(ld)
        icir(:) = C_Ind3_Rev(ipd,ld)
        id1 = icir(iCar)
        id2 = icir(jCar)
        id3 = icir(kCar)
        do ipc=1,nTri_Elem1(lc)
          icir(:) = C_Ind3_Rev(ipc,lc)
          ic1 = icir(iCar)
          ic2 = icir(jCar)
          ic3 = icir(kCar)
          do ipb=1,nTri_Elem1(lb)
            icir(:) = C_Ind3_Rev(ipb,lb)
            ib1 = icir(iCar)
            ib2 = icir(jCar)
            ib3 = icir(kCar)
            do ipa=1,nTri_Elem1(la)
              icir(:) = C_Ind3_Rev(ipa,la)
              ia1 = icir(iCar)
              ia2 = icir(jCar)
              ia3 = icir(kCar)

              tmp = Zero
              do it=1,nt
                do iRys=1,nRys
                  tmp = tmp+PAO(iT,ipa,ipb,ipc,ipd)*xyz2D0(iRys,iT,ia2,ib2,ic2,id2,jCar)*xyz2D0(iRys,iT,ia3,ib3,ic3,id3,kCar)* &
                        xyz2D2(iRys,iT,ia1,ib1,ic1,id1,ix1)
                end do
              end do
              g2(I) = g2(I)+tmp

            end do
          end do
        end do
      end do
    end if
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Assg2: g2 full',' ',g2,1,78)
#endif

return

end subroutine Assg2
