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

implicit real*8(A-H,O-Z)
!#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "iavec.fh"
real*8 g2(78), PAO(nt,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,(lc+1)*(lc+2)/2,(ld+1)*(ld+2)/2), &
       xyz2D0(nRys,nT,0:la+2,0:lb+2,0:lc+2,0:ld+2,3), xyz2D1(nRys,nT,0:la,0:lb,0:lc,0:ld,9), xyz2D2(nRys,nT,0:la,0:lb,0:lc,0:ld,18)
logical IfHss(4,3,4,3)
integer Index2(2,6,3), Index1(3,3), ng(3), nh(3)
! Statement functions
nElem(i) = (i+1)*(i+2)/2
Ind(Icent,Icar,Jcent,jCar) = ((iCent-1)*3+iCar-1)*((iCent-1)*3+iCar)/2+(jCent-1)*3+jCar

!define _DEBUGPRINT_
!iRout = 248
!iPrint = nPrint(iRout)
call dcopy_(78,[Zero],0,g2,1)
#ifdef _DEBUGPRINT_
call RecPrt('Assg2: g2(0)',' ',g2,1,78)
#endif
ii = la*(la+1)*(la+2)/6
jj = lb*(lb+1)*(lb+2)/6
kk = lc*(lc+1)*(lc+2)/6
ll = ld*(ld+1)*(ld+2)/6
kcar = 0 ! dummy initialize

! First we construct the non diagonal derivatives

do iCar=1,3
  do jCar=1,3

    ! Determine the permutation of the cartesian indexes

    if (jCar == iCar) cycle
    if (iCar*jCar == 2) Kcar = 3
    if (iCar*jCar == 6) kCar = 1
    if (iCar*jCar == 3) KCar = 2

    ! Loop over the atomic centre, in order
    ! that the intgrals are calculated in.

    do iDer=1,ng(iCar)
      do jDer=1,ng(jCar)
        iCent = Index1(iDer,iCar)
        jCent = Index1(jDer,jCar)

        if (IfHss(iCent,iCar,jCent,jCar)) then
          I = Ind(iCent,iCar,jCent,jCar)
          ix1 = (iDer-1)*3+iCar
          ix2 = (jDer-1)*3+jCar

          ! Loop over angular momentas

          do ipd=1,nelem(ld)
            id1 = ixyz(iCar,ll+ipd)
            id2 = ixyz(jCar,ll+ipd)
            id3 = ixyz(kCar,ll+ipd)
            do ipc=1,nelem(lc)
              ic1 = ixyz(iCar,kk+ipc)
              ic2 = ixyz(jCar,kk+ipc)
              ic3 = ixyz(kCar,kk+ipc)
              do ipb=1,nelem(lb)
                ib1 = ixyz(iCar,jj+ipb)
                ib2 = ixyz(jCar,jj+ipb)
                ib3 = ixyz(kCar,jj+ipb)
                do ipa=1,nelem(la)
                  ia1 = ixyz(iCar,ii+ipa)
                  ia2 = ixyz(jCar,ii+ipa)
                  ia3 = ixyz(kCar,ii+ipa)

                  ! Loop over Rys-polynomia and exponents of the basis set!

                  tmp = 0.0d0
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
      I = Ind(iCent,iCar,jCent,iCar)

      do ipd=1,nelem(ld)
        id1 = ixyz(iCar,ll+ipd)
        id2 = ixyz(jCar,ll+ipd)
        id3 = ixyz(kCar,ll+ipd)

        do ipc=1,nelem(lc)
          ic1 = ixyz(iCar,kk+ipc)
          ic2 = ixyz(jCar,kk+ipc)
          ic3 = ixyz(kCar,kk+ipc)

          do ipb=1,nelem(lb)
            ib1 = ixyz(iCar,jj+ipb)
            ib2 = ixyz(jCar,jj+ipb)
            ib3 = ixyz(kCar,jj+ipb)
            do ipa=1,nelem(la)
              ia1 = ixyz(iCar,ii+ipa)
              ia2 = ixyz(jCar,ii+ipa)
              ia3 = ixyz(kCar,ii+ipa)

              tmp = 0.0d0
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
