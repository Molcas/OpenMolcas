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
! Copyright (C) Kurt Pfingst                                           *
!***********************************************************************

subroutine CmbnKEr(Rnr,qC,Di,nZeta,la,lb,Zeta,final,nComp,Alpha,nAlpha,Beta,nBeta)
!***********************************************************************
!     Author: Kurt Pfingst                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "nrmf.fh"
#include "rmat.fh"
#include "gam.fh"
real*8 final(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2), Zeta(nZeta), Rnr(nZeta,0:la+lb+2), qC(nZeta,0:la+lb), &
       Di(nZeta,-1:la+lb-1), Alpha(nAlpha), Beta(nBeta)
character*80 Label
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

iRout = 134
iPrint = nPrint(iRout)

iComp = 1
do ixa=0,la
  do ixb=0,lb
    rx1 = dble(ixb*(ixb-1))
    n = ixa+ixb
    do iya=0,la-ixa
      iza = la-ixa-iya
      ipa = Ind(la,ixa,iza)
      do iyb=0,lb-ixb
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)
        ry1 = dble(iyb*(iyb-1))
        m = iya+iyb
        rz1 = dble(izb*(izb-1))
        k = iza+izb

        ! Combine integrals
        ! define various factors
        !***************************************************
        ck1 = 2d0*dble(ixb+iyb+izb)+3d0
        !***************************************************
        !***************************************************
        const1 = rx1*gammath(n+m-2,k)*gammaph(m,n-2)
        !***************************************************
        const2 = ry1*gammath(n+m-2,k)*gammaph(m-2,n)
        !***************************************************
        const3 = rz1*gammath(n+m,k-2)*gammaph(m,n)
        !***************************************************
        CConst1 = const1+const2+const3
        !***************************************************
        !***************************************************
        CConst2 = ck1*gammath(n+m,k)*gammaph(m,n)
        !***************************************************
        CConst3 = gammath(n+m,k)*gammaph(m,n)
        !***************************************************
        ! Constants for Bloch term b1/b2/b3
        na = ixa+iya+iza
        nb = ixb+iyb+izb
        b1 = 0.5d0*(dble(nb)+1.d0)*rmatr**(na+nb+1)
        b1a = 0.5d0*(dble(na)+1.d0)*rmatr**(na+nb+1)
        W = gammath(n+m,k)*gammaph(m,n)

        ibeta = 1
        ialpha = 1
        kc = 1
        do iZeta=1,nZeta
          ralpha = alpha(ialpha)
          rbeta = beta(ibeta)
          b2 = rbeta*rmatr**(na+nb+3)
          b2a = ralpha*rmatr**(na+nb+3)
          b3 = exp(-Zeta(iZeta)*rmatr*rmatr)
          BBLoch = W*b3*((b1-b2)-bParm*(b1-b2)*(b1a-b2a))
          c0 = 0.5d0
          c1 = -rbeta
          c2 = 2d0*rbeta*rbeta
          final(iZeta,iComp,ipa,ipb) = BBloch-(c0*CConst1*Rnr(iZeta,n+m+k-2)+c1*CConst2*Rnr(iZeta,n+m+k)+ &
                                               c2*CConst3*Rnr(iZeta,n+m+k+2))
          if (iZeta == kc*nalpha) then
            ibeta = ibeta+1
            ialpha = 0
            kc = kc+1
          end if
          ialpha = ialpha+1
        end do

      end do
    end do
  end do
end do

!***********************************************************************

if (iPrint >= 99) then
  write(6,*) ' Result in Cmbnker1'
  do ia=1,(la+1)*(la+2)/2
    do ib=1,(lb+1)*(lb+2)/2
      write(Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
      call RecPrt(Label,' ',final(1,1,ia,ib),nZeta,nComp)
    end do
  end do
end if

!***********************************************************************

! Add Coulomb contributions for photoionization calculations

if (abs(qCoul) > Epsq) then
  do ixa=0,la
    do ixb=0,lb
      do iya=0,la-ixa
        iza = la-ixa-iya
        ipa = Ind(la,ixa,iza)
        do iyb=0,lb-ixb
          izb = lb-ixb-iyb
          ipb = Ind(lb,ixb,izb)
          lrs = ixa+ixb+iya+iyb+iza+izb
          lcost = iza+izb
          lsint = ixa+ixb+iya+iyb
          lsinf = iya+iyb
          lcosf = ixa+ixb
          Fact = gammath(lsint,lcost)*gammaph(lsinf,lcosf)
          do iZeta=1,nZeta
            final(iZeta,iComp,ipa,ipb) = final(iZeta,iComp,ipa,ipb)+Fact*qCoul*qC(iZeta,lrs)
          end do
!
        end do
      end do
    end do
  end do
end if

!***********************************************************************

if (iPrint >= 99) then
  write(6,*) ' Result in Cmbnker2'
  do ia=1,(la+1)*(la+2)/2
    do ib=1,(lb+1)*(lb+2)/2
      write(Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
      call RecPrt(Label,' ',final(1,1,ia,ib),nZeta,nComp)
    end do
  end do
end if

!***********************************************************************

! Add DIPOL contributions for photoionization calculations

if (abs(Dipol1) > Epsq) then
  do ixa=0,la
    do ixb=0,lb
      do iya=0,la-ixa
        iza = la-ixa-iya
        ipa = Ind(la,ixa,iza)
        do iyb=0,lb-ixb
          izb = lb-ixb-iyb
          ipb = Ind(lb,ixb,izb)
          lrs = ixa+ixb+iya+iyb+iza+izb
          ! Beitrag der x-Komponente
          lcost = iza+izb
          lsint = ixa+ixb+iya+iyb+1
          lsinf = iya+iyb
          lcosf = ixa+ixb+1
          Fact1 = Dipol(1)*gammath(lsint,lcost)*gammaph(lsinf,lcosf)
          ! Beitrag der y-Komponente
          lcost = iza+izb
          lsint = ixa+ixb+iya+iyb+1
          lsinf = iya+iyb+1
          lcosf = ixa+ixb
          Fact2 = Dipol(2)*gammath(lsint,lcost)*gammaph(lsinf,lcosf)
          ! Beitrag der z-Komponente
          lcost = iza+izb+1
          lsint = ixa+ixb+iya+iyb
          lsinf = iya+iyb
          lcosf = ixa+ixb
          Fact3 = Dipol(3)*gammath(lsint,lcost)*gammaph(lsinf,lcosf)
          ! Summe
          Fact = Fact1+Fact2+Fact3
          do iZeta=1,nZeta
            final(iZeta,iComp,ipa,ipb) = final(iZeta,iComp,ipa,ipb)+Fact*Di(iZeta,lrs)
          end do

        end do
      end do
    end do
  end do
end if

!***********************************************************************

return

end subroutine CmbnKEr
