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

!#define _DEBUGPRINT_
subroutine CmbnKEr(Rnr,qC,Di,nZeta,la,lb,Zeta,rFinal,nComp,Alpha,nAlpha,Beta,nBeta)
!***********************************************************************
!     Author: Kurt Pfingst                                             *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use rmat, only: bParm, Dipol, Dipol1, EpsQ, GammaPh, GammaTh, lCosf, lCost, lSInf, lSInt, QCoul, RMatR
use Constants, only: Two, Three, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp, nAlpha, nBeta
real(kind=wp), intent(in) :: Rnr(nZeta,0:la+lb+2), qC(nZeta,0:la+lb), Di(nZeta,-1:la+lb-1), Zeta(nZeta), Alpha(nAlpha), Beta(nBeta)
real(kind=wp), intent(out) :: rFinal(nZeta,nComp,nTri_Elem1(la),nTri_Elem1(lb))
integer(kind=iwp) :: ialpha, ibeta, iComp, ipa, ipb, ixa, ixb, iya, iyb, iza, izb, iZeta, k, kc, lrs, m, n, na, nb
real(kind=wp) :: b1, b1a, b2, b2a, b3, BBLoch, CConst1, CConst2, CConst3, ck1, const1, const2, const3, Fact, Fact1, Fact2, Fact3, &
                 ralpha, rbeta, rx1, ry1, rz1, W
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ia, ib
character(len=80) :: Label
#endif

iComp = 1
do ixa=0,la
  do ixb=0,lb
    rx1 = real(ixb*(ixb-1),kind=wp)
    n = ixa+ixb
    do iya=0,la-ixa
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)
      do iyb=0,lb-ixb
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)
        ry1 = real(iyb*(iyb-1),kind=wp)
        m = iya+iyb
        rz1 = real(izb*(izb-1),kind=wp)
        k = iza+izb

        ! Combine integrals
        ! define various factors
        !***************************************************
        ck1 = Two*real(ixb+iyb+izb,kind=wp)+Three
        !***************************************************
        !***************************************************
        const1 = rx1*gammath(n+m-2,k)*gammaph(m,n-2)
        !***************************************************
        const2 = ry1*gammath(n+m-2,k)*gammaph(m-2,n)
        !***************************************************
        const3 = rz1*gammath(n+m,k-2)*gammaph(m,n)
        !***************************************************
        !***************************************************
        CConst1 = const1+const2+const3
        !***************************************************
        CConst2 = ck1*gammath(n+m,k)*gammaph(m,n)
        !***************************************************
        CConst3 = gammath(n+m,k)*gammaph(m,n)
        !***************************************************
        ! Constants for Bloch term b1/b2/b3
        na = ixa+iya+iza
        nb = ixb+iyb+izb
        b1 = Half*(real(nb+1,kind=wp))*rmatr**(na+nb+1)
        b1a = Half*(real(na+1,kind=wp))*rmatr**(na+nb+1)
        W = gammath(n+m,k)*gammaph(m,n)

        ibeta = 1
        ialpha = 1
        kc = 1
        do iZeta=1,nZeta
          ralpha = Alpha(ialpha)
          rbeta = Beta(ibeta)
          b2 = rbeta*rmatr**(na+nb+3)
          b2a = ralpha*rmatr**(na+nb+3)
          b3 = exp(-Zeta(iZeta)*rmatr**2)
          BBLoch = W*b3*((b1-b2)-bParm*(b1-b2)*(b1a-b2a))
          rFinal(iZeta,iComp,ipa,ipb) = BBloch-(Half*CConst1*Rnr(iZeta,n+m+k-2)-rbeta*CConst2*Rnr(iZeta,n+m+k)+ &
                                                Two*rbeta**2*CConst3*Rnr(iZeta,n+m+k+2))
          if (iZeta == kc*nAlpha) then
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

#ifdef _DEBUGPRINT_
write(u6,*) ' Result in Cmbnker1'
do ia=1,(la+1)*(la+2)/2
  do ib=1,(lb+1)*(lb+2)/2
    write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
    call RecPrt(Label,' ',rFinal(:,:,ia,ib),nZeta,nComp)
  end do
end do
#endif

!***********************************************************************

! Add Coulomb contributions for photoionization calculations

if (abs(qCoul) > Epsq) then
  do ixa=0,la
    do ixb=0,lb
      do iya=0,la-ixa
        iza = la-ixa-iya
        ipa = C_Ind(la,ixa,iza)
        do iyb=0,lb-ixb
          izb = lb-ixb-iyb
          ipb = C_Ind(lb,ixb,izb)
          lrs = ixa+ixb+iya+iyb+iza+izb
          lcost = iza+izb
          lsint = ixa+ixb+iya+iyb
          lsinf = iya+iyb
          lcosf = ixa+ixb
          Fact = gammath(lsint,lcost)*gammaph(lsinf,lcosf)
          rFinal(:,iComp,ipa,ipb) = rFinal(:,iComp,ipa,ipb)+Fact*qCoul*qC(:,lrs)

        end do
      end do
    end do
  end do
end if

!***********************************************************************

#ifdef _DEBUGPRINT_
write(u6,*) ' Result in Cmbnker2'
do ia=1,(la+1)*(la+2)/2
  do ib=1,(lb+1)*(lb+2)/2
    write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
    call RecPrt(Label,' ',rFinal(:,:,ia,ib),nZeta,nComp)
  end do
end do
#endif

!***********************************************************************

! Add DIPOL contributions for photoionization calculations

if (abs(Dipol1) > Epsq) then
  do ixa=0,la
    do ixb=0,lb
      do iya=0,la-ixa
        iza = la-ixa-iya
        ipa = C_Ind(la,ixa,iza)
        do iyb=0,lb-ixb
          izb = lb-ixb-iyb
          ipb = C_Ind(lb,ixb,izb)
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
          rFinal(:,iComp,ipa,ipb) = rFinal(:,iComp,ipa,ipb)+(Fact1+Fact2+Fact3)*Di(:,lrs)

        end do
      end do
    end do
  end do
end if

!***********************************************************************

return

end subroutine CmbnKEr
