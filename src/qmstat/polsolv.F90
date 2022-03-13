!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine PolSolv(iDT,iFI,iFP,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,FFp,iCNum,r2inv,difac,nSize)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"
dimension xx(nSize,nSize), yy(nSize,nSize), zz(nSize,nSize)
dimension xxi(nSize,nSize), yyi(nSize,nSize), zzi(nSize,nSize)
dimension rr3(nSize,nSize), Gri(nSize,nSize)
dimension FFp(nSize,3)
dimension iDT(3), iFI(3), iFP(3)

! The fields between the solvent molecules are computed as well as
! the image charge contribution. Rather basic formulas, but keeping
! track on indices may be difficult.

do j=1,nPol
  IndCor = j+(iCnum-1)*Ncent
  Inddt = j+(iCnum-1)*nPol
  do i=iCnum+1,nPart
    IndCor = IndCor+nCent
    Inddt = Inddt+nPol
    Agr = Sqrs(IndCor) !Sqrs come from geogen.
    Skal = Work(iDT(1)+Inddt-1)*Cordst(IndCor,1)+Work(iDT(2)+Inddt-1)*Cordst(IndCor,2)+Work(iDT(3)+Inddt-1)*Cordst(IndCor,3)
    Ta = Skal*Agr**2*R2inv
    Tal = -Difac*Ta
    Qimp(Inddt) = tal*agr
    ! Image dipoles: Reflect dipole vector in radial vector. This is vector geometry at its best.
    dim(IndDt,1) = (Tal*Cordst(IndCor,1)*2+DiFac*Work(iDt(1)+IndDt-1))*Agr**3
    dim(IndDt,2) = (Tal*Cordst(IndCor,2)*2+DiFac*Work(iDt(2)+IndDt-1))*Agr**3
    dim(IndDt,3) = (Tal*Cordst(IndCor,3)*2+DiFac*Work(iDt(3)+IndDt-1))*Agr**3
  end do
end do
do j=1,3
  do i=1+(nPol*iCnum),nSize
    Work(iFi(j)+i-1) = 0
  end do
end do

! Here the actual fields are computed, both from the explicit solvent and from its image.

do i=1+(nPol*iCNum),nPart*nPol !The real part.
  do j=1+(nPol*iCnum),nPart*nPol
    idel1 = (i-1)/nPol
    idel2 = (j-1)/nPol
    if (idel1 == idel2) cycle
    Skal = xx(i,j)*Work(iDt(1)+i-1)+yy(i,j)*Work(iDt(2)+i-1)+zz(i,j)*Work(iDt(3)+i-1)
    Skal = Skal*3
    Work(iFi(1)+j-1) = Work(iFi(1)+j-1)-(Work(iDt(1)+i-1)-Skal*xx(i,j))*rr3(i,j)
    Work(iFi(2)+j-1) = Work(iFi(2)+j-1)-(Work(iDt(2)+i-1)-Skal*yy(i,j))*rr3(i,j)
    Work(iFi(3)+j-1) = Work(iFi(3)+j-1)-(Work(iDt(3)+i-1)-Skal*zz(i,j))*rr3(i,j)
  end do
end do
do i=1+(nPol*iCnum),nPart*nPol !The image part.
  do j=1+(nPol*iCnum),nPart*nPol
    Skal = (xxi(i,j)*dim(i,1)+yyi(i,j)*dim(i,2)+zzi(i,j)*dim(i,3))
    Skal = Skal*3
    Work(iFi(1)+j-1) = Work(iFi(1)+j-1)-(dim(i,1)-Skal*xxi(i,j))*Gri(i,j)**3-Qimp(i)*xxi(i,j)*Gri(i,j)**2
    Work(iFi(2)+j-1) = Work(iFi(2)+j-1)-(dim(i,2)-Skal*yyi(i,j))*Gri(i,j)**3-Qimp(i)*yyi(i,j)*Gri(i,j)**2
    Work(iFi(3)+j-1) = Work(iFi(3)+j-1)-(dim(i,3)-Skal*zzi(i,j))*Gri(i,j)**3-Qimp(i)*zzi(i,j)*Gri(i,j)**2
  end do
end do

! Add up the things, and add also the part from the permanent charges.

do i=1+(nPol*iCnum),nSize
  FFp(i,1) = Work(iFi(1)+i-1)+Work(iFp(1)+i-1)
  FFp(i,2) = Work(iFi(2)+i-1)+Work(iFp(2)+i-1)
  FFp(i,3) = Work(iFi(3)+i-1)+Work(iFp(3)+i-1)
end do

! The circle is now complete. I left you when I was about to learn and I meet you as the master.

return

end subroutine PolSolv
