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

subroutine PolSolv(DT,FI,FP,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,FFp,iCNum,r2inv,difac,nSize)

use qmstat_global, only: Cordst, DipIm, nCent, nPart, nPol, Sqrs, Qimp
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iCNum, nSize
real(kind=wp), intent(in) :: DT(3,nPol*nPart), FP(3,nPol*nPart), xx(nSize,nSize), yy(nSize,nSize), zz(nSize,nSize), &
                             rr3(nSize,nSize), xxi(nSize,nSize), yyi(nSize,nSize), zzi(nSize,nSize), Gri(nSize,nSize), r2inv, difac
real(kind=wp), intent(out) :: FI(3,nPol*nPart), FFp(nSize,3)
integer(kind=iwp) :: i, idel1, idel2, IndCor, Inddt, j
real(kind=wp) :: Agr, Skal, Ta, Tal

! The fields between the solvent molecules are computed as well as
! the image charge contribution. Rather basic formulas, but keeping
! track on indices may be difficult.

do j=1,nPol
  IndCor = j+(iCnum-1)*nCent
  Inddt = j+(iCnum-1)*nPol
  do i=iCnum+1,nPart
    IndCor = IndCor+nCent
    Inddt = Inddt+nPol
    Agr = Sqrs(IndCor) !Sqrs come from geogen.
    Skal = DT(1,Inddt)*Cordst(1,IndCor)+DT(2,Inddt)*Cordst(2,IndCor)+DT(3,Inddt)*Cordst(3,IndCor)
    Ta = Skal*Agr**2*R2inv
    Tal = -Difac*Ta
    Qimp(Inddt) = tal*agr
    ! Image dipoles: Reflect dipole vector in radial vector. This is vector geometry at its best.
    DipIm(:,IndDt) = (Tal*Cordst(:,IndCor)*2+DiFac*DT(:,IndDt))*Agr**3
  end do
end do
FI(:,nPol*iCnum+1:nSize) = Zero

! Here the actual fields are computed, both from the explicit solvent and from its image.

do i=1+(nPol*iCNum),nPart*nPol !The real part.
  do j=1+(nPol*iCnum),nPart*nPol
    idel1 = (i-1)/nPol
    idel2 = (j-1)/nPol
    if (idel1 == idel2) cycle
    Skal = xx(i,j)*DT(1,i)+yy(i,j)*DT(2,i)+zz(i,j)*DT(3,i)
    Skal = Skal*3
    FI(1,j) = FI(1,j)-(DT(1,i)-Skal*xx(i,j))*rr3(i,j)
    FI(2,j) = FI(2,j)-(DT(2,i)-Skal*yy(i,j))*rr3(i,j)
    FI(3,j) = FI(3,j)-(DT(3,i)-Skal*zz(i,j))*rr3(i,j)
  end do
end do
do i=1+(nPol*iCnum),nPart*nPol !The image part.
  do j=1+(nPol*iCnum),nPart*nPol
    Skal = (xxi(i,j)*DipIm(1,i)+yyi(i,j)*DipIm(2,i)+zzi(i,j)*DipIm(3,i))
    Skal = Skal*3
    FI(1,j) = FI(1,j)-(DipIm(1,i)-Skal*xxi(i,j))*Gri(i,j)**3-Qimp(i)*xxi(i,j)*Gri(i,j)**2
    FI(2,j) = FI(2,j)-(DipIm(2,i)-Skal*yyi(i,j))*Gri(i,j)**3-Qimp(i)*yyi(i,j)*Gri(i,j)**2
    FI(3,j) = FI(3,j)-(DipIm(3,i)-Skal*zzi(i,j))*Gri(i,j)**3-Qimp(i)*zzi(i,j)*Gri(i,j)**2
  end do
end do

! Add up the things, and add also the part from the permanent charges.

do i=1+(nPol*iCnum),nSize
  FFp(i,:) = Fi(:,i)+FP(:,i)
end do

! The circle is now complete. I left you when I was about to learn and I meet you as the master.

return

end subroutine PolSolv
