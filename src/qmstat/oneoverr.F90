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

subroutine OneOverR(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyz)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Two, Three, Five
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iQ_Atoms, iFil(nTri_Elem(iQ_Atoms),10), iCNum
real(kind=wp) :: Ax, Ay, Az, BoMaH(iQ_Atoms), BoMaO(iQ_Atoms), EEDisp, Eint(nTri_Elem(iQ_Atoms),10), outxyz(MxQCen,3)
integer(kind=iwp) :: i, ijhr, ip, j, jjhr, k
real(kind=wp) :: Gx, Gy, Gz, R2(5), Rab3i(5), Rab(3,5), Rg(5), Se(5), U(3,5)

EEdisp = Zero
!----------------------------------------------------------------------*
! Compute some distances and inverted distances etc. The potential,    *
! the field and etc. and when we already have the numbers, we also do  *
! the dispersion interaction.                                          *
!----------------------------------------------------------------------*
do k=1,nTri_Elem(iQ_Atoms)
  Gx = outxyz(k,1)+Ax
  Gy = outxyz(k,2)+Ay
  Gz = outxyz(k,3)+Az
  do j=iCnum+1,nPart
    i = 1+(j-1)*nCent
    ip = 1+(j-1)*nPol
    !Below follow a lot of distances to and fro.
    Rab(1,1) = Cordst(i,1)-Gx
    Rab(2,1) = Cordst(i,2)-Gy
    Rab(3,1) = Cordst(i,3)-Gz
    Rab(1,2) = Cordst(i+1,1)-Gx
    Rab(2,2) = Cordst(i+1,2)-Gy
    Rab(3,2) = Cordst(i+1,3)-Gz
    Rab(1,3) = Cordst(i+2,1)-Gx
    Rab(2,3) = Cordst(i+2,2)-Gy
    Rab(3,3) = Cordst(i+2,3)-Gz
    Rab(1,4) = Cordst(i+3,1)-Gx
    Rab(2,4) = Cordst(i+3,2)-Gy
    Rab(3,4) = Cordst(i+3,3)-Gz
    Rab(1,5) = Cordst(i+4,1)-Gx
    Rab(2,5) = Cordst(i+4,2)-Gy
    Rab(3,5) = Cordst(i+4,3)-Gz
    R2(1) = Rab(1,1)**2+Rab(2,1)**2+Rab(3,1)**2
    R2(2) = Rab(1,2)**2+Rab(2,2)**2+Rab(3,2)**2
    R2(3) = Rab(1,3)**2+Rab(2,3)**2+Rab(3,3)**2
    R2(4) = Rab(1,4)**2+Rab(2,4)**2+Rab(3,4)**2
    R2(5) = Rab(1,5)**2+Rab(2,5)**2+Rab(3,5)**2
    Rg(1) = sqrt(R2(1))
    Rg(2) = sqrt(R2(2))
    Rg(3) = sqrt(R2(3))
    Rg(4) = sqrt(R2(4))
    Rg(5) = sqrt(R2(5))
    Se(1) = One/Rg(1)
    Se(2) = One/Rg(2)
    Se(3) = One/Rg(3)
    Se(4) = One/Rg(4)
    Se(5) = One/Rg(5)
    ! This term will below turn into the interaction between
    ! charges on water and the MME-charges on the QM-molecule.
    Eint(k,1) = Eint(k,1)-Qsta(1)*Se(2)-Qsta(2)*Se(3)-Qsta(3)*Se(4)-Qsta(4)*Se(5)
    Rab3i(1) = Se(1)/R2(1)
    Rab3i(2) = Se(2)/R2(2)
    Rab3i(3) = Se(3)/R2(3)
    Rab3i(4) = Se(4)/R2(4)
    Rab3i(5) = Se(5)/R2(5)
    !------------------------------------------------------------------*
    ! The dispersion interaction between QM-atoms and solvent is       *
    ! computed, with or without damping. The initial if-clause sees to *
    ! that only atom-centers are included, while bonds and virtual     *
    ! centers are ignored.                                             *
    !------------------------------------------------------------------*
    if (k <= iQ_atoms) then
      call DispEnergy(EEDisp,BoMah(k),BoMaO(k),Rg(1),Rg(2),Rg(3),Rab3i(1),Rab3i(2),Rab3i(3),k)
    end if
    !------------------------------------------------------------------*
    ! Now we wrap up the electrostatics.                               *
    !------------------------------------------------------------------*
    U(1,1) = Rab(1,1)*Se(1)
    U(2,1) = Rab(2,1)*Se(1)
    U(3,1) = Rab(3,1)*Se(1)
    U(1,2) = Rab(1,2)*Se(2)
    U(2,2) = Rab(2,2)*Se(2)
    U(3,2) = Rab(3,2)*Se(2)
    U(1,3) = Rab(1,3)*Se(3)
    U(2,3) = Rab(2,3)*Se(3)
    U(3,3) = Rab(3,3)*Se(3)
    U(1,4) = Rab(1,4)*Se(4)
    U(2,4) = Rab(2,4)*Se(4)
    U(3,4) = Rab(3,4)*Se(4)
    U(1,5) = Rab(1,5)*Se(5)
    U(2,5) = Rab(2,5)*Se(5)
    U(3,5) = Rab(3,5)*Se(5)
    ! These three terms will below turn into the interaction
    ! between water charges and the MME-dipoles on the QM-mol.
    ! Change sign of charge, change sign of vector and then we
    ! should also change sign since when a dipole interacts with
    ! a field we have a minus sign, but this minus sign we have
    ! omitted in hel; therefore this calculation gives the right
    ! number eventually.
    Eint(k,2) = -Qsta(1)*Rab(1,2)*Rab3i(2)-Qsta(2)*Rab(1,3)*Rab3i(3)-Qsta(3)*Rab(1,4)*Rab3i(4)-Qsta(4)*Rab(1,5)*Rab3i(5)+Eint(k,2)
    Eint(k,3) = -Qsta(1)*Rab(2,2)*Rab3i(2)-Qsta(2)*Rab(2,3)*Rab3i(3)-Qsta(3)*Rab(2,4)*Rab3i(4)-Qsta(4)*Rab(2,5)*Rab3i(5)+Eint(k,3)
    Eint(k,4) = -Qsta(1)*Rab(3,2)*Rab3i(2)-Qsta(2)*Rab(3,3)*Rab3i(3)-Qsta(3)*Rab(3,4)*Rab3i(4)-Qsta(4)*Rab(3,5)*Rab3i(5)+Eint(k,4)
    ! And here it is the MME-quadrupoles that are prepared.
    ! Change sign of charges, change sign two times of the
    ! vector (in effect, zero times then) and then in the
    ! energy expression for the interaction between the field
    ! vector from a charge and a quarupole there is a plus
    ! sign, so a minus is the right sign below.
    Eint(k,5) = Eint(k,5)-Qsta(1)*U(1,2)**2*Rab3i(2)-Qsta(2)*U(1,3)**2*Rab3i(3)-Qsta(3)*U(1,4)**2*Rab3i(4)- &
                Qsta(4)*U(1,5)**2*Rab3i(5)
    Eint(k,7) = Eint(k,7)-Qsta(1)*U(2,2)**2*Rab3i(2)-Qsta(2)*U(2,3)**2*Rab3i(3)-Qsta(3)*U(2,4)**2*Rab3i(4)- &
                Qsta(4)*U(2,5)**2*Rab3i(5)
    Eint(k,10) = Eint(k,10)-Qsta(1)*U(3,2)**2*Rab3i(2)-Qsta(2)*U(3,3)**2*Rab3i(3)-Qsta(3)*U(3,4)**2*Rab3i(4)- &
                 Qsta(4)*U(3,5)**2*Rab3i(5)
    Eint(k,6) = Eint(k,6)-Qsta(1)*U(1,2)*U(2,2)*Rab3i(2)-Qsta(2)*U(1,3)*U(2,3)*Rab3i(3)-Qsta(3)*U(1,4)*Rab3i(4)*U(2,4)- &
                Qsta(4)*U(1,5)*Rab3i(5)*U(2,5)
    Eint(k,8) = Eint(k,8)-Qsta(1)*U(1,2)*U(3,2)*Rab3i(2)-Qsta(2)*U(1,3)*U(3,3)*Rab3i(3)-Qsta(3)*U(1,4)*Rab3i(4)*U(3,4)- &
                Qsta(4)*U(1,5)*Rab3i(5)*U(3,5)
    Eint(k,9) = Eint(k,9)-Qsta(1)*U(3,2)*U(2,2)*Rab3i(2)-Qsta(2)*U(3,3)*U(2,3)*Rab3i(3)-Qsta(3)*U(3,4)*Rab3i(4)*U(2,4)- &
                Qsta(4)*U(3,5)*Rab3i(5)*U(2,5)

    !------------------------------------------------------------------*
    ! And now a whole lot of grad(1/r) and higher...                   *
    !------------------------------------------------------------------*
    ! Monopoles.
    Work(iFil(k,1)-1+ip) = Rab(1,1)*Rab3i(1)
    Work(iFil(k,1)-1+ip+1) = Rab(1,2)*Rab3i(2)
    Work(iFil(k,1)-1+ip+2) = Rab(1,3)*Rab3i(3)
    Work(iFil(k,1)-1+nPart*nPol+ip) = Rab(2,1)*Rab3i(1)
    Work(iFil(k,1)-1+nPart*nPol+ip+1) = Rab(2,2)*Rab3i(2)
    Work(iFil(k,1)-1+nPart*nPol+ip+2) = Rab(2,3)*Rab3i(3)
    Work(iFil(k,1)-1+2*nPart*nPol+ip) = Rab(3,1)*Rab3i(1)
    Work(iFil(k,1)-1+2*nPart*nPol+ip+1) = Rab(3,2)*Rab3i(2)
    Work(iFil(k,1)-1+2*nPart*nPol+ip+2) = Rab(3,3)*Rab3i(3)
    ! Dipole -- x-component.
    Work(iFil(k,2)-1+ip) = (Three*U(1,1)**2-One)*Rab3i(1)
    Work(iFil(k,2)-1+ip+1) = (Three*U(1,2)**2-One)*Rab3i(2)
    Work(iFil(k,2)-1+ip+2) = (Three*U(1,3)**2-One)*Rab3i(3)
    Work(iFil(k,2)-1+nPart*nPol+ip) = U(2,1)*U(1,1)*Rab3i(1)*Three
    Work(iFil(k,2)-1+nPart*nPol+ip+1) = U(2,2)*U(1,2)*Rab3i(2)*Three
    Work(iFil(k,2)-1+nPart*nPol+ip+2) = U(2,3)*U(1,3)*Rab3i(3)*Three
    Work(iFil(k,2)-1+2*nPart*nPol+ip) = U(3,1)*U(1,1)*Rab3i(1)*Three
    Work(iFil(k,2)-1+2*nPart*nPol+ip+1) = U(3,2)*U(1,2)*Rab3i(2)*Three
    Work(iFil(k,2)-1+2*nPart*nPol+ip+2) = U(3,3)*U(1,3)*Rab3i(3)*Three
    ! Dipole -- y-component.
    Work(iFil(k,3)-1+ip) = U(2,1)*U(1,1)*Rab3i(1)*Three
    Work(iFil(k,3)-1+ip+1) = U(2,2)*U(1,2)*Rab3i(2)*Three
    Work(iFil(k,3)-1+ip+2) = U(2,3)*U(1,3)*Rab3i(3)*Three
    Work(iFil(k,3)-1+nPart*nPol+ip) = (Three*U(2,1)**2-One)*Rab3i(1)
    Work(iFil(k,3)-1+nPart*nPol+ip+1) = (Three*U(2,2)**2-One)*Rab3i(2)
    Work(iFil(k,3)-1+nPart*nPol+ip+2) = (Three*U(2,3)**2-One)*Rab3i(3)
    Work(iFil(k,3)-1+2*nPart*nPol+ip) = U(3,1)*U(2,1)*Rab3i(1)*Three
    Work(iFil(k,3)-1+2*nPart*nPol+ip+1) = U(3,2)*U(2,2)*Rab3i(2)*Three
    Work(iFil(k,3)-1+2*nPart*nPol+ip+2) = U(3,3)*U(2,3)*Rab3i(3)*Three
    ! Dipole -- z-component.
    Work(iFil(k,4)-1+ip) = U(3,1)*U(1,1)*Rab3i(1)*Three
    Work(iFil(k,4)-1+ip+1) = U(3,2)*U(1,2)*Rab3i(2)*Three
    Work(iFil(k,4)-1+ip+2) = U(3,3)*U(1,3)*Rab3i(3)*Three
    Work(iFil(k,4)-1+nPart*nPol+ip) = U(3,1)*U(2,1)*Rab3i(1)*Three
    Work(iFil(k,4)-1+nPart*nPol+ip+1) = U(3,2)*U(2,2)*Rab3i(2)*Three
    Work(iFil(k,4)-1+nPart*nPol+ip+2) = U(3,3)*U(2,3)*Rab3i(3)*Three
    Work(iFil(k,4)-1+2*nPart*nPol+ip) = (Three*U(3,1)**2-One)*Rab3i(1)
    Work(iFil(k,4)-1+2*nPart*nPol+ip+1) = (Three*U(3,2)**2-One)*Rab3i(2)
    Work(iFil(k,4)-1+2*nPart*nPol+ip+2) = (Three*U(3,3)**2-One)*Rab3i(3)
    ! Quadrupole -- xx-component.
    Work(iFil(k,5)-1+ip) = U(1,1)*(Five*U(1,1)**2-Two)*Rab3i(1)*Se(1)
    Work(iFil(k,5)-1+ip+1) = U(1,2)*(Five*U(1,2)**2-Two)*Rab3i(2)*Se(2)
    Work(iFil(k,5)-1+ip+2) = U(1,3)*(Five*U(1,3)**2-Two)*Rab3i(3)*Se(3)
    Work(iFil(k,5)-1+nPart*nPol+ip) = Five*U(2,1)*U(1,1)**2*Rab3i(1)*Se(1)
    Work(iFil(k,5)-1+nPart*nPol+ip+1) = Five*U(2,2)*U(1,2)**2*Rab3i(2)*Se(2)
    Work(iFil(k,5)-1+nPart*nPol+ip+2) = Five*U(2,3)*U(1,3)**2*Rab3i(3)*Se(3)
    Work(iFil(k,5)-1+2*nPart*nPol+ip) = Five*U(3,1)*U(1,1)**2*Rab3i(1)*Se(1)
    Work(iFil(k,5)-1+2*nPart*nPol+ip+1) = Five*U(3,2)*U(1,2)**2*Rab3i(2)*Se(2)
    Work(iFil(k,5)-1+2*nPart*nPol+ip+2) = Five*U(3,3)*U(1,3)**2*Rab3i(3)*Se(3)
    ! Quadrupole -- yy-component.
    Work(iFil(k,7)-1+ip) = Five*U(2,1)**2*U(1,1)*Rab3i(1)*Se(1)
    Work(iFil(k,7)-1+ip+1) = Five*U(2,2)**2*U(1,2)*Rab3i(2)*Se(2)
    Work(iFil(k,7)-1+ip+2) = Five*U(2,3)**2*U(1,3)*Rab3i(3)*Se(3)
    Work(iFil(k,7)-1+nPart*nPol+ip) = U(2,1)*(Five*U(2,1)**2-Two)*Rab3i(1)*Se(1)
    Work(iFil(k,7)-1+nPart*nPol+ip+1) = U(2,2)*(Five*U(2,2)**2-Two)*Rab3i(2)*Se(2)
    Work(iFil(k,7)-1+nPart*nPol+ip+2) = U(2,3)*(Five*U(2,3)**2-Two)*Rab3i(3)*Se(3)
    Work(iFil(k,7)-1+2*nPart*nPol+ip) = Five*U(3,1)*U(2,1)**2*Rab3i(1)*Se(1)
    Work(iFil(k,7)-1+2*nPart*nPol+ip+1) = Five*U(3,2)*U(2,2)**2*Rab3i(2)*Se(2)
    Work(iFil(k,7)-1+2*nPart*nPol+ip+2) = Five*U(3,3)*U(2,3)**2*Rab3i(3)*Se(3)
    ! Quadrupole -- zz-component.
    Work(iFil(k,10)-1+ip) = Five*U(3,1)**2*U(1,1)*Rab3i(1)*Se(1)
    Work(iFil(k,10)-1+ip+1) = Five*U(3,2)**2*U(1,2)*Rab3i(2)*Se(2)
    Work(iFil(k,10)-1+ip+2) = Five*U(3,3)**2*U(1,3)*Rab3i(3)*Se(3)
    Work(iFil(k,10)-1+nPart*nPol+ip) = Five*U(3,1)**2*U(2,1)*Rab3i(1)*Se(1)
    Work(iFil(k,10)-1+nPart*nPol+ip+1) = Five*U(3,2)**2*U(2,2)*Rab3i(2)*Se(2)
    Work(iFil(k,10)-1+nPart*nPol+ip+2) = Five*U(3,3)**2*U(2,3)*Rab3i(3)*Se(3)
    Work(iFil(k,10)-1+2*nPart*nPol+ip) = U(3,1)*(Five*U(3,1)**2-Two)*Rab3i(1)*Se(1)
    Work(iFil(k,10)-1+2*nPart*nPol+ip+1) = U(3,2)*(Five*U(3,2)**2-Two)*Rab3i(2)*Se(2)
    Work(iFil(k,10)-1+2*nPart*nPol+ip+2) = U(3,3)*(Five*U(3,3)**2-Two)*Rab3i(3)*Se(3)
    ! Quadrupole -- xy-component.
    Work(iFil(k,6)-1+ip) = U(2,1)*(Five*U(1,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,6)-1+ip+1) = U(2,2)*(Five*U(1,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,6)-1+ip+2) = U(2,3)*(Five*U(1,3)**2-One)*Rab3i(3)*Se(3)
    Work(iFil(k,6)-1+nPart*nPol+ip) = U(1,1)*(Five*U(2,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,6)-1+nPart*nPol+ip+1) = U(1,2)*(Five*U(2,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,6)-1+nPart*nPol+ip+2) = U(1,3)*(Five*U(2,3)**2-One)*Rab3i(3)*Se(3)
    Work(iFil(k,6)-1+2*nPart*nPol+ip) = Five*U(3,1)*U(2,1)*U(1,1)*Rab3i(1)*Se(1)
    Work(iFil(k,6)-1+2*nPart*nPol+ip+1) = Five*U(3,2)*U(2,2)*U(1,2)*Rab3i(2)*Se(2)
    Work(iFil(k,6)-1+2*nPart*nPol+ip+2) = Five*U(3,3)*U(2,3)*U(1,3)*Rab3i(3)*Se(3)
    ! Quadrupole -- xz-component.
    Work(iFil(k,8)-1+ip) = U(3,1)*(Five*U(1,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,8)-1+ip+1) = U(3,2)*(Five*U(1,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,8)-1+ip+2) = U(3,3)*(Five*U(1,3)**2-One)*Rab3i(3)*Se(3)
    Work(iFil(k,8)-1+nPart*nPol+ip) = Five*U(3,1)*U(2,1)*U(1,1)*Rab3i(1)*Se(1)
    Work(iFil(k,8)-1+nPart*nPol+ip+1) = Five*U(3,2)*U(2,2)*U(1,2)*Rab3i(2)*Se(2)
    Work(iFil(k,8)-1+nPart*nPol+ip+2) = Five*U(3,3)*U(2,3)*U(1,3)*Rab3i(3)*Se(3)
    Work(iFil(k,8)-1+2*nPart*nPol+ip) = U(1,1)*(Five*U(3,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,8)-1+2*nPart*nPol+ip+1) = U(1,2)*(Five*U(3,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,8)-1+2*nPart*nPol+ip+2) = U(1,3)*(Five*U(3,3)**2-One)*Rab3i(3)*Se(3)
    ! Quadrupole -- yz-component.
    Work(iFil(k,9)-1+ip) = Five*U(3,1)*U(2,1)*U(1,1)*Rab3i(1)*Se(1)
    Work(iFil(k,9)-1+ip+1) = Five*U(3,2)*U(2,2)*U(1,2)*Rab3i(2)*Se(2)
    Work(iFil(k,9)-1+ip+2) = Five*U(3,3)*U(2,3)*U(1,3)*Rab3i(3)*Se(3)
    Work(iFil(k,9)-1+nPart*nPol+ip) = U(3,1)*(Five*U(2,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,9)-1+nPart*nPol+ip+1) = U(3,2)*(Five*U(2,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,9)-1+nPart*nPol+ip+2) = U(3,3)*(Five*U(2,3)**2-One)*Rab3i(3)*Se(3)
    Work(iFil(k,9)-1+2*nPart*nPol+ip) = U(2,1)*(Five*U(3,1)**2-One)*Rab3i(1)*Se(1)
    Work(iFil(k,9)-1+2*nPart*nPol+ip+1) = U(2,2)*(Five*U(3,2)**2-One)*Rab3i(2)*Se(2)
    Work(iFil(k,9)-1+2*nPart*nPol+ip+2) = U(2,3)*(Five*U(3,3)**2-One)*Rab3i(3)*Se(3)
    !------------------------------------------------------------------*
    ! If damping of the field is requested, then do it.                *
    !------------------------------------------------------------------*
    if (FieldDamp) then
      do ijhr=1,10
        do jjhr=0,2
          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip) = Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip)*(One-exp(CAFieldG*Rg(1)))**CFexp
          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+1) = Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+1)*(One-exp(CBFieldG*Rg(2)))**CFexp
          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+2) = Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+2)*(One-exp(CBFieldG*Rg(3)))**CFexp
        end do
      end do
    end if
  end do
end do

return

end subroutine OneOverR
