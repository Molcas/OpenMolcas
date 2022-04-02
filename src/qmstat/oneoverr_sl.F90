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
! Copyright (C) 2011, Jose Manuel Hermida Ramon                        *
!***********************************************************************

subroutine OneOverR_Sl(Fil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyz,Eint_Nuc)
!----------------------------------------------------------------------*
!     Jose. Modification of the OneOverR subroutine to include the     *
!     electrostatic penetration of the charge density in the           *
!     electrostatic operator of the Hamiltonian. 2011-05-30            *
!----------------------------------------------------------------------*

use qmstat_global, only: CAFieldG, CBFieldG, CFexp, Cordst, Cut_Elc, DifSlExp, FieldDamp, lMltSlC, lQuad, MxMltp, nCent, nMlt, &
                         nPart, nPol, nSlSiteC, Qsta, SlExpC, SlExpQ, SlFactC, SlPC
use Index_Functions, only: nTri3_Elem1, nTri_Elem
use Constants, only: Zero, One, Two, Three, Five
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, iCNum
real(kind=wp), intent(out) :: Fil(nPol*nPart,3,nTri_Elem(iQ_Atoms),10), EEDisp
real(kind=wp), intent(in) :: Ax, Ay, Az, BoMaH(iQ_Atoms), BoMaO(iQ_Atoms), outxyz(3,nTri_Elem(iQ_Atoms))
real(kind=wp), intent(inout) :: Eint(nTri_Elem(iQ_Atoms),10), Eint_Nuc(iQ_Atoms)
integer(kind=iwp) :: i, ijhr, ip, j, k, nMltTemp
real(kind=wp) :: distMin, EintSl(nTri3_Elem1(MxMltp)), EintSl_Nuc, EintSlTemp, G(3), R2(5), Rab3i(5), Rab(3,5), Rg(5), Se(5), U(3,5)
logical(kind=iwp) :: lAtom, Skip

EEdisp = Zero
!----------------------------------------------------------------------*
! Compute some distances and inverted distances etc. The potential,    *
! the field and etc. and when we already have the numbers, we also do  *
! the dispersion interaction.                                          *
!----------------------------------------------------------------------*
do k=1,nTri_Elem(iQ_Atoms)
  lAtom = (k <= iQ_atoms)
  G(1) = outxyz(1,k)+Ax
  G(2) = outxyz(2,k)+Ay
  G(3) = outxyz(3,k)+Az
  do j=iCnum+1,nPart
    i = (j-1)*nCent
    ip = (j-1)*nPol
    !Below follows a lot of distances to and fro.
    Rab(:,1) = Cordst(:,i+1)-G
    Rab(:,2) = Cordst(:,i+2)-G
    Rab(:,3) = Cordst(:,i+3)-G
    Rab(:,4) = Cordst(:,i+4)-G
    Rab(:,5) = Cordst(:,i+5)-G
    R2(:) = Rab(1,:)**2+Rab(2,:)**2+Rab(3,:)**2
    Rg(:) = sqrt(R2)
    Se(:) = One/Rg
    Rab3i(:) = Se/R2
    !------------------------------------------------------------------*
    ! The dispersion interaction between QM-atoms and solvent is       *
    ! computed, with or without damping. The initial if-clause sees to *
    ! that only atom-centers are included, while bonds and virtual     *
    ! centers are ignored.                                             *
    !------------------------------------------------------------------*
    if (lAtom) call DispEnergy(EEDisp,BoMaH(k),BoMaO(k),Rg(1),Rg(2),Rg(3),Rab3i(1),Rab3i(2),Rab3i(3),k)
    U(1,:) = Rab(1,:)*Se
    U(2,:) = Rab(2,:)*Se
    U(3,:) = Rab(3,:)*Se
    !------------------------------------------------------------------*
    ! Now we wrap up the electrostatics.                               *
    !------------------------------------------------------------------*
    ! Jose. First we check if distance with at least one of the centers*
    ! of the clasical molecule is inside the Cut-off                   *
    !------------------------------------------------------------------*
    distMin = minval(Rg)
    Skip = .false.
    if (distMin <= Cut_Elc) then
      if (lQuad) then
        ! this is done because for Sl_Grad charge is L=0, dipole is L=1 and
        ! quadrupole is L=2. In QmStat they are 1, 2 and 3, respectively.
        nMltTemp = nMlt-1

        call Sl_Grad(nSlSiteC,lMltSlC,Rab,Rg,Se,SlExpC,SlFactC,SlPC,nMltTemp,SlExpQ(:,k),DifSlExp,EintSl,EintSl_Nuc,lAtom)

        !--------------------------------------------------------------*
        ! Change in the order of field gradients because subroutine    *
        ! Sl_Grad has the same order as Molcas: xx=5, xy=6, xz=7, yy=8,*
        ! yz=9 and zz=10 and we need the QmStat order: xx=5, xy=6,     *
        ! yy=7, xz=8, yz=9 and zz=10                                   *
        !--------------------------------------------------------------*
        EintSlTemp = EintSl(7)
        EintSl(7) = EintSl(8)
        EintSl(8) = EintSlTemp

        Eint(k,:) = Eint(k,:)-EintSl(:) !Check below why it is a subtraction and not a sum
        if (lAtom) Eint_Nuc(k) = Eint_Nuc(k)-EintSl_Nuc
        Skip = .true.
      else
        ! this is done because for Sl_Grad charge is L=0, dipole is L=1 and
        ! quadrupole is L=2. In QmStat they are 1, 2 and 3, respectively.
        nMltTemp = nMlt-1

        ijhr = min(nMltTemp,1)
        call Sl_Grad(nSlSiteC,lMltSlC,Rab,Rg,Se,SlExpC,SlFactC,SlPC,ijhr,SlExpQ(:,k),DifSlExp,EintSl,EintSl_Nuc,lAtom)

        Eint(k,1:4) = Eint(k,1:4)-EintSl(1:4) ! Check below why it is a subtraction and not a sum
        if (lAtom) Eint_Nuc(k) = Eint_Nuc(k)-EintSl_Nuc
      end if
    else
      !----------------------------------------------------------------*
      ! The Eint(k,1) term will below turn into the interaction between*
      ! charges on water and the MME-charges on the QM-molecule. Plus  *
      ! the interaction between the charge densities in Water and the  *
      ! charge densities in QM-Molecule when the Penetration is        *
      ! evaluated                                                      *
      !----------------------------------------------------------------*
      Eint(k,1) = -Qsta(1)*Se(2)-Qsta(2)*Se(3)-Qsta(3)*Se(4)-Qsta(4)*Se(5)+Eint(k,1)
      if (lAtom) Eint_Nuc(k) = -Qsta(1)*Se(2)-Qsta(2)*Se(3)-Qsta(3)*Se(4)-Qsta(4)*Se(5)+Eint_Nuc(k)

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
    end if

    if (.not. Skip) then
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
    end if

    !------------------------------------------------------------------*
    ! And now a whole lot of grad(1/r) and higher...                   *
    !------------------------------------------------------------------*
    ! Monopoles.
    Fil(ip+1:ip+3,1,k,1) = Rab(1,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,2,k,1) = Rab(2,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,3,k,1) = Rab(3,1:3)*Rab3i(1:3)
    ! Dipole -- x-component.
    Fil(ip+1:ip+3,1,k,2) = (Three*U(1,1:3)**2-One)*Rab3i(1:3)
    Fil(ip+1:ip+3,2,k,2) = Three*U(2,1:3)*U(1,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,3,k,2) = Three*U(3,1:3)*U(1,1:3)*Rab3i(1:3)
    ! Dipole -- y-component.
    Fil(ip+1:ip+3,1,k,3) = Three*U(2,1:3)*U(1,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,2,k,3) = (Three*U(2,1:3)**2-One)*Rab3i(1:3)
    Fil(ip+1:ip+3,3,k,3) = Three*U(3,1:3)*U(2,1:3)*Rab3i(1:3)
    ! Dipole -- z-component.
    Fil(ip+1:ip+3,1,k,4) = Three*U(3,1:3)*U(1,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,2,k,4) = Three*U(3,1:3)*U(2,1:3)*Rab3i(1:3)
    Fil(ip+1:ip+3,3,k,4) = (Three*U(3,1:3)**2-One)*Rab3i(1:3)
    ! Quadrupole -- xx-component.
    Fil(ip+1:ip+3,1,k,5) = U(1,1:3)*(Five*U(1,1:3)**2-Two)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,5) = Five*U(2,1:3)*U(1,1:3)**2*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,5) = Five*U(3,1:3)*U(1,1:3)**2*Rab3i(1:3)*Se(1:3)
    ! Quadrupole -- yy-component.
    Fil(ip+1:ip+3,1,k,7) = Five*U(2,1:3)**2*U(1,1:3)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,7) = U(2,1:3)*(Five*U(2,1:3)**2-Two)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,7) = Five*U(3,1:3)*U(2,1:3)**2*Rab3i(1:3)*Se(1:3)
    ! Quadrupole -- zz-component.
    Fil(ip+1:ip+3,1,k,10) = Five*U(3,1:3)**2*U(1,1:3)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,10) = Five*U(3,1:3)**2*U(2,1:3)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,10) = U(3,1:3)*(Five*U(3,1:3)**2-Two)*Rab3i(1:3)*Se(1:3)
    ! Quadrupole -- xy-component.
    Fil(ip+1:ip+3,1,k,6) = U(2,1:3)*(Five*U(1,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,6) = U(1,1:3)*(Five*U(2,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,6) = Five*U(3,1:3)*U(2,1:3)*U(1,1:3)*Rab3i(1:3)*Se(1:3)
    ! Quadrupole -- xz-component.
    Fil(ip+1:ip+3,1,k,8) = U(3,1:3)*(Five*U(1,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,8) = Five*U(3,1:3)*U(2,1:3)*U(1,1:3)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,8) = U(1,1:3)*(Five*U(3,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    ! Quadrupole -- yz-component.
    Fil(ip+1:ip+3,1,k,9) = Five*U(3,1:3)*U(2,1:3)*U(1,1:3)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,2,k,9) = U(3,1:3)*(Five*U(2,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    Fil(ip+1:ip+3,3,k,9) = U(2,1:3)*(Five*U(3,1:3)**2-One)*Rab3i(1:3)*Se(1:3)
    !------------------------------------------------------------------*
    ! If damping of the field is requested, then do it.                *
    !------------------------------------------------------------------*
    if (FieldDamp) then
      Fil(ip+1,:,k,:) = Fil(ip+1,:,k,:)*(One-exp(CAFieldG*Rg(1)))**CFexp
      Fil(ip+2,:,k,:) = Fil(ip+2,:,k,:)*(One-exp(CBFieldG*Rg(2)))**CFexp
      Fil(ip+3,:,k,:) = Fil(ip+3,:,k,:)*(One-exp(CBFieldG*Rg(3)))**CFexp
    end if
  end do
end do

return

end subroutine OneOverR_Sl
