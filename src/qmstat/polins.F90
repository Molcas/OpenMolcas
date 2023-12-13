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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  Polins
!
!> @brief
!>   Add the field from the QM-region onto the solvent.
!>   Include the field from the polarizabilities in the solvent onto the QM-region.
!>   (The effect of the static field is taken care of in helstate)
!> @author A. Ohrn
!>
!> @details
!> First we compute the field from the QM-region onto the solvent. The
!> central quantity is the density matrix, which gives us how the
!> QM-region polarizes. The reaction field due to the QM-region is also
!> accounted for. Then we allow the polarized field from the solvent
!> to interact with the QM-region. The static field from the solvent,
!> in other word that from the charges, is already coupled in
!> ::helstate.
!>
!> @param[in,out] Energy   The energy of the electrostatic interaction
!> @param[in]     iAtom2   Number of particles in the solvent, times number of polarizabilities per solvent molecule
!> @param[in]     iCi      Number of centers in QM-molecule
!> @param[in]     Fil      The static field from the solvent
!> @param[out]    VpolMat  The polarization matrix
!> @param[in,out] FFp      The field from the induced dipoles in the solvent
!> @param[in]     polfac   A factor for the computation of the image
!> @param[out]    poli     The solvent polarized field on QM
!> @param[out]    xyzmyq   Total dipole of QM-region
!> @param[in]     xyzmyi   Total induced dipole of solvent
!> @param[in]     xyzmyp   Total permanent dipole of solvent
!> @param[in]     iCstart  Number to keep track of solvent molecules
!> @param[in]     iQ_Atoms
!> @param[in]     qtot     Total charge of QM-region
!> @param[in]     ChaNuc
!> @param[in]     RoMatSt
!> @param[out]    xyzQuQ
!> @param[in]     CT
!***********************************************************************

subroutine Polins(Energy,iAtom2,iCi,Fil,VpolMat,FFp,polfac,poli,xyzmyq,xyzmyi,xyzmyp,iCstart,iQ_Atoms,qtot,ChaNuc,RoMatSt,xyzQuQ,CT)

use qmstat_global, only: ChargedQM, Cordst, nCent, nPart, nPol, nState, outxyzRAS, Pol, RasCha, RasDip, RasQua
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAtom2, iCI, iCstart, iQ_Atoms
real(kind=wp), intent(inout) :: Energy, FFp(nPart*nPol,3)
real(kind=wp), intent(in) :: Fil(nPart*nPol,3,iCi,10), polfac, xyzmyi(3), xyzmyp(3), qtot, ChaNuc(iQ_Atoms), &
                             RoMatSt(nTri_Elem(nState)), CT(3)
real(kind=wp), intent(out) :: VpolMat(nTri_Elem(nState)), Poli(iCi,10), xyzmyq(3), xyzQuQ(6)
integer(kind=iwp) :: i, iCnum, iS, Iu, j, jS, k, kaunt, l
real(kind=wp) :: CofC(3), Gunnar(10), G(3), qD(6), qK(6), qs, qQ(6), Trace1, Trace2, xyzMyC(3)
real(kind=wp), allocatable :: Dm(:,:), Eil(:,:), Qm(:), QQm(:,:)

!----------------------------------------------------------------------*
! Begin with some zeros.                                               *
!----------------------------------------------------------------------*
call mma_allocate(Qm,iCi,label='Qm')
call mma_allocate(Dm,iCi,3,label='Dm')
call mma_allocate(QQm,iCi,6,label='QQm')
call mma_allocate(Eil,iAtom2,3,label='Eil')
iCnum = iCStart/nCent
Qm(1:iQ_Atoms) = -ChaNuc
Qm(iQ_Atoms+1:) = Zero
Dm(:,:) = Zero
QQm(:,:) = Zero
Eil(:,:) = Zero
!----------------------------------------------------------------------*
! Below we compute how the MME of the QM-molecule changes with the new *
! density matrix Romat. In this step we connect a state density with   *
! various multipoles, which will go on and interact with the solvent.  *
! This is one of the pivotal steps in coupling the QM-region with the  *
! classical region.                                                    *
!----------------------------------------------------------------------*
kaunt = 0
do i=1,nState
  do j=1,i
    kaunt = kaunt+1
    Qm(:) = Qm(:)+RasCha(kaunt,:)*RomatSt(kaunt)
    Dm(:,1) = Dm(:,1)+RasDip(kaunt,1,:)*RomatSt(kaunt)
    Dm(:,2) = Dm(:,2)+RasDip(kaunt,2,:)*RomatSt(kaunt)
    Dm(:,3) = Dm(:,3)+RasDip(kaunt,3,:)*RomatSt(kaunt)
    QQm(:,1) = QQm(:,1)+RasQua(kaunt,1,:)*RomatSt(kaunt)
    QQm(:,3) = QQm(:,3)+RasQua(kaunt,3,:)*RomatSt(kaunt)
    QQm(:,6) = QQm(:,6)+RasQua(kaunt,6,:)*RomatSt(kaunt)
    QQm(:,2) = QQm(:,2)+RasQua(kaunt,2,:)*RomatSt(kaunt)
    QQm(:,4) = QQm(:,4)+RasQua(kaunt,4,:)*RomatSt(kaunt)
    QQm(:,5) = QQm(:,5)+RasQua(kaunt,5,:)*RomatSt(kaunt)
  end do
end do
xyzMyQ(:) = Zero
xyzMyC(:) = Zero
xyzQuQ(:) = Zero
CofC(:) = Zero
qQ(:) = Zero
qD(:) = Zero
qK(:) = Zero
do i=1,iCi
  xyzMyQ(:) = xyzMyQ(:)+Dm(i,:)+Qm(i)*outxyzRAS(:,i)
  qQ(1) = qQ(1)+Qm(i)*(outxyzRAS(1,i)-CT(1))*(outxyzRAS(1,i)-CT(1))
  qQ(2) = qQ(2)+Qm(i)*(outxyzRAS(1,i)-CT(1))*(outxyzRAS(2,i)-CT(2))
  qQ(3) = qQ(3)+Qm(i)*(outxyzRAS(1,i)-CT(1))*(outxyzRAS(3,i)-CT(3))
  qQ(4) = qQ(4)+Qm(i)*(outxyzRAS(2,i)-CT(2))*(outxyzRAS(2,i)-CT(2))
  qQ(5) = qQ(5)+Qm(i)*(outxyzRAS(2,i)-CT(2))*(outxyzRAS(3,i)-CT(3))
  qQ(6) = qQ(6)+Qm(i)*(outxyzRAS(3,i)-CT(3))*(outxyzRAS(3,i)-CT(3))
  qD(1) = qD(1)+Two*Dm(i,1)*(outxyzRAS(1,i)-CT(1))
  qD(2) = qD(2)+Dm(i,1)*(outxyzRAS(2,i)-CT(2))+Dm(i,2)*(outxyzRAS(1,i)-CT(1))
  qD(3) = qD(3)+Dm(i,1)*(outxyzRAS(3,i)-CT(3))+Dm(i,3)*(outxyzRAS(1,i)-CT(1))
  qD(4) = qD(4)+Two*Dm(i,2)*(outxyzRAS(2,i)-CT(2))
  qD(5) = qD(5)+Dm(i,2)*(outxyzRAS(3,i)-CT(3))+Dm(i,3)*(outxyzRAS(2,i)-CT(2))
  qD(6) = qD(6)+Two*Dm(i,3)*(outxyzRAS(3,i)-CT(3))
  qK(1) = qK(1)+QQm(i,1)
  qK(2) = qK(2)+QQm(i,2)
  qK(3) = qK(3)+QQm(i,4)
  qK(4) = qK(4)+QQm(i,3)
  qK(5) = qK(5)+QQm(i,5)
  qK(6) = qK(6)+QQm(i,6)
end do
Trace1 = qQ(1)+qQ(4)+qQ(6)
Trace2 = qD(1)+qD(4)+qD(6)
Trace1 = Trace1
Trace2 = Trace2
xyzQuQ(:) = OneHalf*(qQ+qD)+qK
xyzQuQ(1) = xyzQuQ(1)-Half*(Trace1+Trace2)
xyzQuQ(4) = xyzQuQ(4)-Half*(Trace1+Trace2)
xyzQuQ(6) = xyzQuQ(6)-Half*(Trace1+Trace2)
if (ChargedQM) then  !If charged system, then do...
  qs = Zero
  do i=1,iCi
    CofC(:) = CofC+abs(qm(i))*outxyzRAS(:,i) !Center of charge
    qs = qs+abs(qm(i))
  end do
  CofC(:) = CofC/qs
  G(:) = CofC-outxyzRAS(:,1)+Cordst(:,1) !Where C-of-C is globally
  xyzMyC(:) = xyzMyC+qtot*G !Dipole
end if
! The energy of the induced dipole in its reaction field. It is ok since polfac*(xyzMyQ+xyzMyi) is
! the field from the induced dipole according to the image approximation. And the sought energy
! is -0.5*my_perm*R_ind, which is the thing below, see Boethcer p. 145.
do i=1,3
  Energy = Energy+Polfac*xyzMyQ(i)*(xyzMyQ(i)+xyzMyi(i))
end do
!----------------------------------------------------------------------*
! The multipoles of the QM-region, modified due to the polarization,   *
! now interact with each polarizability in the solvent.                *
!----------------------------------------------------------------------*
do i=1,iCi
  do j=1+(nPol*iCnum),iAtom2
    Eil(j,:) = Eil(j,:)+Fil(j,:,i,1)*Qm(i)+ &
               Fil(j,:,i,2)*Dm(i,1)+Fil(j,:,i,3)*Dm(i,2)+Fil(j,:,i,4)*Dm(i,3)+ &
               Fil(j,:,i,5)*QQm(i,1)+Fil(j,:,i,7)*QQm(i,3)+Fil(j,:,i,10)*QQm(i,6)+ &
               Fil(j,:,i,6)*QQm(i,2)*Two+Fil(j,:,i,8)*QQm(i,4)*Two+Fil(j,:,i,9)*QQm(i,5)*Two
  end do
end do
! THIS IS LEBENSGEFAHRLICH (original comment says it all!)
!----------------------------------------------------------------------*
! We add up the field from the QM-region to the field on all the       *
! solvent polarizabilities.                                            *
!----------------------------------------------------------------------*
do i=1+(nPol*iCNum),iAtom2
  Iu = i-((i-1)/nPol)*nPol
  ! Here we add the QM-molecule image to the solvent polarizabilities.
  ! Good old classical dielectric cavity model!
  ! FFp now contains the field on each solvent polarizability from all different sources.
  do j=1,3
    FFp(i,j) = FFp(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
    Energy = Energy+FFp(i,j)*Eil(i,j)*Pol(iu)
  end do
end do
call mma_deallocate(Qm)
call mma_deallocate(Dm)
call mma_deallocate(QQm)
call mma_deallocate(Eil)
!----------------------------------------------------------------------*
! Now we wish to make the induced field from the solvent interact with *
! the QM-region. The static field has already interacted in helstate.  *
! The reaction field of the QM-region in the dielectric cavity is      *
! also included, excluding the quadrupoles and higher; they are        *
! small anyway, so this is not a major restriction.                    *
!----------------------------------------------------------------------*
Gunnar(:) = Zero
Gunnar(2:4) = PolFac*(xyzMyP+xyzMyQ+xyzMyI+xyzMyC)
do l=1,iCi
  ! The potential from the dipole (the 1/r*r*r is in Fil(...).
  Gunnar(1) = Gunnar(2)*outxyzRAS(1,l)+Gunnar(3)*outxyzRAS(2,l)+Gunnar(4)*outxyzRAS(3,l)
  do i=1,10
    Poli(l,i) = Gunnar(i)
    do j=1+(nPol*iCnum),iAtom2
      Iu = j-((j-1)/nPol)*nPol
      do k=1,3
        !Poli is the polarized field of the solvent. Remember that
        !FFp() is the total field on the polarizabilities.
        Poli(l,i) = Poli(l,i)-FFp(j,k)*Pol(iu)*Fil(j,k,l,i)
      end do
    end do
  end do
end do
VpolMat(:) = Zero
kaunt = 0
! Attention! The reason we use RasCha etc. and not the computed Qm, Dm etc. from above is
! that the density we want to describe is the density of the basis-functions. Compare with
! ordinary <psi_i|V_el|psi_j>.
do iS=1,nState
  do jS=1,iS
    kaunt = kaunt+1
    do j=1,iCi
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,1)*RasCha(kaunt,j)+ &
                       Poli(j,2)*RasDip(kaunt,1,j)+Poli(j,3)*RasDip(kaunt,2,j)+Poli(j,4)*RasDip(kaunt,3,j)+ &
                       Poli(j,5)*RasQua(kaunt,1,j)+Poli(j,7)*RasQua(kaunt,3,j)+Poli(j,10)*RasQua(kaunt,6,j)+ &
                       Poli(j,6)*RasQua(kaunt,2,j)*Two+Poli(j,8)*RasQua(kaunt,4,j)*Two+Poli(j,9)*RasQua(kaunt,5,j)*Two
    end do
  end do
end do
! This is how the nuclei interact with the induced field (in equil2 exists a corresponding
! term for the interaction with the static field).
! This way interaction between a charged molecule and the induced/permanent potential is included.
do i=1,iQ_Atoms
  Energy = Energy-Two*Poli(i,1)*ChaNuc(i)
end do

return

end subroutine Polins
