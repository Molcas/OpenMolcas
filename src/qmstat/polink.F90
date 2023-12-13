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
!  Polink
!
!> @brief
!>   Add the field from the QM-region onto the solvent. Include the field from the
!>   polarizabilities in the solvent onto the QM-region.
!>   (The effect of the static field is taken care of in hel)
!> @author A. Ohrn
!>
!> @details
!> To begin with we obtain the charge distribution of the QM-region
!> as it exists due to the pressent density matrix (recall that it
!> is the changes in the density matrix that causes the QM-region to
!> be polarized). The field from these new multipoles are added on to
!> solvent. We also include the reaction field from the QM-region.
!> Then, with the new field from the QM-region included, we compute
!> the field from the polarizabiolities in the solvent onto the QM-region,
!> which is done just like in ::hel.
!>
!> @param[in,out] Energy  The energy of the electrostatic interaction
!> @param[in]     iAtom2  Number of particles in the solvent, times number of polarizabilities per solvent molecule
!> @param[in]     iCi     Number of centers in QM-molecule
!> @param[in]     Fil     The static field from the solvent
!> @param[out]    VpolMat The matrix due to polarization
!> @param[in,out] FFp     The field from the induced dipoles in the solvent
!> @param[in]     polfac  A factor for the computation of the image
!> @param[out]    poli    The solvent polarized field on QM-region
!> @param[in]     iCstart Number to keep track of solvent molecules
!> @param[in]     iTri    ``iOrb(1)*(iOrb(1)+1)/2``
!> @param[in]     iQ_Atoms
!> @param[in]     qTot
!> @param[in]     ChaNuc
!> @param[out]    xyzMyQ
!> @param[in]     xyzMyI
!> @param[in]     xyzMyP
!> @param[in]     RoMat
!> @param[out]    xyzQuQ
!> @param[in]     CT
!***********************************************************************

subroutine Polink(Energy,iAtom2,iCi,Fil,VpolMat,FFp,polfac,poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP,RoMat, &
                  xyzQuQ,CT)

use qmstat_global, only: Cha, ChargedQM, Cordst, DipMy, nCent, nPart, nPol, outxyz, Pol, Quad
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAtom2, iCi, iCstart, iTri, iQ_Atoms
real(kind=wp), intent(inout) :: Energy, FFp(nPart*nPol,3)
real(kind=wp), intent(out) :: VpolMat(iTri), Poli(iCi,10), xyzMyQ(3), xyzQuQ(6)
real(kind=wp), intent(in) :: Fil(nPart*nPol,3,iCi,10), polfac, qTot, ChaNuc(iQ_Atoms), xyzMyI(3), xyzMyP(3), RoMat(iTri), CT(3)
integer(kind=iwp) :: i, iCnum, Iu, j, l
real(kind=wp) :: CofC(3), Gunnar(10), G(3), qD(6), qK(6), qQ(6), qs, Trace1, Trace2, xyzMyC(3)
real(kind=wp), allocatable :: Dm(:,:), Eil(:,:), Qm(:), QQm(:,:)

!----------------------------------------------------------------------*
! Begin with some zeros.                                               *
!----------------------------------------------------------------------*
call mma_allocate(Qm,iCi,label='Qm')
call mma_allocate(Dm,iCi,3,label='Dm')
call mma_allocate(QQm,iCi,6,label='QQm')
call mma_allocate(Eil,iAtom2,3,label='Eil')
iCnum = iCStart/nCent
Qm(1:iQ_Atoms) = -ChaNuc !Here is nuclear contribution added to atoms.
Qm(iQ_Atoms+1:) = Zero
Dm(:,:) = Zero
QQm(:,:) = Zero
Eil(:,:) = Zero
!----------------------------------------------------------------------*
! Below we compute how the MME of the QM-molecule changes with the new *
! density matrix Romat. What we actually do is a HF-SCF procedure with *
! a MME-expanded density. A change in the density has the effect that  *
! the set of multipoles in the MME are slightly perturbed.             *
!----------------------------------------------------------------------*
do i=1,iTri
  Qm(:) = Cha(i,:)*Romat(i)+Qm(:)
  Dm(:,1) = Dm(:,1)+DipMy(i,1,:)*Romat(i)
  Dm(:,2) = Dm(:,2)+DipMy(i,2,:)*Romat(i)
  Dm(:,3) = Dm(:,3)+DipMy(i,3,:)*Romat(i)
  QQm(:,1) = QQm(:,1)+Quad(i,1,:)*Romat(i)
  QQm(:,3) = QQm(:,3)+Quad(i,3,:)*Romat(i)
  QQm(:,6) = QQm(:,6)+Quad(i,6,:)*Romat(i)
  QQm(:,2) = QQm(:,2)+Quad(i,2,:)*Romat(i)
  QQm(:,4) = QQm(:,4)+Quad(i,4,:)*Romat(i)
  QQm(:,5) = QQm(:,5)+Quad(i,5,:)*Romat(i)
end do
xyzMyQ(:) = Zero
xyzMyC(:) = Zero
xyzQuQ(:) = Zero
CofC(:) = Zero
qQ(:) = Zero
qD(:) = Zero
qK(:) = Zero
! Observe one tricky thing about xyzmyq: the electric multipoles
! we use above are actually of opposite sign, so how can xyzmyq be
! the dipole in the qm-region unless we change sign (which we do
! not)? The reason is that the density matrix elements will also
! have opposite sign, which in turn has no physical meaning.
! We also compute the quadupole moment -- a messy formula.
do i=1,iCi
  xyzMyQ(:) = xyzMyQ(:)+Dm(i,:)+Qm(i)*outxyz(:,i)
  qQ(1) = qQ(1)+Qm(i)*(outxyz(1,i)-CT(1))*(outxyz(1,i)-CT(1))
  qQ(2) = qQ(2)+Qm(i)*(outxyz(1,i)-CT(1))*(outxyz(2,i)-CT(2))
  qQ(3) = qQ(3)+Qm(i)*(outxyz(1,i)-CT(1))*(outxyz(3,i)-CT(3))
  qQ(4) = qQ(4)+Qm(i)*(outxyz(2,i)-CT(2))*(outxyz(2,i)-CT(2))
  qQ(5) = qQ(5)+Qm(i)*(outxyz(2,i)-CT(2))*(outxyz(3,i)-CT(3))
  qQ(6) = qQ(6)+Qm(i)*(outxyz(3,i)-CT(3))*(outxyz(3,i)-CT(3))
  qD(1) = qD(1)+Two*Dm(i,1)*(outxyz(1,i)-CT(1))
  qD(2) = qD(2)+Dm(i,1)*(outxyz(2,i)-CT(2))+Dm(i,2)*(outxyz(1,i)-CT(1))
  qD(3) = qD(3)+Dm(i,1)*(outxyz(3,i)-CT(3))+Dm(i,3)*(outxyz(1,i)-CT(1))
  qD(4) = qD(4)+Two*Dm(i,2)*(outxyz(2,i)-CT(2))
  qD(5) = qD(5)+Dm(i,2)*(outxyz(3,i)-CT(3))+Dm(i,3)*(outxyz(2,i)-CT(2))
  qD(6) = qD(6)+Two*Dm(i,3)*(outxyz(3,i)-CT(3))
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
if (ChargedQM) then !If charged system, then do...
  qs = Zero
  do i=1,iCi
    CofC(:) = CofC+abs(Qm(i))*outxyz(:,i) !Center of charge
    qs = qs+abs(Qm(i))
  end do
  CofC(:) = CofC/qs
  G(:) = CofC-outxyz(:,1)+Cordst(:,1) !Where C-of-C is globally
  xyzMyC(:) = xyzMyC+qtot*G !Dipole
end if
do i=1,3
  ! Change sign on both the dipoles, which in effect gives
  ! no sign change, all in order with Boettcher, p.145.
  Energy = Energy+Polfac*xyzMyQ(i)*(xyzMyQ(i)+xyzMyi(i))
end do
!----------------------------------------------------------------------*
! The multipoles of the QM-region, modified due to the polarization,   *
! now interacts with each polarizability in the solvent.               *
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
  do j=1,3
    FFp(i,j) = FFp(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
    ! How much the induced dipoles in solvent interacts with the field from the QM-region.
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
  ! Potential from the apparent surface charge, see Boettcher (4.22).
  Gunnar(1) = Gunnar(2)*outxyz(1,l)+Gunnar(3)*outxyz(2,l)+Gunnar(4)*outxyz(3,l)
  do i=1,10
    Poli(l,i) = Gunnar(i)
    do j=1+(nPol*iCnum),iAtom2
      Iu = j-((j-1)/nPol)*nPol
      ! Compute the generalized field from induced dipoles in solvent on the QM-region cites.
      Poli(l,i) = Poli(l,i)-FFp(j,1)*Pol(iu)*Fil(j,1,l,i)-FFp(j,2)*Pol(iu)*Fil(j,2,l,i)-FFp(j,3)*Pol(iu)*Fil(j,3,l,i)
    end do
  end do
end do
VpolMat(:) = Zero
do i=1,iTri
  do j=1,iCi
    Vpolmat(i) = Vpolmat(i)+Poli(j,1)*Cha(i,j)+ &
                 Poli(j,2)*DipMy(i,1,j)+Poli(j,3)*DipMy(i,2,j)+Poli(j,4)*DipMy(i,3,j)+ &
                 Poli(j,5)*Quad(i,1,j)+Poli(j,7)*Quad(i,3,j)+Poli(j,10)*Quad(i,6,j)+ &
                 Poli(j,6)*Quad(i,2,j)*Two+Poli(j,8)*Quad(i,4,j)*Two+Poli(j,9)*Quad(i,5,j)*Two
  end do
end do
! This is how the nuclei interact with the induced field (in equil2 exists a corresponding
! term for the interaction with the static field).
! This way interaction between a charged molecule and the induced/permanent potential is included.
do i=1,iQ_Atoms
  Energy = Energy-2*Poli(i,1)*ChaNuc(i)
end do

return

end subroutine Polink
