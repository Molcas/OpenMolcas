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
!>   (The effect of the static field is taken care of in hel.f)
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
!> which is done just like in hel.f.
!>
!> @param[out]    Energy  The energy of the electrostatic interaction
!> @param[in,out] iCall   An integer that tells if this is the first call in the iteration.
!>                        Necessary for the copy of the one-particle Hamiltonian
!> @param[in]     iAtom2  Number of particles in the solvent, times number of polarizabilities per solvent molecule
!> @param[in]     iCi     Number of centers in QM-molecule
!> @param[in]     iFil    Pointer to the static field from the solvent
!> @param[out]    VpolMat The matrix due to polarization
!> @param[in,out] fil     The field from the induced dipoles in the solvent
!> @param[in]     polfac  A factor for the computation of the image
!> @param[out]    poli    The solvent polarized field on QM-region
!> @param[in]     iCstart Number to keep track of solvent molecules
!> @param[in]     iTri    ``iOrb(1)*(iOrb(1)+1)/2``
!***********************************************************************

subroutine Polink(Energy,iCall,iAtom2,iCi,iFil,VpolMat,fil,polfac,poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP, &
                  RoMat,xyzQuQ,CT)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, OneHalf
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qm1.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iCall, iAtom2, iCi, iFil(iCi,10), iCstart, iTri, iQ_Atoms
real(kind=wp) :: Energy, VpolMat(iTri), Fil(npart*npol,3), polfac, Poli(iCi,10), qTot, ChaNuc(MxAt), xyzMyQ(3), xyzMyI(3), &
                 xyzMyP(3), RoMat(iTri), xyzQuQ(6), CT(3)
integer(kind=iwp) :: i, iCnum, Iu, j, k, kk, l
real(kind=wp) :: CofC(3), Gunnar(10), Gx, Gy, Gz, qD(6), qK(6), qQ(6), qs, Trace1, Trace2, xyzMyC(3)
real(kind=wp), allocatable :: Dm(:,:), Eil(:,:), Qm(:), QQm(:,:)

!----------------------------------------------------------------------*
! Begin with some zeros.                                               *
!----------------------------------------------------------------------*
call mma_allocate(Qm,iCi,label='Qm')
call mma_allocate(Dm,iCi,3,label='Dm')
call mma_allocate(QQm,iCi,6,label='QQm')
call mma_allocate(Eil,iAtom2,3,label='Eil')
iCnum = iCStart/Ncent
do i=1,iCi
  Qm(i) = Zero
  if (i <= iQ_Atoms) Qm(i) = -ChaNuc(i) !Here is nuclear contribution added to atoms.
  do j=1,3
    Dm(i,j) = Zero
    QQm(i,j) = Zero
    QQm(i,j+3) = Zero
  end do
end do
do i=1,iAtom2
  do j=1,3
    Eil(i,j) = Zero
  end do
end do
!----------------------------------------------------------------------*
! Below we compute how the MME of the QM-molecule changes with the new *
! density matrix Romat. What we actually do is a HF-SCF procedure with *
! a MME-expanded density. A change in the density has the effect that  *
! the set of multipoles in the MME are slightly perturbed.             *
!----------------------------------------------------------------------*
do i=1,iTri
  do j=1,iCi
    Qm(j) = Cha(i,j)*Romat(i)+Qm(j)
    Dm(j,1) = Dm(j,1)+DipMy(i,1,j)*Romat(i)
    Dm(j,2) = Dm(j,2)+DipMy(i,2,j)*Romat(i)
    Dm(j,3) = Dm(j,3)+DipMy(i,3,j)*Romat(i)
    QQm(j,1) = QQm(j,1)+Quad(i,1,j)*Romat(i)
    QQm(j,3) = QQm(j,3)+Quad(i,3,j)*Romat(i)
    QQm(j,6) = QQm(j,6)+Quad(i,6,j)*Romat(i)
    QQm(j,2) = QQm(j,2)+Quad(i,2,j)*Romat(i)
    QQm(j,4) = QQm(j,4)+Quad(i,4,j)*Romat(i)
    QQm(j,5) = QQm(j,5)+Quad(i,5,j)*Romat(i)
  end do
end do
do kk=1,3
  xyzMyQ(kk) = Zero
  xyzMyC(kk) = Zero
  CofC(kk) = Zero
end do
do kk=1,6
  xyzQuQ(kk) = Zero
  qQ(kk) = Zero
  qD(kk) = Zero
  qK(kk) = Zero
end do
! Observe one tricky thing about xyzmyq: the electric multipoles
! we use above are actually of opposite sign, so how can xyzmyq be
! the dipole in the qm-region unless we change sign (which we do
! not)? The reason is that the density matrix elements will also
! have opposite sign, which in turn has no physical meaning.
! We also compute the quadupole moment -- a messy formula.
do i=1,iCi
  xyzMyQ(1) = xyzMyQ(1)+Dm(i,1)+Qm(i)*outxyz(i,1)
  xyzMyQ(2) = xyzMyQ(2)+Dm(i,2)+Qm(i)*outxyz(i,2)
  xyzMyQ(3) = xyzMyQ(3)+Dm(i,3)+Qm(i)*outxyz(i,3)
  qQ(1) = qQ(1)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,1)-CT(1))
  qQ(2) = qQ(2)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,2)-CT(2))
  qQ(3) = qQ(3)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,3)-CT(3))
  qQ(4) = qQ(4)+Qm(i)*(outxyz(i,2)-CT(2))*(outxyz(i,2)-CT(2))
  qQ(5) = qQ(5)+Qm(i)*(outxyz(i,2)-CT(2))*(outxyz(i,3)-CT(3))
  qQ(6) = qQ(6)+Qm(i)*(outxyz(i,3)-CT(3))*(outxyz(i,3)-CT(3))
  qD(1) = qD(1)+Two*Dm(i,1)*(outxyz(i,1)-CT(1))
  qD(2) = qD(2)+Dm(i,1)*(outxyz(i,2)-CT(2))+Dm(i,2)*(outxyz(i,1)-CT(1))
  qD(3) = qD(3)+Dm(i,1)*(outxyz(i,3)-CT(3))+Dm(i,3)*(outxyz(i,1)-CT(1))
  qD(4) = qD(4)+Two*Dm(i,2)*(outxyz(i,2)-CT(2))
  qD(5) = qD(5)+Dm(i,2)*(outxyz(i,3)-CT(3))+Dm(i,3)*(outxyz(i,2)-CT(2))
  qD(6) = qD(6)+Two*Dm(i,3)*(outxyz(i,3)-CT(3))
  qK(1) = qK(1)+QQm(i,1)
  qK(2) = qK(2)+QQm(i,2)
  qK(3) = qK(3)+QQm(i,4)
  qK(4) = qK(4)+QQm(i,3)
  qK(5) = qK(5)+QQm(i,5)
  qK(6) = qK(6)+QQm(i,6)
end do
Trace1 = qQ(1)+qQ(4)+qQ(6)
Trace2 = qD(1)+qD(4)+qD(6)
Trace1 = Trace1/Three
Trace2 = Trace2/Three
xyzQuQ(1) = OneHalf*(qQ(1)+qD(1)-Trace1-Trace2)+qK(1)
xyzQuQ(2) = OneHalf*(qQ(2)+qD(2))+qK(2)
xyzQuQ(3) = OneHalf*(qQ(3)+qD(3))+qK(3)
xyzQuQ(4) = OneHalf*(qQ(4)+qD(4)-Trace1-Trace2)+qK(4)
xyzQuQ(5) = OneHalf*(qQ(5)+qD(5))+qK(5)
xyzQuQ(6) = OneHalf*(qQ(6)+qD(6)-Trace1-Trace2)+qK(6)
if (ChargedQM) then !If charged system, then do...
  qs = Zero
  do i=1,iCi
    CofC(1) = CofC(1)+abs(Qm(i))*outxyz(i,1) !Center of charge
    CofC(2) = CofC(2)+abs(Qm(i))*outxyz(i,2)
    CofC(3) = CofC(3)+abs(Qm(i))*outxyz(i,3)
    qs = qs+abs(Qm(i))
  end do
  CofC(1) = CofC(1)/qs
  CofC(2) = CofC(2)/qs
  CofC(3) = CofC(3)/qs
  Gx = CofC(1)-outxyz(1,1)+Cordst(1,1) !Where C-of-C is globally
  Gy = CofC(2)-outxyz(1,2)+Cordst(1,2)
  Gz = CofC(3)-outxyz(1,3)+Cordst(1,3)
  xyzMyC(1) = xyzMyC(1)+qtot*Gx !Dipole
  xyzMyC(2) = xyzMyC(2)+qtot*Gy
  xyzMyC(3) = xyzMyC(3)+qtot*Gz
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
    do k=1,3
      Eil(j,k) = Eil(j,k)+Work(iFil(i,1)-1+j+(k-1)*nPart*nPol)*Qm(i)
      do l=1,3
        Eil(j,k) = Eil(j,k)+Work(iFil(i,l+1)-1+j+(k-1)*nPart*nPol)*Dm(i,l)
      end do
      Eil(j,k) = Eil(j,k)+Work(iFil(i,5)-1+j+(k-1)*nPart*nPol)*QQm(i,1)
      Eil(j,k) = Eil(j,k)+Work(iFil(i,7)-1+j+(k-1)*nPart*nPol)*QQm(i,3)
      Eil(j,k) = Eil(j,k)+Work(iFil(i,10)-1+j+(k-1)*nPart*nPol)*QQm(i,6)
      Eil(j,k) = Eil(j,k)+Work(iFil(i,6)-1+j+(k-1)*nPart*nPol)*QQm(i,2)*Two
      Eil(j,k) = Eil(j,k)+Work(iFil(i,8)-1+j+(k-1)*nPart*nPol)*QQm(i,4)*Two
      Eil(j,k) = Eil(j,k)+Work(iFil(i,9)-1+j+(k-1)*nPart*nPol)*QQm(i,5)*Two
    end do
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
    Fil(i,j) = Fil(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
    ! How much the induced dipoles in solvent interacts with the field from the QM-region.
    Energy = Energy+Fil(i,j)*Eil(i,j)*Pol(iu)
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
do i=1,10
  Gunnar(i) = 0
end do
Gunnar(2) = PolFac*(xyzMyP(1)+xyzMyQ(1)+xyzMyI(1)+xyzMyC(1))
Gunnar(3) = PolFac*(xyzMyP(2)+xyzMyQ(2)+xyzMyI(2)+xyzMyC(2))
Gunnar(4) = PolFac*(xyzMyP(3)+xyzMyQ(3)+xyzMyI(3)+xyzMyC(3))
do l=1,iCi
  ! Potential from the apparent surface charge, see Boettcher (4.22).
  Gunnar(1) = Gunnar(2)*outxyz(l,1)+Gunnar(3)*outxyz(l,2)+Gunnar(4)*outxyz(l,3)
  do i=1,10
    Poli(l,i) = Gunnar(i)
    do j=1+(nPol*iCnum),iAtom2
      Iu = j-((j-1)/nPol)*nPol
      ! Compute the generalized field from induced dipoles in solvent on the QM-region cites.
      do k=1,3
        Poli(l,i) = Poli(l,i)-Fil(j,k)*Pol(iu)*Work(iFil(l,i)-1+j+(k-1)*nPart*nPol)
      end do
    end do
  end do
end do
do i=1,iTri
  VpolMat(i) = 0
end do
do i=1,iTri
  do j=1,iCi
    Vpolmat(i) = Vpolmat(i)+Poli(j,1)*Cha(i,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,2)*DipMy(i,1,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,3)*DipMy(i,2,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,4)*DipMy(i,3,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,5)*Quad(i,1,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,7)*Quad(i,3,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,10)*Quad(i,6,j)
    Vpolmat(i) = Vpolmat(i)+Poli(j,6)*Quad(i,2,j)*Two
    Vpolmat(i) = Vpolmat(i)+Poli(j,8)*Quad(i,4,j)*Two
    Vpolmat(i) = Vpolmat(i)+Poli(j,9)*Quad(i,5,j)*Two
  end do
end do
! This is how the nuclei interact with the induced field (in equil2 exists a corresponding
! term for the interaction with the static field).
! This way interaction between a charged molecule and the induced/permanent potential is included.
do i=1,iQ_Atoms
  Energy = Energy-2*Poli(i,1)*ChaNuc(i)
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iCall)

end subroutine Polink
