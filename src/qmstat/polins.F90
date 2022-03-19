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
!>   (The effect of the static field is taken care of in helstate.f)
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
!> @param[out]    Energy  The energy of the electrostatic interaction
!> @param[in]     iAtom2  Number of particles in the solvent, times number of polarizabilities per solvent molecule
!> @param[in]     iCi     Number of centers in QM-molecule
!> @param[in]     iFil    Pointer to the static field from the solvent
!> @param[out]    VpolMat The polarization matrix
!> @param[in,out] fil     The field from the induced dipoles in the solvent
!> @param[in]     polfac  A factor for the computation of the image
!> @param[out]    poli    The solvent polarized field on QM
!> @param[in]     xyzmyq  Total dipole of QM-region
!> @param[in]     xyzmyi  Total induced dipole of solvent
!> @param[in]     xyzmyp  Total permanent dipole of solvent
!> @param[in]     qtot    Total charge of QM-region
!> @param[in]     iCstart Number to keep track of solvent molecules
!***********************************************************************

subroutine Polins(Energy,iAtom2,iCi,iFil,VpolMat,fil,polfac,poli,xyzmyq,xyzmyi,xyzmyp,iCstart,iQ_Atoms,qtot,ChaNuc,RoMatSt,xyzQuQ, &
                  CT)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, OneHalf
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iAtom2, iCI, iFil(iCi,10), iCstart, iQ_Atoms
real(kind=wp) :: Energy, VpolMat(nTri_Elem(nState)), Fil(npart*npol,3), polfac, Poli(iCi,10), xyzmyq(3), xyzmyi(3), xyzmyp(3), &
                 qtot, ChaNuc(MxAt), RoMatSt(nTri_Elem(nState)), xyzQuQ(6), CT(3)
integer(kind=iwp) :: i, iCnum, iS, Iu, j, jS, k, kaunt, kk, l
real(kind=wp) :: CofC(3), Gunnar(10), Gx, Gy, Gz, qD(6), qK(6), qs, qQ(6), Trace1, Trace2, xyzMyC(3)
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
  if (i <= iQ_Atoms) Qm(i) = -ChaNuc(i)
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
! density matrix Romat. In this step we connect a state density with   *
! various multipoles, which will go on and interact with the solvent.  *
! This is one of the pivotal steps in coupling the QM-region with the  *
! classical region.                                                    *
!----------------------------------------------------------------------*
kaunt = 0
do i=1,nState
  do j=1,i
    kaunt = kaunt+1
    do k=1,iCi
      Qm(k) = Qm(k)+RasCha(kaunt,k)*RomatSt(kaunt)
      Dm(k,1) = Dm(k,1)+RasDip(kaunt,1,k)*RomatSt(kaunt)
      Dm(k,2) = Dm(k,2)+RasDip(kaunt,2,k)*RomatSt(kaunt)
      Dm(k,3) = Dm(k,3)+RasDip(kaunt,3,k)*RomatSt(kaunt)
      QQm(k,1) = QQm(k,1)+RasQua(kaunt,1,k)*RomatSt(kaunt)
      QQm(k,3) = QQm(k,3)+RasQua(kaunt,3,k)*RomatSt(kaunt)
      QQm(k,6) = QQm(k,6)+RasQua(kaunt,6,k)*RomatSt(kaunt)
      QQm(k,2) = QQm(k,2)+RasQua(kaunt,2,k)*RomatSt(kaunt)
      QQm(k,4) = QQm(k,4)+RasQua(kaunt,4,k)*RomatSt(kaunt)
      QQm(k,5) = QQm(k,5)+RasQua(kaunt,5,k)*RomatSt(kaunt)
    end do
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
do i=1,iCi
  do kk=1,3
    xyzMyQ(kk) = xyzMyQ(kk)+Dm(i,kk)+Qm(i)*outxyzRAS(i,kk)
  end do
  qQ(1) = qQ(1)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,1)-CT(1))
  qQ(2) = qQ(2)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,2)-CT(2))
  qQ(3) = qQ(3)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,3)-CT(3))
  qQ(4) = qQ(4)+Qm(i)*(outxyzRAS(i,2)-CT(2))*(outxyzRAS(i,2)-CT(2))
  qQ(5) = qQ(5)+Qm(i)*(outxyzRAS(i,2)-CT(2))*(outxyzRAS(i,3)-CT(3))
  qQ(6) = qQ(6)+Qm(i)*(outxyzRAS(i,3)-CT(3))*(outxyzRAS(i,3)-CT(3))
  qD(1) = qD(1)+Two*Dm(i,1)*(outxyzRAS(i,1)-CT(1))
  qD(2) = qD(2)+Dm(i,1)*(outxyzRAS(i,2)-CT(2))+Dm(i,2)*(outxyzRAS(i,1)-CT(1))
  qD(3) = qD(3)+Dm(i,1)*(outxyzRAS(i,3)-CT(3))+Dm(i,3)*(outxyzRAS(i,1)-CT(1))
  qD(4) = qD(4)+Two*Dm(i,2)*(outxyzRAS(i,2)-CT(2))
  qD(5) = qD(5)+Dm(i,2)*(outxyzRAS(i,3)-CT(3))+Dm(i,3)*(outxyzRAS(i,2)-CT(2))
  qD(6) = qD(6)+Two*Dm(i,3)*(outxyzRAS(i,3)-CT(3))
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
if (ChargedQM) then  !If charged system, then do...
  qs = Zero
  do i=1,iCi
    CofC(1) = CofC(1)+abs(qm(i))*outxyzRAS(i,1) !Center of charge
    CofC(2) = CofC(2)+abs(qm(i))*outxyzRAS(i,2)
    CofC(3) = CofC(3)+abs(qm(i))*outxyzRAS(i,3)
    qs = qs+abs(qm(i))
  end do
  CofC(1) = CofC(1)/qs
  CofC(2) = CofC(2)/qs
  CofC(3) = CofC(3)/qs
  Gx = CofC(1)-outxyzRAS(1,1)+Cordst(1,1) !Where C-of-C is globally
  Gy = CofC(2)-outxyzRAS(1,2)+Cordst(1,2)
  Gz = CofC(3)-outxyzRAS(1,3)+Cordst(1,3)
  xyzMyC(1) = xyzMyC(1)+qtot*Gx !Dipole
  xyzMyC(2) = xyzMyC(2)+qtot*Gy
  xyzMyC(3) = xyzMyC(3)+qtot*Gz
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
  ! Fil now contains the field on each solvent polarizability from all different sources.
  do j=1,3
    Fil(i,j) = Fil(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
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
  Gunnar(i) = Zero
end do
Gunnar(2) = PolFac*(xyzMyP(1)+xyzMyQ(1)+xyzMyI(1)+xyzMyC(1))
Gunnar(3) = PolFac*(xyzMyP(2)+xyzMyQ(2)+xyzMyI(2)+xyzMyC(2))
Gunnar(4) = PolFac*(xyzMyP(3)+xyzMyQ(3)+xyzMyI(3)+xyzMyC(3))
do l=1,iCi
  ! The potential from the dipole (the 1/r*r*r is in Work(iFil...).
  Gunnar(1) = Gunnar(2)*outxyzRAS(l,1)+Gunnar(3)*outxyzRAS(l,2)+Gunnar(4)*outxyzRAS(l,3)
  do i=1,10
    Poli(l,i) = Gunnar(i)
    do j=1+(nPol*iCnum),iAtom2
      Iu = j-((j-1)/nPol)*nPol
      do k=1,3
        !Poli is the polarized field of the solvent. Remember that
        !Fil() is the total field on the polarizabilities.
        Poli(l,i) = Poli(l,i)-Fil(j,k)*Pol(iu)*Work(iFil(l,i)-1+j+(k-1)*nPart*nPol)
      end do
    end do
  end do
end do
do i=1,nTri_Elem(nState)
  VpolMat(i) = Zero
end do
kaunt = 0
! Attention! The reason we use RasCha etc. and not the computed Qm, Dm etc. from above is
! that the density we want to describe is the density of the basis-functions. Compare with
! ordinary <psi_i|V_el|psi_j>.
do iS=1,nState
  do jS=1,iS
    kaunt = kaunt+1
    do j=1,iCi
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,1)*RasCha(kaunt,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,2)*RasDip(kaunt,1,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,3)*RasDip(kaunt,2,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,4)*RasDip(kaunt,3,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,5)*RasQua(kaunt,1,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,7)*RasQua(kaunt,3,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,10)*RasQua(kaunt,6,j)
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,6)*RasQua(kaunt,2,j)*Two
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,8)*RasQua(kaunt,4,j)*Two
      Vpolmat(kaunt) = Vpolmat(kaunt)+Poli(j,9)*RasQua(kaunt,5,j)*Two
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
