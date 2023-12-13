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
!  Geogen
!
!> @brief
!>   Generate the new configuration in typical random manner
!> @author A. Ohrn
!>
!> @details
!> The creator of random geometries in typical Monte-Carlo fashion.
!> The changes made are:
!>
!> 1. The dielectric radius is modified
!> 2. Each coordinate *except* the \c iSta-1 -th molecule is changed
!> 3. All molecules (with the above exception) are rotated around the oxygen (which approximately equals the CM)
!> 4. Every molecule except the fixed ones are rotated slightly around one of the global \f$ x \f$-, \f$ y \f$- or \f$ z \f$-axes;
!>    the purpose of this is to emulate a rotation of the central molecule and therefore make the system more dynamic.
!>
!> @param[in,out] Ract     The dielectric radius on input and the slightly perturbed radius on output
!> @param[out]    Rold     Stores the input dielectric radius
!> @param[in]     iCNum    How many solvent places that are taken up by the QM-molecule
!> @param[in]     iQ_Atoms
!***********************************************************************

subroutine Geogen(Ract,Rold,iCNum,iQ_Atoms)

use qmstat_global, only: CordIm, Cordst, DelFi, DelR, delX, Diel, DipIm, iSeed, iSta, nCent, nCha, nPart, nPol, OldGeo, Sqrs, Qim, &
                         Qimp, Qmeq, QmProd, Qsta, xyzMyp
use Constants, only: Zero, One, Two, Three, Ten, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Ract
real(kind=wp), intent(out) :: Rold
integer(kind=iwp), intent(in) :: iCNum, iQ_Atoms
integer(kind=iwp) :: i, ii, iImage, ij, Ind, iQsta, j, k
real(kind=wp) :: A, A2, B, CB, Cx, Cy, Cz, DiFac, Dx, q, qq, S2, SB, Sqrts2, x, xNy, y, yNy, z, zNy
real(kind=wp), external :: Random_Molcas

!----------------------------------------------------------------------*
! Store old configuration.                                             *
!----------------------------------------------------------------------*
OldGeo(:,:) = Cordst
!----------------------------------------------------------------------*
! Query which type of simulation this is, and if quantum then change   *
! the quantum molecule.                                                *
!----------------------------------------------------------------------*
if (Qmeq .or. QmProd) then !Which coordinates to keep fixed.
  iSta = iCNum+1  !This sees to that the QM-molecule is excluded from the moves below.
  do j=1,3
    Dx = delX*(Random_Molcas(iSeed)-Half)
    Cordst(j,1:iQ_Atoms) = Cordst(j,1:iQ_Atoms)+Dx !Move QM-mol.
  end do
end if
!----------------------------------------------------------------------*
! Obtain the random-stuff and make small geometry change.              *
!----------------------------------------------------------------------*
Rold = Ract
Ract = Ract+(Random_Molcas(iSeed)-Half)*DelR !Change in cavity radius
do i=iSta,nPart !Which molecules to give new coordinates.
  ij = (i-1)*nCent
  do j=1,3
    Dx = DelX*(Random_Molcas(iSeed)-Half)
    Cordst(j,ij+1:ij+nCent) = Cordst(j,ij+1:ij+nCent)+Dx !Make translation
  end do
  Cx = Cordst(1,ij+1) !The oxygen, around which we rotate
  Cy = Cordst(2,ij+1)
  Cz = Cordst(3,ij+1)
  B = (Random_Molcas(iSeed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent !Rotate around the oxygen in yz-plane, i.e. around x-axis.
    y = Cordst(2,ij+k)-Cy
    z = Cordst(3,ij+k)-Cz
    yNy = y*CB+z*SB !This is a rotation matrix
    zNy = z*CB-y*SB
    Cordst(2,ij+k) = yNy+Cy
    Cordst(3,ij+k) = zNy+Cz
  end do
  B = (Random_Molcas(iSeed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent !And now rotate in xz-plane
    x = Cordst(1,ij+k)-Cx
    z = Cordst(3,ij+k)-Cz
    xNy = x*CB+z*SB
    zNy = z*CB-x*SB
    Cordst(1,ij+k) = xNy+Cx
    Cordst(3,ij+k) = zNy+Cz
  end do
  B = (Random_Molcas(iSeed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent  !To your surprise, here we rotate in the xy-plane
    x = Cordst(1,ij+k)-Cx
    y = Cordst(2,ij+k)-Cy
    xNy = x*CB+y*SB
    yNy = y*CB-x*SB
    Cordst(1,ij+k) = xNy+Cx
    Cordst(2,ij+k) = yNy+Cy
  end do
end do
! Here all other water molecules rotate around one of the three axes,
! except the ones we fix, which in a QM-simulation is the quantum particle.
A = Random_Molcas(iSeed)
B = (Random_Molcas(iSeed)-Half)*DelFi/Ten
CB = cos(B)
SB = sin(B)
if (A*Three <= One) then !make it random whether we rotate around x, y or z.
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cy = Cordst(2,ii)
      Cz = Cordst(3,ii)
      Cordst(2,ii) = Cy*CB+Cz*SB
      Cordst(3,ii) = Cz*CB-Cy*SB
    end do
  end do
else if (A*Three <= Two) then
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cx = Cordst(1,ii)
      Cz = Cordst(3,ii)
      Cordst(1,ii) = Cx*CB+Cz*SB
      Cordst(3,ii) = CB*Cz-SB*Cx
    end do
  end do
else
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cx = Cordst(1,ii)
      Cy = Cordst(2,ii)
      Cordst(1,ii) = Cx*CB+Cy*SB
      Cordst(2,ii) = CB*Cy-SB*Cx
    end do
  end do
end if
!----------------------------------------------------------------------*
! Generate the image points that correspond with the new coordinates.  *
! We follow Friedman. Since no image is created here for the qm-       *
! molecule, start with query.                                          *
!----------------------------------------------------------------------*
Ind = 0 !The SM-defaults.
iImage = 1
if (Qmeq .or. QmProd) then
  Ind = iCNum*nCent !Makes sure that the first slots in CordIm are empty
  iImage = iSta
end if
iQsta = nCent-nCha+1
A2 = Ract**2
DiFac = -(Diel-One)/(Diel+One)
xyzMyp(:) = Zero
do i=iImage,nPart
  do j=1,nCent
    Ind = Ind+1
    S2 = Zero
    do k=1,3
      S2 = S2+Cordst(k,Ind)**2
    end do
    S2 = A2/S2
    Sqrts2 = sqrt(S2)
    Sqrs(Ind) = Sqrts2
    if (j <= nPol) then
      QImp(Ind) = Zero
      DipIm(:,Ind) = Zero
    end if
    if (j >= iQsta) then
      qq = Qsta(j-nCent+nCha)
      q = DiFac*Sqrts2*qq
      Qim(Ind) = q
    else
      qq = Zero
      Qim(Ind) = Zero
    end if
    xyzMyp(:) = xyzMyp-qq*Cordst(:,Ind) !Total dipole of the cavity; used in polink.
    CordIm(:,Ind) = Cordst(:,Ind)*S2
  end do
end do

!----------------------------------------------------------------------*
! Exit.                                                                *
!----------------------------------------------------------------------*
return

end subroutine Geogen
