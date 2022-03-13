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
!> @param[in,out] Ract  The dielectric radius on input and the slightly perturbed radius on output
!> @param[out]    Rold  Stores the input dielectric radius
!> @param[in]     iCNum How many solvent places that are taken up by the QM-molecule
!***********************************************************************

subroutine Geogen(Ract,Rold,iCNum,iQ_Atoms)

use Constants, only: Zero, One, Two, Three, Ten, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Ract, Rold
integer(kind=iwp) :: iCNum, iQ_Atoms
#include "maxi.fh"
#include "qmcom.fh"
#include "qminp.fh"
integer(kind=iwp) :: i, iAt, ii, iImage, ij, Ind, iQsta, j, k
real(kind=wp) :: A, A2, B, CB, Cx, Cy, Cz, DiFac, Dq(3), Dx, q, qq, S2, SB, Sqrts2, x, xNy, y, yNy, z, zNy
real(kind=wp), external :: Ranf

!----------------------------------------------------------------------*
! Store old configuration.                                             *
!----------------------------------------------------------------------*
do i=1,nPart*nCent
  do j=1,3
    OldGeo(i,j) = Cordst(i,j)
  end do
end do
!----------------------------------------------------------------------*
! Query which type of simulation this is, and if quantum then change   *
! the quantum molecule.                                                *
!----------------------------------------------------------------------*
if (Qmeq .or. QmProd) then !Which coordinates to keep fixed.
  iSta = iCNum+1  !This sees to that the QM-molecule is excluded from the moves below.
  Dq(1) = delX*(ranf(iseed)-Half)
  Dq(2) = delX*(ranf(iseed)-Half)
  Dq(3) = delX*(ranf(iseed)-Half)
  do iAt=1,iQ_Atoms
    do ii=1,3
      Cordst(iAt,ii) = Cordst(iAt,ii)+Dq(ii) !Move QM-mol.
    end do
  end do
end if
!----------------------------------------------------------------------*
! Obtain the random-stuff and make small geometry change.              *
!----------------------------------------------------------------------*
Rold = Ract
Ract = Ract+(ranf(iseed)-Half)*DelR !Change in cavity radius
do i=iSta,nPart !Which molecules to give new coordinates.
  ij = (i-1)*nCent
  do j=1,3
    Dx = DelX*(ranf(iseed)-Half)
    do k=1,nCent
      ii = ij+k
      Cordst(ii,j) = Cordst(ii,j)+Dx !Make translation
    end do
  end do
  Cx = Cordst(ij+1,1) !The oxygen, around which we rotate
  Cy = Cordst(ij+1,2)
  Cz = Cordst(ij+1,3)
  B = (ranf(iseed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent !Rotate around the oxygen in yz-plane, i.e. around x-axis.
    y = Cordst(ij+k,2)-Cy
    z = Cordst(ij+k,3)-Cz
    yNy = y*CB+z*SB !This is a rotation matrix
    zNy = z*CB-y*SB
    Cordst(ij+k,2) = yNy+Cy
    Cordst(ij+k,3) = zNy+Cz
  end do
  B = (ranf(iseed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent !And now rotate in xz-plane
    x = Cordst(ij+k,1)-Cx
    z = Cordst(ij+k,3)-Cz
    xNy = x*CB+z*SB
    zNy = z*CB-x*SB
    Cordst(ij+k,1) = xNy+Cx
    Cordst(ij+k,3) = zNy+Cz
  end do
  B = (ranf(iseed)-Half)*DelFi
  CB = cos(B)
  SB = sin(B)
  do k=2,nCent  !To your surprise, here we rotate in the xy-plane
    x = Cordst(ij+k,1)-Cx
    y = Cordst(ij+k,2)-Cy
    xNy = x*CB+y*SB
    yNy = y*CB-x*SB
    Cordst(ij+k,1) = xNy+Cx
    Cordst(ij+k,2) = yNy+Cy
  end do
end do
! Here all other water molecules rotate around one of the three axes,
! except the ones we fix, which in a QM-simulation is the quantum particle.
A = ranf(iseed)
B = (ranf(iseed)-Half)*DelFi/Ten
CB = cos(B)
SB = sin(B)
if (A*Three <= One) then !make it random whether we rotate around x, y or z.
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cy = Cordst(ii,2)
      Cz = Cordst(ii,3)
      Cordst(ii,2) = Cy*CB+Cz*SB
      Cordst(ii,3) = Cz*CB-Cy*SB
    end do
  end do
else if (A*Three <= Two) then
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cx = Cordst(ii,1)
      Cz = Cordst(ii,3)
      Cordst(ii,1) = Cx*CB+Cz*SB
      Cordst(ii,3) = CB*Cz-SB*Cx
    end do
  end do
else
  do i=iSta,nPart
    ij = (i-1)*nCent
    do k=1,nCent
      ii = ij+k
      Cx = Cordst(ii,1)
      Cy = Cordst(ii,2)
      Cordst(ii,1) = Cx*CB+Cy*SB
      Cordst(ii,2) = CB*Cy-SB*Cx
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
DiFac = -(DiEl-One)/(DiEl+One)
do i=1,3
  xyzMyp(i) = Zero
end do
do i=iImage,nPart
  do j=1,nCent
    Ind = Ind+1
    S2 = Zero
    do k=1,3
      S2 = S2+Cordst(Ind,k)**2
    end do
    S2 = A2/S2
    Sqrts2 = sqrt(S2)
    Sqrs(Ind) = Sqrts2
    if (j <= nPol) then
      QImp(Ind) = Zero
      do k=1,3
        dim(Ind,k) = Zero
      end do
    end if
    if (j >= iQsta) then
      qq = Qsta(j-nCent+nCha)
      q = DiFac*Sqrts2*qq
      Qim(Ind) = q
    else
      qq = Zero
      Qim(Ind) = Zero
    end if
    do k=1,3
      xyzMyp(k) = xyzMyp(k)-qq*Cordst(Ind,k) !Total dipole of the cavity; used in polink.
      CordIm(Ind,k) = Cordst(Ind,k)*S2
    end do
  end do
end do

!----------------------------------------------------------------------*
! Exit.                                                                *
!----------------------------------------------------------------------*
return

end subroutine Geogen
