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
!  OrbRot2
!
!> @brief
!>   Rotate orbitals of the solvent
!> @author A. Ohrn
!>
!> @details
!> Given a rotation matrix \f$ R \f$ it is easy to rotate the p-orbitals since
!> they behave just like the three axes. The d-orbitals are more
!> intricate.
!>
!> To obtain a representation
!> of the d-orbital \f$ \mathrm{d}_{xy} \f$ we perform the outer product
!> \f$ \mathrm{d}_x \mathrm{d}_y + \mathrm{d}_y \mathrm{d}_x \f$,
!> which generates a two-dimensional matrix. Now rotate each px and py
!> and perform this multiplication. We obtain a symmetric matrix from
!> whose elements we can compute how the \f$ \mathrm{d}_{xy} \f$ transforms when
!> rotated. The same is done for the other d-orbitals and we can
!> construct a transformation matrix, which we apply to the MO-coeff.
!> Other ways to look at it are available, but the present one can be
!> seen as a generation of table 3 in \cite Iva1996-JPC-100-6342.
!> The present method is efficient
!> with no trigonometric functions and ample use of BLAS_UTIL.
!>
!> @note
!> ::Qfread as well as ::transrot must precede.
!>
!> @param[in]     Rot  The rotation matrix
!> @param[in,out] Cmo  The MO-coefficients
!> @param[in]     iQ   The angular type of the \f$ i \f$ -th basis
!>                     (observe, not the \f$ i \f$ -th basis function, see givemeinfo)
!> @param[in]     iOrb Number of orbitals
!> @param[in]     nBas Number of basis functions
!> @param[in]     lMax Number of bases (not basis functions), see ::qfread
!> @param[in]     nCnC Number of contracted basis functions of same type as the \f$ i \f$ -th basis.
!>                     For example 7s4p will have vector ``7,7,7,7,7,7,7,4,4,4,4``.
!***********************************************************************

subroutine OrbRot2(Rot,Cmo,iQ,iOrb,nBas,lMax,nCnC)

use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Rot(3,3)
integer(kind=iwp), intent(in) :: lMax, iQ(lMax), iOrb, nBas, nCnC(lMax)
real(kind=wp), intent(inout) :: Cmo(nBas,iOrb)
integer(kind=iwp) :: i, iIn, IqSuckOut, ISkutt, j
real(kind=wp) :: Ctemp1, Ctemp2, Ctemp3, Ctemp4, Ctemp5, DBlock(5,5), Dxmin(6), Dxy(6), Dxz(6), Dyz(6), Dzz(6), PBlock(3,3), r1, &
                 r2, r3, r4, Resx(3), Resy(3), Resz(3)
logical(kind=iwp) :: NewIq
real(kind=wp), parameter :: dUnit(6) = [-One,Zero,Zero,-One,Zero,-One], Px(3) = [One,Zero,Zero], Py(3) = [Zero,One,Zero], &
                            Pz(3) = [Zero,Zero,One]
#include "warnings.h"

! Dgemv multiply the matrix Rot on the vector Px. For details about
! various parameters see the source code for the routine; it is very
! detailed. Here we get the p-orbitals which also are used to obtain
! the transformation of the other orbitals.

call dGeMV_('N',3,3,One,Rot,3,Px,1,Zero,Resx,1)
call dGeMV_('N',3,3,One,Rot,3,Py,1,Zero,Resy,1)
call dGeMV_('N',3,3,One,Rot,3,Pz,1,Zero,Resz,1)

! Construct the p-block. Easy since they are linear polynomials in the cartesian R-matrix.

PBlock(:,1) = Resx
PBlock(:,2) = Resy
PBlock(:,3) = Resz

! Generation of quadratic polynomial of R-matrix elements for d_xy.

Dxy(:) = Zero
call Dspr2('L',3,One,Resx,1,Resy,1,Dxy)

! Generation of quadratic polynomial of R-matrix elements for d_xz.

Dxz(:) = Zero
call Dspr2('L',3,One,Resx,1,Resz,1,Dxz)

! Generation of quadratic polynomial of R-matrix elements for d_yz.

Dyz(:) = Zero
call Dspr2('L',3,One,Resy,1,Resz,1,Dyz)

! Generation of quadratic polynomial of R-matrix elements for d_xx-yy.

Dxmin(:) = Zero
call Dspr('L',3,One,Resx,1,Dxmin)
call Dspr('L',3,-One,Resy,1,Dxmin)

! Generation of quadratic polynomial of R-matrix elements for d_3zz-1.

Dzz(:) = dUnit
call Dspr('L',3,Three,Resz,1,Dzz)

! Construct the d-block. This requires knowledge of how Molcas
! handles d-orbitals, especially their normalization, which is the
! reason for the constants r1,r2,r3,r4.

r1 = One
r2 = sqrt(Three)
r3 = One
r4 = One/sqrt(Three)
DBlock(1,1) = r1*Dxy(2)                    !dxy-->dxy
DBlock(2,1) = r1*Dxy(5)                    !dxy-->dyz
DBlock(3,1) = (r1/r4)*Half*Dxy(6)          !dxy-->d3zz-1
DBlock(4,1) = r1*Dxy(3)                    !dxy-->dxz
DBlock(5,1) = (r1/r3)*(Dxy(1)+Half*Dxy(6)) !dxy-->dxx-yy
DBlock(1,2) = r1*Dyz(2)                    !dyz-->dxy
DBlock(2,2) = r1*Dyz(5)                    !dyz-->dyz
DBlock(3,2) = (r1/r4)*Half*Dyz(6)          !dyz-->d3zz-1
DBlock(4,2) = r1*Dyz(3)                    !dyz-->dxz
DBlock(5,2) = (r1/r3)*(Dyz(1)+Half*Dyz(6)) !dyz-->dxx-yy
DBlock(1,3) = r4*Dzz(2)                    !d3zz-->dxy
DBlock(2,3) = r4*Dzz(5)                    !d3zz-->dyz
DBlock(3,3) = r1*Half*Dzz(6)               !d3zz-->d3zz
DBlock(4,3) = r4*Dzz(3)                    !d3zz-->dxz
DBlock(5,3) = (r1/r2)*(Dzz(1)+Half*Dzz(6)) !d3zz-->dxx-yy
DBlock(1,4) = r1*Dxz(2)                    !dxz-->dxy
DBlock(2,4) = r1*Dxz(5)                    !dxz-->dyz
DBlock(3,4) = (r1/r4)*Half*Dxz(6)          !dxz-->d3zz-1
DBlock(4,4) = r1*Dxz(3)                    !dxz-->dxz
DBlock(5,4) = (r1/r3)*(Dxz(1)+Half*Dxz(6)) !dxz-->dxx-yy
DBlock(1,5) = r3*Dxmin(2)                  !dxx-yy-->dxy
DBlock(2,5) = r3*Dxmin(5)                  !dxx-yy-->dyz
DBlock(3,5) = r2*Half*Dxmin(6)             !dxx-yy-->d3zz
DBlock(4,5) = r3*Dxmin(3)                  !dxx-yy-->dxz
DBlock(5,5) = r1*(Dxmin(1)+Half*Dxmin(6))  !dxx-yy-->dxx-yy

! With the proper number of blocks at hand, we make transformation.
! The brute force way is to construct the entire transformation matrix
! and make a matrix multiplication. Since so many elements in that
! matrix are zero, more effecient ways to transform are available.
! That is what is used below. The formulas follow from considering the
! transformation matrix multiplied with the CMO-matrix. Multio
! importante is it to observe how the basis functions in Molcas are
! ordered!

do i=1,iOrb
  ! OBSERVE!!! WE ARE ASSUMING THAT THE FIRST BASIS IS OF S-TYPE!!! IF YOU ARE IMPLEMENTING SYMMETRY, THIS
  ! MIGHT NOT BE A VALID ASSUMPTION SO THEN THE NEWIQ-CONSTRUCT BELOW MUST BE ALTERED!!!
  iIn = 1
  do j=2,lMax
    iIn = iIn+1
    IqSuckOut = iQ(j)
    NewIq = iQ(j) /= iQ(j-1)
    !This if-clause controls the jumping when new angular basis function appears.
    if (Newiq) iIn = iIn+(2*iQ(j-1)-2)*nCnC(j-1)
    select case (iqSuckOut)
      case (1) !This is s-function
      case (2) !This is p-function
        iSkutt = nCnC(j)
        Ctemp1 = PBlock(1,1)*Cmo(iIn,i)+PBlock(1,2)*Cmo(iIn+iSkutt,i)+PBlock(1,3)*Cmo(iIn+2*iSkutt,i)
        Ctemp2 = PBlock(2,1)*Cmo(iIn,i)+PBlock(2,2)*Cmo(iIn+iSkutt,i)+PBlock(2,3)*Cmo(iIn+2*iSkutt,i)
        Ctemp3 = PBlock(3,1)*Cmo(iIn,i)+PBlock(3,2)*Cmo(iIn+iSkutt,i)+PBlock(3,3)*Cmo(iIn+2*iSkutt,i)
        Cmo(iIn,i) = Ctemp1
        Cmo(iIn+iSkutt,i) = Ctemp2
        Cmo(iIn+2*iSkutt,i) = Ctemp3
      case (3) !This is d-function
        iSkutt = nCnC(j)
        Ctemp1 = DBlock(1,1)*Cmo(iIn,i)+DBlock(1,2)*Cmo(iIn+iSkutt,i)+DBlock(1,3)*Cmo(iIn+2*iSkutt,i)+ &
                 DBlock(1,4)*Cmo(iIn+3*iSkutt,i)+DBlock(1,5)*Cmo(iIn+4*iSkutt,i)
        Ctemp2 = DBlock(2,1)*Cmo(iIn,i)+DBlock(2,2)*Cmo(iIn+iSkutt,i)+DBlock(2,3)*Cmo(iIn+2*iSkutt,i)+ &
                 DBlock(2,4)*Cmo(iIn+3*iSkutt,i)+DBlock(2,5)*Cmo(iIn+4*iSkutt,i)
        Ctemp3 = DBlock(3,1)*Cmo(iIn,i)+DBlock(3,2)*Cmo(iIn+iSkutt,i)+DBlock(3,3)*Cmo(iIn+2*iSkutt,i)+ &
                 DBlock(3,4)*Cmo(iIn+3*iSkutt,i)+DBlock(3,5)*Cmo(iIn+4*iSkutt,i)
        Ctemp4 = DBlock(4,1)*Cmo(iIn,i)+DBlock(4,2)*Cmo(iIn+iSkutt,i)+DBlock(4,3)*Cmo(iIn+2*iSkutt,i)+ &
                 DBlock(4,4)*Cmo(iIn+3*iSkutt,i)+DBlock(4,5)*Cmo(iIn+4*iSkutt,i)
        Ctemp5 = DBlock(5,1)*Cmo(iIn,i)+DBlock(5,2)*Cmo(iIn+iSkutt,i)+DBlock(5,3)*Cmo(iIn+2*iSkutt,i)+ &
                 DBlock(5,4)*Cmo(iIn+3*iSkutt,i)+DBlock(5,5)*Cmo(iIn+4*iSkutt,i)
        Cmo(iIn,i) = Ctemp1
        Cmo(iIn+iSkutt,i) = Ctemp2
        Cmo(iIn+2*iSkutt,i) = Ctemp3
        Cmo(iIn+3*iSkutt,i) = Ctemp4
        Cmo(iIn+4*iSkutt,i) = Ctemp5
      case default !Here we go if non-implemented angular quantum number appears.
        write(u6,*)
        write(u6,*) ' ERROR in OrbRot2. Not ready for f-orbitals'
        call Quit(_RC_GENERAL_ERROR_)
    end select
  end do
end do

return

end subroutine OrbRot2
