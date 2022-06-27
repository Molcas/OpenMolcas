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

subroutine FndBnd(IOut,AlBond,MxBond,NAtoms,IAn,C,Nbond,IBond,IBType,PBO)
! Generate connectivity based on bond distances alone.  The criteria
! are contained in routine IPBO.
! Cartesian coords. are in Angstroms

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut, MxBond, NAtoms, IAn(NAtoms)
logical(kind=iwp), intent(in) :: AlBond
integer(kind=iwp), intent(out) :: Nbond(NAtoms), IBond(MxBond,NAtoms), IBType(MxBond,NAtoms)
real(kind=wp), intent(in) :: C(3,NAtoms)
real(kind=wp), intent(out) :: PBO(MxBond,NAtoms)
integer(kind=iwp) :: I, IbondO, J
real(kind=wp) :: BondOr, RIJ
integer(kind=iwp), external :: IPBO

IBond(:,:) = 0
IBType(:,:) = 0
BondOr = Zero
NBond(:) = 0
do I=1,NAtoms
  do J=1,NAtoms
    if (J == I) cycle
    RIJ = sqrt((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2+(C(3,I)-C(3,J))**2)
    IBondO = IPBO(IAn(I),IAn(J),RIJ,BondOr)
    if ((IBondO > 0) .or. AlBond) then
      NBond(I) = NBond(I)+1
      if (NBond(I) > MxBond) then
        write(IOut,1000) MxBond,I
        call Abend()
      end if
      IBond(NBond(I),I) = J
      IBType(NBond(I),I) = IBondO
      PBO(NBond(I),I) = BondOr
    end if
  end do
end do

return

1000 format(' Maximum number of bonds=',I3,' exceeded for atom',I4,'.')

end subroutine FndBnd
