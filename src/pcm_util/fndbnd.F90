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

subroutine FndBnd(IOut,IPrint,AlBond,ToAng,MxBond,NAtoms,IAn,C,Nbond,IBond,IBType,PBO,Re)
! Generate connectivity based on bond distances alone.  The criteria
! are contained in routine IPBO.
! Cartesian coords. are in Angstroms

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut, IPrint, MxBond, NAtoms, IAn(NAtoms)
logical(kind=iwp), intent(in) :: AlBond
integer(kind=iwp), intent(out) :: Nbond(NAtoms), IBond(MxBond,NAtoms), IBType(MxBond,NAtoms)
real(kind=wp), intent(in) :: ToAng, C(3,NAtoms), Re(*)
real(kind=wp), intent(out) :: PBO(MxBond,NAtoms)
integer(kind=iwp) :: I, IbondO, J
real(kind=wp) :: BondOr, RIJ
integer(kind=iwp), external :: IPBO

do I=1,MxBond
  do J=1,NAtoms
    IBond(I,J) = 0
    IBType(I,J) = 0
  end do
end do
BondOr = Zero
do I=1,NAtoms
  NBond(I) = 0
  do J=1,NAtoms
    if (J == I) goto 90
    RIJ = sqrt((C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2+(C(3,I)-C(3,J))**2)
    IBondO = IPBO(ToAng,IAn(I),IAn(J),RIJ,BondOr)
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
90  continue
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(IPrint)
  call Unused_real_array(Re)
end if

1000 format(' Maximum number of bonds=',I3,' exceeded for atom',I4,'.')

end subroutine FndBnd
