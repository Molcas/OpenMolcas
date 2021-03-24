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

subroutine Fluctuating(AInv,nAtoms,Lambda,dQ,nij,nPert,iANr,rMP,nElem,EC,Alpha)

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nij, nPert, iANr(nAtoms), nElem
real(kind=wp), intent(in) :: AInv(nAtoms,nAtoms), EC(3,nij), Alpha
real(kind=wp), intent(out) :: Lambda(nAtoms), dQ(nAtoms)
real(kind=wp), intent(inout) :: rMP(nij,0:nElem-1,0:nPert-1)
integer(kind=iwp) :: iAtom, ii, ij, iPert, jAtom, jj
real(kind=wp) :: A(3), B(3), R_BS_i, R_BS_j, ri, rij02, rij2, rj
real(kind=wp), external :: Bragg_Slater

do iPert=1,6
  do iAtom=1,nAtoms
    ii = iAtom*(iAtom+1)/2
    dQ(iAtom) = rMP(ii,0,0)-rMP(ii,0,iPert)
  end do
  !call RecPrt('dQ',' ',dQ,1,nAtoms)
  call DGEMM_('N','N',nAtoms,1,nAtoms,One,AInv,nAtoms,dQ,nAtoms,Zero,Lambda,nAtoms)
  !call RecPrt('Lambda',' ',Lambda,1,nAtoms)

  do iAtom=1,nAtoms
    R_BS_i = Bragg_Slater(iANr(iAtom))
    ii = iAtom*(iAtom+1)/2
    A(:) = EC(:,ii)

    do jAtom=1,iAtom-1
      R_BS_j = Bragg_Slater(iANr(jAtom))
      jj = jAtom*(jAtom+1)/2
      B(:) = EC(:,jj)
      rij2 = (A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
      ri = Lambda(iAtom)
      rj = Lambda(jAtom)
      rij02 = ((R_BS_i+R_BS_j))**2
      ij = iAtom*(iAtom-1)/2+jAtom
      !write(u6,*) ij,ri,rj,rij02
      rMP(ij,0,iPert) = -(ri-rj)*exp(-Alpha*(rij2/rij02))/Two
      !write(u6,*) 'RMP',iAtom,jAtom,rMP(ij,0,iPert)
    end do

  end do

end do
!call RecPrt('rMP',' ',rMP,nij,nElem*nPert)

return

end subroutine Fluctuating
