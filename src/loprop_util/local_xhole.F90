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

subroutine Local_XHole(XHole2,nAtoms,nBas1,nTemp,iCenter,Ttot,Ttot_Inv,nij,EC,iANr,Bond_Threshold,XHLoc2)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas1, nTemp, iCenter(nBas1), nij, iANr(nAtoms)
real(kind=wp), intent(inout) :: XHole2(nBas1,nBas1)
real(kind=wp), intent(in) :: Ttot(nTemp), Ttot_Inv(nTemp), EC(3,nij), Bond_Threshold
real(kind=wp), intent(out) :: XHLoc2(nij)
integer(kind=iwp) :: i, iAtom, ij, j, jAtom, nDenno
real(kind=wp) :: Acc, d2Loc(nij)
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Dens(:), Sq_Temp(:,:), Temp(:,:)

! Binomial stuff

call Set_Binom()

! Transform the density matrix to the LoProp basis

call Qpg_dArray('D1ao',Found,nDenno)
if (Found .and. (nDenno /= 0)) then
  call mma_allocate(Dens,nDenno,Label='Dens')
else
  write(u6,*) 'Local XHole: D1ao not found.'
  call Abend()
end if
call Get_D1ao(Dens,nDenno)
call mma_allocate(Sq_Temp,nBas1,nBas1,label='Sq_Temp')
call DSq(Dens,Sq_Temp,1,nBas1,nBas1)
call mma_deallocate(Dens)

call mma_allocate(Temp,nBas1,nBas1,label='Temp')
call DGEMM_('N','T',nBas1,nBas1,nBas1,One,Sq_Temp,nBas1,Ttot_Inv,nBas1,Zero,Temp,nBas1)
call DGEMM_('N','N',nBas1,nBas1,nBas1,One,Ttot_Inv,nBas1,Temp,nBas1,Zero,Sq_Temp,nBas1)

! Transform the exchange-hole stuff to LoProp basis

call Transmu(XHole2,nBas1,Ttot,Temp)
call mma_deallocate(Temp)

! Now localize

do iAtom=1,nAtoms
  do jAtom=1,iAtom

    ! Sum up contributions to the domain ij

    Acc = Zero
    do j=1,nBas1
      do i=1,nBas1
        if (((iCenter(i) == iAtom) .and. (iCenter(j) == jAtom)) .or. ((iCenter(i) == jAtom) .and. (iCenter(j) == iAtom))) then
          Acc = Acc+Sq_Temp(i,j)*XHole2(i,j)
        end if
      end do
    end do
    ij = iAtom*(iAtom-1)/2+jAtom
    d2Loc(ij) = Acc
  end do   ! jAtom
end do   ! iAtom

! Distributes the contributions from the bonds that doesn't fulfill the requirement
! Bond_Length <= Bond_Threshold*(Bragg_Slater(iAtom)+Bragg_Slater(jAtom) to the
! two atoms involved in the bond.

call Move_XHole(d2Loc,EC,nAtoms,nij,iANr,Bond_Threshold)
XHLoc2(:) = d2Loc(:)

return

end subroutine Local_XHole
