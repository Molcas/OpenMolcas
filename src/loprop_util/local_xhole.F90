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

subroutine Local_XHole(ipXHole2,dMolExpec,nAtoms,nBas1,nTemp,iCenter,Ttot,Ttot_Inv,Coor,nij,EC,iANr,Bond_Threshold,iPrint,XHLoc2)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
dimension Coor(3,nAtoms), Sq_Temp(nTemp), Ttot_Inv(nTemp)
dimension Temp(nTemp), A(3), B(3), d2Loc(nij), EC(3,nij)
dimension Ttot(nTemp), XHLoc2(nij)
dimension iCenter(nBas1), iANr(nAtoms)
logical Found
real*8, allocatable :: Dens(:)

! Binomal stuff

call Set_Binom()

! Transform the density matrix to the LoProp basis

call Qpg_dArray('D1ao',Found,nDenno)
if (Found .and. (nDenno /= 0)) then
  call mma_Allocate(Dens,nDenno,Label='Dens')
else
  write(6,*) 'Local XHole: D1ao not found.'
  call Abend()
end if
call Get_D1ao(Dens,nDenno)
call DSq(Dens,Sq_Temp,1,nBas1,nBas1)
call mma_deallocate(Dens)

call DGEMM_('N','T',nBas1,nBas1,nBas1,1.0d0,Sq_Temp,nBas1,Ttot_Inv,nBas1,0.0d0,Temp,nBas1)
call DGEMM_('N','N',nBas1,nBas1,nBas1,1.0d0,Ttot_Inv,nBas1,Temp,nBas1,0.0d0,Sq_Temp,nBas1)

! Transform the exchange-hole stuff to LoProp basis

call Transmu(Work(ipXHole2),nBas1,Ttot,Temp)

! Now localize

do iAtom=1,nAtoms
  call dcopy_(3,Coor(1,iAtom),1,A,1)
  do jAtom=1,iAtom
    call dcopy_(3,Coor(1,jAtom),1,B,1)

    ! Sum up contributions to the domain ij

    Acc = Zero
    iOffO = ipXHole2-1
    do j=1,nBas1
      do i=1,nBas1
        if (((iCenter(i) == iAtom) .and. (iCenter(j) == jAtom)) .or. ((iCenter(i) == jAtom) .and. (iCenter(j) == iAtom))) then
          ij = (j-1)*nBas1+i
          Acc = Acc+Sq_Temp(ij)*Work(ij+iOffO)
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
call dcopy_(nij,d2Loc,1,XHLoc2,1)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(dMolExpec)
  call Unused_integer(iPrint)
end if

end subroutine Local_XHole
