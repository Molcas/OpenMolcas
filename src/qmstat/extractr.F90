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

subroutine ExtractR(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nState,HMatOld,xyzQu,lExtr,iExtr_Eig,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iLu, i9, iC, iDt(3), nState, iExtr_Eig, iExtr_Atm(*), ip_ExpVal, ip_ExpCento
real(kind=wp) :: Etot, xyzMy(3), Hmat(*), HMatOld(*), xyzQu(6), ENR, ENP
logical(kind=iwp) :: lExtr(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ind, j, k, nDim

write(iLu,*) '<<<<<<<Configuration ',i9,'>>>>>>>'
if (lExtr(1)) then
  write(iLu,*) 'Total Energy'
  write(iLu,'(F15.8)') Etot
end if
if (lExtr(2)) then
  write(iLu,*) 'QM-Dipole'
  write(iLu,'(3(F12.5))') (xyzMy(k),k=1,3)
end if
if (lExtr(3)) then
  write(iLu,*) 'QM-Quadrupole'
  write(iLu,'(6(F12.5))') (xyzQu(k),k=1,6)
end if
if (lExtr(4)) then
  write(iLu,*) 'Eigenvalues of RASSI-matrix'
  do i=1,iExtr_Eig
    ind = nTri_Elem(i)
    write(iLu,'(F15.8)') Hmat(ind)
  end do
end if
if (lExtr(5)) then
  write(iLu,*) 'Corresponding eigenvectors'
  do j=0,iExtr_Eig-1
    write(iLu,'(5(F15.8))') (Work(iC+j*nState+k),k=0,nState-1)
  end do
end if
if (lExtr(6)) then
  write(iLu,*) 'Expectation values (H_0,V_el,V_pol,V_pp)'
  write(iLu,*) '  Nuc cont:',ENR
  if (lExtr(4)) nDim = iExtr_Eig
  if (.not. lExtr(4)) nDim = nState
  do i=1,nDim
    write(iLu,'(4(F15.8))') (Work(ip_ExpVal+4*(i-1)+k),k=0,3)
  end do
  call GetMem('ExpVals','Free','Real',ip_ExpVal,4*nDim)
end if
if (lExtr(7)) then
  write(iLu,*) 'Expectation values partial V_el, V_pol'
  write(iLu,*) '  Nuc cont:',ENP
  if (lExtr(4)) nDim = iExtr_Eig
  if (.not. lExtr(4)) nDim = nState
  do j=1,nDim
    write(iLu,'(2(F15.8))') (Work(ip_ExpCento+4*(j-1)+k),k=1,2)
  end do
  call GetMem('ExpVals','Free','Real',ip_ExpCento,4*nDim)
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(iDt)
  call Unused_real_array(HMatOld)
  call Unused_integer_array(iExtr_Atm)
end if

end subroutine ExtractR
