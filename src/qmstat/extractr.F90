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

subroutine ExtractR(iLu,i9,Etot,xyzMy,Hmat,C,nState,xyzQu,lExtr,iExtr_Eig,ExpVal,ExpCento,ENR,ENP)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iLu, i9, nState, iExtr_Eig
real(kind=wp) :: Etot, xyzMy(3), Hmat(nTri_Elem(iExtr_Eig)), C(nState,nState), xyzQu(6), ExpVal(4,nState), ExpCento(4,nState), &
                 ENR, ENP
logical(kind=iwp) :: lExtr(*)
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
  do j=1,iExtr_Eig
    write(iLu,'(5(F15.8))') C(:,j)
  end do
end if
if (lExtr(6)) then
  write(iLu,*) 'Expectation values (H_0,V_el,V_pol,V_pp)'
  write(iLu,*) '  Nuc cont:',ENR
  if (lExtr(4)) then
    nDim = iExtr_Eig
  else
    nDim = nState
  end if
  do i=1,nDim
    write(iLu,'(4(F15.8))') ExpVal(:,i)
  end do
end if
if (lExtr(7)) then
  write(iLu,*) 'Expectation values partial V_el, V_pol'
  write(iLu,*) '  Nuc cont:',ENP
  if (lExtr(4)) then
    nDim = iExtr_Eig
  else
    nDim = nState
  end if
  do j=1,nDim
    write(iLu,'(2(F15.8))') ExpCento(:,i)
  end do
end if

return

end subroutine ExtractR
