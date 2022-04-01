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

subroutine Extract(iLu,i9,Etot,xyzMy,Hmat,C,nMatBas,xyzQuQ,ExpVal,ExpCento,ENR,ENP)

use qmstat_global, only: iExtr_Eig, lExtr, QmType
use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iLu, i9, nMatBas
real(kind=wp), intent(in) :: Etot, xyzMY(3), Hmat(nTri_Elem(iExtr_Eig)), C(nMatBas,nMatBas), xyzQuQ(6), ExpVal(4,*), &
                             ExpCento(4,*), ENR, ENP

! Just pass on the numbers according to QM-method.

if (QmType(1:4) == 'RASS') then
  call ExtractR(iLu,i9,Etot,xyzMy,Hmat,C,nMatBas,xyzQuQ,lExtr,iExtr_Eig,ExpVal,ExpCento,ENR,ENP)
else if (QmType(1:3) == 'SCF') then
  call ExtractS(iLu,i9,Etot,xyzMy,xyzQuQ,lExtr,ExpVal,ExpCento,ENR,ENP)
end if

return

end subroutine Extract
