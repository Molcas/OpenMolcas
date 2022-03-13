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

subroutine Extract(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld,xyzQuQ,ip_ExpVal,ip_ExpCento,ENR,ENP)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
dimension xyzMy(3), Hmat(*), HMatOld(*), xyzQuQ(6)
dimension iDt(3)

! Just pass on the numbers according to QM-method.

if (QmType(1:4) == 'RASS') then
  call ExtractR(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld,xyzQuQ,lExtr,iExtr_Eig,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)
else if (QmType(1:3) == 'SCF') then
  call ExtractS(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,xyzQuQ,lExtr,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)
end if

return

end subroutine Extract
