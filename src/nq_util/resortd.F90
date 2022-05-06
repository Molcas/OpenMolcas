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

subroutine ResortD(D_Old,D_New,iBas,iCmp,jBas,jCmp)

implicit real*8(a-h,o-z)
real*8 D_Old(iBas,jBas,iCmp,jCmp), D_New(iBas,iCmp,jBas,jCmp)

do jC=1,jCmp
  do jB=1,jBas
    do iC=1,iCmp
      do iB=1,iBas
        D_New(iB,iC,jB,jC) = D_Old(iB,jB,iC,jC)
      end do
    end do
  end do
end do

return

end subroutine ResortD
