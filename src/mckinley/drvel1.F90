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

subroutine Drvel1(Grad)

use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
logical, external :: TF
real*8 Grad(*)

idisp = 0
do jIrrep=0,nirrep-1
  do Jcar=1,3
    iirrep = irrfnc(2**(jcar-1))
    if (jirrep == iirrep) then
      mdc = 0
      do iCnttp=1,nCnttp
        ZA = dbsc(iCnttp)%Charge
        do iCnt=1,dbsc(iCnttp)%nCntr
          mdc = mdc+1
          do iCar=1,3
            iComp = 2**(iCar-1)
            if (TF(mdc,jIrrep,iComp)) then
              idisp = idisp+1
              if (icar == jcar) Grad(idisp) = ZA
            end if
          end do
        end do
      end do
    end if
  end do
end do

return

end subroutine Drvel1
