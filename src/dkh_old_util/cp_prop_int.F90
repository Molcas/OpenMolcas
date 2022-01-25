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

subroutine Cp_Prop_Int(AInt,nAInt,BInt,nBInt,nBas,nIrrep,Label)
!***********************************************************************
!                                                                      *
! Object: replace the diagonal blocks of the property integrals.       *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 aint(nAInt+4), BInt(nBInt+4)
integer nBas(0:nIrrep-1)

iCmp = 1
iExp = 1
do iIrrep=0,nIrrep-1
  do jIrrep=0,iIrrep
    ij = ieor(iIrrep,jIrrep)
    if (iand(Label,2**ij) == 0) Go To 20
    if (iIrrep == jIrrep) then

      Len = nBas(iIrrep)*(nBas(iIrrep)+1)/2
      do iLen=0,Len-1
        !write(6,*) AInt(iExp+iLen), BInt(iCmp+iLen)
        aint(iExp+iLen) = BInt(iCmp+iLen)
      end do
      iCmp = iCmp+Len
      iExp = iExp+Len
    else
      Len = nBas(iIrrep)*nBas(jIrrep)
      iExp = iExp+Len
    end if
20  continue
  end do
end do

return

end subroutine Cp_Prop_Int
