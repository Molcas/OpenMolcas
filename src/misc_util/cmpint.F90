!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine CmpInt(XInt,nInt,nBas,nIrrep,Label)
!***********************************************************************
!                                                                      *
! Object: to remove the offdiagonal nonzero blocks of matrix elements  *
!         for an operator.                                             *
!                                                                      *
!         XInt(1:nInt):array with nonzero elements                     *
!                                                                      *
!         nBas(0:nIrrep-1):number of basis functions in each irrep     *
!                                                                      *
!         Label: symmetry label of the operator for which the          *
!                matrix elements where computed.                       *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             March 1991                                               *
!***********************************************************************

implicit real*8(A-H,O-Z)
real*8 XInt(nInt+4)
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
        XInt(iLen+iCmp) = Xint(iLen+iExp)
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
do iadd=0,3
  XInt(iCmp+iadd) = XInt(iExp+iadd)
end do
nInt = iCmp-1

return

end subroutine CmpInt
