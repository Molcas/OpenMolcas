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
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine PrMx2(Label,iComp,lOper,result,Mem)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!                                                                      *
!     Modified by AB 950620                                            *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep, iOper

implicit real*8(A-H,O-Z)
! Local arrays
real*8 result(Mem)
character Label*(*), Line*80
integer lOper
logical type

ip1 = 1
type = .true.
do iIrrep=0,nIrrep-1
  if (nBas(iIrrep) <= 0) cycle
  do jIrrep=0,iIrrep
    lop = ieor(iOper(iIrrep),iOper(jIrrep))
    if (lop /= loper) cycle
    if (nBas(jIrrep) <= 0) cycle
    if (type) then
      type = .false.
      write(6,*)
      write(6,*)
      write(6,'(A,A,A,I2)') ' SO Integral gradients of the ',Label,' Component ',iComp
    end if
    Line = ''
    if (iIrrep == jIrrep) then
      write(Line,'(1X,A,I1)') ' Diagonal Symmetry Block ',iIrrep+1
      !call TriPrt(Line,' ',Result(ip1),nBas(iIrrep))
      ip1 = ip1+nBas(iIrrep)*(nBas(iIrrep)+1)/2
    else
      write(Line,'(1X,A,I1,A,I1)') ' Off-diagonal Symmetry Block ',iIrrep+1,',',jIrrep+1
      call RecPrt(Line,' ',result(ip1),nBas(iIrrep),nBas(jIrrep))
      ip1 = ip1+nBas(iIrrep)*nBas(jIrrep)
    end if
  end do
end do

return

end subroutine PrMx2
