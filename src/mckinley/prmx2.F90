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

subroutine PrMx2(Label,iComp,lOper,Rslt,Mem)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!                                                                      *
!     Modified by AB 950620                                            *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep, iOper
use Definitions, only: wp, iwp, u6

implicit none
character(len=*) :: Label
integer(kind=iwp) :: iComp, lOper, Mem
real(kind=wp) :: Rslt(Mem)
integer(kind=iwp) :: iIrrep, ip1, jIrrep, lop
character(len=80) :: Line
logical(kind=iwp) :: ltype

ip1 = 1
ltype = .true.
do iIrrep=0,nIrrep-1
  if (nBas(iIrrep) <= 0) cycle
  do jIrrep=0,iIrrep
    lop = ieor(iOper(iIrrep),iOper(jIrrep))
    if (lop /= loper) cycle
    if (nBas(jIrrep) <= 0) cycle
    if (ltype) then
      ltype = .false.
      write(u6,*)
      write(u6,*)
      write(u6,'(A,A,A,I2)') ' SO Integral gradients of the ',Label,' Component ',iComp
    end if
    Line = ''
    if (iIrrep == jIrrep) then
      write(Line,'(1X,A,I1)') ' Diagonal Symmetry Block ',iIrrep+1
      !call TriPrt(Line,' ',Rslt(ip1),nBas(iIrrep))
      ip1 = ip1+nBas(iIrrep)*(nBas(iIrrep)+1)/2
    else
      write(Line,'(1X,A,I1,A,I1)') ' Off-diagonal Symmetry Block ',iIrrep+1,',',jIrrep+1
      call RecPrt(Line,' ',Rslt(ip1),nBas(iIrrep),nBas(jIrrep))
      ip1 = ip1+nBas(iIrrep)*nBas(jIrrep)
    end if
  end do
end do

return

end subroutine PrMx2
