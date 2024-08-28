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

subroutine PrMtrx(Label,lOper,nComp,ip,Matrix)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden, January '91                  *
!***********************************************************************

use Basis_Info, only: nBas
use Gateway_global, only: PrPrt
use Symmetry_Info, only: nIrrep
use Definitions, only: u6

implicit none
integer nComp
character(len=*) Label
real*8 Matrix(*)
integer ip(nComp), lOper(nComp)
character(len=80) Line
logical type
integer iComp, ip1, iSmLbl, iIrrep, jIrrep

do iComp=1,nComp
  ip1 = ip(iComp)
  iSmLbl = lOper(iComp)
  if (Prprt) iSmLbl = iand(1,iSmLbl)
  type = .true.
  do iIrrep=0,nIrrep-1
    if (nBas(iIrrep) <= 0) cycle
    do jIrrep=0,iIrrep
      if (nBas(jIrrep) <= 0) cycle
      if (iand(iSmLbl,2**ieor(iIrrep,jIrrep)) == 0) cycle
      if (type) then
        type = .false.
        write(u6,*)
        write(u6,*)
        write(u6,'(A,A,A,I2)') ' SO Integrals of type ',Label,' Component ',iComp
      end if
      Line = ''
      if (iIrrep == jIrrep) then
        write(Line,'(1X,A,I1)') ' Diagonal Symmetry Block ',iIrrep+1
        call TriPrt(Line,' ',Matrix(ip1),nBas(iIrrep))
        ip1 = ip1+nBas(iIrrep)*(nBas(iIrrep)+1)/2
      else
        write(Line,'(1X,A,I1,A,I1)') ' Off-diagonal Symmetry Block ',iIrrep+1,',',jIrrep+1
        call RecPrt(Line,' ',Matrix(ip1),nBas(iIrrep),nBas(jIrrep))
        ip1 = ip1+nBas(iIrrep)*nBas(jIrrep)
      end if
    end do
  end do
end do

return

end subroutine PrMtrx
