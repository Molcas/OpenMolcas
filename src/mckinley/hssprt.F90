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

subroutine HssPrt(Hess,nHess)

use Symmetry_Info, only: nIrrep, lIrrep

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "real.fh"
integer nDisp(0:7)
character Label*39
real*8 Hess(nHess)
real*8, allocatable :: Temp(:)
! Statement function
Ind(idisp,jdisp) = idisp*(idisp-1)/2+jdisp

!                                                                      *
!***********************************************************************
!                                                                      *
iDisp = 0
do iIrrep=0,nIrrep-1
  nDisp(iIrrep) = iDisp
  iDisp = iDisp+lDisp(iIrrep)
end do

if (nirrep == 1) then
  write(Label,100) 'Hessian in Irrep ',lIrrep(0)
  call TriPrt(Label,' ',Hess,ldisp(0))
else
  call mma_allocate(Temp,nHess,Label='Temp')
  do iIrrep=0,nIrrep-1
    write(Label,100) 'Hessian in Irrep ',lIrrep(iIrrep)
    do i=1,lDisp(iirrep)
      do j=1,i
        ii = ind(i,j)
        jj = ind(ndisp(iirrep)+i,ndisp(iirrep)+j)
        Temp(ii) = Hess(jj)
      end do
    end do
    call TriPrt(Label,' ',Temp,ldisp(iirrep))
  end do
  call mma_deallocate(Temp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *

return

100 format(A,A)

end subroutine HssPrt
