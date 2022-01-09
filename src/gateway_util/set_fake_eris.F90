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

subroutine Set_Fake_ERIs

use Basis_Info, only: nBas
use RICD_Info, only: Do_RI, Cholesky
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
#include "cholesky.fh"
character(LEN=16) NamRfil
integer, dimension(:), allocatable :: iSOShl
integer nVec_RI(8)

write(6,*)
write(6,*) '   *** Skipping anything related to ERIs ***'
write(6,*)

if (.not. (Cholesky .or. Do_RI)) return

call Get_NameRun(NamRfil)
call NameRun('AUXRFIL')

call Get_iScalar('ChoVec Address',CHO_ADRVEC)
nBasT = nBas(0)
do i=1,nIrrep-1
  nBasT = nBasT+nBas(i)
end do
call mma_allocate(iSOShl,nBasT)
call Get_dScalar('Cholesky Threshold',THRCOM)
call Get_iArray('NumCho',NumCho,nIrrep)
call Get_iArray('nVec_RI',nVec_RI,nIrrep)
call Get_iArray('iSOShl',ISOSHL,NBAST)

call NameRun(NamRfil)
call Put_iArray('iSOShl',ISOSHL,NBAST)
call mma_deallocate(iSOShl)
call Put_iArray('NumCho',NumCho,nIrrep)
call Put_iArray('nVec_RI',nVec_RI,nIrrep)
call Put_iScalar('ChoVec Address',CHO_ADRVEC)
call Put_dScalar('Cholesky Threshold',THRCOM)

return

end subroutine Set_Fake_ERIs
