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

subroutine Set_Fake_ERIs()

use Basis_Info, only: nBas
use RICD_Info, only: Chol => Cholesky, Do_RI
use Cholesky, only: CHO_ADRVEC, NumCho, ThrCom
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, nBasT, nVec_RI(8)
integer(kind=iwp), allocatable :: iSOShl(:)

write(u6,*)
write(u6,*) '   *** Skipping anything related to ERIs ***'
write(u6,*)

if (.not. (Chol .or. Do_RI)) return

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

call NameRun('#Pop')

call Put_iArray('iSOShl',ISOSHL,NBAST)
call mma_deallocate(iSOShl)
call Put_iArray('NumCho',NumCho,nIrrep)
call Put_iArray('nVec_RI',nVec_RI,nIrrep)
call Put_iScalar('ChoVec Address',CHO_ADRVEC)
call Put_dScalar('Cholesky Threshold',THRCOM)

return

end subroutine Set_Fake_ERIs
