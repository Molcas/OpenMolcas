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

subroutine read_dmrg_parameter_for_mclr()

use input_mclr, only: ERASSCF
use dmrginfo, only: DoDMRG, DoMCLR, LRRAS2, MS2_RGLR, nEle_RGLR, nStates_RGLR, RGRAS2
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, ierr

open(unit=100,file='dmrg_for_mclr.parameters',status='OLD',action='READ',iostat=ierr)
if (ierr /= 0) then
  doDMRG = .false.
  doMCLR = .false.
else
  read(100,'(11X,L1,4X)') doDMRG
  read(100,'(4X,I8,4X)') nele_RGLR
  read(100,'(4X,I8,4X)') ms2_RGLR
  !write(u6,*) doDMRG,dmrg_state%nactel,dmrg_state%ms2
  do i=1,8
    read(100,'(4X,I3)',advance='no') RGras2(i)
  end do
  read(100,*)
  do i=1,8
    read(100,'(4X,I3)',advance='no') LRras2(i)
  end do
  !write(u6,*) RGras2
  !write(u6,*) LRras2
  read(100,*)
  read(100,'(4X,I8,4X)') nstates_RGLR
  !call mma_allocate(checkpoint,nstates_RGLR,label='checkpoint')
  !checkpoint = ''
  do i=1,nstates_RGLR
    read(100,*)
    read(100,'(G20.12)') ERASSCF(i)
    write(6,*) 'RASSCF energy',ERASSCF(i)
  end do
  ! It is redundant
  doMCLR = .true.
end if
close(100)

write(u6,*) 'doDMRG, nele_dmrg, ms2_dmrg'
write(u6,*) doDMRG,nele_rglr,ms2_rglr
call xflush(u6)

end subroutine read_dmrg_parameter_for_mclr
