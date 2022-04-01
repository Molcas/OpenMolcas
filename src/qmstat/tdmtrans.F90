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

subroutine TdmTrans(nBas)

use qmstat_global, only: iPrint, MoAveRed, MxSymQ, nRedMO, NrFiles, NrStates, nState, RassiM, EigV
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas(MxSymQ)
integer(kind=iwp) :: i, iDisk, j, kaunt, Lu, nStatePrim
character(len=6) :: TDMchar
logical(kind=iwp) :: Exists
real(kind=wp), allocatable :: Eig2(:,:), NonH(:), NonS(:)
#include "warnings.h"

! Sag hej till publiken.

write(u6,*)
write(u6,*) '     Transforming the transition density matrices.'

! Inquire if the ToFile is in WorkDir.

call f_Inquire(RassiM,Exists)
if (.not. Exists) then
  write(u6,*)
  write(u6,*) 'No Transition density matrix file found.'
  write(u6,*) 'Did you use the TOFIle keyword in RASSI?'
  call Quit(_RC_IO_ERROR_READ_)
end if
call f_Inquire(EigV,Exists)
if (.not. Exists) then
  write(u6,*)
  write(u6,*) 'No Rassi eigenvectors found.'
  write(u6,*) 'Did you use the TOFIle keyword in RASSI?'
  call Quit(_RC_IO_ERROR_READ_)
end if

! Compute number of 'primitive' states.

nStatePrim = 0
do i=1,NrFiles
  nStatePrim = nStatePrim+NrStates(i)
end do

! Open EigV file and read information.

Lu = 92
call DaName(Lu,EigV)
iDisk = 0

! Read RASSCF overlap and H-matrix.

call mma_allocate(NonH,nTri_Elem(nStatePrim),label='NonOrtH')
call mma_allocate(NonS,nTri_Elem(nStatePrim),label='NonOrtS')
kaunt = 0
do i=1,nStatePrim
  do j=1,i
    kaunt = kaunt+1
    call dDaFile(Lu,2,NonH(kaunt),1,iDisk)
  end do
end do
kaunt = 0
do i=1,nStatePrim
  do j=1,i
    kaunt = kaunt+1
    call dDaFile(Lu,2,NonS(kaunt),1,iDisk)
  end do
end do
if (iPrint >= 10) then
  call TriPrt('RASSCF Hamiltonian',' ',NonH,nStatePrim)
  call TriPrt('RASSCF Overlaps',' ',NonS,nStatePrim)
end if
call DaClos(Lu)

! Construct CASSI state basis.

call mma_allocate(Eig2,nStatePrim,nStatePrim,label='RedEigV1')
call ContRASBas(nStatePrim,NonH,NonS,Eig2)
call mma_deallocate(NonH)
call mma_deallocate(NonS)

! Now transform from 'primitive' RASSCF to 'contracted' RASSI states.

call RasRasTrans(nBas(1),nStatePrim,Eig2,iPrint)
call mma_deallocate(Eig2)

! If requested, obtain reduced MO-basis, otherwise just go as usual.

if (MoAveRed) then
  call MoReduce(nBas,nRedMO)
  write(TDMchar,'(A)') 'TDMSCR'
  call FetchTDM(nRedMO,nState,TDMchar)
else
  write(u6,*) '     ----- Use AO-representation of the transition density matrix.'
  nRedMO = 0 !Only a dummy.
end if

! Finished!

write(u6,*) '     ...Done!'

return

end subroutine TdmTrans
