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

use qmstat_global, only: iBigT, ipAvRed, iPrint, MoAveRed, nRedMO, NrFiles, NrStates, nState, RassiM, EigV
use Index_Functions, only: nTri_Elem
use Definitions, only: iwp, u6

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
#include "warnings.h"
integer(kind=iwp) :: nBas(MxSym)
integer(kind=iwp) :: i, iDisk, iEig2, iNonH, iNonS, j, kaunt, Lu, nSize, nStatePrim
character(len=6) :: TDMchar
logical(kind=iwp) :: Exists

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

nSize = nTri_Elem(nStatePrim)
call GetMem('NonOrtH','Allo','Real',iNonH,nSize)
call GetMem('NonOrtS','Allo','Real',iNonS,nSize)
kaunt = 0
do i=1,nStatePrim
  do j=1,i
    call dDaFile(Lu,2,Work(iNonH+kaunt),1,iDisk)
    kaunt = kaunt+1
  end do
end do
kaunt = 0
do i=1,nStatePrim
  do j=1,i
    call dDaFile(Lu,2,Work(iNonS+kaunt),1,iDisk)
    kaunt = kaunt+1
  end do
end do
if (iPrint >= 10) then
  call TriPrt('RASSCF Hamiltonian',' ',Work(iNonH),nStatePrim)
  call TriPrt('RASSCF Overlaps',' ',Work(iNonS),nStatePrim)
end if
call DaClos(Lu)

! Construct CASSI state basis.

call ContRASBas(nStatePrim,iNonH,iNonS,iEig2)
call GetMem('NonOrtH','Free','Real',iNonH,nSize)
call GetMem('NonOrtS','Free','Real',iNonS,nSize)

! Now transform from 'primitive' RASSCF to 'contracted' RASSI states.

call RasRasTrans(nBas(1),nStatePrim,iEig2,iPrint)

! If requested, obtain reduced MO-basis, otherwise just go as usual.

if (MoAveRed) then
  call MoReduce(nBas,nRedMO,ipAvRed)
  write(TDMchar,'(A)') 'TDMSCR'
  call FetchTDM(nRedMO,nState,iBigT,TDMchar)
else
  write(u6,*) '     ----- Use AO-representation of the transition density matrix.'
  nRedMO = 0 !Only a dummy.
end if

! Finished!

write(u6,*) '     ...Done!'

return

end subroutine TdmTrans
