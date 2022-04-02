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

! GET NUMBERS FROM SAMPFILE.
subroutine Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)

use qmstat_global, only: Cordst, delFi, delR, delX, iLuSaIn, iPrint, iTcSim, nCent, nMacro, nMicro, nPart
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Ract
integer(kind=iwp), intent(in) :: iQ_Atoms, info_atom(iQ_Atoms)
real(kind=wp), intent(in) :: Coord(3,iQ_Atoms)
integer(kind=iwp), intent(inout) :: iDiskSa
integer(kind=iwp) :: i
real(kind=wp) :: Esub, Etot, Gamold, GaOld
character(len=200) :: Head
real(kind=wp), allocatable :: CT(:)

call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Gamold,GaOld,Esub)
iDiskSa = iTcSim(1)
call mma_allocate(CT,nPart*nCent,label='CTemp')
do i=1,3
  call dDaFile(iLuSaIn,2,CT,nPart*nCent,iDiskSa)
  Cordst(i,:) = CT
  iDiskSa = iTcSim(i+1)
end do
call mma_deallocate(CT)

! We dummy-read the induced dipoles from the sampfile.

!call mma_allocate(Dum,nPart*nPol,label='Dummy')
!do i=1,3
!  call dDaFile(iLuSaIn,2,Dum,nPol*nPart,iDiskSa)
!end do
!call mma_deallocate(Dum)

! And now we place the QM-molecule in proper place and set some
! numbers to zero or one so we only collect configurations from
! the sampfile.

call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
delX = Zero
delFi = Zero
delR = Zero
nMacro = 1
nMicro = 1

! Some printing if requested.

if (iPrint >= 15) then
  write(Head,*) 'Coordinates after substitution in configuration read from sampfile.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine Get9
