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

subroutine Freq_Molden(Freq,nFreq,Vectors,nVectors,nSym,Intens,mDisp,RedMas)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nFreq, nVectors, nSym, mDisp(nSym)
real(kind=wp), intent(in) :: Freq(nFreq), Vectors(nVectors), Intens(nFreq), RedMas(nFreq)
integer(kind=iwp) :: iCoor, iCoord, iFreq, Lu_9, nAll_Atoms, nCoord, nUnique_Atoms
real(kind=wp), allocatable :: Coord(:,:), NMode(:,:,:)
character(len=2), allocatable :: Element(:)
integer(kind=iwp), external :: isFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('Freq',' ',Freq,1,nFreq)
call RecPrt('Intens',' ',Intens,1,nFreq)
call RecPrt('Vectors',' ',Vectors,1,nVectors)
write(u6,*) 'mDisp=',mDisp
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Open input file for MOLDEN

Lu_9 = 9
Lu_9 = isFreeUnit(Lu_9)
call molcas_open(Lu_9,'MD_FREQ')
!                                                                      *
!***********************************************************************
!                                                                      *
write(Lu_9,*) '[Molden Format]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write frequecnies to Molden input file

write(Lu_9,*) '[N_FREQ]'
write(Lu_9,*) nFreq
write(Lu_9,*) '[FREQ]'
do iFreq=1,nFreq
  write(Lu_9,*) Freq(iFreq)
end do
write(Lu_9,*) '[INT]'
do iFreq=1,nFreq
  write(Lu_9,*) Intens(iFreq)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Write coordinates of all centers

call Get_nAtoms_All(nCoord)
call mma_allocate(Coord,3,nCoord,label='Coord')
call Get_Coord_All(Coord,nCoord)
#ifdef _DEBUGPRINT_
call RecPrt('Coord(all)',' ',Coord,3,nCoord)
#endif
call mma_allocate(Element,nCoord,label='Element')
call Get_Name_All(Element)
write(Lu_9,*) '[NATOM]'
write(Lu_9,*) nCoord

write(Lu_9,*) '[FR-COORD]'
do iCoord=1,nCoord
  write(Lu_9,*) Element(iCoord),Coord(:,iCoord)
end do
call mma_deallocate(Coord)
call mma_deallocate(Element)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write normal modes, observe that the order here is the same as
! for the coordinates.

call Get_iScalar('Unique atoms',nUnique_Atoms)
call Get_nAtoms_All(nAll_Atoms)
call mma_allocate(NMode,3,nAll_Atoms,nFreq,label='NMode')
NMode(:,:,:) = Zero
call Get_NMode_All(Vectors,nVectors,nFreq,nUnique_Atoms,NMode,nAll_Atoms,mDisp)
write(Lu_9,*) '[FR-NORM-COORD]'
#ifdef _DEBUGPRINT_
call RecPrt('Normal Modes',' ',NMode,3*nAll_Atoms,nFreq)
#endif
do iFreq=1,nFreq
  write(Lu_9,*) 'vibration ',iFreq
  do iCoor=1,nAll_Atoms
    write(Lu_9,*) NMode(:,iCoor,iFreq)
  end do
end do
call mma_deallocate(NMode)
!                                                                      *
!***********************************************************************
!                                                                      *
! Alessio Valentini 2018 - add reduced masses to freq.molden file
! in order to use this file for computation of initial conditions
! in semiclassical molecular dynamics

write(Lu_9,*) '[RMASS]'
do iFreq=1,nFreq
  write(Lu_9,*) RedMas(iFreq)
end do

close(Lu_9)

return

end subroutine Freq_Molden
