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
! Copyright (C) Thomas Dresselhaus                                     *
!***********************************************************************

subroutine embPotOutput(nAtoms,dens)
!***********************************************************************
!                                                                      *
! Object: Calculates and writes out the electrostatic potential on a   *
!         grid.                                                        *
!                                                                      *
! Called from: OneEl                                                   *
!              RASSCF (version below)                                  *
!                                                                      *
! Calling    : IniSewM                                                 *
!              Drv1_Pot                                                *
!              molcas_open                                             *
!              embpot routines                                         *
!                                                                      *
!     Author: Thomas Dresselhaus                                       *
!                                                                      *
!***********************************************************************

use Embedding_Global, only: embOutEspPath, embWriteEsp, nEmbGridPoints, embGridCoord
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: dens(*)

integer(kind=iwp) :: i, iUnitEmb, j, nOrdOp
real(kind=wp) :: distGpAtom
real(kind=wp), allocatable :: coordsEmb(:,:), chargesEmb(:), espGrid(:)
integer(kind=iwp), external :: isFreeUnit

call EmbPotInit(.true.)
! I dunno whether the coordinates/charges are still stored
! somewhere, but I'll just read them in again here.
call mma_allocate(coordsEmb,3,nAtoms,label='Coor')
call Get_dArray('Unique Coordinates',coordsEmb,3*nAtoms)
call mma_allocate(chargesEmb,nAtoms,label='Charge')
call Get_dArray('Nuclear charge',chargesEmb,nAtoms)
!write(u6,*) 'dens(1:3)=',dens(1:3)
!write(u6,*) 'Coords of gridpt 1:',embGridCoord(:,1)
!write(u6,*) 'Coords of gridpt 3:',embGridCoord(:,3)
! Get the electrostatic potential on the grid
!write(u6,*) 'nEmbGridPoints=', nEmbGridPoints
if (embWriteEsp) then
  call mma_allocate(espGrid,nEmbGridPoints,label='ESP')
  !write(u6,*) 'espGrid(1:3)=',espGrid(1:3)
  call inisewm('mltpl',0)
  !write(u6,*) 'inisewm done!'
  nordop = 0
  call Drv1_Pot(dens,embGridCoord,espGrid,nEmbGridPoints,1,nordop)
  !write(u6,*) 'Coords:'
  !write(u6,*) coordsEmb(:,1)
  !write(u6,*) coordsEmb(:,2)
  !write(u6,*) coordsEmb(:,3)
  !write(u6,*) '----------------------------------------------------'
  !write(u6,*) 'Charges:'
  !write(u6,*) chargesEmb(1:3)
  iUnitEmb = isFreeUnit(1)
  call molcas_open(iUnitEmb,embOutEspPath)
end if
do i=1,nEmbGridPoints
  do j=1,nAtoms
    distGpAtom = sqrt((coordsEmb(1,j)-embGridCoord(1,i))**2+ &
                      (coordsEmb(2,j)-embGridCoord(2,i))**2+ &
                      (coordsEmb(3,j)-embGridCoord(3,i))**2)
    if (embWriteEsp) espGrid(i) = espGrid(i)+chargesEmb(j)/distGpAtom
  end do
  if (embWriteEsp) then
    !write(u6,*) 'ESP on grid point ',i,': ',espGrid(i)
    ! Write to file
    write(iUnitEmb,'(es18.10)') -espGrid(i)
    !write(iUnitEmb,'(es18.10)') espGrid(i)
  end if
end do
if (embWriteEsp) then
  close(iUnitEmb)
  call mma_deallocate(espGrid)
end if
call mma_deallocate(coordsEmb)
call mma_deallocate(chargesEmb)
call embpotfreemem()

end subroutine embPotOutput

! Actually the incoming densities are in an AO basis and treated
! as such...
subroutine embPotOutputMODensities(nAtoms,nSym,densInact,densAct,nBasPerSym,nBasTotSquare)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nSym, nBasPerSym(nSym), nBasTotSquare
real(kind=wp), intent(in) :: densInact(nBasTotSquare), densAct(nBasTotSquare)
integer(kind=iwp) :: i, iCnt, iCnt2, iCol, iDensDimPck, iRow
real(kind=wp), allocatable :: totDens(:), totDensP(:)

iDensDimPck = 0
do i=1,nSym
  iDensDimPck = iDensDimPck+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
end do

call mma_allocate(totDens,nBasTotSquare,label='TotD')
call mma_allocate(totDensP,iDensDimPck,label='TotDp')

call dcopy_(nBasTotSquare,densInact,1,totDens,1)
call daxpy_(nBasTotSquare,One,densAct,1,totDens,1)
!Pack density
iCnt = 0
iCnt2 = 0
do i=1,nSym
  do iRow=1,nBasPerSym(i)
    do iCol=1,iRow
      iCnt = iCnt+1
      if (iRow == iCol) then
        totDensP(iCnt) = totDens((iRow-1)*nBasPerSym(i)+iCol+iCnt2)
      else
        totDensP(iCnt) = Two*totDens((iRow-1)*nBasPerSym(i)+iCol+iCnt2)
      end if
    end do
  end do
  iCnt2 = iCnt2+nBasPerSym(i)*nBasPerSym(i)
end do

call embPotOutput(nAtoms,totDensP)

call mma_deallocate(totDens)
call mma_deallocate(totDensP)

end subroutine embPotOutputMODensities
