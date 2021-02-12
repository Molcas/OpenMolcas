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

subroutine embPotOutput(nAtoms,mAdDns)
!***********************************************************************
!                                                                      *
! Object: Calculates and writes out the electrostatic potential on a   *
!         grid.                                                        *
!                                                                      *
! Called from: OneEl                                                   *
!              RASSCF (version below)                                  *
!                                                                      *
! Calling    : GetMem                                                  *
!              IniSewM                                                 *
!              Drv1_Pot                                                *
!              molcas_open                                             *
!              embpot routines                                         *
!                                                                      *
!     Author: Thomas Dresselhaus                                       *
!                                                                      *
!***********************************************************************

use Embedding_Global, only: embOutEspPath, embWriteEsp, nEmbGridPoints, posEmbGridCoord
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, mAdDns

#include "WrkSpc.fh"
integer(kind=iwp) :: i, ipCoordsEmb, ipChargesEmb, ipEspGrid, iUnitEmb, j, nOrdOp
real(kind=wp) :: distGpAtom
integer(kind=iwp), external :: isFreeUnit

call EmbPotInit(.true.)
! I dunno whether the coordinates/charges are still stored
! somewhere, but I'll just read them in again here.
call GetMem('Coor','ALLO','REAL',ipCoordsEmb,3*nAtoms)
call Get_dArray('Unique Coordinates',Work(ipCoordsEmb),3*nAtoms)
call GetMem('Charge','ALLO','REAL',ipChargesEmb,nAtoms)
call Get_dArray('Nuclear charge',Work(ipChargesEmb),nAtoms)
!write(u6,*) 'Work(mAdDns:+2)=',Work(mAdDns),Work(mAdDns+1),Work(mAdDns+2)
!write(u6,*) 'Coords of gridpt 1:',Work(posEmbGridCoord),Work(posEmbGridCoord+1),Work(posEmbGridCoord+2)
!write(u6,*) 'Coords of gridpt 3:',Work(posEmbGridCoord+6),Work(posEmbGridCoord+7),Work(posEmbGridCoord+8)
! Get the electrostatic potential on the grid
!write(u6,*) 'nEmbGridPoints=', nEmbGridPoints
if (embWriteEsp) then
  call GetMem('ESP','ALLO','REAL',ipEspGrid,nEmbGridPoints)
  !write(u6,*) 'Work(ipEspGrid:+2)=',Work(ipEspGrid),Work(ipEspGrid+1), Work(ipEspGrid+2)
  call inisewm('mltpl',0)
  !write(u6,*) 'inisewm done!'
  nordop = 0
  call Drv1_Pot(Work(mAdDns),Work(posEmbGridCoord),Work(ipEspGrid),nEmbGridPoints,1,nordop)
  !write(u6,*) 'Coords:'
  !write(u6,*) Work(ipCoordsEmb),Work(ipCoordsEmb+1),Work(ipCoordsEmb+2)
  !write(u6,*) Work(ipCoordsEmb+3),Work(ipCoordsEmb+4),Work(ipCoordsEmb+5)
  !write(u6,*) Work(ipCoordsEmb+u6),Work(ipCoordsEmb+7),Work(ipCoordsEmb+8)
  !write(u6,*) '----------------------------------------------------'
  !write(u6,*) 'Charges:'
  !write(u6,*) Work(ipChargesEmb),Work(ipChargesEmb+1),Work(ipChargesEmb+2)
  iUnitEmb = isFreeUnit(1)
  call molcas_open(iUnitEmb,embOutEspPath)
end if
do i=0,nEmbGridPoints-1
  do j=0,nAtoms-1
    distGpAtom = sqrt((Work(ipCoordsEmb+3*j)-Work(posEmbGridCoord+3*i))*(Work(ipCoordsEmb+3*j)-Work(posEmbGridCoord+3*i))+ &
                      (Work(ipCoordsEmb+3*j+1)-Work(posEmbGridCoord+3*i+1))*(Work(ipCoordsEmb+3*j+1)-Work(posEmbGridCoord+3*i+1))+ &
                      (Work(ipCoordsEmb+3*j+2)-Work(posEmbGridCoord+3*i+2))*(Work(ipCoordsEmb+3*j+2)-Work(posEmbGridCoord+3*i+2)))
    if (embWriteEsp) Work(ipEspGrid+i) = Work(ipEspGrid+i)+Work(ipChargesEmb+j)/distGpAtom
  end do
  if (embWriteEsp) then
    !write(u6,*) 'ESP on grid point ',i,': ',Work(ipEspGrid+i)
    ! Write to file
    write(iUnitEmb,'(es18.10)') (-1.0*Work(ipEspGrid+i))
    !write(iUnitEmb,'(es18.10)') Work(ipEspGrid+i)
  end if
end do
if (embWriteEsp) then
  close(iUnitEmb)
  call GetMem('ESP','FREE','REAL',ipEspGrid,nEmbGridPoints)
end if
call GetMem('Coor','FREE','REAL',ipCoordsEmb,3*nAtoms)
call GetMem('Charge','FREE','REAL',ipChargesEmb,nAtoms)
call embpotfreemem()

end subroutine embPotOutput

! Actually the incoming densities are in an AO basis and treated
! as such...
subroutine embPotOutputMODensities(nAtoms,nSym,ipDensInact,ipDensAct,nBasPerSym,nBasTotSquare)

use Constants, only: One
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nSym, ipDensInact, ipDensAct, nBasPerSym(nSym), nBasTotSquare
integer(kind=iwp) :: i, iCnt, iCnt2, iCol, iDensDim, iDensDimPck, ipTotDens, ipTotDensP, iRow

#include "WrkSpc.fh"

iDensDimPck = 0
do i=1,nSym
  iDensDimPck = iDensDimPck+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
end do

call GetMem('TotD','ALLO','REAL',ipTotDens,nBasTotSquare)
call GetMem('TotDp','ALLO','REAL',ipTotDensP,iDensDimPck)

call dcopy_(nBasTotSquare,Work(ipDensInact),1,Work(ipTotDens),1)
call daxpy_(nBasTotSquare,One,Work(ipDensAct),1,Work(ipTotDens),1)
!Pack density
iCnt = 0
iCnt2 = 0
do i=1,nSym
  do iRow=1,nBasPerSym(i)
    do iCol=1,iRow
      iCnt = iCnt+1
      if (iRow == iCol) then
        Work(ipTotDensP+iCnt-1) = Work(ipTotDens+(iRow-1)*nBasPerSym(i)+iCol-1+iCnt2)
      else
        Work(ipTotDensP+iCnt-1) = 2*Work(ipTotDens+(iRow-1)*nBasPerSym(i)+iCol-1+iCnt2)
      end if
    end do
  end do
  iCnt2 = iCnt2+nBasPerSym(i)*nBasPerSym(i)
end do

call embPotOutput(nAtoms,ipTotDensP)

call GetMem('TotD','FREE','REAL',ipTotDens,iDensDim)
call GetMem('TotDp','FREE','REAL',ipTotDensP,iDensDimPck)

end subroutine embPotOutputMODensities
