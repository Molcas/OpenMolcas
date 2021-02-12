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

#include "WrkSpc.fh"

#include "embpotdata.fh"
real*8 distGpAtom

call EmbPotInit(.true.)
! I dunno whether the coordinates/charges are still stored
! somewhere, but I'll just read them in again here.
call GetMem('Coor','ALLO','REAL',ipCoordsEmb,3*nAtoms)
call Get_dArray('Unique Coordinates',Work(ipCoordsEmb),3*nAtoms)
call GetMem('Charge','ALLO','REAL',ipChargesEmb,nAtoms)
call Get_dArray('Nuclear charge',Work(ipChargesEmb),nAtoms)
!write(6,*) 'Work(mAdDns:+2)=',Work(mAdDns),Work(mAdDns+1),Work(mAdDns+2)
!write(6,*) 'Coords of gridpt 1:',Work(posEmbGridCoord),Work(posEmbGridCoord+1),Work(posEmbGridCoord+2)
!write(6,*) 'Coords of gridpt 3:',Work(posEmbGridCoord+6),Work(posEmbGridCoord+7),Work(posEmbGridCoord+8)
! Get the electrostatic potential on the grid
!write(6,*) 'nEmbGridPoints=', nEmbGridPoints
if (embWriteEsp) then
  call GetMem('ESP','ALLO','REAL',ipEspGrid,nEmbGridPoints)
  !write(6,*) 'Work(ipEspGrid:+2)=',Work(ipEspGrid),Work(ipEspGrid+1), Work(ipEspGrid+2)
  call inisewm('mltpl',0)
  !write(6,*) 'inisewm done!'
  nordop = 0
  call Drv1_Pot(Work(mAdDns),Work(posEmbGridCoord),Work(ipEspGrid),nEmbGridPoints,1,nordop)
  !write(6,*) 'Coords:'
  !write(6,*) Work(ipCoordsEmb),Work(ipCoordsEmb+1),Work(ipCoordsEmb+2)
  !write(6,*) Work(ipCoordsEmb+3),Work(ipCoordsEmb+4),Work(ipCoordsEmb+5)
  !write(6,*) Work(ipCoordsEmb+6),Work(ipCoordsEmb+7),Work(ipCoordsEmb+8)
  !write(6,*) '----------------------------------------------------'
  !write(6,*) 'Charges:'
  !write(6,*) Work(ipChargesEmb),Work(ipChargesEmb+1),Work(ipChargesEmb+2)
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
    !write(6,*) 'ESP on grid point ',i,': ',Work(ipEspGrid+i)
    ! Write to file
    write(iUnitEmb,'(e18.10)') (-1.0*Work(ipEspGrid+i))
    !write(iUnitEmb,'(e18.10)') Work(ipEspGrid+i)
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

#include "WrkSpc.fh"

  integer nBasTotSquare, nSym
  integer nBasPerSym(nSym)

  iDensDimPck = 0
  do i=1,nSym
    iDensDimPck = iDensDimPck+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
  end do

  call GetMem('TotD','ALLO','REAL',ipTotDens,nBasTotSquare)
  call GetMem('TotDp','ALLO','REAL',ipTotDensP,iDensDimPck)

  call dcopy_(nBasTotSquare,Work(ipDensInact),1,Work(ipTotDens),1)
  call daxpy_(nBasTotSquare,1.0d0,Work(ipDensAct),1,Work(ipTotDens),1)
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
