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

subroutine Analyze_Q(iQ_Atoms)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "warnings.h"
parameter(iHUltraMax=1000)
dimension iCo(3)
dimension gR(MxAt,3,iHUltraMax)
data Dum/0.0d0/
dimension iDum(1)

!----------------------------------------------------------------------*
! Some numbers and defaults.                                           *
!----------------------------------------------------------------------*
iHMax = 0
iCStart = (((iQ_Atoms-1)/nAtom)+1)*nCent+1
iCNum = iCStart/nCent
dR = 0.1d0
!----------------------------------------------------------------------*
! Just say what we are doing.                                          *
!----------------------------------------------------------------------*
call NiceOutPut('AAA',Dum,Dum,Dum)
!----------------------------------------------------------------------*
! Open sampfile. Get some numbers about how many sampled etc.          *
!----------------------------------------------------------------------*
call DaName(iLuSaIn,SaFilIn)
iDiskSa = 0
call iDaFile(iLuSaIn,2,iDum,1,iDiskSa)
iHowMSamp = iDum(1)
iDiskTemp = iDiskSa
call WrRdSim(iLuSaIn,2,iDiskSa,iTCSim,64,Etot,Ract,nPart,Dum,Dum,Dum)
iDiskSa = iDiskTemp
!----------------------------------------------------------------------*
! Say something about these numbers to the user.                       *
!----------------------------------------------------------------------*
write(6,*)
write(6,*) 'The sampfile ',SaFilIn,' contains ',iHowMSamp,' sampled configurations.'
write(6,*) 'Total number of particles:',nPart
!----------------------------------------------------------------------*
! BEGIN ANALYZING!                                                     *
!----------------------------------------------------------------------*
do iSamp=1,iHowMSamp
!----------------------------------------------------------------------*
! Begin by getting the coordinates for this configuration. They are    *
! stored in Work(iCo(i)) where i=1 means x-coordinate, i=2 y-coordinate*
! and i=3 z-coordinate.                                                *
!----------------------------------------------------------------------*
  call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Dum,Dum,Dum)
  iDiskSa = iTcSim(1)
  do i=1,3
    call GetMem('Coordinates','Allo','Real',iCo(i),nPart*nCent)
    call dDafile(iLuSaIn,2,Work(iCo(i)),nPart*nCent,iDiskSa)
  end do
  !--------------------------------------------------------------------*
  ! Once we have coordinates, lets compute some distances and start    *
  ! building various distribution functions.                           *
  !--------------------------------------------------------------------*
  do i=1,iQ_Atoms
    do j=1,nAtom
      do k=1,nPart-iCNum
        dist2 = 0.0d0
        do l=1,3
          ind = iCStart+(j-1)+(k-1)*nCent
          dist2 = dist2+(Work(iCo(l)+i-1)-Work(iCo(l)+ind-1))**2
        end do
        dist = sqrt(dist2)
        iH = int((dist+dR*0.5d0)/dR)
        if (iH > iHMax) then
          iHMax = iH
          if (iH > iHUltraMax) then
            write(6,*)
            write(6,*) 'Too fine sections for g(r). Increase section size or allocate more memory.'
            call Quit(_RC_INTERNAL_ERROR_)
          end if
        end if
        gR(i,j,iH) = gR(i,j,iH)+1/dist2
      end do
    end do
  end do
  !--------------------------------------------------------------------*
  ! End loop over sampled configurations.                              *
  !--------------------------------------------------------------------*
  do i=1,3
    call GetMem('Coordinates','Free','Real',iCo(i),nPart*nCent)
  end do
end do
!----------------------------------------------------------------------*
! Time to generate a nice output.                                      *
!----------------------------------------------------------------------*
write(6,*)
write(6,*) 'SUMMARY OF RESULTS FOR SAMPFILE ANALYSIS.'
write(6,*)
do i=1,iQ_Atoms
  write(6,*)
  write(6,*) 'Quantum atom ',i
  write(6,'(5X,A,5X,5(A,I2,1X))') 'Separation',('Solvent atom',k,k=1,nAtom)
  do iH=1,iHMax
    write(6,'(F15.7,5(F15.7))') dR*iH,(gR(i,j,iH),j=1,nAtom)
  end do
end do

call DaClos(iLuSaIn)

return

end subroutine Analyze_Q
