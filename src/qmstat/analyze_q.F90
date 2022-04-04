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

use qmstat_global, only: iLuSaIn, iTcSim, nAtom, nCent, nPart, SaFilIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms
integer(kind=iwp) :: i, iCNum, iCStart, iDiskSa, iDiskTemp, iDum(1), iH, iHMax, iHowMSamp, ind, iSamp, j, k, l
real(kind=wp) :: dist, dist2, dR, Dum, Etot, Ract
real(kind=wp), allocatable :: Co(:,:), gR(:,:,:)
integer(kind=iwp), parameter :: iHUltraMax = 1000
#include "warnings.h"

Dum = Zero
!----------------------------------------------------------------------*
! Some numbers and defaults.                                           *
!----------------------------------------------------------------------*
iHMax = 0
iCStart = (((iQ_Atoms-1)/nAtom)+1)*nCent+1
iCNum = iCStart/nCent
dR = 0.1_wp
!----------------------------------------------------------------------*
! Just say what we are doing.                                          *
!----------------------------------------------------------------------*
call NiceOutPut('AAA')
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
write(u6,*)
write(u6,*) 'The sampfile ',SaFilIn,' contains ',iHowMSamp,' sampled configurations.'
write(u6,*) 'Total number of particles:',nPart
!----------------------------------------------------------------------*
! BEGIN ANALYZING!                                                     *
!----------------------------------------------------------------------*
call mma_allocate(gR,iQ_Atoms,nAtom,iHUltraMax,label='gR')
gR(:,:,:) = Zero
do iSamp=1,iHowMSamp
  !--------------------------------------------------------------------*
  ! Begin by getting the coordinates for this configuration. They are  *
  ! stored in Co(:,i) where i=1 means x-coordinate, i=2 y-coordinate   *
  ! and i=3 z-coordinate.                                              *
  !--------------------------------------------------------------------*
  call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Dum,Dum,Dum)
  iDiskSa = iTcSim(1)
  call mma_allocate(Co,nPart*nCent,3,label='Coordinates')
  call dDafile(iLuSaIn,2,Co,3*nPart*nCent,iDiskSa)
  !--------------------------------------------------------------------*
  ! Once we have coordinates, lets compute some distances and start    *
  ! building various distribution functions.                           *
  !--------------------------------------------------------------------*
  do i=1,iQ_Atoms
    do j=1,nAtom
      do k=1,nPart-iCNum
        dist2 = Zero
        do l=1,3
          ind = iCStart+(j-1)+(k-1)*nCent
          dist2 = dist2+(Co(i,l)-Co(ind,l))**2
        end do
        dist = sqrt(dist2)
        iH = int((dist+dR*Half)/dR)
        if (iH > iHMax) then
          iHMax = iH
          if (iH > iHUltraMax) then
            write(u6,*)
            write(u6,*) 'Too fine sections for g(r). Increase section size or allocate more memory.'
            call Quit(_RC_INTERNAL_ERROR_)
          end if
        end if
        gR(i,j,iH) = gR(i,j,iH)+One/dist2
      end do
    end do
  end do
  !--------------------------------------------------------------------*
  ! End loop over sampled configurations.                              *
  !--------------------------------------------------------------------*
  call mma_deallocate(Co)
end do
!----------------------------------------------------------------------*
! Time to generate a nice output.                                      *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,*) 'SUMMARY OF RESULTS FOR SAMPFILE ANALYSIS.'
write(u6,*)
do i=1,iQ_Atoms
  write(u6,*)
  write(u6,*) 'Quantum atom ',i
  write(u6,'(5X,A,5X,5(A,I2,1X))') 'Separation',('Solvent atom',k,k=1,nAtom)
  do iH=1,iHMax
    write(u6,'(F15.7,5(F15.7))') dR*iH,gR(i,:,iH)
  end do
end do

call mma_deallocate(gR)

call DaClos(iLuSaIn)

return

end subroutine Analyze_Q
