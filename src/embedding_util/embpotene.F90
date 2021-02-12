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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
! Object: calculate the energy due to an embedding potential           !
!                                                                      !
! Called from: RASSCF                                                  !
!                                                                      !
! Calling    : ddot_                                                   !
!              dcopy                                                   !
!              daxpy                                                   !
!                                                                      !
!     Author: Thomas Dresselhaus                                       !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!real*8 function embPotEne1(density,length)
!
!implicit none
!real*8, intent(in) :: density(length)
!integer, intent(in) :: length
!
!! Local variables
!integer :: mAdEmbInts, readCheck
!real*8 :: embpotenescf
!
!! Includes
!#include "WrkSpc.fh"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! Allocate memory for the integrals
!call GetMem('EmbPotInts','Allo','Real',mAdEmbInts,length)
!
!! Read in the one-electron integrals due to the embedding potential
!readCheck = -1
!call RdOne(readCheck,6,'embpot  ',1,Work(mAdEmbInts),1)
!if (readCheck /= 0) then
!  write (6,*) 'R1Inta: Error readin ONEINT'
!  write (6,'(a,a)') 'Label=embpot  '
!  call Abend()
!end if
!
!! Calculate
!embPotEne1 = embPotEneSCF(density,Work(mAdEmbInts),length)
!
!! Clean Memory
!call GetMem('EmbPotInts','Free','Real',mAdEmbInts,length)
!
!end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function embPotEneSCF(density,embeddingInts,length)

implicit none
integer, intent(in) :: length
real*8, intent(in) :: density(length)
real*8, intent(in) :: embeddingInts(length)
real*8, external :: ddot_

! In this case the energy is just the dot-product of the vectors which resemble
! the packed matrices.
embPotEneSCF = DDot_(length,density,1,embeddingInts,1)

end function embPotEneSCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function embPotEneMODensities(densityInactive,densityActive,embeddingInts,nBasPerSym,nBasTotSquare,nFrozenOrbs,nSym)

implicit none
real*8, intent(in) :: densityInactive(*), densityActive(*)
real*8, intent(in) :: embeddingInts(*)
!real*8, intent(in) :: coefficientMatrix(nBasFunc**2)
integer, intent(in) :: nBasTotSquare, nFrozenOrbs
integer, intent(in) :: nBasPerSym(*)

! Local variables
real*8, allocatable :: totalDensity(:)
real*8, allocatable :: totalDensityPacked(:)
!real*8, allocatable :: transformedEmbInts(:)
real*8 :: embpotenescf
integer :: i, dens_dim_packed, iRow, iCol, counter, nSym, counter2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO
!if (nFrozenOrbs > 0) then
!  write(6,*) 'ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!'
!  call Abend()
!end if

allocate(totalDensity(nBasTotSquare))

dens_dim_packed = 0
do i=1,nSym
  dens_dim_packed = dens_dim_packed+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
end do
allocate(totalDensityPacked(dens_dim_packed))

call dcopy_(nBasTotSquare,densityInactive,1,totalDensity,1)
call daxpy_(nBasTotSquare,1.0d0,densityActive,1,totalDensity,1)

!write(6,*) 'Unpacked total density'
!do i=1, dens_dim
!  write(6,*) totalDensity(i)
!end do

!Pack density
counter = 0
counter2 = 0
do i=1,nSym
  do iRow=1,nBasPerSym(i)
    do iCol=1,iRow
      counter = counter+1
      if (iRow == iCol) then
        totalDensityPacked(counter) = totalDensity((iRow-1)*nBasPerSym(i)+iCol+counter2)
      else
        totalDensityPacked(counter) = 2*totalDensity((iRow-1)*nBasPerSym(i)+iCol+counter2)
      end if
    end do
  end do
  counter2 = counter2+nBasPerSym(i)*nBasPerSym(i)
end do
!write(6,*) 'nSym=', nSym
!write(6,*) 'nTot2=', nBasTotSquare
!write(6,*) 'dens_dim_packed=', dens_dim_packed
!write(6,*) 'FOOOO Total dens    |  pot'
!do i=1, dens_dim_packed
!  write(6,*) totalDensityPacked(i), '  |  ', embeddingInts(i)
!end do

!! Transformation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Well, a TraFo does not seem to be necessary.                               !
! Justification:                                                             !
! 1. The totalDensity matrix (sum of the incoming ones) is exactly symmetric !
! 2. With the currently activated scheme the integrated density on the grid  !
!    produces the correct number of electrons.                               !
! 3. In a CASSCF calculation with 0 active electrons/orbitals (i.e. actually !
!    a HF calculation) the totalDensity matrix should be diagonal with only  !
!    '2.0's and '0.0's on the diagonal if it is expressed in MO basis. This  !
!    is not observed though.
! Conclusion: the incoming matrices actually are expressed in the AO basis;  !
!             they are just not packed.                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call transAOtoMO(coefficientMatrix,embeddingInts,transformedEmbInts,nBasFunc,nFrozenOrbs)
!
!embPotEneMODensities = embPotEneSCF(totalDensity,transformedEmbInts,dens_dim)
embPotEneMODensities = embPotEneSCF(totalDensityPacked,embeddingInts,dens_dim_packed)

deallocate(totalDensity)
deallocate(totalDensityPacked)
!deallocate(transformedEmbInts)

! Avoid unused argument warnings
return
if (.false.) call Unused_integer(nFrozenOrbs)

end function embPotEneMODensities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Just to be clean this is in a seperate subroutine          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine transAOtoMO (coefficientMatrix,inMatrix,outMatrix,nBasFunc,nFrozenOrbs)
!
!implicit none
!
!real*8, intent(in) :: coefficientMatrix(nBasFunc**2), inMatrix(nBasFunc*(nBasFunc+1)/2)
!real*8, intent(out) :: outMatrix(nBasFunc*(nBasFunc+1)/2)
!rnteger, intent(in) :: nBasFunc, nFrozenOrbs
!
!! Local variables
!Real*8 :: unpackedMatrix(nBasFunc**2), halfTrafoMat(nBasFunc**2)
!
!! TODO
!if (nFrozenOrbs > 0) then
!  write(6,*) 'ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!'
!  call Abend()
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call Square(inMatrix,unpackedMatrix,1,nBasFunc,nBasFunc)
!!call MXMA(unpackedMatrix,1,nBasFunc, &
!           coefficientMatrix(nFrozenOrbs*nBasFunc),1, &
!           nBasFunc,halfTrafoMat,1, &
!           nBasFunc,nBasFunc,nBasFunc,nBasFunc)
!! stefan: check T (transpose) whether it is correct...
!call DGEMM_('T','N',nBasFunc,nBasFunc,nBasFunc, &
!            1.0d0,unpackedMatrix,nBasFunc, &
!            coefficientMatrix(nFrozenOrbs*nBasFunc),nBasFunc, &
!            0.0d0,halfTrafoMat,nBasFunc)
!call MXMT(halfTrafoMat,nBasFunc,1, &
!          coefficientMatrix(nFrozenOrbs*nBasFunc),1, &
!          nBasFunc,outMatrix,nBasFunc,nBasFunc)
!
!end subroutine
