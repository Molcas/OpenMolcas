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

!function embPotEne1(density,length)
!
!use stdalloc, only: mma_allocate, mma_deallocate
!use Definitions, only: wp, iwp
!
!implicit none
!real(kind=wp) :: embPotEne1
!integer(kind=iwp), intent(in) :: length
!real(kind=wp), intent(in) :: density(length)
!
!! Local variables
!integer(kind=iwp) :: readCheck
!real(kind=wp), allocatable :: embInts(:)
!real(kind=wp), external :: embpotenescf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! Allocate memory for the integrals
!call mma_allocate(embInts,length,label='EmbPotInts')
!
!! Read in the one-electron integrals due to the embedding potential
!readCheck = -1
!iOpt = ibset(ibset(0,sNoOri),sNoNuc)
!call RdOne(readCheck,iOpt,'embpot  ',1,embInts,1)
!if (readCheck /= 0) then
!  write(u6,*) 'R1Inta: Error readin ONEINT'
!  write(u6,'(a,a)') 'Label=embpot  '
!  call Abend()
!end if
!
!! Calculate
!embPotEne1 = embPotEneSCF(density,embInts,length)
!
!! Clean Memory
!call mma_deallocate(embInts)
!
!end function embPotEne1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function embPotEneSCF(density,embeddingInts,length)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: embPotEneSCF
integer(kind=iwp), intent(in) :: length
real(kind=wp), intent(in) :: density(length), embeddingInts(length)
real(kind=wp), external :: ddot_

! In this case the energy is just the dot-product of the vectors which resemble
! the packed matrices.
embPotEneSCF = DDot_(length,density,1,embeddingInts,1)

end function embPotEneSCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function embPotEneMODensities(densityInactive,densityActive,embeddingInts,nBasPerSym,nBasTotSquare,nSym)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: embPotEneMODensities
integer(kind=iwp), intent(in) :: nSym, nBasPerSym(nSym), nBasTotSquare
real(kind=wp), intent(in) :: densityInactive(nBasTotSquare), densityActive(nBasTotSquare), embeddingInts(*)
!real(kind=wp), intent(in) :: coefficientMatrix(nBasFunc**2)

! Local variables
real(kind=wp), allocatable :: totalDensity(:), totalDensityPacked(:)
!real(kind=wp), allocatable :: transformedEmbInts(:)
real(kind=wp) :: embpotenescf
integer(kind=iwp) :: i, dens_dim_packed, iRow, iCol, counter, counter2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO
!if (nFrozenOrbs > 0) then
!  write(u6,*) 'ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!'
!  call Abend()
!end if

call mma_allocate(totalDensity,nBasTotSquare)

dens_dim_packed = 0
do i=1,nSym
  dens_dim_packed = dens_dim_packed+((nBasPerSym(i)*(nBasPerSym(i)+1))/2)
end do
call mma_allocate(totalDensityPacked,dens_dim_packed)

call dcopy_(nBasTotSquare,densityInactive,1,totalDensity,1)
call daxpy_(nBasTotSquare,One,densityActive,1,totalDensity,1)

!write(u6,*) 'Unpacked total density'
!do i=1, dens_dim
!  write(u6,*) totalDensity(i)
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
!write(u6,*) 'nSym=', nSym
!write(u6,*) 'nTot2=', nBasTotSquare
!write(u6,*) 'dens_dim_packed=', dens_dim_packed
!write(u6,*) 'FOOOO Total dens    |  pot'
!do i=1, dens_dim_packed
!  write(u6,*) totalDensityPacked(i), '  |  ', embeddingInts(i)
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

call mma_deallocate(totalDensity)
call mma_deallocate(totalDensityPacked)
!call mma_deallocate(transformedEmbInts)

return

end function embPotEneMODensities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Just to be clean this is in a seperate subroutine          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine transAOtoMO (coefficientMatrix,inMatrix,outMatrix,nBasFunc,nFrozenOrbs)
!
!use stdalloc, only: mma_allocate, mma_deallocate
!use Definitions, only: wp, iwp
!
!implicit none
!integer(kind=iwp), intent(in) :: nBasFunc, nFrozenOrbs
!real(kind=wp), intent(in) :: coefficientMatrix(nBasFunc**2), inMatrix(nBasFunc*(nBasFunc+1)/2)
!real(kind=wp), intent(out) :: outMatrix(nBasFunc*(nBasFunc+1)/2)
!
!! Local variables
!real(kind=wp), allocatable :: unpackedMatrix(:), halfTrafoMat(:)
!
!call mma_allocate(unpackedMatrix,nBasFunc**2,label='unpackedMatrix')
!call mma_allocate(halfTrafoMat,nBasFunc**2,label='halfTrafoMat')
!
!! TODO
!if (nFrozenOrbs > 0) then
!  write(u6,*) 'ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!'
!  call Abend()
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call Square(inMatrix,unpackedMatrix,1,nBasFunc,nBasFunc)
!! stefan: check T (transpose) whether it is correct...
!call DGEMM_('T','N',nBasFunc,nBasFunc,nBasFunc,One,unpackedMatrix,nBasFunc,coefficientMatrix(nFrozenOrbs*nBasFunc),nBasFunc,Zero, &
!            halfTrafoMat,nBasFunc)
!call DGEMM_Tri('T','N',nBasFunc,nBasFunc,nBasFunc,One,halfTrafoMat,nBasFunc,coefficientMatrix(nFrozenOrbs*nBasFunc),nBasFunc, &
!               Zero,outMatrix,nBasFunc)
!
!call mma_deallocate(unpackedMatrix)
!call mma_deallocate(halfTrafoMat)
!
!end subroutine transAOtoMO
