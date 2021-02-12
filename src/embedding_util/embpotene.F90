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

!Real*8 function embPotEne1(density, length)
! Implicit none
!
! Real*8,  Intent(in) :: density(length)
! Integer, Intent(in) :: length
!
! ! Local variables
! Integer :: mAdEmbInts, readCheck
! Real*8  :: embpotenescf
!
! ! Includes
!#include "WrkSpc.fh"
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ! Allocate memory for the integrals
! Call GetMem('EmbPotInts','Allo','Real',mAdEmbInts,length)
!
! ! Read in the one-electron integrals due to the embedding potential
! readCheck=-1
! Call RdOne(readCheck,6,'embpot  ',1,Work(mAdEmbInts),1)
! If (readCheck.ne.0) Then
!    Write (6,*) 'R1Inta: Error readin ONEINT'
!    Write (6,'(A,A)') 'Label=embpot  '
!    Call Abend()
! End If
!
! ! Calculate
! embPotEne1=embPotEneSCF(density, Work(mAdEmbInts), length)
!
! ! Clean Memory
! Call GetMem('EmbPotInts','Free','Real',mAdEmbInts,length)
!
!end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real*8 function embPotEneSCF(density, embeddingInts, length)
 Implicit none

 Integer, Intent(in) :: length
 Real*8,  Intent(in) :: density(length)
 Real*8,  Intent(in) :: embeddingInts(length)
 real*8, external    :: ddot_

 ! In this case the energy is just the dot-product of the vectors which resemble
 ! the packed matrices.
 embPotEneSCF = DDot_(length,density,1,embeddingInts,1)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real*8 function embPotEneMODensities( &
        densityInactive, densityActive, embeddingInts, &
        nBasPerSym, nBasTotSquare, nFrozenOrbs, nSym)
 Implicit none

 Real*8,  Intent(in) :: densityInactive(*), densityActive(*)
 Real*8,  Intent(in) :: embeddingInts(*)
! Real*8,  Intent(in) :: coefficientMatrix(nBasFunc**2)
 Integer, Intent(in) :: nBasTotSquare, nFrozenOrbs
 Integer, Intent(in) :: nBasPerSym(*)

 ! Local variables
 Real*8, allocatable :: totalDensity(:)
 Real*8, allocatable :: totalDensityPacked(:)
! Real*8, allocatable :: transformedEmbInts(:)
 Real*8              :: embpotenescf

 Integer :: i, dens_dim_packed, iRow, iCol, counter, nSym, counter2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! TODO
! if (nFrozenOrbs .gt. 0) then
!  write(6,*) "ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!"
!  call Abend()
! end if

 allocate(totalDensity( nBasTotSquare ))

 dens_dim_packed = 0
 do i=1, nSym
   dens_dim_packed = dens_dim_packed + ((nBasPerSym(i) * (nBasPerSym(i)+1)) / 2)
 end do
 allocate(totalDensityPacked( dens_dim_packed ))

 call dcopy_(nBasTotSquare,densityInactive,1,totalDensity,1)
 call daxpy_(nBasTotSquare,1.0d0,densityActive,1,totalDensity,1)

! write(*,*) "Unpacked total density"
! do i=1, dens_dim
!  write(*,*) totalDensity(i)
! end do

 !Pack density
 counter = 0
 counter2 = 0
 do i=1, nSym
  do iRow=1, nBasPerSym(i)
   do iCol=1, iRow
    counter = counter+1
    if (iRow .eq. iCol) then
     totalDensityPacked(counter) = totalDensity((iRow-1)*nBasPerSym(i)+iCol+counter2)
    else
     totalDensityPacked(counter) = 2*totalDensity((iRow-1)*nBasPerSym(i)+iCol+counter2)
    end if
   end do
  end do
  counter2 = counter2 + nBasPerSym(i)*nBasPerSym(i)
 end do
! write(*,*) "nSym=", nSym
! write(*,*) "nTot2=", nBasTotSquare
! write(*,*) "dens_dim_packed=", dens_dim_packed
! write(*,*) "FOOOO Total dens    |  pot"
! do i=1, dens_dim_packed
!  write(*,*) totalDensityPacked(i), "  |  ", embeddingInts(i)
! end do


! ! Transformation
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
! call transAOtoMO(coefficientMatrix, embeddingInts, transformedEmbInts,&
!                  nBasFunc, nFrozenOrbs)
!
! embPotEneMODensities = embPotEneSCF( totalDensity,&
!                                      transformedEmbInts,&
!                                      dens_dim)
 embPotEneMODensities = embPotEneSCF( totalDensityPacked,&
                                      embeddingInts,&
                                      dens_dim_packed)

 deallocate(totalDensity)
 deallocate(totalDensityPacked)
! deallocate(transformedEmbInts)

! Avoid unused argument warnings
 return
 if (.false.) call Unused_integer(nFrozenOrbs)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Just to be clean this is in a seperate subroutine          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine transAOtoMO (coefficientMatrix, inMatrix, outMatrix,&
!                        nBasFunc, nFrozenOrbs)
!
! implicit none
!
! Real*8,  intent(in)  :: coefficientMatrix(nBasFunc**2)
! Real*8,  intent(in)  :: inMatrix(&
!                          nBasFunc*(nBasFunc+1)/2 )
! Real*8,  intent(out) :: outMatrix(&
!                          nBasFunc*(nBasFunc+1)/2 )
! Integer, intent(in)  :: nBasFunc, nFrozenOrbs
!
! ! Local variables
! Real*8 :: unpackedMatrix(nBasFunc**2), halfTrafoMat(nBasFunc**2)
!
! ! TODO
! if (nFrozenOrbs .gt. 0) then
!  write(6,*) "ERROR! Usage of an embedding potential is not yet able to treat frozen orbitals!"
!  call Abend()
! end if
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Call Square(inMatrix,unpackedMatrix,1,nBasFunc,nBasFunc)
!!Call MXMA(unpackedMatrix,1,nBasFunc,&
!!          coefficientMatrix(nFrozenOrbs*nBasFunc),1,&
!!          nBasFunc,halfTrafoMat,1,nBasFunc,&
!!          nBasFunc,nBasFunc,nBasFunc)
!! stefan: check T (transpose) whether it is correct...
! CALL DGEMM_('T','N',nBasFunc,nBasFunc,nBasFunc, &
!             1.0d0,unpackedMatrix,nBasFunc, &
!             coefficientMatrix(nFrozenOrbs*nBasFunc),nBasFunc,&
!             0.0d0,halfTrafoMat,nBasFunc)
! Call MXMT(halfTrafoMat,nBasFunc,1, &
!           coefficientMatrix(nFrozenOrbs*nBasFunc),1,&
!           nBasFunc,outMatrix,nBasFunc,nBasFunc)
!
!end subroutine
