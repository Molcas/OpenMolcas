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

subroutine Tully(CIBigArray,NSTATE,NCI)

use Surfacehop_globals, only: decoherence, tullySubVerb, fixedrandL, iseedL, DECO, Ethreshold, RandThreshold, FixedRand, &
                              NSUBSTEPS, InitSeed, rassi_ovlp, Run_rassi, firststep
#ifdef _HDF5_
use Surfacehop_globals, only: lH5Restart
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, cZero, cOne, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSTATE, NCI
real(kind=wp), intent(inout) :: CIBigArray(NCI,NSTATE)
integer(kind=iwp) :: cayley_info, i, ii, irlxroot, iseed, ISTATE2, j, jjj, k, maxhop, nciquery, nhop, nsatom, nstatesq, &
                     RASSI_time_run, root_ovlp_el, stateRi, temproot, values(8)
real(kind=wp) :: DT, ediffcheck, EKIN, Etot, hstep, LO, populOS, prod, root_ovlp, SumProb, temp, tloc
complex(kind=wp) :: ArelaxPrev
logical(kind=iwp) :: found, HOPPED, lmaxHop, lnhop, normalTully
character(len=10) :: time
character(len=8) :: date
character(len=5) :: zone
integer(kind=iwp), allocatable :: decVec(:), stateORDER(:)
real(kind=wp), allocatable :: Bmatrix(:,:), CIBigArrayP(:,:), CIBigArrayPP(:,:), currOVLP(:,:), currOVLP_ras(:,:), currPHASE(:), &
                              D12matrix(:,:), D32matrix(:,:), Dmatrix(:,:), ExtrInter(:,:), ExtrSlope(:,:), Gprobab(:), Popul(:), &
                              prevOVLP(:,:), prevPHASE(:), readOVLP(:,:), saveOVLP(:,:), saveOVLP_bk(:,:), sp(:,:), TAU(:), &
                              tempVector(:), tempVector2(:), V(:), Venergy(:), VenergyInter(:), VenergyP(:), VenergySlope(:)
complex(kind=wp), allocatable :: Amatrix(:,:), Mminus(:,:), Uprop(:,:)
real(kind=wp), external :: Random_Molcas

#include "warnings.h"

call mma_allocate(CIBigArrayP,NCI,NSTATE,Label='CIBIgArrayP')
call mma_allocate(CIBigArrayPP,NCI,NSTATE,Label='CIBIgArrayPP')
CIBigArrayP(:,:) = Zero
CIBigArrayPP(:,:) = Zero
call mma_allocate(Venergy,NSTATE,Label='Venergy')

write(u6,*)
write(u6,*) '------------------------------------------'
write(u6,*) '            TULLY ALGORITHM'
write(u6,*) '------------------------------------------'
write(u6,*)

if (rassi_ovlp) then
  write(u6,*) 'Using RASSI for WF overlap'
  write(u6,*)
else
  write(u6,*) 'Using CI vector product for WF overlap'
  write(u6,*)
end if

call get_darray('Last energies',Venergy,NSTATE)
call Get_iScalar('Relax CASSCF root',iRlxRoot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!             Watch if dt exists.  then if Amatrix exists              !
!        if exists, it gets it, if not, it creates a new one           !
!                                                                      !

call qpg_dscalar('Timestep',Found)

#ifdef _HDF5_
if (.not. Found .and. lH5Restart) then
  call restart_surfacehop()
  Found = .true.
  call Get_iScalar('Relax CASSCF root',iRlxRoot)
  call Put_iScalar('NumGradRoot',iRlxRoot)
  call Put_iScalar('Relax Original root',iRlxRoot)
  call Put_dScalar('Last energy',Venergy(iRlxRoot))
end if
#endif

if (Found) call Get_dScalar('Timestep',DT)

call mma_allocate(Amatrix,NSTATE,NSTATE,Label='Amatrix')
call Qpg_zArray('AmatrixV',Found,nStateSq)
write(u6,*) 'Did the density matrix exists? ',Found
if (.not. Found) then
  Amatrix(:,:) = cZero
  Amatrix(iRlxRoot,iRlxRoot) = cOne

  call Put_zArray('AmatrixV',Amatrix,NSTATE**2)
else
  call get_zarray('AmatrixV',Amatrix,NSTATE**2)
end if

call mma_allocate(Popul,NSTATE,Label='Popul')
call mma_allocate(VenergyP,NSTATE,Label='VenergyP')

call Qpg_dArray('AllCIP',Found,nCiQuery)
write(u6,*) 'Did the Pre-coefficients array exists? ',Found
if (.not. Found) then
  call put_darray('AllCIP',CIBigArray,NCI*NSTATE)
  call put_dArray('VenergyP',Venergy,NSTATE)
  do i=1,NSTATE
    Popul(i) = real(Amatrix(i,i))
  end do
  write(u6,*) 'Gnuplot:',Popul(:),Venergy(:),Venergy(iRlxRoot)
  write(u6,*) 'Cannot do deltas at first step, see you later! '
  firststep = .true.
  call cleanup()
  return
else
  firststep = .false.
  call Get_dArray('AllCIP',CIBigArrayP,NCI*NSTATE)
  call Get_dArray('VenergyP',VenergyP,NSTATE)
end if

! Check if overlap exists and if RASSI is run yet for this timestep

if (rassi_ovlp) then
  call Qpg_iscalar('SH RASSI run',Found)
  !write(u6,*) 'Has RASSI ever been run?', Found
  if (.not. Found) then ! RASSI never run before (step 2) or restart
    Run_rassi = .true.
  else
    call get_iscalar('SH RASSI run',RASSI_time_run)
    !write(u6,*) 'RASSI_time_run variable = ',RASSI_time_run
    if (RASSI_time_run == 0) then ! Need to run RASSI for this timestep still
      Run_rassi = .true.
    else if (RASSI_time_run == 1) then
      Run_rassi = .false. ! RASSI already run for this timestep
    else
      write(u6,*) 'Problem checking if RASSI previously run'
      call Abend()
    end if
  end if

  ! return to call RASSI in surfacehop.f90 if not yet run

  if (Run_rassi) then
    write(u6,*) 'Calling RASSI...'
    RASSI_time_run = 1
    call put_iscalar('SH RASSI run',RASSI_time_run)
    call cleanup()
    return
  else
    write(u6,*) 'RASSI already called, continuing...'
    RASSI_time_run = 0 ! Reset for next iteration
    call put_iscalar('SH RASSI run',RASSI_time_run)
    write(u6,*)
    call mma_allocate(readOVLP,2*NSTATE,2*NSTATE,Label='readOVPL')
    call get_dArray('State Overlaps',readOVLP,4*NSTATE**2)
    !do i=1,2*NSTATE
    !  write(u6,*) readOVLP(i,:)
    !end do
    call mma_allocate(currOVLP_ras,NSTATE,NSTATE,Label='currOVLP_ras')
    currOVLP_ras(:,:) = readOVLP(1:NSTATE,NSTATE+1:2*NSTATE)  ! Transpose to match RASSI printed version
    call mma_deallocate(readOVLP)
  end if
end if

! now check for the CI coefficients at Pre-Pre-Step (PP)

call mma_allocate(prevOVLP,NSTATE,NSTATE,Label='prevOVLP')
call mma_allocate(prevPHASE,NSTATE,Label='prevPHASE')

call Qpg_dArray('AllCIPP',Found,nCiQuery)
write(u6,*) 'Did the Pre-Pre-coefficients array exists? ',Found
if (.not. Found) then
  normalTully = .true.
  prevPHASE(:) = One
  call put_darray('AllCIPP',CIBigArrayP,NCI*NSTATE)
  call put_darray('AllCIP',CIBigArray,NCI*NSTATE)
  write(u6,*) 'At second step we will use normal Tully Algorithm:'
else
  normalTully = .false.
  call Get_dArray('AllCIPP',CIBigArrayPP,NCI*NSTATE)
  if (rassi_ovlp) then
    write(u6,*) 'Grabbing previous overlap <t-2dt|t-dt>'
    call Get_dArray('SH_Ovlp_save',prevOVLP,NSTATE**2)
    call Get_dArray('Old_Phase',prevPHASE,NSTATE)
    write(u6,*)
    write(u6,*) '<t-2dt|t-dt> RASSI Overlap'
    do i=1,NSTATE
      write(u6,*) prevOVLP(i,:)
    end do
    write(u6,*)
  end if
end if

call Get_dScalar('MD_Etot',Etot)
call Get_iScalar('Unique atoms',nsAtom)

write(u6,*) 'Density Matrix elements (i=1..#states):'
do i=1,NSTATE
  write(u6,*) Amatrix(i,:)
end do

! Timestep:                      DT
! Total Energy                   Etot
! Coefficients:                  CIBigArray(i,j)    length = NCI,NSTATE
! Prev step coefficients:        CIBigArrayP(i,j)   length = NCI,NSTATE
! Prev-Prev coefficients:        CIBigArrayPP(i,j)  length = NCI,NSTATE
! V energy:                      Venergy(i)         length = NSTATE
! V Prev step energy:            VenergyP(i)        length = NSTATE
! MatriX A:                      Amatrix(i,j)       length = NSTATE,NSTATE
! IF RASSI
! Overlap <t-dt|t> uncorrected   currOVLP_ras       length = NSTATE,NSTATE
! Overlap <t-dt|t> corrected     currOVLP           length = NSTATE,NSTATE
! Overlap <t-dt|t> corr to save  saveOVLP           length = NSTATE,NSTATE
! Overlap <t-2dt|t-dt>           prevOVLP           length = NSTATE,NSTATE
! Last timestep phase diff       prevPHASE          length = NSTATE
! Curr timestep phase diff       currPHASE          length = NSTATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                       Sign Corrector                                 !
!                                                                      !

call mma_allocate(stateORDER,NSTATE,Label='stateORDER')

if (.not. rassi_ovlp) then
  ! so first of all I create 2 temp vectors that store the absolute value
  ! of the energy difference

  call mma_allocate(tempVector,NSTATE,Label='tempVector')
  call mma_allocate(tempVector2,NSTATE,Label='tempVector')

  ! Original sign corrector/root reordering using CI vector product
  write(u6,*) 'Using CI vector products for sign correction/root ordering'
  tempVector(:) = abs(Venergy(:)-Venergy(irlxRoot))
  tempVector2(:) = tempVector(:)

  ! then I sort one of them, (relaxroot becomes first, it's zero)

  do j=1,NSTATE-1
    do k=j+1,NSTATE
      if (tempVector(j) > tempVector(k)) then
        temp = tempVector(j)
        tempVector(j) = tempVector(k)
        tempVector(k) = temp
      end if
    end do
  end do

  ! I get the right order I need to process roots

  call mma_allocate(decVec,NSTATE,Label='decVec')
  call mma_allocate(sp,NSTATE,NSTATE,Label='sp')

  do i=1,NSTATE
    do j=1,NSTATE
      if (tempVector(i) == tempVector2(j)) stateORDER(i) = j
    end do
  end do
  ! ii counter on CURRENT STEP
  do ii=1,NSTATE
    do i=1,NSTATE
      sp(ii,i) = sum(CIBigArray(:,ii)*CIBigArrayP(:,i))
    end do
  end do
  decVec(:) = 1

  do ii=1,NSTATE
    stateRi = stateORDER(ii)
    prod = Zero
    jjj = 0
    do i=1,NSTATE
      if (decVec(i) == 1) then
        if (abs(sp(stateRi,i)) > abs(prod)) then
          prod = sp(stateRi,i)
          jjj = i
        end if
      end if
    end do
    decVec(jjj) = 0
    if (prod < Zero) CIBigArray(:,stateRi) = -CIBigArray(:,stateRi)
  end do

  call mma_deallocate(decVec)
  call mma_deallocate(sp)
  call mma_deallocate(tempVector)
  call mma_deallocate(tempVector2)

else

  call mma_allocate(currOVLP,NSTATE,NSTATE,Label='currOVLP')
  call mma_allocate(saveOVLP,NSTATE,NSTATE,Label='saveOVLP')

  ! Sign correction and root reordering using RASSI overlap matrix <t-dt|t>
  ! Product of all prev phase diff. applied to columns (t-dt)
  ! Product of all prev phase diff. including current step applied to rows (t)
  write(u6,*)
  write(u6,*) 'Using RASSI overlap matrix for sign correction/root ordering'

  ! Root flipping

  do i=1,NSTATE
    stateORDER(i) = i
    root_ovlp = abs(currOVLP_ras(i,i)) ! Diagonal element
    root_ovlp_el = i
    do j=i,NSTATE
      if (abs(currOVLP_ras(j,i)) > root_ovlp) then
        root_ovlp_el = j
        root_ovlp = abs(currOVLP_ras(j,i))
      end if
    end do
    if (root_ovlp < 0.4_wp) write(u6,*) 'WARNING: No overlap greater than 0.4 for root:',i
    if (root_ovlp_el /= i) then
      write(u6,*) 'Root rotation detected, swapping roots ',i,root_ovlp_el
      stateORDER(i) = root_ovlp_el
      currOVLP(i,:) = currOVLP_ras(root_ovlp_el,:)
      currOVLP_ras(root_ovlp_el,:) = currOVLP_ras(i,:)
      currOVLP_ras(i,:) = currOVLP(i,:)
    else
      currOVLP(i,:) = currOVLP_ras(i,:)
    end if
  end do

  call mma_deallocate(currOVLP_ras)

  write(u6,*) '<t-dt|t> RASSI Overlap Uncorrected'
  do i=1,NSTATE
    write(u6,*) currOVLP(i,:)
  end do

  write(u6,*) 'State ordering',stateORDER(:)

  call mma_allocate(currPHASE,NSTATE,Label='currPHASE')

  ! Sign correction
  ! Get current phase - phase read from RASSI * previous phase array
  do i=1,NSTATE
    if (currOVLP(i,i) < Zero) then
      currPHASE(i) = -prevPHASE(i)
    else
      currPHASE(i) = prevPHASE(i)
    end if
  end do

  write(u6,*) 'Current Phase'
  write(u6,*) currPHASE(:)
  write(u6,*) 'Previous Phase'
  write(u6,*) prevPHASE(:)

  ! Correct columns of current RASSI using prevPHASE

  do i=1,NSTATE
    currOVLP(:,i) = currOVLP(:,i)*prevPHASE(i)
  end do

  if (tullySubVerb) then
    write(u6,*) '<t-dt|t> RASSI Overlap Columns Corrected'
    do i=1,NSTATE
      write(u6,*) currOVLP(i,:)
    end do
  end if

  ! Correct rows of current RASSI using currPHASE

  do i=1,NSTATE
    currOVLP(i,:) = currOVLP(i,:)*currPHASE(i)
  end do

  write(u6,*) '<t-dt|t> RASSI Overlap phase Corrected'
  do i=1,NSTATE
    write(u6,*) currOVLP(i,:)
  end do

  call mma_allocate(saveOVLP_bk,NSTATE,NSTATE,Label='saveOVLP_bk')

  ! Swap back rows if root flipping occured

  do i=1,NSTATE
    saveOVLP_bk(i,:) = currOVLP(stateORDER(i),:)
  end do

  ! Swap columns if root flipping occured to give final 'to save' RASSI

  do i=1,NSTATE
    saveOVLP(:,i) = saveOVLP_bk(:,stateORDER(i))
  end do

  call mma_deallocate(saveOVLP_bk)

  write(u6,*) '<t-dt|t> RASSI Overlap to save (flipping incl.)'
  do i=1,NSTATE
    write(u6,*) currOVLP(i,:)
  end do
end if

call mma_deallocate(prevPHASE)
call mma_deallocate(stateORDER)

!                                                                      !
!                        end of sign corrector                         !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_allocate(Dmatrix,NSTATE,NSTATE,Label='Dmatrix')
call mma_allocate(ExtrInter,NSTATE,NSTATE,Label='ExtrInter')
call mma_allocate(ExtrSlope,NSTATE,NSTATE,Label='ExtrSlope')
call mma_allocate(VenergyInter,NSTATE,Label='VenergyInter')
call mma_allocate(VenergySlope,NSTATE,Label='VenergySlope')

if_normaltully: if (normalTully) then

  ! create D matrix normal tully

  write(u6,*)
  write(u6,*) 'Executing Normal Tully !!'
  write(u6,*)

  if (.not. rassi_ovlp) then
    write(u6,*) 'Using CI vector product to calculate D matrix'
    do i=1,NSTATE
      do j=1,NSTATE
        if (i == j) then
          Dmatrix(i,i) = Zero
        else
          Dmatrix(i,j) = -sum(CIBigArray(:,i)*CIBigArrayP(:,j))/DT
        end if
      end do
    end do
  else
    write(u6,*) 'Using RASSI overlap to calculate D matrix'
    Dmatrix(:,:) = currOVLP(:,:)/DT
    do i=1,NSTATE
      Dmatrix(i,i) = Zero
    end do
  end if

  normalTully = .false.

  ExtrInter(:,:) = Dmatrix(:,:)
  ExtrSlope(:,:) = Zero
  VenergyInter(:) = VenergyP(:)
  VenergySlope(:) = (Venergy(:)-VenergyP(:))/DT

else if_normaltully

  call mma_allocate(D12matrix,NSTATE,NSTATE,Label='D12matrix')
  call mma_allocate(D32matrix,NSTATE,NSTATE,Label='D32matrix')

  ! Create D matrix according to Hammes-Schiffer-Tully (interpolating extrapolating)

  if (.not. rassi_ovlp) then
    write(u6,*) 'Using CI vector product to calculate D Matrix'

    ! D32matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i == j) then
          D32matrix(i,i) = Zero
        else
          D32matrix(i,j) = sum(CIBigArrayPP(:,i)*CIBigArrayP(:,j)-CIBigArrayP(:,i)*CIBigArrayPP(:,j))/(Two*DT)
        end if
      end do
    end do

    ! D12matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i == j) then
          D12matrix(i,i) = Zero
        else
          D12matrix(i,j) = sum(CIBigArrayP(:,i)*CIBigArray(:,j)-CIBigArray(:,i)*CIBigArrayP(:,j))/(Two*DT)
        end if
      end do
    end do

  else ! RASSI overlaps

    write(u6,*) 'Using RASSI overlap to calculate D matrix'

    ! D32matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i == j) then
          D32matrix(i,j) = Zero
        else
          D32matrix(i,j) = (prevOVLP(i,j)-prevOVLP(j,i))/(Two*DT)
        end if
      end do
    end do

    ! D12matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i == j) then
          D12matrix(i,j) = Zero
        else
          D12matrix(i,j) = (currOVLP(i,j)-currOVLP(j,i))/(Two*DT)
        end if
      end do
    end do

  end if

  ! definition of Y intercept (ExtrInter) and slope (ExtrSlope) for EXTRapolation line

  ExtrSlope(:,:) = (D12matrix(:,:)-D32matrix(:,:))/DT
  ExtrInter(:,:) = D12matrix(:,:)
  VenergyInter(:) = VenergyP(:)
  VenergySlope(:) = (Venergy(:)-VenergyP(:))/DT

  call mma_deallocate(D12matrix)
  call mma_deallocate(D32matrix)

end if if_normaltully

call mma_deallocate(prevOVLP)

if (rassi_ovlp) call mma_deallocate(currOVLP)

! UNCOMMENT to print coefficients !!!
! Just a few coefficients
! write(u6,*) 'WaveFunctionsCoefficients are: ',NCI,'*',NSTATE
! write(u6,*) '       This step        Previous step:        PP:'

! write(u6,*) CIBigArray(1,1),CIBigArrayP(1,1),CIBigArrayPP(1,1)
! write(u6,*) CIBigArray(2,1),CIBigArrayP(2,1),CIBigArrayPP(2,1)
! write(u6,*) CIBigArray(3,1),CIBigArrayP(3,1),CIBigArrayPP(3,1)
! write(u6,*) CIBigArray(4,1),CIBigArrayP(4,1),CIBigArrayPP(4,1)

! All coefficients
!do i=1,NSTATE
!  do j=1,NCI
!    write(u6,*) CIBigArray(j,i),CIBigArrayP(j,i),CIBigArrayPP(j,i)
!  end do
!end do
! UNCOMMENT to print coefficients !!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                        GET RANDOM number LO                          !
!                                                                      !

! Check if Random number is fixed by user
if (fixedrandL) then
  write(u6,*) 'From input Random is Fixed at:',FixedRand
  LO = FixedRand/real(NSUBSTEPS,kind=wp)
  write(u6,*) 'Normalized by Substep number:',LO
else
  ! check if seed number is read from input / RunFile
  if (iseedL) then
    call qpg_iscalar('Seed',Found)
    if (Found) then
      call get_iscalar('Seed',iseed)
      write(u6,*) 'Seed number read from the RunFile: ',iseed
    else
      iseed = InitSeed ! initial seed read from input
      write(u6,*) 'Seed number read from input file: ',iseed
    end if
  else
    ! or generate a new random seed number
    call date_and_time(date,time,zone,values)
    ! Just milliseconds multiplied by seconds
    iseed = ((values(7)+1)*values(8)+1)
  end if
  LO = max(RandThreshold,Random_Molcas(iseed))
  call put_iscalar('Seed',iseed)
end if
!                                                                      !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                             INTEGRATOR                               !
!                   NSUBSTEPS is the substep number                    !
!                                                                      !

write(u6,*) 'executing the integration:'
write(u6,*) 'Time Step is:    ',DT
write(u6,*) 'Substeps are: ',NSUBSTEPS

temproot = iRlxRoot

call mma_allocate(Bmatrix,NSTATE,NSTATE,Label='Bmatrix')
call mma_allocate(Gprobab,NSTATE,Label='Gprobab')
call mma_allocate(Mminus,NSTATE,NSTATE,Label='Mminus')
call mma_allocate(TAU,NSTATE,Label='TAU')
call mma_allocate(Uprop,NSTATE,NSTATE,Label='Uprop')
call mma_allocate(V,NSTATE,Label='V')

substeps: do ii=1,NSUBSTEPS
  if (.not. fixedrandL) then
    LO = max(RandThreshold,Random_Molcas(iseed))
    call put_iscalar('Seed',iseed)
  end if

  hstep = DT/real(NSUBSTEPS,kind=wp)
  tloc = (real(ii,kind=wp)-Half)*hstep
  Dmatrix(:,:) = ExtrInter(:,:)+ExtrSlope(:,:)*(tloc-DT*Half)
  V(:) = VenergyInter(:)+VenergySlope(:)*tloc

  ! Enforce strict antisymmetry of the time-derivative coupling matrix.
  do i=1,NSTATE
    Dmatrix(i,i) = Zero
    do j=i+1,NSTATE
      temp = Half*(Dmatrix(i,j)-Dmatrix(j,i))
      Dmatrix(i,j) = temp
      Dmatrix(j,i) = -temp
    end do
  end do

  ! Unitary Cayley propagation of the electronic density matrix.
  ! The generator G = D + iV gives dA/dt = A G - G A.
  ! The diagonal energies are shifted by the current active-state energy
  ! before building G: a constant shift only adds a global phase to the
  ! propagator, which cancels exactly in U^dagger*A*U, but it is essential
  ! for the accuracy of the Cayley approach, which requires
  ! |h*G/2| << 1.

  Uprop(:,:) = cmplx(Dmatrix(:,:),kind=wp)
  do i=1,NSTATE
    Uprop(i,i) = Uprop(i,i)+(V(i)-V(temproot))*Onei
  end do
  Uprop(:,:) = Half*hstep*Uprop(:,:)
  Mminus(:,:) = -Uprop(:,:)
  do i=1,NSTATE
    Uprop(i,i) = Uprop(i,i)+cOne
    Mminus(i,i) = Mminus(i,i)+cOne
  end do

  ! Uprop is overwritten by U = inv(I - hG/2) * (I + hG/2).
  call solve_cayley_system()
  if (cayley_info /= 0) then
    write(u6,*) 'Cayley propagation failed in Tully. Singular pivot = ',cayley_info
    call Abend()
  end if

  ! Mminus used for scratch
  call zgemm_('N','N',NSTATE,NSTATE,NSTATE,cOne,Amatrix,NSTATE,Uprop,NSTATE,cZero,Mminus,NSTATE)
  call zgemm_('C','N',NSTATE,NSTATE,NSTATE,cOne,Uprop,NSTATE,Mminus,NSTATE,cZero,Amatrix,NSTATE)

  Bmatrix(:,:) = -Two*real(conjg(Amatrix(:,:))*Dmatrix(:,:))
  do i=1,NSTATE
    !B(i,i) not used
    Bmatrix(i,i) = Two*aimag(conjg(Amatrix(i,i))*V(i))
  end do

  do i=1,NSTATE
    if (i == temproot) then
      Gprobab(i) = Zero
    else
      Gprobab(i) = max(Zero,Bmatrix(i,temproot)*DT/(real(Amatrix(temproot,temproot))*NSUBSTEPS))
    end if
  end do

  SumProb = Zero
  do i=1,NSTATE
    if (i == temproot) cycle
    SumProb = SumProb+Gprobab(i)
    if (LO <= SumProb) then
      write(u6,*) 'Following Tully, should hop from',temproot,'to',i
      if (V(i) < Etot) then
        write(u6,*) 'this root has an energy lower than the total',Etot,' (thus, permitted)'
        Ediffcheck = V(i)-V(temproot)
        write(u6,*) 'Ediffcheck is:',Ediffcheck
        Ediffcheck = abs(Ediffcheck)
        if (Ediffcheck < Ethreshold) then
          write(u6,*) 'lower than the threshold:',Ethreshold
          write(u6,*) 'temproot set to:',i
          temproot = i
        end if
      end if
      exit
    end if
  end do

  ! Persico-Granucci all in a THEN branch

  ! Note: the 1/2 factor in the exponents corrects for an error in eq. 17 of doi:10.1063/1.2715585
  ! (it should correct population instead of amplitudes)

  pg: if (decoherence) then
    EKIN = Etot-V(temproot)
    if (EKIN <= Zero) then
      write(u6,*) 'WARNING! Negative Kinetic Energy. Ekin= ',EKIN,' a.u.'
      write(u6,*) 'Kinetic energy rescaled to 10 e-5.'
      EKIN = 0.00001_wp
    end if
    do i=1,NSTATE
      if (i /= temproot) TAU(i) = abs(One/(V(temproot)-V(i)))*(One+DECO/EKIN)
    end do

    do i=1,NSTATE
      if (i == temproot) cycle
      do j=1,NSTATE
        if (j /= temproot) Amatrix(i,j) = Amatrix(i,j)*exp(-Half*(DT/NSUBSTEPS)/TAU(i))*exp(-Half*(DT/NSUBSTEPS)/TAU(j))
      end do
    end do

    populOS = Zero
    do i=1,NSTATE
      if (i /= temproot) populOS = populOS+real(Amatrix(i,i))
    end do

    ArelaxPrev = Amatrix(temproot,temproot)
    Amatrix(temproot,temproot) = One-populOS

    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= temproot .and. j == temproot) then
          Amatrix(i,j) = Amatrix(i,j)*exp(-Half*(DT/NSUBSTEPS)/TAU(i))*sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
        else if (i == temproot .and. j /= temproot) then
          Amatrix(i,j) = Amatrix(i,j)*exp(-Half*(DT/NSUBSTEPS)/TAU(j))*sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
        end if
      end do
    end do
  end if pg

  do i=1,NSTATE
    Popul(i) = real(Amatrix(i,i))
  end do

  if (tullySubVerb) then
    write(u6,*) 'Substep:',ii
    write(u6,*) 'Temproot:',temproot
    write(u6,*) 'MATRIX D:'
    do i=1,NSTATE
      write(u6,*) Dmatrix(i,:)
    end do
    write(u6,*) 'Density Matrix elements (i=1..#states):'
    do i=1,NSTATE
      write(u6,*) Amatrix(i,:)
    end do
    write(u6,*) 'B Matrix:'
    do i=1,NSTATE
      write(u6,*) Bmatrix(i,:)
    end do
    write(u6,*) 'Probabilities:'
    write(u6,*) Gprobab(:)
    write(u6,*) 'Random Number is: ',LO
    write(u6,*) 'Populations Energies:',Popul(:),V(:),V(temproot)
    write(u6,*)
  end if

end do substeps

call mma_deallocate(Bmatrix)
call mma_deallocate(Dmatrix)
call mma_deallocate(ExtrInter)
call mma_deallocate(ExtrSlope)
call mma_deallocate(Gprobab)
call mma_deallocate(Mminus)
call mma_deallocate(TAU)
call mma_deallocate(Uprop)
call mma_deallocate(V)
call mma_deallocate(VenergyInter)
call mma_deallocate(VenergySlope)

!                                                                      !
!                        END OF INTEGRATOR                             !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(u6,*) 'Gnuplot:',Popul(:),Venergy(:),Venergy(temproot)
!write(u6,*) 'Gnuplot:',Popul(:),V(:),V(temproot)

call Add_Info('Pop',Popul,NSTATE,5)

!write(u6,*) 'CASSCF rlxrt', iRlxRoot

if (temproot == iRlxRoot) then
  HOPPED = .false.
else
  ! Is the keyword MAXHOP set?
  call Qpg_iScalar('MaxHopsTully',lmaxHop)
  if (lmaxHop) then
    call get_iScalar('MaxHopsTully',maxHop)
    call qpg_iScalar('Number of Hops',lnHop)
    if (.not. lnHop) then
      nHop = 0
    else
      call get_iScalar('Number of Hops',nHop)
    end if
    write(u6,*) 'User set a max of',maxHop
    write(u6,*) 'nhop is:',nHop
    if (maxHop <= nHop) then
      write(u6,*) 'This surface HOP is not allowed'
      HOPPED = .false.
    else
      nHop = nHop+1
      call put_iScalar('Number of Hops',nHop)
      HOPPED = .true.
    end if
  else
    ! In case temproot is different then iRlxRoot, and no lmaxHop has been set
    HOPPED = .true.
  end if
end if

! give the "HOP TO:" state name ISTATE2 !!!!
! start the HOPPING procedure
if (HOPPED) then
  ISTATE2 = temproot
  write(u6,*) 'HOP ALLOWED'
  write(u6,'(6X,A,A,A)') '+',repeat('-',118),'+'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,A,2(47X,A))') '|','A HOP event is detected!|'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,A,44X,2(A,I3,4X),40X,A)') '|','From state:',iRlxRoot,'To state:',ISTATE2,'|'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,A,A,A,//)') '+',repeat('-',118),'+'

  call Put_iScalar('NumGradRoot',ISTATE2)
  call Put_iScalar('Relax CASSCF root',ISTATE2)
  call Put_iScalar('Relax Original root',ISTATE2)
  call Put_dScalar('Last energy',Venergy(ISTATE2))
end if

! scale velocities
!
!call get_dArray('Velocities',vel,nsAtom*3)
!
!write(u6,*) 'Velocities before Hop:'
!do i=1,nsAtom
!  write(u6,*) vel(i*3-2),vel(i*3-1),vel(i*3)
!end do
!EKIN = Etot-Venergy(iRlxRoot)
!EKIN_target = Etot-Venergy(temproot)
!scalfac = sqrt(Ekin_target/Ekin)
!write(u6,*) Etot,Venergy(iRlxRoot),Venergy(temproot),EKIN,EKIN_target,scalfac
!vel(:,:) = scalfac*vel(:,:)
!write(u6,*) 'Velocities after Hop:'
!do i=1,nsAtom
!  write(u6,*) vel(:,i)
!end do
!
!call put_dArray('Velocities',vel,nsAtom*3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!                             SAVING                                   !
!                                                                      !

call put_lscalar('hopped',HOPPED)
call Put_dArray('VenergyP',Venergy,NSTATE)
call Put_dArray('AllCIPP',CIBigArrayP,NCI*NSTATE)
call Put_dArray('AllCIP',CIBigArray,NCI*NSTATE)
call put_zarray('AmatrixV',Amatrix,NSTATE**2)

if (rassi_ovlp) then
  call Put_dArray('SH_Ovlp_Save',saveOVLP,NSTATE**2)
  call Put_dArray('Old_Phase',currPHASE,NSTATE)
  call mma_deallocate(currPHASE)
  call mma_deallocate(saveOVLP)
end if
!                                                                      !
!                           END SAVING                                 !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call cleanup()

contains

subroutine solve_cayley_system()

  integer(kind=iwp) :: i, k, piv
  real(kind=wp) :: piv_abs, test_abs
  complex(kind=wp) :: factor
  complex(kind=wp) :: tmp(NSTATE)

  cayley_info = 0

  do k=1,NSTATE
    piv = k
    piv_abs = abs(Mminus(k,k))
    do i=k+1,NSTATE
      test_abs = abs(Mminus(i,k))
      if (test_abs > piv_abs) then
        piv = i
        piv_abs = test_abs
      end if
    end do

    if (piv_abs <= Zero) exit

    if (piv /= k) then
      tmp(:) = Mminus(k,:)
      Mminus(k,:) = Mminus(piv,:)
      Mminus(piv,:) = tmp(:)
      tmp(:) = Uprop(k,:)
      Uprop(k,:) = Uprop(piv,:)
      Uprop(piv,:) = tmp(:)
    end if

    factor = One/Mminus(k,k)
    if (abs(factor) > Zero) then
      Mminus(k,:) = factor*Mminus(k,:)
      Uprop(k,:) = factor*Uprop(k,:)
    end if

    do i=1,NSTATE
      if (i == k) cycle
      factor = Mminus(i,k)
      if (abs(factor) > Zero) then
        Mminus(i,:) = Mminus(i,:)-factor*Mminus(k,:)
        Uprop(i,:) = Uprop(i,:)-factor*Uprop(k,:)
      end if
    end do
  end do

  if (k <= NSTATE) cayley_info = k

end subroutine solve_cayley_system

subroutine cleanup()

  call mma_deallocate(Amatrix)
  call mma_deallocate(CIBigArrayP)
  call mma_deallocate(CIBigArrayPP)
  call mma_deallocate(Popul)
  call mma_deallocate(Venergy)
  call mma_deallocate(VenergyP,safe='*')

end subroutine cleanup

end subroutine Tully
