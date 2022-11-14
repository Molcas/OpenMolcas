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

use Tully_variables, only: decoherence, tullySubVerb, fixedrandL, iseedL, DECO, Ethreshold, RandThreshold, FixedRand, NSUBSTEPS, &
                           InitSeed, rassi_ovlp, Run_rassi, firststep
#ifdef _HDF5_
use Surfacehop_globals, only: lH5Restart
#endif
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSTATE, NCI
real(kind=wp), intent(inout) :: CIBigArray(NCI*NSTATE)
#include "warnings.h"

integer(kind=iwp) :: values(8)
character(len=8) :: date
character(len=10) :: time
character(len=5) :: zone
logical(kind=iwp) :: HOPPED, normalTully, found, lmaxHop, lnhop
integer(kind=iwp) :: maxhop, nhop, RASSI_time_run
real(kind=wp) :: DT, LO, EKIN, TAU(NSTATE)
real(kind=wp) :: CIBigArrayP(NCI*NSTATE), readOVLP(NSTATE*2*NSTATE*2)
real(kind=wp) :: CIBigArrayPP(NCI*NSTATE), Etot, ediffcheck, currOVLP_ras(NSTATE,NSTATE)
real(kind=wp) :: currOVLP(NSTATE,NSTATE), prevOVLP(NSTATE,NSTATE), saveOVLP(NSTATE,NSTATE), saveOVLP_bk(NSTATE,NSTATE)
real(kind=wp) :: currPHASE(NSTATE), prevPHASE(NSTATE)
real(kind=wp) :: Dmatrix(NSTATE,NSTATE), sp(NSTATE,NSTATE)
real(kind=wp) :: D32matrix(NSTATE,NSTATE), D12matrix(NSTATE,NSTATE)
real(kind=wp) :: ExtrSlope(NSTATE,NSTATE), ExtrInter(NSTATE,NSTATE)
real(kind=wp) :: VenergySlope(NSTATE), tempVector(NSTATE)
real(kind=wp) :: tempVector2(NSTATE), VenergyInter(NSTATE)
real(kind=wp) :: V(NSTATE,NSTATE), Bmatrix(NSTATE,NSTATE)
real(kind=wp) :: Gprobab(NSTATE), Popul(NSTATE)
real(kind=wp) :: VenergyP(NSTATE), Venergy(NSTATE), temp, root_ovlp
real(kind=wp) :: SumProb, scalarprod, prod, populOS
integer(kind=iwp) :: k, l, j, i, ii, jjj, t, tt, o, root_ovlp_el
integer(kind=iwp) :: rightOrder(NSTATE), decVec(NSTATE), stateORDER(NSTATE)
integer(kind=iwp) :: nstatesq, nciquery, stateRi, temproot, nsatom
integer(kind=iwp) :: ISTATE2, iseed, irlxroot
complex(kind=wp) :: Amatrix(NSTATE,NSTATE), AmatrixDT(NSTATE,NSTATE), ArelaxPrev
real(kind=wp), external :: Random_Molcas

CIBigArrayP(:) = Zero
CIBigArrayPP(:) = Zero
V(:,:) = Zero
TAU(:) = Zero
write(u6,*) ''
write(u6,*) '------------------------------------------'
write(u6,*) '            TULLY ALGORITHM'
write(u6,*) '------------------------------------------'
write(u6,*) ''

if (rassi_ovlp) then
  write(u6,*) 'Using RASSI for WF overlap'
  write(u6,*) ''
else
  write(u6,*) 'Using CI vector product for WF overlap'
  write(u6,*) ''
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

if (Found) then
  call Get_dScalar('Timestep',DT)
end if

call Qpg_zArray('AmatrixV',Found,nStateSq)
write(u6,*) 'Did the density matrix exists? ',Found
if (.not. Found) then
  Amatrix(:,:) = (Zero,Zero)
  Amatrix(iRlxRoot,iRlxRoot) = (One,Zero)

  call Put_zArray('AmatrixV',Amatrix,NSTATE*NSTATE)
else
  call get_zarray('AmatrixV',Amatrix,NSTATE*NSTATE)
end if

call Qpg_dArray('AllCIP',Found,nCiQuery)
write(u6,*) 'Did the Pre-coefficients array exists? ',Found
if (.not. Found) then
  call put_darray('AllCIP',CIBigArray,NCI*NSTATE)
  call put_dArray('VenergyP',Venergy,NSTATE)
  do i=1,NSTATE
    Popul(i) = real(Amatrix(i,i))
  end do
  write(u6,*) 'Gnuplot:',(Popul(j),j=1,NSTATE,1),(Venergy(j),j=1,NSTATE,1),Venergy(iRlxRoot)
  write(u6,*) 'Cannot do deltas at first step, see you later! '
  firststep = .true.
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
    return
  else
    write(u6,*) 'RASSI already called, continuing...'
    RASSI_time_run = 0 ! Reset for next iteration
    call put_iscalar('SH RASSI run',RASSI_time_run)
    write(u6,*) ''
    call get_dArray('State Overlaps',readOVLP,NSTATE*2*NSTATE*2)
    !do i=1,NSTATE**4
    !  write(u6,*) readOVLP(i)
    !end do
    do t=1,NSTATE
      do tt=1,NSTATE
        o = ((2*t-1)*NSTATE)+tt
        currOVLP_ras(tt,t) = readOVLP(o)  ! Transpose to match RASSI printed version
      end do
    end do
  end if
end if

! now check for the CI coefficients at Pre-Pre-Step (PP)

call Qpg_dArray('AllCIPP',Found,nCiQuery)
write(u6,*) 'Did the Pre-Pre-coefficients array exists? ',Found
if (.not. Found) then
  normalTully = .true.
  do i=1,NSTATE
    prevPHASE(i) = 1
  enddo
  call put_darray('AllCIPP',CIBigArrayP,NCI*NSTATE)
  call put_darray('AllCIP',CIBigArray,NCI*NSTATE)
  write(u6,*) 'At second step we will use normal Tully Algorithm:'
else
  normalTully = .false.
  call Get_dArray('AllCIPP',CIBigArrayPP,NCI*NSTATE)
  if (rassi_ovlp) then
    write(u6,*) 'Grabbing previous overlap <t-2dt|t-dt>'
    call Get_dArray('SH_Ovlp_save',prevOVLP,NSTATE*NSTATE)
    call Get_dArray('Old_Phase',prevPHASE,NSTATE)
    write(u6,*) ''
    write(u6,*) '<t-2dt|t-dt> RASSI Overlap'
    do i=1,NSTATE
      write(u6,*) (prevOVLP(i,j),j=1,NSTATE,1)
    end do
    write(u6,*) ''
  end if
end if

call Get_dScalar('MD_Etot',Etot)
call Get_iScalar('Unique atoms',nsAtom)

write(u6,*) 'Density Matrix elements (i=1..#states):'
do i=1,NSTATE
  write(u6,*) (Amatrix(i,j),j=1,NSTATE,1)
end do

! Timestep:                      DT
! Total Energy                   Etot
! Coefficients:                  CIBigArray(i)      length = NCI*NSTATE
! Prev step coefficients:        CIBigArrayP(i)     length = NCI*NSTATE
! Prev-Prev coefficients:        CIBigArrayPP(i)    length = NCI*NSTATE
! V energy:                      Venergy(i)         length = NSTATE
! V Prev step energy:            VenergyP(i)        length = NSTATE
! MatriX A:                      Amatrix(i,j)       length = NSTATE*NSTATE
! IF RASSI
! Overlap <t-dt|t> uncorrected   currOVLP_ras       length = NSTATE*NSTATE
! Overlap <t-dt|t> corrected     currOVLP           length = NSTATE*NSTATE
! Overlap <t-dt|t> corr to save  saveOVLP           length = NSTATE*NSTATE
! Overlap <t-2dt|t-dt>           prevOVLP           length = NSTATE*NSTATE
! Last timestep phase diff       prevPHASE          length = NSTATE
! Curr timestep phase diff       currPHASE          length = NSTATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!                         Sign Corrector                                 !
!                                                                        !

! so first of all I create 2 temp vectors that store the absolute value
! of the energy difference

if (.not. rassi_ovlp) then
  ! Original sign corrector/root reordering using CI vector product
  write(u6,*) 'Using CI vector products for sign correction/root ordering'
  do i=1,NSTATE
    tempVector(i) = abs(Venergy(i)-Venergy(irlxRoot))
    tempVector2(i) = abs(Venergy(i)-Venergy(irlxRoot))
  end do

  ! then I sort one of them, (relaxroot becomes first, it's zero)

  do j=1,(NSTATE-1)
    do k=(j+1),NSTATE
      if (tempVector(j) > tempVector(k)) then
        temp = tempVector(j)
        tempVector(j) = tempVector(k)
        tempVector(k) = temp
      end if
    end do
  end do

  ! I get the right order I need to process roots

  do i=1,NSTATE
    do j=1,NSTATE
      if (tempVector(i) == tempVector2(j)) then
        rightOrder(i) = j
      end if
    end do
  end do
  ! ii counter on CURRENT STEP
  do ii=1,NSTATE
    do i=1,NSTATE
      scalarprod = Zero
      do j=1,NCI
        scalarprod = scalarprod+CIBigArray(NCI*(ii-1)+j)*CIBigArrayP(NCI*(i-1)+j)
      end do
      sp(ii,i) = scalarprod
    end do
    decVec(ii) = 1
  end do

  do ii=1,NSTATE
    stateRi = rightOrder(ii)
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
    if (prod < 0) then
      do k=1,NCI
        CIBigArray(NCI*(stateRi-1)+k) = -CIBigArray(NCI*(stateRi-1)+k)
      end do
    end if
  end do

else

  ! Sign correction and root reordering using RASSI overlap matrix <t-dt|t>
  ! Product of all prev phase diff. applied to columns (t-dt)
  ! Product of all prev phase diff. including current step applied to rows (t)
  write(u6,*) ''
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
    if (root_ovlp < 0.4) then
      write(u6,*) 'WARNING: No overlap greater than 0.4 for root:',i
    end if
    if (root_ovlp_el /= i) then
      write(u6,*) 'Root rotation detected, swapping roots ', i, root_ovlp_el
      do ii=1,NSTATE
        stateORDER(i) = root_ovlp_el
        currOVLP(i,ii) = currOVLP_ras(root_ovlp_el,ii)
        currOVLP_ras(root_ovlp_el,ii) = currOVLP_ras(i,ii)
        currOVLP_ras(i,ii) = currOVLP(i,ii)
      end do
    else
      do ii=1,NSTATE
        currOVLP(i,ii) = currOVLP_ras(i,ii)
      end do
    end if
  end do

  write(u6,*) '<t-dt|t> RASSI Overlap Uncorrected'
  do i=1,NSTATE
    write(u6,*) (currOVLP(i,j),j=1,NSTATE,1)
  end do

  write(u6,*) 'State ordering', (stateORDER(i),i=1,NSTATE,1)

  ! Sign correction
  ! Get current phase - phase read from RASSI * previous phase array
  do i=1,NSTATE
    currPHASE(i) = prevPHASE(i)
    if (currOVLP(i,i) < 0) then
      currPHASE(i) = -1*prevPHASE(i)
    end if
  enddo

  write(u6,*) 'Current Phase'
  write(u6,*) (currPHASE(i),i=1,NSTATE,1)
  write(u6,*) 'Previous Phase'
  write(u6,*) (prevPHASE(i),i=1,NSTATE,1)

  ! Correct columns of current RASSI using prevPHASE

  do i=1,NSTATE
    do ii=1, NSTATE
      currOVLP(ii,i) = currOVLP(ii,i)*prevPHASE(i)
    enddo
  enddo

  if (tullySubVerb) then
    write(u6,*) '<t-dt|t> RASSI Overlap Columns Corrected'
    do i=1,NSTATE
      write(u6,*) (currOVLP(i,j),j=1,NSTATE,1)
    end do
  end if

  ! Correct rows of current RASSI using currPHASE

  do i=1,NSTATE
    do ii=1, NSTATE
      currOVLP(i,ii) = currOVLP(i,ii)*currPHASE(i)
    enddo
  enddo

  write(u6,*) '<t-dt|t> RASSI Overlap phase Corrected'
  do i=1,NSTATE
    write(u6,*) (currOVLP(i,j),j=1,NSTATE,1)
  end do

  ! Swap back rows if root flipping occured

  do i=1,NSTATE
    do ii=1,NSTATE
      saveOVLP_bk(i,ii) = currOVLP(stateorder(i),ii)
    enddo
  enddo

  ! Swap columns if root flipping occured to give final 'to save' RASSI

  do i=1,NSTATE
    do ii=1,NSTATE
      saveOVLP(ii,i) = saveOVLP_bk(ii,stateorder(i))
    enddo
  enddo

  write(u6,*) '<t-dt|t> RASSI Overlap to save (flipping incl.)'
  do i=1,NSTATE
    write(u6,*) (currOVLP(i,j),j=1,NSTATE,1)
  end do
end if
!                                                                       !
!                         end of sign corrector                         !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! create D matrix normal tully

if_normaltully: if (normalTully) then

  write(u6,*) ''
  write(u6,*) 'Executing Normal Tully !!'
  write(u6,*) ''

  if (.not. rassi_ovlp) then
    write(u6,*) 'Using CI vector product to calculate D matrix'
    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          Dmatrix(i,j) = Zero
          do ii=1,NCI
            Dmatrix(i,j) = Dmatrix(i,j)+CIBigArray(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
          end do
          Dmatrix(i,j) = -Dmatrix(i,j)/DT
        else
          Dmatrix(i,i) = Zero
        end if
      end do
    end do
  else
    write(u6,*) 'Using RASSI overlap to calculate D matrix'
    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          Dmatrix(i,j) = currOVLP(i,j)/DT
        else
          Dmatrix(i,i) = Zero
        end if
      end do
    end do
  end if

  normalTully = .false.

  do i=1,NSTATE
    do j=1,NSTATE
      ExtrInter(i,j) = Dmatrix(i,j)
      ExtrSlope(i,j) = Zero
    end do
    VenergyInter(i) = VenergyP(i)
    VenergySlope(i) = (Venergy(i)-VenergyP(i))/DT
  end do

else if_normaltully

  ! Create D matrix according to Hammes-Schiffer-Tully (interpolating extrapolating)

  ! D32matrix
  if (.not. rassi_ovlp) then
    write(u6,*) 'Using CI vector product to calculate D Matrix'
    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          D32matrix(i,j) = Zero
          do ii=1,NCI
            D32matrix(i,j) = D32matrix(i,j)+CIBigArrayPP(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
            D32matrix(i,j) = D32matrix(i,j)-CIBigArrayP(NCI*(i-1)+ii)*CIBigArrayPP(NCI*(j-1)+ii)
          end do
          D32matrix(i,j) = D32matrix(i,j)/(2*DT)
        else
          D32matrix(i,i) = Zero
        end if
      end do
    end do

    ! D12matrix

    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          D12matrix(i,j) = Zero
          do ii=1,NCI
            D12matrix(i,j) = D12matrix(i,j)+CIBigArrayP(NCI*(i-1)+ii)*CIBigArray(NCI*(j-1)+ii)
            D12matrix(i,j) = D12matrix(i,j)-CIBigArray(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
          end do
          D12matrix(i,j) = D12matrix(i,j)/(2*DT)
        else
          D12matrix(i,i) = Zero
        end if
      end do
    end do

  else ! RASSI overlaps

    write(u6,*) 'Using RASSI overlap to calculate D matrix'

    ! D32matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          D32matrix(i,j) = (prevOVLP(i,j)-prevOVLP(j,i))/(2*DT)
        else
          D32matrix(i,j) = Zero
        end if
      end do
    end do

    ! D12matrix
    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= j) then
          D12matrix(i,j) = (currOVLP(i,j)-currOVLP(j,i))/(2*DT)
        else
          D12matrix(i,j) = Zero
        end if
      end do
    end do

  end if

  ! definition of Y intercept (ExtrInter) and slope (ExtrSlope) for EXTRapolation line

  do i=1,NSTATE
    do j=1,NSTATE
      ExtrSlope(i,j) = (D12matrix(i,j)-D32matrix(i,j))/DT
      ExtrInter(i,j) = D12matrix(i,j)
    end do
    VenergyInter(i) = VenergyP(i)
    VenergySlope(i) = (Venergy(i)-VenergyP(i))/DT
  end do

end if if_normaltully

! UNCOMMENT to print coefficients !!!
! Just a few coefficients
! write(u6,*) 'WaveFunctionsCoefficients are: ',NCI,'*',NSTATE
! write(u6,*) '       This step        Previous step:        PP:'

! write(u6,*) CIBigArray(1),CIBigArrayP(1),CIBigArrayPP(1)
! write(u6,*) CIBigArray(2),CIBigArrayP(2),CIBigArrayPP(2)
! write(u6,*) CIBigArray(3),CIBigArrayP(3),CIBigArrayPP(3)
! write(u6,*) CIBigArray(4),CIBigArrayP(4),CIBigArrayPP(4)

! All coefficients
!do i=1,NCI*NSTATE
!  write(u6,*) CIBigArray(i),CIBigArrayP(i),CIBigArrayPP(i)
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
  LO = Random_Molcas(iseed)
  call put_iscalar('Seed',iseed)
  if (LO < RandThreshold) then
    LO = RandThreshold
  end if
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

substeps: do ii=1,NSUBSTEPS
  if (.not. fixedrandL) then
    LO = Random_Molcas(iseed)
    call put_iscalar('Seed',iseed)
    if (LO < RandThreshold) then
      LO = RandThreshold
    end if
  end if

  do i=1,NSTATE
    do j=1,NSTATE
      Dmatrix(i,j) = ExtrInter(i,j)+ExtrSlope(i,j)*(((ii-1)*DT/NSUBSTEPS)-DT/2)
      if (i /= j) then
        V(i,j) = Zero
      else
        V(i,j) = VenergyInter(i)+VenergySlope(i)*(ii-1)*DT/NSUBSTEPS
      end if
    end do
  end do

  do i=1,NSTATE
    do j=1,NSTATE
      AmatrixDT(i,j) = (Zero,Zero)
      do l=1,NSTATE
        AmatrixDT(i,j) = AmatrixDT(i,j)+Amatrix(i,l)*Dmatrix(l,j)-Amatrix(l,j)*Dmatrix(i,l)
        AmatrixDT(i,j) = AmatrixDT(i,j)+(Zero,One)*Amatrix(i,l)*V(l,j)-(Zero,One)*Amatrix(l,j)*V(i,l)
      end do
    end do
  end do

  do i=1,NSTATE
    do j=1,NSTATE
      Amatrix(i,j) = Amatrix(i,j)+AmatrixDT(i,j)*DT/NSUBSTEPS
    end do
  end do

  do i=1,NSTATE
    do j=1,NSTATE
      if (i == j) then
        !B(i,i) not used
        Bmatrix(i,j) = 2*aimag(conjg(Amatrix(i,j))*V(i,j))
      else
        Bmatrix(i,j) = -2*real(conjg(Amatrix(i,j))*Dmatrix(i,j))
      end if
    end do
  end do

  do i=1,NSTATE
    if (i /= temproot) then
      Gprobab(i) = Bmatrix(i,temproot)*DT/(real(Amatrix(temproot,temproot))*NSUBSTEPS)
      if (Gprobab(i) < Zero) then
        Gprobab(i) = Zero
      end if
    else
      Gprobab(i) = Zero
    end if
  end do

  SumProb = Zero
  do i=1,NSTATE
    if (i /= temproot) then
      SumProb = SumProb+Gprobab(i)
      if (LO <= SumProb) then
        write(u6,*) 'Following Tully, should hop from',temproot,'to',i
        if (V(i,i) < Etot) then
          write(u6,*) 'this root has an energy lower than the total',Etot,' (thus, permitted)'
          Ediffcheck = (V(i,i)-V(temproot,temproot))
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
    end if
  end do

  ! Persico-Granucci all in a THEN branch

  pg: if (decoherence) then
    EKIN = Etot-V(temproot,temproot)
    if (EKIN <= Zero) then
      write(u6,*) 'WARNING! Negative Kinetic Energy. Ekin= ',EKIN,' a.u.'
      write(u6,*) 'Kinetic energy rescaled to 10 e-5.'
      EKIN = 0.00001_wp
    end if
    do i=1,NSTATE
      if (i /= temproot) then
        TAU(i) = abs(One/(V(temproot,temproot)-V(i,i)))*(One+DECO/EKIN)
      end if
    end do

    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= temproot .and. j /= temproot) then
          Amatrix(i,j) = Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(i))*exp(-(DT/NSUBSTEPS)/TAU(j))
        end if
      end do
    end do

    populOS = Zero
    do i=1,NSTATE
      if (i /= temproot) then
        populOS = populOS+real(Amatrix(i,i))
      end if
    end do

    ArelaxPrev = Amatrix(temproot,temproot)
    Amatrix(temproot,temproot) = One-populOS

    do i=1,NSTATE
      do j=1,NSTATE
        if (i /= temproot .and. j == temproot) then
          Amatrix(i,j) = Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(i))*sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
        else if (i == temproot .and. j /= temproot) then
          Amatrix(i,j) = Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(j))*sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
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
      write(u6,*) (Dmatrix(i,j),j=1,NSTATE,1)
    end do
    write(u6,*) 'Density Matrix elements (i=1..#states):'
    do i=1,NSTATE
      write(u6,*) (Amatrix(i,j),j=1,NSTATE,1)
    end do
    write(u6,*) 'B Matrix:'
    do i=1,NSTATE
      write(u6,*) (Bmatrix(i,j),j=1,NSTATE,1)
    end do
    write(u6,*) 'Probabilities:'
    write(u6,*) (Gprobab(j),j=1,NSTATE,1)
    write(u6,*) 'Random Number is: ',LO
    write(u6,*) 'Populations Energies:',(Popul(j),j=1,NSTATE,1),(V(j,j),j=1,NSTATE,1),V(temproot,temproot)
    write(u6,*) ''
  end if

end do substeps

!                                                                      !
!                        END OF INTEGRATOR                             !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(u6,*) 'Gnuplot:',(Popul(j),j=1,NSTATE,1),(Venergy(j),j=1,NSTATE,1),Venergy(temproot)
!write(u6,*) 'Gnuplot:',(Popul(j),j=1,NSTATE,1),(V(j,j),j=1,NSTATE,1),V(temproot,temproot)

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
  write(u6,'(6X,120A1)') '+',('-',i=1,118),'+'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,A,2(47X,A))') '|','A HOP event is detected!|'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,A,44X,2(A,I3,4X),40X,A)') '|','From state:',iRlxRoot,'To state:',ISTATE2,'|'
  write(u6,'(6X,A,118X,A)') '|','|'
  write(u6,'(6X,120A1,//)') '+',('-',i=1,118),'+'

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
!do i=1,nsAtom
!  do j=1,3
!    vel(3*(i-1)+j) = scalfac*vel(3*(i-1)+j)
!  end do
!end do
!write(u6,*) 'Velocities after Hop:'
!do i=1,nsAtom
!  write(u6,*) vel(i*3-2),vel(i*3-1),vel(i*3)
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
call put_zarray('AmatrixV',Amatrix,NSTATE*NSTATE)

if (rassi_ovlp) then
  call Put_dArray('SH_Ovlp_Save',saveOVLP,NSTATE*NSTATE)
  call Put_dArray('Old_Phase',currPHASE,NSTATE)
end if
!                                                                      !
!                           END SAVING                                 !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

return

end subroutine Tully
