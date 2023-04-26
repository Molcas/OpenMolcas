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

subroutine Alaska_Super_Driver(iRC)

use Alaska_Info, only: Auto, DefRoot, ForceNAC, iRlxRoot
use Para_Info, only: nProcs
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iRC
#include "warnings.h"
#include "nac.fh"
integer(kind=iwp) :: Columbus, iGo, iMp2Prpt, iPL, iReturn, istatus, LuInput, LuSpool, LuSpool2, nGrad, nsAtom, nSym
logical(kind=iwp) :: Do_Cholesky, Numerical, Do_DF, Do_ESPF, StandAlone, Exists, Do_Numerical_Cholesky, Do_1CCD, MCLR_Ready
character(len=128) :: FileName
character(len=180) :: Line
character(len=80) :: KSDFT
character(len=16) :: mstate1, mstate2, StdIn
character(len=8) :: Method
real(kind=wp), allocatable :: Grad(:)
integer(kind=iwp), external :: iPrintLevel, isFreeUnit
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
!
! Check the input for the numerical keyword

! copy STDINP to LuSpool

LuSpool = 37
call SpoolInp(LuSpool)
call Chk_Numerical(LuSpool,Numerical)

! Call the numerical procedure if numerical option is available.
! Otherwise hope that the analytic code know how to handle the
! case.

call Get_cArray('Relax Method',Method,8)
call Get_iScalar('Columbus',Columbus)
!                                                                      *
!***********************************************************************
!                                                                      *
! Default for Cholesky or RI/DF is to do the numerical procedure.
! However, for pure DFT we have analytic gradients.

call DecideOnCholesky(Do_Cholesky)
call DecideOnDF(Do_DF)
call DecideOn1CCD(Do_1CCD)
call Get_iScalar('NSYM',nSym)

! Default for MBPT2 is the numerical procedure but if variational
! densities are calculated analytical gradients shall be used.
iMp2Prpt = 0
if (Method == 'MBPT2   ') then
  call Get_iScalar('mp2prpt',iMp2Prpt)

  ! Make sure that the analytic procedure is used if possible.

  if ((nSym == 1) .and. (iMp2Prpt /= 2) .and. (.not. Numerical)) then
    call WarningMessage(2,'Error in Alaska_Super_Driver')
    write(u6,*) 'Alaska: the MBPT2 module was run without the Grdt option!'
    write(u6,*) '   Correct the input and restart the calculation!'
    call Abend()
  end if
else if (Method == 'CASPT2') then
  call Get_iScalar('mp2prpt',iMp2Prpt)
end if

Do_Numerical_Cholesky = Do_Cholesky .or. Do_DF

if ((Method == 'KS-DFT  ') .and. Do_Numerical_Cholesky) then
  call Get_cArray('DFT functional',KSDFT,80)

  !   RI/DF                         1C-CD
  if (Do_DF .or. (Do_Cholesky .and. Do_1CCD .and. (nSym == 1))) then
    Do_Numerical_Cholesky = .false.
  end if

end if

if ((Do_DF .or. (Do_Cholesky .and. Do_1CCD .and. (nSym == 1)))) then

  if ((Method == 'KS-DFT  ') .or. (Method == 'UHF-SCF ') .or. (Method == 'RHF-SCF ') .or. (Method == 'CASSCF  ') .or. &
      (Method == 'RASSCF  ') .or. (Method == 'GASSCF  ') .or. (Method == 'DMRGSCF ') .or. (Method == 'CASSCFSA') .or. &
      (Method == 'RASSCFSA') .or. (Method == 'CASPT2  ') .or. (Method == 'MCPDFT  ') .or. (Method == 'MSPDFT  ')) then
    Do_Numerical_Cholesky = .false.
  else if ((Method == 'MBPT2   ') .and. (nSym == 1)) then
    Do_Numerical_Cholesky = .false.
    !write(u6,*) 'Do Numerical', Do_Numerical_Cholesky
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call DecideOnESPF(Do_ESPF)
! No ESPF for NAC, right?
!if (isNAC) Do_ESPF = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
if (Method == 'DMRGSCFS') then
  call Get_iScalar('SA ready',iGo)
end if

if (Numerical .or. Do_Numerical_Cholesky .or. (Method == 'GASSCFSA') .or. ((Method == 'DMRGSCFS') .and. (iGo /= 2)) .or. &
    ((Method == 'CASPT2') .and. (iMp2Prpt /= 2)) .or. ((Method == 'MBPT2') .and. (iMp2Prpt /= 2)) .or. &
    (Method == 'CCSDT') .or. (Method == 'EXTERNAL')) then
  if (isNAC) then
    call Store_Not_Grad(0,NACstates(1),NACstates(2))
    call WarningMessage(2,'Numerical nonadiabatic coupling not implemented')
    if (Auto) then
      iRC = 0
      return
    else
      call Abend()
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Numerical gradients to be used!
  ! Alaska will automatically generate the input for CASPT2_Gradient
  ! and signal to AUTO (iRC=2) to run the input file Stdin.x.

  if (iPL >= 3) then
    write(u6,*)
    write(u6,*) ' Alaska requests the Numerical_Gradient module to be executed!'
    write(u6,*)
  end if

  LuInput = 11
  LuInput = IsFreeUnit(LuInput)
  call StdIn_Name(StdIn)
  call Molcas_Open(LuInput,StdIn)

  write(LuInput,'(A)') '>ECHO OFF'
  write(LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
  write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

  write(LuInput,'(A)') ' &NUMERICAL_GRADIENT &End'
  write(LuInput,'(A)') 'End of Input'
  write(LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
  write(LuInput,'(A)') '>ECHO ON'
  close(LuInput)
  call Finish(_RC_INVOKED_OTHER_MODULE_)
!                                                                      *
!***********************************************************************
!                                                                      *
!   These do not work in parallel. Warn and stop early, better than
!   crash or give wrong results
!
!else if (Do_Cholesky .and. (Method == 'CASSCFSA') .and. (nProcs > 1)) then
!  call WarningMessage(2,'Error in Alaska_Super_Driver')
!  write(u6,*) 'RI SA-CASSCF analytical gradients do not work correctly in parallel (yet).'
!  call Abend()
else if (Do_Cholesky .and. (Method == 'MBPT2') .and. (nProcs > 1)) then
  call WarningMessage(2,'Error in Alaska_Super_Driver')
  write(u6,*) 'RI MBPT2 analytical gradients do not work correctly in parallel (yet).'
  call Abend()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((Method == 'CASSCFSA') .or. (Method == 'RASSCFSA') .or. ((Method == 'DMRGSCFS') .and. (iGo /= 2)) .or. &
         (Method == 'CASPT2  ')) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! State-Average CASSCF / DMRGSCF

  call Get_iScalar('SA ready',iGo)
  call Get_iScalar('Relax CASSCF root',iRlxRoot)

  if (iRlxRoot == 0) iRlxRoot = 1
  if (isNAC) then
    write(mstate1,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
  else
    write(mstate1,'(I16)') iRlxRoot
  end if

  ! iGo = -1 non-equivalent multi state SA-CASSCF
  ! iGo = 0  equivalent multi state SA-CASSCF
  ! iGo = 2  single root SA-CASSCF
  mstate2 = ''
  if (iGo /= 2) then
    call Get_cArray('MCLR Root',mstate2,16)
  end if

  ! If an explicit root was requested in MCLR and none in ALASKA,
  ! go for it
  if (DefRoot) then
    if (mstate2(1:1) == '+') then
      mstate1 = mstate2
      if (index(mstate2,',') /= 0) then
        read(mstate2,'(1X,I7,1X,I7)') NACStates(1),NACStates(2)
        ForceNAC = .true.
      else
        read(mstate2,'(1X,I15)') iRlxRoot
      end if
    end if
  end if
  mstate1(1:1) = mstate2(1:1)
  MCLR_Ready = (iGO == 1) .and. (mstate1 == mstate2)

  if (MCLR_Ready .or. (iGO > 1)) then
    call Alaska(LuSpool,iRC)

    ! Add ESPF contribution

    if (Do_ESPF) then
      StandAlone = .false.
      call ESPF(iReturn,StandAlone)
      if (iReturn /= 0) then
        call WarningMessage(2,'Error in Alaska_Super_Driver')
        write(u6,*) 'Alaska: ESPF finish with non-zero return code!'
        call Abend()
      end if
    end if
    ! Reset iGO to 0 to allow for new MCLR/ALASKA calculations
    if (iGo == 1) iGo = 0
    call Put_iScalar('SA ready',iGo)
  else if (iGO == -1) then
    call WarningMessage(2,'Error in Alaska_Super_Driver')
    write(u6,*) 'Gradients not implemented for SA-CASSCF with non-equivalent weights!'
    call Abend()
  else
    if (iPL >= 3) then
      write(u6,*)
      write(u6,*) ' Alaska requests MCLR to be run before it starts again!'
      write(u6,*)
    end if

    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_open(LuInput,StdIn)

    write(LuInput,'(A)') '>ECHO OFF'
    write(LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
    write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

    write(LuInput,'(A)') ' &MCLR &End'
    if (isNAC) then
      write(LuInput,'(A)') 'NAC'
      write(LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
    end if
    write(LuInput,'(A)') 'End of Input'
    write(LuInput,'(A)') ' '

    FileName = 'ALASKINP'
    call f_inquire(Filename,Exists)

    if (Exists) then
      LuSpool2 = 77
      LuSpool2 = IsFreeUnit(LuSpool2)
      call Molcas_Open(LuSpool2,Filename)

      do
        read(LuSpool2,'(A)',iostat=istatus) Line
        if (istatus > 0) call Abend()
        if (istatus < 0) exit
        write(LuInput,'(A)') Line
      end do

      close(LuSpool2)

    else

      write(LuInput,'(A)') ' &Alaska &End'
      !write(LuInput,'(A)') 'Show'
      write(LuInput,'(A)') 'CutOff'
      write(LuInput,'(A)') '1.0D-7'
      write(LuInput,'(A)') 'End of Input'

    end if

    write(LuInput,'(A)') '>RM -FORCE $Project.MckInt'
    write(LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
    write(LuInput,'(A)') '>ECHO ON'
    close(LuInput)
    call Finish(_RC_INVOKED_OTHER_MODULE_)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((Method == 'MCPDFT') .or. (Method == 'MSPDFT')) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! MC-PDFT calculation

  Do_ESPF = .false.
  call Get_iScalar('SA ready',iGo)
  call Get_iScalar('Relax CASSCF root',iRlxRoot)

  ! Andrew - I need to identify the root and make sure it is not a
  ! state averaged calculation.  iGo=1 means do MCLR
  ! iGo=99 means the potentials were not calculated during the
  ! MCPDFT step, which is required for analytic gradients.
  if (iGO == 99) then
    call WarningMessage(2,'Error in Alaska_Super_Driver')
    write(u6,*) 'MC-PDFT was run without the GRADient keyword.  Analytic gradients require this keyword.  Please use the '// &
                'GRADient keyword in the preceeding MC-PDFT step.'
    call Abend()
  end if

  if (iRlxRoot == 0) iRlxRoot = 1
  if (isNAC) then
    write(mstate1,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
  else
    write(mstate1,'(I16)') iRlxRoot
  end if

  ! iGo = -1 non-equivalent multi state SA-CASSCF
  ! iGo = 0  equivalent multi state SA-CASSCF
  ! iGo = 2  single root SA-CASSCF
  mstate2 = ''
  if (iGo /= 2) then
    call Get_cArray('MCLR Root',mstate2,16)
  end if

  ! If an explicit root was requested in MCLR and none in ALASKA,
  ! go for it
  if (DefRoot) then
    if (mstate2(1:1) == '+') then
      mstate1 = mstate2
      if (index(mstate2,',') /= 0) then
        read(mstate2,'(1X,I7,1X,I7)') NACStates(1),NACStates(2)
        ForceNAC = .true.
      else
        read(mstate2,'(1X,I15)') iRlxRoot
      end if
    end if
  end if
  mstate1(1:1) = mstate2(1:1)
  MCLR_Ready = (iGO == 1) .and. (mstate1 == mstate2)

  if (MCLR_Ready .or. (iGO > 1)) then
    call Alaska(LuSpool,iRC)

    ! Add ESPF contribution

    !if (Do_ESPF) Then
    !  StandAlone = .false.
    !  call ESPF(iReturn,StandAlone)
    !  if (iReturn /= 0) then
    !     call WarningMessage(2,'Error in Alaska_Super_Driver')
    !     write(u6,*) 'Alaska: ESPF finish with non-zero return code!'
    !     call Abend()
    !  end if
    !end if
    ! Reset iGO to 0 to allow for new MCLR/ALASKA calculations
    if (iGo == 1) iGo = 0
    call Put_iScalar('SA ready',iGo)
  else if (iGO == -1) then
    call WarningMessage(2,'Error in Alaska_Super_Driver')
    write(u6,*) 'Gradients not implemented for SA-CASSCF with non-equivalent weights!'
    call Abend()
  else
    if (iPL >= 3) then
      write(u6,*)
      write(u6,*) ' Alaska requests MCLR to be run before it starts again!'
      write(u6,*)
    end if

    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_open(LuInput,StdIn)

    write(LuInput,'(A)') '>ECHO OFF'
    write(LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
    write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

    write(LuInput,'(A)') ' &MCLR &End'
    write(LuInput,'(A)') ' PRINT = 100'
    if (isNAC) then
      write(LuInput,'(A)') 'NAC'
      write(LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
    end if
    write(LuInput,'(A)') 'End of Input'
    write(LuInput,'(A)') ' '

    FileName = 'ALASKINP'
    call f_inquire(Filename,Exists)

    if (Exists) then
      LuSpool2 = 77
      LuSpool2 = IsFreeUnit(LuSpool2)
      call Molcas_Open(LuSpool2,Filename)

      do
        read(LuSpool2,'(A)',iostat=istatus) Line
        if (istatus > 0) call Abend()
        if (istatus < 0) exit
        write(LuInput,'(A)') Line
      end do

      close(LuSpool2)

    else

      write(LuInput,'(A)') ' &Alaska &End'
      !write(LuInput,'(A)') 'Show'
      write(LuInput,'(A)') 'CutOff'
      write(LuInput,'(A)') '1.0D-7'
      write(LuInput,'(A)') 'End of Input'

    end if

    write(LuInput,'(A)') '>RM -FORCE $Project.MckInt'
    write(LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
    write(LuInput,'(A)') '>ECHO ON'
    close(LuInput)
    call Finish(_RC_INVOKED_OTHER_MODULE_)

  end if

!*********** columbus interface ****************************************

else if ((method == 'MR-CISD ') .and. (Columbus == 1)) then

  ! COLUMBUS MR-CI gradient:
  !     effective density matrix on RUNFILE 'densao_var'
  !     effective fock matrix  on RUNFILE   'FockOcc'
  !     effective  D2          on GAMMA
  !     G_TOC indices          on binary file gtoc
  !     RUNFILE read in drvh1
  !     GAMMA read in pget0
  !     gtoc read in drvg1

  call Get_iScalar('Relax CASSCF root',iRlxRoot)
  call Alaska(LuSpool,iRC)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Read the root, as it could be a CASSCF state-specific excited
  ! state calculation

  iRlxRoot = 0
  call Qpg_iScalar('Relax CASSCF root',Exists)
  if (Exists) call Get_iScalar('Relax CASSCF root',iRlxRoot)
  if (iRlxRoot == 0) iRlxRoot = 1

  ! Go ahead and compute the gradients

  call Alaska(LuSpool,iRC)

  ! Add ESPF contribution

  if (Do_ESPF) then
    StandAlone = .false.
    call ESPF(iReturn,StandAlone)
    if (iReturn /= 0) then
      call WarningMessage(2,'Error in Alaska_Super_Driver')
      write(u6,*) 'Alaska: ESPF finish with non-zero return code!'
      call Abend()
    end if
  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the gradient and store in the GRADS file
! It is done here because the gradient could have been modified by ESPF
! and we do not want to pass the root to ESPF (yet)

call Get_iScalar('Unique atoms',nsAtom)
call mma_Allocate(Grad,3*nsAtom,Label='Grad')
nGrad = 3*nsAtom
call Get_dArray_chk('GRAD',Grad,nGrad)
if (isNAC) then
  call Store_Grad(Grad,nGrad,0,NACstates(1),NACstates(2))
else
  call Store_Grad(Grad,nGrad,iRlxRoot,0,0)
end if
call mma_deallocate(Grad)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Alaska_Super_Driver
