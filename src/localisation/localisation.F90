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
! Copyright (C) Yannick Carissan                                       *
!               2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Localisation(iReturn)
! Author: Yannick Carissan
!
! Modifications:
!    - October 10, 2005 (Thomas Bondo Pedersen):
!      completely restructed; introduce Boys and Cholesky
!      localisations.
!    - December 2005 / January 2006 (Thomas Bondo Pedersen):
!      Edmiston-Ruedenberg, PAO, and pair domain analysis included.

use Localisation_globals, only: AnaAtom, Analysis, AnaPAO, AnaPAO_Save, BName, CMO, DoCNOs, DoDomain, EOrb, EvalER, Ind, iWave, &
                                LC_FileOrb, LocCanOrb, LocModel, LocNatOrb, LocPAO, LuSpool, MOrig, NamAct, nBas, nCMO, nFro, &
                                nOrb, nOrb2Loc, nSym, Occ, Order, PrintMOs, Silent, Skip, Test_Localisation, Timing, Wave
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
#include "debug.fh"
integer(kind=iwp) :: ibo, iCheck, icHour, icMin, IndT(7,8), iOff, iPRway, irc, iSym, iTol, iUHF, iwHour, iwMin, j, jbo, jInd, &
                     jPrt, jTyp, k, kCMO, kEor, kIndT, kOcc, lMOrig, lOff, LU_, nbo
real(kind=wp) :: AddInfoVal, C1, C1_Loc, C2, C2_Loc, CPUtot, cSec, ERFun(2), Functional, W1, W1_Loc, W2, W2_Loc, WLLtot, wSec, &
                 xnr0(8), xnr1, xNrm
character(len=180) :: Line
character(len=80) :: Title, Txt
character(len=20) :: NameFile
character(len=15) :: AddInfoString
character(len=7) :: matname
character(len=6) :: Filename
character(len=4) :: Model
character(len=2) :: PreFix
real(kind=wp), allocatable :: CMO2(:), CMO3(:), jXarray(:)
character(len=*), parameter :: SecNam = 'Localisation'
integer(kind=iwp), external :: isFreeUnit !vv , LocUtil_Models
real(kind=wp), external :: ddot_
character(len=180), external :: Get_Ln

! Start timing.
! -------------

call CWTime(C1,W1)

! Print banner.
! -------------

!call Banner_Localisation()

! Read basic info from runfile and INPORB.
! ----------------------------------------

! Quick and dirty read of the FileOrb name before INPORB is opened
LuSpool = 17
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'LOCALISATION')
LC_FileOrb = ' '
do
  Line = Get_Ln(LuSpool)
  call UpCase(Line)
  if (Line(1:4) == 'FILE') then
    Line = Get_Ln(LuSpool)
    call FileOrb(Line,LC_FileOrb)
    exit
  else if (Line(1:4) == 'END ') then
    exit
  end if
end do
call Close_LuSpool(LuSpool)

call GetInfo_Localisation_0()

! Read and process input.
! -----------------------

call ReadInp_Localisation()

! Check that all is OK.
! ---------------------

irc = 0
call Chk_Input(irc)
if (irc > 0) then
  call SysAbendMsg(SecNam,'Inconsistent input',' ')
else if (irc < 0) then
  iReturn = 0
  call Error(0) ! nothing to do; exit...
  return
end if

! If test option was specified, we need to keep a copy of the
! original MOs.
! -----------------------------------------------------------

if (Test_Localisation .or. Analysis .or. AnaPAO .or. AnaPAO_Save .or. LocNatOrb .or. LocCanOrb) then
  lMOrig = nCMO
  call mma_allocate(MOrig,lMOrig,label='MOrig')
  MOrig(:) = CMO(:)
end if

! Print initial orbitals.
! -----------------------

if ((.not. Silent) .and. PrintMOs) then
  write(Title,'(80x)')
  write(Title,'(a)') 'Initial MO''s'
  iPrWay = 2 ! short
  call PriMO_Localisation(Title,.true.,.true.,-One,1.0e5_wp,nSym,nBas,nOrb,BName,EOrb,Occ,CMO,iPrWay,Ind)
end if

if (Debug) then
  write(u6,'(A,A,I2)') SecNam,': debug info at start:'
  write(u6,'(A,8I8)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
  write(u6,'(A,8I8)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
  write(u6,'(A,8I8)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
  write(u6,'(A,8I8)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
end if

! Evaluate ER functional for initial orbitals.
! --------------------------------------------

if (EvalER) then
  ERFun(1) = Zero
  call ComputeFuncER(ERFun(1),CMO,nBas,nOrb2Loc,nFro,nSym,Timing)
end if

! Localise orbitals (if the user did not request us to skip it).
! --------------------------------------------------------------

if (Skip) then
  write(u6,'(//,5X,A,//)') 'NOTICE: LOCALISATION PROCEDURE SKIPPED BY USER REQUEST!'
  AddInfoString = 'SKIPPED LOCALI '
  AddInfoVal = Zero
  iTol = 15
else
  if (Timing) call CWTime(C1_Loc,W1_Loc)
  if ((LocModel == 1) .or. (LocModel == 2) .or. (LocModel == 4)) then
    if (LocModel == 1) then
      Model = 'Pipe'
      AddInfoString = 'PIPEKFUNCTIONAL'
    else if (LocModel == 2) then
      Model = 'Boys'
      AddInfoString = 'BOYSFUNCTIONAL '
    else if (LocModel == 4) then
      Model = 'Edmi'
      AddInfoString = 'ERFUNCTIONAL   '
    end if
    irc = 0
    call Localise_Iterative(irc,Model,Functional)
    if (irc /= 0) then
      write(Txt,'(A,I3)') 'Return code from Localise_Iterative:',irc
      call SysAbendMsg(SecNam,'Localisation failed!',Txt)
    end if
    AddInfoVal = Functional
    iTol = 4
  else if (LocModel == 3) then
    if (LocPAO) then
      Model = 'PAO '
      AddInfoString = 'CHOLESKY PAO   '
    else
      Model = 'Chol'
      AddInfoString = 'CHOLESKY NORM  '
    end if
    irc = 0
    call Localise_Noniterative(irc,Model,xNrm)
    if (irc /= 0) then
      write(Txt,'(A,I3)') 'Return code from Localise_Noniterative:',irc
      call SysAbendMsg(SecNam,'Localisation failed!',Txt)
    end if
    AddInfoVal = xNrm
    iTol = 4
  else if (Wave) then ! wavelet transform
    irc = 0
    call Wavelet_Transform(irc,CMO,nSym,nBas,nFro,nOrb2Loc,iWave,.false.,xNrm)
    if (irc /= 0) then
      write(Txt,'(A,I3)') 'Return code from Wavelet_Transform:',irc
      call SysAbendMsg(SecNam,'Localisation failed!',Txt)
    end if
    AddInfoVal = xNrm
    iTol = 4
  else if (DoCNOs) then ! Constrained NOs (analysis)
    irc = 0
    call Get_CNOs(irc,nFro,nOrb2Loc,xNrm)
    if (irc /= 0) then
      write(Txt,'(A,I3)') 'Return code from Get_CNOs:',irc
      call SysAbendMsg(SecNam,'Localisation failed!',Txt)
    end if
    AddInfoVal = xNrm
    iTol = 4
  else
    call SysAbendMsg(SecNam,'Logical bug','(LocModel)')
    AddInfoString = '?!?!?!?!?!?!?!?'
    AddInfoVal = -huge(AddInfoVal)
    iTol = 15
  end if
  if (Timing) then
    call CWTime(C2_Loc,W2_Loc)
    write(u6,'(/,1X,A,F10.2,A)') 'CPU  time for localisation procedure:',C2_Loc-C1_Loc,' seconds'
    write(u6,'(1X,A,F10.2,A,/)') 'Wall time for localisation procedure:',W2_Loc-W1_Loc,' seconds'
  end if
end if

! Info for check system.
! ----------------------

call Add_Info(AddInfoString,[AddInfoVal],1,iTol)

! Order local orbitals according to Cholesky ordering.
! ----------------------------------------------------

if (Order) then
  write(u6,'(/,1X,A)') 'Sorting local orbitals according to Cholesky ordering. (Based on overlap U=X^TSC.)'
  call Sort_Localisation(CMO,nBas,nOrb2Loc,nFro,nSym)
end if

! Evaluate ER functional for local orbitals.
! ------------------------------------------

if (EvalER) then
  ERFun(2) = Zero
  call ComputeFuncER(ERFun(2),CMO,nBas,nOrb2Loc,nFro,nSym,Timing)
  write(u6,'(/,1X,A,1P,D15.8,/,1X,A,D15.8,/)') 'ER functional for initial orbitals: ',ERFun(1), &
                                               'ER functional for local   orbitals: ',ERFun(2)
end if

! Test section.
! -------------

if (Test_Localisation) then
  write(u6,'(//,1X,A)') 'Testing orbital localisation...'
  irc = -1
  call TestLoc(irc)
  if (irc == 0) then
    write(u6,*) '...OK!'
  else
    write(u6,*) SecNam,': localisation error detected!'
    write(u6,*) ' TestLoc returned ',irc
    call Error(1) ! exit
    return
  end if
end if

! Analysis.
! ---------

if (Analysis) then
  PreFix = 'L_'
  if (AnaAtom) then
    call BitMap_Localisation_Atom(PreFix)
  else
    call BitMap_Localisation(PreFix)
  end if
end if

! Orbital domains.
! ----------------

if (DoDomain) then
  irc = 0
  call Domain_Localisation(irc)
  if (irc /= 0) then
    write(u6,*) SecNam,': Domain error!'
    write(u6,*) ' Domain_Localisation returned ',irc
    write(u6,*) ' Program continues nevertheless...'
  end if
end if

! Local Natural Orbitals
! ----------------------

if (LocNatOrb .or. LocCanOrb) then
  nbo = 0
  do iSym=1,nSym
    nbo = nbo+nBas(iSym)*nOrb2Loc(iSym)
  end do
  call mma_allocate(CMO2,nbo,label='XCMO')
  call mma_allocate(CMO3,nbo,label='XCMO')

  ibo = 0
  jbo = 1
  do iSym=1,nSym
    kCMO = ibo+nBas(iSym)*nFro(iSym)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),MOrig(kCMO+1),1,CMO2(jbo),1)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),CMO(kCMO+1),1,CMO3(jbo),1)
    ibo = ibo+nBas(iSym)**2
    jbo = jbo+nBas(iSym)*nOrb2Loc(iSym)
  end do

  if (LocNatOrb) then
    matname = 'Density'
    call mma_allocate(jXarray,size(Occ),label='jXarray')
    jXarray(:) = Occ(:)
  else
    matname = 'Fock'
    call mma_allocate(jXarray,size(EOrb),label='jXarray')
    jXarray(:) = EOrb(:)
  end if
  lOff = 1
  do iSym=1,nSym
    xnr0(iSym) = ddot_(nOrb2Loc(iSym),[One],0,jXarray(lOff+nFro(iSym)),1)
    lOff = lOff+nBas(iSym)
  end do

  call Loc_Nat_Orb(irc,CMO2,CMO3,jXarray,nOrb2Loc)
  if (irc /= 0) then
    write(u6,*) SecNam,': localisation error detected!'
    write(u6,*) ' Loc_Nat_Orb returned ',irc
    call Error(1) ! exit
    return
  end if

  write(u6,*)
  write(u6,*) ' ------------------------------------------------- '
  write(u6,*) ' Sum of the eigenvalues of the partial ',matname,' '
  write(u6,*) ' matrix of the orbitals that have been localised   '
  write(u6,*) ' ------------------------------------------------- '
  write(u6,*) '    Symm.        before     / after localisation   '
  write(u6,*) ' ------------------------------------------------- '
  lOff = 1
  do iSym=1,nSym
    xnr1 = ddot_(nOrb2Loc(iSym),[One],0,jXarray(lOff+nFro(iSym)),1)
    lOff = lOff+nBas(iSym)
    write(u6,'(3X,I4,8X,F11.5,4X,F11.5)') iSym,xnr0(iSym),xnr1
  end do
  write(u6,*) ' ------------------------------------------------- '
  write(u6,*)

  if (LocNatOrb) then
    Occ(:) = jXarray(:)
  else
    EOrb(:) = jXarray(:)
  end if
  call mma_deallocate(jXarray)

  ibo = 0
  jbo = 1
  do iSym=1,nSym
    kCMO = ibo+nBas(iSym)*nFro(iSym)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),CMO3(jbo),1,CMO(kCMO+1),1)
    ibo = ibo+nBas(iSym)**2
    jbo = jbo+nBas(iSym)*nOrb2Loc(iSym)
  end do
  call mma_deallocate(CMO2)
  call mma_deallocate(CMO3)
end if

!-TBP, July 2010 (in connection with fixing deleted orbitals bug, patch
! 7.7.073_Localisation):
! Moved the following section outside print section, which is only
! executed if PrintMOs=.True., implying that the info on LOCORB
! would differ depending on this print flag!!

! Print a warning if localisation is done in a subset of orbitals
! belonging to more than one subset (e.g. mixing occupied and
! virtual orbitals and thus breaking the variational principle).
! Set orbital energies to zero (unless localisation was skipped).
! Also zero occupations if localisation mixed different subsets.
! --------------------------------------------------------------

iCheck = 0
iOff = 0
iSym = 1
jTyp = min(6,max(2,Ind(nFro(1)+1))) ! Fro=Ina and Del=Vir
do while ((iSym <= nSym) .and. (iCheck == 0))
  jInd = iOff+nFro(iSym)+1
  j = 1
  do while ((j < nOrb2Loc(iSym)) .and. (iCheck == 0))
    iCheck = min(6,max(2,Ind(jInd+j)))-jTyp
    j = j+1
  end do
  iOff = iOff+nBas(iSym)
  iSym = iSym+1
end do
if ((.not. LocCanOrb) .and. (.not. Skip)) then
  iOff = 0
  do iSym=1,nSym
    kEor = iOff+nFro(iSym)+1
    call FZero(EOrb(kEor),nOrb2Loc(iSym))
    iOff = iOff+nBas(iSym)
  end do
end if
if (iCheck == 0) then
  jPrt = 0
  if ((jTyp > 3) .and. (jTyp < 6)) jPrt = 1
else
  write(u6,*) '****  WARNING  ***  WARNING  ***  WARNING  ****'
  write(u6,*) ' The orbitals for which localisation is requested belong to more than one of the subspaces:'
  write(u6,*) ' (Froz+Inac|RAS1|RAS2|RAS3|Sec+Del). '
  write(u6,*) '***********************************************'
  jPrt = 1
end if
if ((.not. LocNatOrb) .and. (.not. Skip) .and. (.not. DoCNOs)) then
  if (jPrt == 1) then
    iOff = 0
    do iSym=1,nSym
      kOcc = iOff+nFro(iSym)+1
      call FZero(Occ(kOcc),nOrb2Loc(iSym))
      iOff = iOff+nBas(iSym)
    end do
  end if
end if

! Print final localised orbitals.
! -------------------------------

if (PrintMOs) then
  write(Title,'(80x)')
  write(Title,'(a)') 'Final localised MO''s'
  iPrWay = 2 ! short
  call PriMO_Localisation(Title,.true.,.true.,-One,1.0e5_wp,nSym,nBas,nOrb,BName,EOrb,Occ,CMO,iPrWay,Ind)
end if

! Write LOCORB file.
! ------------------

write(Namefile,'(A)') 'LOCORB'
write(Title,'(80X)')
write(Title,'(A)') 'Localised orbitals'
LU_ = isFreeUnit(11)
j = 0
IndT(:,:) = 0
do iSym=1,nSym
  do k=1,nOrb(iSym)
    kIndT = Ind(j+k)
    if ((kIndT > 0) .and. (kIndT <= 7)) then
      IndT(kIndT,iSym) = IndT(kIndT,iSym)+1
    else
      call WarningMessage(2,'Localisation: Illegal orbital type')
      write(u6,'(A,I6,A,I2,A,I9)') 'Orbital',k,' of sym.',iSym,' has illegal type:',kIndT
      call Abend()
    end if
  end do
  j = j+nBas(iSym)
end do
call WrVec_Localisation(Namefile,LU_,'COEI',nSym,nBas,nBas,CMO,Occ,EOrb,IndT,Title)
if (.not. Silent) then
  write(u6,'(1X,A)') 'The LOCORB file has been written.'
end if

! Write MOLDEN file.
! ------------------

iUHF = 0
Filename = 'MD_LOC'
call Molden_Interface(iUHF,Namefile,Filename)
if (.not. Silent) then
  write(u6,'(1X,A)') 'The MOLDEN file has been written.'
end if

! Set return code.
! ----------------

iReturn = 0

! Free memory (by flushing).
! --------------------------

call Error(0)

! That's it!
! ----------

contains

subroutine Error(code)
  integer(kind=iwp), intent(in) :: code
  if (code /= 0) iReturn = code

  call mma_deallocate(CMO)
  call mma_deallocate(Occ)
  call mma_deallocate(EOrb)
  call mma_deallocate(Ind)
  call mma_deallocate(BName)
  if (allocated(MOrig)) call mma_deallocate(MOrig)
  if (allocated(NamAct)) call mma_deallocate(NamAct)

  ! Print timing.
  ! -------------

  if (.not. Silent) then
    call CWTime(C2,W2)
    CPUtot = C2-C1
    WLLtot = W2-W1
    call Cho_CnvTim(CPUtot,icHour,icMin,cSec)
    call Cho_CnvTim(WLLtot,iwHour,iwMin,wSec)
    write(u6,'(/,1X,A,I8,A,I2,A,F6.2,A)') '*** Total localisation time (CPU) : ',icHour,' hours ',icMin,' minutes ',cSec, &
                                          ' seconds ***'
    write(u6,'(1X,A,I8,A,I2,A,F6.2,A,/)') '*** Total localisation time (Wall): ',iwHour,' hours ',iwMin,' minutes ',wSec, &
                                          ' seconds ***'
  end if
end subroutine Error

end subroutine Localisation
