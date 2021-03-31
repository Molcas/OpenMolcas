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

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "inflocal.fh"
#include "debug.fh"
#include "real.fh"
character*80 Title, Txt
character NameFile*20, Filename*6
character*180 Line, Get_Ln
external Get_Ln

!vv integer  LocUtil_Models
!vv external LocUtil_Models

character*7 matname
character*2 PreFix
character*4 Model
character*12 SecNam
character*15 AddInfoString
parameter(SecNam='Localisation')

real*8 ERFun(2), xnr0(8)
integer IndT(7,8)

! Start timing.
! -------------

call CWTime(C1,W1)

! Dummy allocation used to flush memory at the end.
! -------------------------------------------------

l_Dum = 1
call GetMem('LocDum','Allo','Real',ip_Dum,l_Dum)

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
100 Line = Get_Ln(LuSpool)
call UpCase(Line)
if (Line(1:4) == 'FILE') goto 200
if (Line(1:4) == 'END ') goto 999
goto 100
200 Line = Get_Ln(LuSpool)
call FileOrb(Line,LC_FileOrb)
999 call Close_LuSpool(LuSpool)

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
  Go To 1 ! nothing to do; exit...
end if

! If test option was specified, we need to keep a copy of the
! original MOs.
! -----------------------------------------------------------

if (Test_Localisation .or. Analysis .or. AnaPAO .or. AnaPAO_Save .or. LocNatOrb .or. LocCanOrb) then
  lMOrig = nCMO
  call GetMem('MOrig','Allo','Real',ipMOrig,lMOrig)
  call dCopy_(lMOrig,Work(ipCMO),1,Work(ipMOrig),1)
else
  ipMOrig = -99999999
  lMOrig = 0
end if

! Print initial orbitals.
! -----------------------

if ((.not. Silent) .and. PrintMOs) then
  write(Title,'(80x)')
  write(Title,'(a)') 'Initial MO''s'
  iPrWay = 2 ! short
  call PriMO_Localisation(Title,.true.,.true.,-1.0d0,1.0d5,nSym,nBas,nOrb,Name,Work(ipEor),Work(ipOcc),Work(ipCMO),iPrWay, &
                          iWork(ipInd))
end if

if (Debug) then
  write(6,'(A,A,I2)') SecNam,': debug info at start:'
  write(6,'(A,8I8)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
  write(6,'(A,8I8)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
  write(6,'(A,8I8)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
  write(6,'(A,8I8)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
end if

! Evaluate ER functional for initial orbitals.
! --------------------------------------------

if (EvalER) then
  ERFun(1) = 0.0d0
  call ComputeFuncER(ERFun(1),Work(ipCMO),nBas,nOrb2Loc,nFro,nSym,Timing)
end if

! Localise orbitals (if the user did not request us to skip it).
! --------------------------------------------------------------

if (Skip) then
  write(6,'(//,5X,A,//)') 'NOTICE: LOCALISATION PROCEDURE SKIPPED BY USER REQUEST!'
  AddInfoString = 'SKIPPED LOCALI '
  AddInfoVal = 0.0d0
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
    call Wavelet_Transform(irc,ipCMO,nSym,nBas,nFro,nOrb2Loc,iWave,.false.,xNrm)
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
    AddInfoVal = -9.9d15
    iTol = 15
  end if
  if (Timing) then
    call CWTime(C2_Loc,W2_Loc)
    write(6,'(/,1X,A,F10.2,A)') 'CPU  time for localisation procedure:',C2_Loc-C1_Loc,' seconds'
    write(6,'(1X,A,F10.2,A,/)') 'Wall time for localisation procedure:',W2_Loc-W1_Loc,' seconds'
  end if
end if

! Info for check system.
! ----------------------

call Add_Info(AddInfoString,[AddInfoVal],1,iTol)

! Order local orbitals according to Cholesky ordering.
! ----------------------------------------------------

if (Order) then
  write(6,'(/,1X,A)') 'Sorting local orbitals according to Cholesky ordering. (Based on overlap U=X^TSC.)'
  call Sort_Localisation(Work(ipCMO),nBas,nOrb2Loc,nFro,nSym)
end if

! Evaluate ER functional for local orbitals.
! ------------------------------------------

if (EvalER) then
  ERFun(2) = 0.0d0
  call ComputeFuncER(ERFun(2),Work(ipCMO),nBas,nOrb2Loc,nFro,nSym,Timing)
  write(6,'(/,1X,A,1P,D15.8,/,1X,A,D15.8,/)') 'ER functional for initial orbitals: ',ERFun(1), &
                                              'ER functional for local   orbitals: ',ERFun(2)
end if

! Test section.
! -------------

if (Test_Localisation) then
  write(6,'(//,1X,A)') 'Testing orbital localisation...'
  irc = -1
  call TestLoc(irc)
  if (irc == 0) then
    write(6,*) '...OK!'
  else
    write(6,*) SecNam,': localisation error detected!'
    write(6,*) ' TestLoc returned ',irc
    iReturn = 1
    Go To 1 ! exit
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
    write(6,*) SecNam,': Domain error!'
    write(6,*) ' Domain_Localisation returned ',irc
    write(6,*) ' Program continues nevertheless...'
  end if
end if

! Local Natural Orbitals
! ----------------------

if (LocNatOrb .or. LocCanOrb) then
  nbo = 0
  do iSym=1,nSym
    nbo = nbo+nBas(iSym)*nOrb2Loc(iSym)
  end do
  call GetMem('XCMO','Allo','Real',jCMO,2*nbo)
  jXMO = jCMO+nbo
!
  ibo = 0
  jbo = 0
  do iSym=1,nSym
    kCMO = ibo+nBas(iSym)*nFro(iSym)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(ipMOrig+kCMO),1,Work(jbo+jCMO),1)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(ipCMO+kCMO),1,Work(jbo+jXMO),1)
    ibo = ibo+nBas(iSym)**2
    jbo = jbo+nBas(iSym)*nOrb2Loc(iSym)
  end do
!
  if (LocNatOrb) then
    matname = 'Density'
    jXarray = ipOcc
  else
    matname = 'Fock'
    jXarray = ipEor
  end if
  lOff = 0
  do iSym=1,nSym
    xnr0(iSym) = ddot_(nOrb2Loc(iSym),[1.0d0],0,Work(jXarray+lOff+nFro(iSym)),1)
    lOff = lOff+nBas(iSym)
  end do
!
  call Loc_Nat_Orb(irc,Work(jCMO),Work(jXMO),Work(jXarray),nOrb2Loc)
  if (irc /= 0) then
    write(6,*) SecNam,': localisation error detected!'
    write(6,*) ' Loc_Nat_Orb returned ',irc
    iReturn = 1
    Go To 1 ! exit
  end if
!
  write(6,*)
  write(6,*) ' ------------------------------------------------- '
  write(6,*) ' Sum of the eigenvalues of the partial ',matname,' '
  write(6,*) ' matrix of the orbitals that have been localised   '
  write(6,*) ' ------------------------------------------------- '
  write(6,*) '    Symm.        before     / after localisation   '
  write(6,*) ' ------------------------------------------------- '
  lOff = 0
  do iSym=1,nSym
    xnr1 = ddot_(nOrb2Loc(iSym),[1.0d0],0,Work(jXarray+lOff+nFro(iSym)),1)
    lOff = lOff+nBas(iSym)
    write(6,'(3X,I4,8X,F11.5,4X,F11.5)') iSym,xnr0(iSym),xnr1
  end do
  write(6,*) ' ------------------------------------------------- '
  write(6,*)
!
  ibo = 0
  jbo = 0
  do iSym=1,nSym
    kCMO = ibo+nBas(iSym)*nFro(iSym)
    call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(jbo+jXMO),1,Work(ipCMO+kCMO),1)
    ibo = ibo+nBas(iSym)**2
    jbo = jbo+nBas(iSym)*nOrb2Loc(iSym)
  end do
  call GetMem('XCMO','Free','Real',jCMO,2*nbo)
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
jTyp = min(6,max(2,iWork(ipInd+nFro(1)))) ! Fro=Ina and Del=Vir
do while ((iSym <= nSym) .and. (iCheck == 0))
  jInd = ipInd+iOff+nFro(iSym)
  j = 1
  do while ((j < nOrb2Loc(iSym)) .and. (iCheck == 0))
    iCheck = min(6,max(2,iWork(jInd+j)))-jTyp
    j = j+1
  end do
  iOff = iOff+nBas(iSym)
  iSym = iSym+1
end do
if ((.not. LocCanOrb) .and. (.not. Skip)) then
  iOff = 0
  do iSym=1,nSym
    kEor = ipEor+iOff+nFro(iSym)
    call FZero(Work(kEor),nOrb2Loc(iSym))
    iOff = iOff+nBas(iSym)
  end do
end if
if (iCheck == 0) then
  jPrt = 0
  if ((jTyp > 3) .and. (jTyp < 6)) jPrt = 1
else
  write(6,*) '****  WARNING  ***  WARNING  ***  WARNING  ****'
  write(6,*) ' The orbitals for which localisation is requested belong to more than one of the subspaces:'
  write(6,*) ' (Froz+Inac|RAS1|RAS2|RAS3|Sec+Del). '
  write(6,*) '***********************************************'
  jPrt = 1
end if
if ((.not. LocNatOrb) .and. (.not. Skip) .and. (.not. DoCNOs)) then
  if (jPrt == 1) then
    iOff = 0
    do iSym=1,nSym
      kOcc = ipOcc+iOff+nFro(iSym)
      call FZero(Work(kOcc),nOrb2Loc(iSym))
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
  call PriMO_Localisation(Title,.true.,.true.,-1.0d0,1.0d5,nSym,nBas,nOrb,Name,Work(ipEor),Work(ipOcc),Work(ipCMO),iPrWay, &
                          iWork(ipInd))
end if

! Write LOCORB file.
! ------------------

write(Namefile,'(A)') 'LOCORB'
write(Title,'(80X)')
write(Title,'(A)') 'Localised orbitals'
LU_ = isFreeUnit(11)
j = ipInd-1
call iZero(IndT,56)
do iSym=1,nSym
  do k=1,nOrb(iSym)
    kIndT = iWork(j+k)
    if ((kIndT > 0) .and. (kIndT <= 7)) then
      IndT(kIndT,iSym) = IndT(kIndT,iSym)+1
    else
      call WarningMessage(2,'Localisation: Illegal orbital type')
      write(6,'(A,I6,A,I2,A,I9)') 'Orbital',k,' of sym.',iSym,' has illegal type:',kIndT
      call Abend()
    end if
  end do
  j = j+nBas(iSym)
end do
call WrVec_Localisation(Namefile,LU_,'COEI',nSym,nBas,nBas,Work(ipCMO),Work(ipOcc),Work(ipEor),IndT,Title)
if (.not. Silent) then
  write(6,'(1X,A)') 'The LOCORB file has been written.'
end if

! Write MOLDEN file.
! ------------------

iUHF = 0
Filename = 'MD_LOC'
call Molden_Interface(iUHF,Namefile,Filename)
if (.not. Silent) then
  write(6,'(1X,A)') 'The MOLDEN file has been written.'
end if

! Set return code.
! ----------------

iReturn = 0

! Free memory (by flushing).
! --------------------------

1 call GetMem('LocDum','Flus','Real',ip_Dum,l_Dum)
call GetMem('LocDum','Free','Real',ip_Dum,l_Dum)

! Print timing.
! -------------

if (.not. Silent) then
  call CWTime(C2,W2)
  CPUtot = C2-C1
  WLLtot = W2-W1
  call Cho_CnvTim(CPUtot,icHour,icMin,cSec)
  call Cho_CnvTim(WLLtot,iwHour,iwMin,wSec)
  write(6,'(/,1X,A,I8,A,I2,A,F6.2,A)') '*** Total localisation time (CPU) : ',icHour,' hours ',icMin,' minutes ',cSec,' seconds ***'
  write(6,'(1X,A,I8,A,I2,A,F6.2,A,/)') '*** Total localisation time (Wall): ',iwHour,' hours ',iwMin,' minutes ',wSec,' seconds ***'
end if

! That's it!
! ----------

end subroutine Localisation
