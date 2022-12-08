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
! Copyright (C) 1992, Markus P. Fuelscher                              *
!               1995, Martin Schuetz                                   *
!               2004,2005, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine RdInp(CMO,Eall,Eocc,Eext,iTst,ESCF)
!***********************************************************************
!                                                                      *
!     Locate input stream and read commands                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     modified by MGS                                                  *
!                 TBP                                                  *
!     University of Lund, Sweden, 1992/95, 2004/05.                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

#include "intent.fh"

use MBPT2_Global, only: DelGhost, DoCholesky, DoDF, DoLDF, iDel, iFro, iPL, NamAct, nBas, nDel1, nDel2, nFro1, nFro2, nTit, &
                        Thr_ghs, Title
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(_OUT_) :: CMO(*), Eall(*), Eocc(*), Eext(*)
integer(kind=iwp), intent(out) :: iTst
real(kind=wp), intent(out) :: ESCF
integer(kind=iwp) :: i, iCom, iCount, iDNG, iDummy(1), iErr, iExt, iLow, iOrb, ip, iostatus, iPrt, iSym, iUpp, j, jCom, jDel, &
                     jFro, jOcc, l_Occup, LC, LEE, LEO, LSQ, Lu_orb, LuSpool, nExtT, nFre, nOccT
logical(kind=iwp) :: FrePrt, ERef_UsrDef, DecoMP2_UsrDef, DNG, NoGrdt, lTit, lFro, lFre, lDel, lSFro, lSDel, lExt, lPrt, LumOrb, &
                     Skip
character(len=4) :: Command
character(len=8) :: emiloop, inGeo
character(len=80) :: VecTitle
character(len=180) :: Line
integer(kind=iwp), allocatable :: SQ(:)
real(kind=wp), allocatable :: C(:), EE(:), EO(:), Occup(:)
character(len=*), parameter :: ComTab(39) = ['TITL','FROZ','DELE','SFRO','SDEL','EXTR','PRIN','TEST','PRPT','LUMO', &
                                             'EREF','VIRA','T1AM','GRDT','LAPL','GRID','BLOC','CHOA','DECO','NODE', &
                                             'THRC','SPAN','MXQU','PRES','CHKI','FORC','VERB','NOVE','FREE','PREC', &
                                             'SOSM','OEDT','OSFA','LOVM','DOMP','FNOM','GHOS','NOGR','END ']
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: ChoMP2_ChkPar, Reduce_Prt
character(len=180), external :: Get_Ln
#include "chomp2_cfg.fh"
#include "corbinf.fh"
#include "warnings.h"
#include "Molcas.fh"

!----------------------------------------------------------------------*
!     Locate "start of input"                                          *
!----------------------------------------------------------------------*
lTit = .false.
lFro = .false.
lFre = .false.
lDel = .false.
lSFro = .false.
lSDel = .false.
lExt = .false.
lPrt = .false.
LovMP2 = .false.
DoMP2 = .false.
FNOMP2 = .false.
LumOrb = .false.
all_vir = .false.
DoT1amp = .false.
ESCF = Zero
Thr_ghs = Half
DelGHost = .false.
ERef_UsrDef = .false.
DecoMP2_UsrDef = .false.
nActa = 0
call Put_iScalar('mp2prpt',0)

! copy input from standard input to a local scratch file

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'MBPT2')
!----------------------------------------------------------------------*
!     Define default values                                            *
!----------------------------------------------------------------------*
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

if (DoCholesky) then
  ChoAlg = 2
else
  ChoAlg = -999999
end if
DecoMP2 = Decom_Def
ThrMP2 = -huge(ThrMP2)
SpanMP2 = Span_Def
MxQualMP2 = MxQual_Def
ChkDecoMP2 = .false.
ForceBatch = .false.
if (iPL >= 3) then
  Verbose = .true.
else
  Verbose = .false.
end if
! Scaled Opposite-Spin MP2
SOS_mp2 = .false.
set_cd_thr = .true.
OED_Thr = 1.0e-8_wp
C_os = 1.3_wp
EOSMP2 = Zero
! Frozen natural orbitals
DoFNO = .false.
! MP2-gradient/1pdm
DoDens = .false.
DoGrdt = .false.
NoGrdt = .false.
NoGamma = .false.
! Laplace transform
Laplace = .false.
Laplace_nGridPoints = 0
Laplace_BlockSize = Laplace_BlockSize_Def
! LDF settings
if (DoLDF) then
  SOS_MP2 = .true.
  Laplace = .true.
end if

nTit = 0
iTst = 0
!nTstP = -1
iPrt = 0
nFro1(:) = 0
nFro2(:) = 0
nDel1(:) = 0
nDel2(:) = 0
call Get_iArray('Non valence orbitals',nFro1,nSym)
nFre = 0
!----------------------------------------------------------------------*
!     Read the input stream line by line and identify key command      *
!----------------------------------------------------------------------*
Skip = .false.
jCom = 0
outer: do
  if (.not. Skip) then
    Line = Get_Ln(LuSpool)
    call StdFmt(Line,Command)
    jCom = 0
    do iCom=1,size(ComTab)
      if (Command == ComTab(iCom)) jCom = iCom
    end do
    if (jCom == 0) then
      write(u6,*) 'RdInp: Illegal keyword!'
      write(u6,'(A,A)') 'Command=',Command
      call Abend()
    end if
  end if
  Skip = .false.
  !--------------------------------------------------------------------*
  !     Branch to the processing of the command sections               *
  !--------------------------------------------------------------------*
  select case (ComTab(jCom))

    case ('TITL')
      !---- Process the "TITL" input card -----------------------------*
      if (lTit) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'Title option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lTit = .true.
      do
        Line = Get_Ln(LuSpool)
        call StdFmt(Line,Command)
        jCom = 0
        do iCom=1,size(ComTab)
          if (Command == ComTab(iCom)) jCom = iCom
        end do
        if (jCom /= 0) then
          Skip = .true.
          cycle outer
        end if
        nTit = nTit+1
        if (nTit > size(Title)) exit
        Title(nTit) = Line(1:len(Title))
      end do

    case ('FROZ')
      !---- Process the "FROZ" input card -----------------------------*
      if (lFre .or. lFro) then
        write(u6,*) 'RdInp: Error while reading input!'
        if (lFro) write(u6,*) 'Frozen option already processed!'
        if (lFre) write(u6,*) 'Freeze option and Frozen option are incompatible!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lFro = .true.
      Line = Get_Ln(LuSpool)

      write(u6,*)
      write(u6,'(A)') 'WARNING!'
      write(u6,'(A)') 'Default frozen orbitals as non valence orbitals is overwritten by user input.'
      write(u6,'(A,8I4)') 'Default values:',(nFro1(iSym),iSym=1,nSym)
      write(u6,*)

      read(Line,*,iostat=iostatus) (nFro1(iSym),iSym=1,nSym)
      if (iostatus > 0) call error()

    case ('DELE')
      !---- Process the "DELE" input card -----------------------------*
      if (lDel) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'Delete option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lDel = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) (nDel1(iSym),iSym=1,nSym)
      if (iostatus > 0) call error()

    case ('SFRO')
      !---- Process the "SFRO" input card -----------------------------*
      if (lSFro) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'SFrozen option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lSFro = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) (nFro2(iSym),iSym=1,nSym)
      if (iostatus > 0) call error()
      call mma_allocate(iFro,8,maxval(nFro2),label='iFro')
      iFro(:,:) = 0
      do iSym=1,nSym
        Line = Get_Ln(LuSpool)
        read(Line,*,iostat=iostatus) (iFro(iSym,iOrb),iOrb=1,nFro2(iSym))
        if (iostatus > 0) call error()
      end do

    case ('SDEL')
      !---- Process the "SDEL" input card -----------------------------*
      if (lSDel) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'SDelete option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lSDel = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) (nDel2(iSym),iSym=1,nSym)
      if (iostatus > 0) call error()
      call mma_allocate(iDel,8,maxval(nDel2),label='iDel')
      iDel(:,:) = 0
      do iSym=1,nSym
        Line = Get_Ln(LuSpool)
        read(Line,*,iostat=iostatus) (iDel(iSym,iOrb),iOrb=1,nDel2(iSym))
        if (iostatus > 0) call error()
      end do

    case ('EXTR')
      !---- Process the "Extract" input card --------------------------*
      if (lExt) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'Extract option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      write(u6,*) 'RdInp: EXTRACT option is redundant and is ignored!'

    case ('PRIN')
      !---- Process the "Print" input card ----------------------------*
      if (lPrt) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'Print option already processed!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if
      lPrt = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) iPrt
      if (iostatus > 0) call error()

    case ('TEST')
      !---- Process the "Test" input card -----------------------------*
      iTst = 1

    !case ('TSTP')
    !  !---- Process the "TSTP" input card -----------------------------*
    !  iTst = 1
    !  Line = Get_Ln(LuSpool)
    !  read(Line,*,iostat=iostatus) nTstP
    !  if (iostatus > 0) call error()

    case ('PRPT')
      !---- Process the "PRPT" input card -----------------------------*
      call Put_iScalar('mp2prpt',1)
      DoDens = .true.
      NoGamma = .true.

    case ('LUMO')
      !---- Process the "LUMO" input card -----------------------------*
      LumOrb = .true.

    case ('EREF')
      !---- Process the "EREF" input card -----------------------------*
      if (.not. LumOrb) then
        write(u6,*) 'RdInp: Error while reading input!'
        write(u6,*) 'EREF can be used only with LumOrb.'
        write(u6,*) '(Note: LumOrb keyword must precede EREF)'
        call Abend()
      end if
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) ESCF
      if (iostatus > 0) call error()
      ERef_UsrDef = .true.

    case ('VIRA')
      !---- Process the "VIRA" input card -----------------------------*
      all_Vir = .true.

    case ('T1AM')
      !---- Process the "T1AM" input card -----------------------------*
      DoT1amp = .true.
      if (.not. DoCholesky) then
        write(u6,*) 'RdInp: T1AM is available only with Cholesky/RI .'
        call Abend()
      end if

    case ('GRDT')
      !---- Process the "GRDT" input card -----------------------------*
      call Put_iScalar('mp2prpt',2)
      DoDens = .true.
      DoGrdt = .true.

    case ('LAPL')
      !---- Process the "LAPLace" input card --------------------------*
      ! Laplace transform with default or previously specified grid.
      Laplace = .true.

    case ('GRID')
      !---- Process the "GRID" input card -----------------------------*
      ! Read number of Laplace grid points (activates Laplace as well)
      Laplace = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) Laplace_nGridPoints
      if (iostatus > 0) call error()
      Laplace_nGridPoints = max(0,Laplace_nGridPoints)
      if (Laplace_nGridPoints > Laplace_mGridPoints) then
        call WarningMessage(2,'Input Error')
        write(u6,'(A,I6)') 'Number of Laplace grid points specified:',Laplace_nGridPoints
        write(u6,'(A,I6)') 'Max allowed:                            ',Laplace_mGridPoints
        call Quit(_RC_INPUT_ERROR_)
      end if

    case ('BLOC')
      !---- Process the "BLOC" input card -----------------------------*
      ! Read vector block size for CD/DF-Laplace-SOS-MP2
      ! (Activates Laplace as well)
      Laplace = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) Laplace_BlockSize
      if (iostatus > 0) call error()
      Laplace_BlockSize = max(0,Laplace_BlockSize)
      if (Laplace_BlockSize == 0) then
        Laplace_BlockSize = Laplace_BlockSize_Def
      end if

    case ('CHOA')
      !---- Process the "CHOAlgorithm" input card ---------------------*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) ChoAlg
      if (iostatus > 0) call error()
      if (ChoAlg < 0) then
        ChoAlg = 0
      else if (ChoAlg > 2) then
        ChoAlg = 2
      end if

    case ('DECO')
      !---- Process the "DECOmpose MP2 integrals" input card ----------*
      DecoMP2 = .true.
      DecoMP2_UsrDef = .true.

    case ('NODE')
      !---- Process the "NODEcomposition of MP2 integrals" input card -*
      DecoMP2 = .false.
      DecoMP2_UsrDef = .true.

    case ('THRC','PREC')
      !---- Process the "THRCholesky" input card ----------------------*
      !---- Process the "PRECision" input card (= "THRCholesky" card) -*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) ThrMP2
      if (iostatus > 0) call error()
      set_cd_thr = .false.

    case ('SPAN')
      !---- Process the "SPAN" input card -----------------------------*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) SpanMP2
      if (iostatus > 0) call error()

    case ('MXQU')
      !---- Process the "MXQUal" input card ---------------------------*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) MxQualMP2
      if (iostatus > 0) call error()
      if (MxQualMP2 < 1) MxQualMP2 = MxQual_Def

    case ('PRES')
      !---- Process the "PRESort input card (OBSOLETE) ----------------*

    case ('CHKI')
      !---- Process the "CHKI" input card -----------------------------*
      ChkDecoMP2 = .true.

    case ('FORC')
      !---- Process the "FORCebatching" input card --------------------*
      ForceBatch = .true.

    case ('VERB')
      !---- Process the "VERBose" input card --------------------------*
      Verbose = .true.

    case ('NOVE')
      !---- Process the "NOVErbose" input card ------------------------*
      Verbose = .false.

    case ('FREE')
      !---- Process the "FREEze" input card  --------------------------*
      if (lFre .or. lFro) then
        write(u6,*) 'RdInp: Error while reading input!'
        if (lFre) write(u6,*) 'Freeze option already processed!'
        if (lFro) write(u6,*) 'Frozen option and Freeze option are incompatible!'
        write(u6,'(A,A)') 'Last read line:',Line
        call Abend()
      end if

      write(u6,*)
      write(u6,'(A)') 'WARNING!'
      write(u6,'(A)') 'Default frozen orbitals as non valence orbitals is overwritten by user input.'
      write(u6,'(A,8I4)') 'Default values:',(nFro1(iSym),iSym=1,nSym)
      write(u6,*)
      nFro1(:) = 0

      lFre = .true.
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) nFre
      if (iostatus > 0) call error()

    case ('SOSM')
      !---- Process the "SOSMp2" input card ---------------------------*
      SOS_MP2 = .true.
      if (.not. DoLDF) then
        DecoMP2 = .true.
        if (ChoMP2_ChkPar()) then
          call WarningMessage(2,'SOS-MP2 is not implemented for parallel runs. !! SORRY !!')
          call Quit(_RC_NOT_AVAILABLE_)
        end if
      end if

    case ('OEDT')
      !---- Process the "OEDThreshold" input card ---------------------*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) OED_Thr
      if (iostatus > 0) call error()

    case ('OSFA')
      !---- Process the "OSFActor" input card -------------------------*
      Line = Get_Ln(LuSpool)
      read(Line,*,iostat=iostatus) C_os
      if (iostatus > 0) call error()

    case ('LOVM')
      !---- Process the "LovMP2" input card ---------------------------*
      if (.not. DoCholesky) then
        write(u6,*)
        write(u6,*) '********************* ERROR ***********************'
        write(u6,*) ' LovMP2 not implemented with conventional ERIs.'
        write(u6,*) ' Please, use Cholesky or RI options.'
        write(u6,*) '***************************************************'
        call Abend()
      else
        LovMP2 = .true.
      end if
      read(LuSpool,*) nActa,ThrLov
      ! nActa = number of active atoms
      ! ThrLov = threshold for orbital selection

      if ((ThrLov < Zero) .or. (ThrLov >= One)) then
        write(u6,*) ' Threshold out of range! Must be in [0,1[ '
        call Abend()
      end if
      ! namAct = names of active atoms (symm. indep. centers)
      call mma_allocate(NamAct,nActa,label='NamAct')
      do
        read(LuSpool,'(A)',iostat=iostatus) Line
        if (iostatus < 0) call error()
        if ((Line(1:1) /= '*') .and. (Line /= '')) exit
      end do
      call UpCase(Line)
      do i=1,nActa
        if (Line == '') call error()
        Line = adjustl(Line)
        j = index(Line,' ')
        namAct(i) = Line(1:j-1)
        Line(1:j-1) = ''
      end do

    case ('DOMP')
      !---- Process the "DoMP2" input card ----------------------------*
      DoMP2 = .true.

    case ('FNOM')
      !---- process the "FNOM" input card -----------------------------*
      if (.not. DoCholesky) then
        write(u6,*)
        write(u6,*) '********************* ERROR ***********************'
        write(u6,*) ' FNO-MP2 not implemented with conventional ERIs.   '
        write(u6,*) ' Please, use Cholesky or RI options.'
        write(u6,*) '***************************************************'
        call Abend()
      else
        FnoMP2 = .true.
      end if
      read(LuSpool,*) vkept

      if ((vkept <= Zero) .or. (vkept > One)) then
        write(u6,*) ' Requested fraction of virtual space out of range! '
        write(u6,*) ' Must be in ]0,1] '
        call Abend()
      end if

    case ('GHOS')
      !---- process the "GHOSt" input card ----------------------------*
      ! Removal of GHOST virtual space
      do
        read(LuSpool,'(A)',iostat=iostatus) Line
        if (iostatus < 0) call error()
        if ((Line(1:1) /= '*') .and. (Line /= '')) exit
      end do
      read(Line,*,iostat=iostatus) Thr_ghs
      if (iostatus > 0) call error()
      if ((thr_ghs < Zero) .or. (thr_ghs >= One)) then
        write(u6,*) ' GHOST threshold out of range! Must be in [0,1[ '
        call Abend()
      end if
      DelGHost = .true.

    case ('NOGR')
      !---  Process the "NOGR" input card -----------------------------*
      NoGrdt = .true.

    case default
      exit

  end select
end do outer
!----------------------------------------------------------------------*
!     "End of input"                                                   *
!----------------------------------------------------------------------*

if (.not. allocated(NamAct)) call mma_allocate(NamAct,nActa,label='NamAct')
if (.not. allocated(iFro)) call mma_allocate(iFro,8,0,label='iFro')
if (.not. allocated(iDel)) call mma_allocate(iDel,8,0,label='iDel')

! Postprocessing for SOS-MP2 and Laplace
if (SOS_MP2) then
  if (.not. (DoCholesky .or. DoDF .or. DoLDF)) then
    call WarningMessage(2,'SOS-MP2 only implemented for CD/DF/LDF')
    call Quit(_RC_INPUT_ERROR_)
  end if
end if
if (Laplace) then
  if (.not. (DoCholesky .or. DoDF .or. DoLDF)) then
    call WarningMessage(2,'Laplace transformation only implemented for CD/DF/LDF')
    call Quit(_RC_INPUT_ERROR_)
  end if
  if (.not. SOS_MP2) then
    call WarningMessage(2,'Laplace transformation only implemented for CD/DF/LDF-SOS-MP2')
    call Quit(_RC_INPUT_ERROR_)
  end if
  if (.not. DecoMP2_UsrDef) then
    ! SOS-MP2 without Laplace automatically turns on DecoMP2.
    ! With Laplace, however, we should use the same default
    ! as for standard MP2. So, change back to default unless
    ! the user explicitly asked for something else.
    DecoMP2 = Decom_Def
  end if
  if (DoDens .or. DoGrdt) then
    call WarningMessage(2,'Laplace transformation is incompatible with PRPT/GRDT')
    call Quit(_RC_INPUT_ERROR_)
  end if
end if

! Check if the calculation is inside a loop and make analytical
! gradients default in this case

! Numerical gradients requested in GATEWAY
call Qpg_iScalar('DNG',DNG)
if (DNG) then
  call Get_iScalar('DNG',iDNG)
  DNG = iDNG == 1
end if
DNG = NoGrdt .or. DNG

! Inside LAST_ENERGY we do not need analytical gradients
if (SuperName(1:11) == 'last_energy') DNG = .true.

! Inside NUMERICAL_GRADIENT override input!
if (SuperName(1:18) == 'numerical_gradient') then
  call Put_iScalar('mp2prpt',0)
  DNG = .true.
  DoDens = .false.
  DoGrdt = .false.
end if

if (nSym == 1) then
  call GetEnvF('EMIL_InLoop',emiloop)
  if (emiloop == ' ') emiloop = '0'
  call GetEnvF('MOLCAS_IN_GEO',inGeo)
  if ((emiloop(1:1) /= '0') .and. (inGeo(1:1) /= 'Y') .and. (.not. DNG)) then
    call Put_iScalar('mp2prpt',2)
    DoDens = .true.
    DoGrdt = .true.
  end if
end if

if (LumOrb) then
  Lu_orb = 7
  call DaName(Lu_orb,'INPORB')
  l_Occup = nOrb(1)
  do iSym=2,nSym
    l_Occup = l_Occup+nOrb(iSym)
  end do
  call mma_allocate(Occup,l_Occup,label='Occup')
  call RDVEC('INPORB',Lu_orb,'COE',nSym,nBas,nOrb,CMO,Occup,Eall,iDummy,VecTitle,0,iErr)
  if (iErr /= 0) then
    write(u6,'(A,I4)') 'ERROR: RdVec returned code',iErr
    call Abend()
  end if
  call DaClos(Lu_orb)
  write(u6,*)
  write(u6,*) ' Input Orbitals read from INPORB: ',VecTitle
  if (.not. ERef_UsrDef) then
    write(u6,*) ' WARNING: reference energy read from RunFile'
    write(u6,*) '          (may not correspond to orbitals)'
    call Get_dScalar('SCF energy',Escf)
  end if
  write(u6,*)
  iErr = 0
  ip = 0
  do iSym=1,nSym
    iCount = 0
    do i=1,nOrb(iSym)
      ip = ip+1
      if (abs(Occup(ip)) > 1.0e-14_wp) iCount = iCount+1
    end do
    if (iCount /= nOcc(iSym)) then
      iErr = iErr+1
      write(u6,'(A,I2,A,I6,A,I6)') 'WARNING: number of occupied orbitals in symmetry',iSym,' is',iCount, &
                                   ' according to INPORB; from RunFile:',nOcc(iSym)
    end if
  end do
  if (iErr /= 0) then
    write(u6,'(A,A)') 'WARNING: occupation mismatch between RunFile and INPORB. ','RunFile occupation will be used:'
    write(u6,'(8I6)') (nOcc(iSym),iSym=1,nSym)
    iErr = 0
  end if
  call mma_deallocate(Occup)
else
  call Get_dScalar('SCF energy',Escf)
end if

if (lFre .and. (nFre > 0)) then ! freeze lowest occupied orbitals
  FrePrt = iPrt > 0
  call Freezer(Eall,nFre,nFro,nFro1,nOcc,nBas,nSym,FrePrt)
end if
do iSym=1,nSym
  nFro(iSym) = nFro(iSym)+nFro1(iSym)
  nOcc(iSym) = nOcc(iSym)-nFro1(iSym)
  nDel(iSym) = nDel(iSym)+nDel1(iSym)
  nExt(iSym) = nExt(iSym)-nDel1(iSym)
  nOrb(iSym) = nOrb(iSym)-nDel1(iSym)
end do
call Put_iArray('nFroPT',NFRO,NSYM)
call Put_iArray('nDelPT',NDEL,NSYM)
jOcc = 0
iExt = 0
iCount = 0
nOccT = 0
nExtT = 0
do i=1,nSym
  nOccT = nOccT+nOcc(i)
  nExtT = nExtT+nExt(i)
end do

do iSym=1,nSym
  iLow = iCount+nFro(iSym)+1
  iUpp = iLow+nOcc(iSym)-1
  do i=iLow,iUpp
    jOcc = jOcc+1
    Eocc(jOcc) = Eall(i)
  end do
  iLow = iUpp+1
  iUpp = iLow+nExt(iSym)-1
  do i=iLow,iUpp
    iExt = iExt+1
    Eext(iExt) = Eall(i)
  end do
  iCount = iCount+nBas(iSym)
end do
iCount = 0

if (DoDens) then
  jFro = 0
  jDel = 0
  do iSym=1,nSym

    iLow = iCount+1
    iUpp = iLow+nFro(iSym)-1
    do i=iLow,iUpp
      jFro = jFro+1
      EOcc(jFro+nOccT) = Eall(i)
    end do
    iLow = iCount+nFro(iSym)+nOcc(iSym)+nExt(iSym)+1
    iUpp = iLow+nDel(iSym)-1
    do i=iLow,iUpp
      jDel = jDel+1
      Eext(jDel+nExtT) = Eall(i)
    end do
    iCount = iCount+nBas(iSym)
  end do
end if
if (lSFro .or. lSDel) then
  LC = 0
  LEO = 0
  LEE = 0
  LSQ = 0
  do iSym=1,nSym
    LC = LC+nBas(iSym)**2
    LSQ = max(LSQ,nBas(iSym))
    LEO = LEO+nOcc(iSym)
    LEE = LEE+nExt(iSym)
  end do
  call mma_allocate(C,LC,label='C')
  call mma_allocate(EO,LEO,label='EO')
  call mma_allocate(EE,LEE,label='EE')
  call mma_allocate(SQ,LSQ,label='SQ')
  C(:) = CMO(1:LC)
  EO(:) = Eocc(1:LEO)
  EE(:) = Eext(1:LEE)
  call FrzDel(nFro2,iFro,Eocc,EO,nDel2,iDel,Eext,EE,CMO,C,SQ)
  call mma_deallocate(C)
  call mma_deallocate(EO)
  call mma_deallocate(EE)
  call mma_deallocate(SQ)
end if
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*

! remove local copy of standard input

call Close_LuSpool(LuSpool)

if (LumOrb .and. DoT1amp) then
  write(u6,'(/,A)') 'ERROR!  Keywords incompatibility.'
  write(u6,'(/,A)') 'Both LUMOrb and T1AM were selected.'
  write(u6,'(/,A)') '***  I must shut down MBPT2 ! ***'
  call Abend()
end if

if (all_Vir .and. DoMP2) then
  write(u6,'(/,A)') 'WARNING!'
  write(u6,'(/,A)') 'Both VirAll and DoMP2 were selected.'
  write(u6,'(/,A)') '***  I turn off DoMP2 ! ***'
  DoMP2 = .false.
end if

! turn off decomposition for parallel runs.

if (DecoMP2 .and. ChoMP2_ChkPar()) then
  write(u6,'(/,A)') 'WARNING!'
  write(u6,'(A)') 'Decomposition of MP2 integrals is not possible. Turning off decomposition.'
  DecoMP2 = .false.
end if

call xFlush(u6)

return

contains

subroutine error()
  !--------------------------------------------------------------------*
  !     Error Exit                                                     *
  !--------------------------------------------------------------------*
  write(u6,*) 'RdInp: Error while reading input!'
  write(u6,'(A,A)') 'Last read line:',Line
  call Abend()
end subroutine error

end subroutine RdInp
