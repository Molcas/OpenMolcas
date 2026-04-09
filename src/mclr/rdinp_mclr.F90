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

subroutine RdInp_MCLR()
!***********************************************************************
!                                                                      *
!     Locate input stream and read commands                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: sOpSiz
use Fock_util_global, only: Deco, dmpk, Estimate, Nscreen, Update
use MCLR_Data, only: DspVec, ESTERR, FANCY_PRECONDITIONER, ISMECIMSPD, ISNAC, ISTATE, lDisp, NACSTATES, NewPre, nexp_max, nGP, &
                     NoFile, NSSA, OVERRIDE, SA, SwLbl
use input_mclr, only: CasInt, Debug, double, Eps, iBreak, IsPop, kPrint, lCalc, lRoots, lSave, mTit, nAtoms, nDisp, NewCho, nIter, &
                      nsRot, nSym, ntPert, nUserPT, Omega, Page, RASSI, SpinPol, StepType, TimeDep, TitleIn, TwoStep, UserP, UserT
use PCM_grad, only: RFPERT
use cgs_mod, only: CGS
use Molcas, only: MxAtom
use RASDim, only: MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp) :: I, ICOM, ICOMP, ID, iDum(1), iMass, IOPT, IP, IPP, IRC, IRRFNC, istatus, ISYLBL, ISYM, ITIT, J, JCOM, nDiff
logical(kind=iwp) :: DoRys, Epsilon_Undef, Skip
character(len=72) :: Line
character(len=8) :: Label, SewLab
character(len=4) :: Command
character(len=2) :: Element(MxAtom)
real(kind=wp), allocatable :: umass(:)
character(len=3), allocatable :: cmass(:)
integer(kind=iwp), parameter :: nCom = 38
character(len=*), parameter :: ComTab(nCom) = ['TITL','DEBU','ROOT','    ','CGS ','RFPE','ITER','THRE','END ','TIME', &
                                               '    ','NOFI','SEWA','NOCO','NOTW','SPIN','PRIN','PCGD','RESI','NOTO', &
                                               'EXPD','NEGP','LOWM','    ','SAVE','RASS','DISO','CASI','SALA','NODE', &
                                               'ESTE','    ','MASS','NAC ','    ','THER','CHOF','TWOS']

!----------------------------------------------------------------------*
!     Locate "start of input"                                          *
!----------------------------------------------------------------------*
call RdNLst(u5,'MCLR')
!----------------------------------------------------------------------*
!     Define default values                                            *
!----------------------------------------------------------------------*
debug = .false.
Epsilon_Undef = .true.
! Calling Basis_Info_Get(), Center_Info_Get(), and Get_info_Static() are replaced with IniSew below
nDiff = 0
DoRys = .true.
call IniSew(DoRys,nDiff)
istate = 1     ! State for which the Lagrangian is calc.
override = .false.
if (debug) write(u6,*) 'Got Basis_Info and Center_Info'
lRoots = -1
kprint = 0
ngp = .false.
NoFile = .false.
mTit = 0
Omega = Zero
TimeDep = .false.
Page = .false.
ibreak = 2
nIter = 200
RASSI = .false.
spinpol = .false.
SA = .false.
esterr = .false.
FANCY_PRECONDITIONER = .true.
lsave = .false.
lCalc(:) = .true.
DspVec(1:nDisp) = [(i,i=1,nDisp)]
CASINT = .true.
NACstates(1) = 0
NACstates(2) = 0
NSSA(1) = 0
NSSA(2) = 0
isNAC = .false.
isMECIMSPD = .false.
NewCho = .false.
! Cholesky. Cannot modify it in the input (yet?)
dmpk = 1.0e-2_wp
Nscreen = 10
Deco = .true.
Update = .true.
Estimate = .false.
TwoStep = .false.
StepType = 'xxxx'
CGS = .false.
RFPERT = .false.
!----------------------------------------------------------------------*
!     Read the input stream line by line and identify key command      *
!----------------------------------------------------------------------*
Skip = .false.
jCom = 0 ! dummy init
outer: do
  if (Skip) then
    Skip = .false.
  else
    read(u5,'(A)',iostat=istatus) Line
    if (istatus /= 0) call Error(istatus)
    Line = adjustl(Line)
    if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) cycle
    call StdFmt(Line,Command)
    jCom = 0
    do iCom=1,nCom
      if (Command == ComTab(iCom)) jCom = iCom
    end do
  end if
  if (jCom == 0) then
    write(u6,'(A,A)') 'RdInp: illegal command:',Command
    call Abend()
  end if
  !--------------------------------------------------------------------*
  !     Branch to the processing of the command sections               *
  !--------------------------------------------------------------------*
  select case (jCom)
    case (1)
      !---- TITL ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if (Line(1:1) == '*') cycle
        call StdFmt(Line,Command)
        jCom = 0
        do iCom=1,nCom
          if (Command == ComTab(iCom)) jCom = iCom
        end do
        if (jCom /= 0) then
          Skip = .true.
          exit
        end if
        mTit = mTit+1
        if (mTit > mxTit) exit
        read(Line,'(18A4)') (TitleIN(iTit),iTit=(mTit-1)*18+1,mTit*18)
      end do

    case (2)
      !---- DEBU ------------------------------------------------------*
      debug = .true.

    case (3)
      !---- ROOT ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if ((Line(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) lRoots
      if (istatus /= 0) call Error(istatus)
      if (debug) write(u6,*) 'LROOT'

    case (5)
      !---- CGS  ------------------------------------------------------*
      CGS = .true.

    case (6)
      !---- RFPErt ----------------------------------------------------*
      RFPERT = .true.

    case (7)
      !---- ITER ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if ((Line(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) nIter
      if (istatus /= 0) call Error(istatus)

    case (8)
      !---- THRE ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if ((Line(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) Eps
      if (istatus /= 0) call Error(istatus)
      Epsilon_Undef = .false.

    case (9)
      !---- END  ------------------------------------------------------*
      exit outer

    case (10)
      !---- TIME ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if ((Line(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) Omega
      if (istatus /= 0) call Error(istatus)
      TimeDep = .true.
      nIter = 100

    case (12)
      !---- NOFI ------------------------------------------------------*
      Nofile = .true.
      if (debug) write(u6,*) 'NOFILE'

    case (13)
      !---- SEWA ------------------------------------------------------*
      do
        do
          read(u5,'(A)',iostat=istatus) Line
          if (istatus /= 0) call Error(istatus)
          Line = adjustl(Line)
          call StdFmt(Line,Command)
          if ((Command(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
        end do
        if (debug) write(u6,*) 'SEWARD INPUT'
        if ((Command(1:4) == 'END ') .or. (Command(1:4) == 'ENDS')) exit
        read(Line,'(A8,I2,I2)',iostat=istatus) SewLab,isym,ip
        if (istatus /= 0) call Error(istatus)
        iRc = -1
        iOpt = ibset(0,sOpSiz)
        iComp = ip
        iSyLbl = ibset(0,isym)
        Label = SewLab
        call iRdOne(iRc,iOpt,Label,iComp,idum,iSyLbl)
        if (iRc /= 0) then
          write(u6,*) 'RdInp: Error reading ONEINT'
          write(u6,'(A,A)') 'Label=',Label
          call Abend()
        end if

        !---- read number of symm. species ----------------------------*

        ipp = sum(lDisp(isym+1:nsym))
        DspVec(nDisp-ipp+2:nDisp+1) = dspVec(nDisp-ipp+1:nDisp)
        ntpert(nDisp-ipp+2:nDisp+1) = ntpert(nDisp-ipp+1:nDisp)
        lcalc(nDisp-ipp+2:nDisp+1) = lcalc(nDisp-ipp+1:nDisp)
        id = ndisp-ipp+1
        DspVec(id) = ip
        ldisp(isym) = ldisp(isym)+1
        ndisp = ndisp+1
        ntpert(id) = 2
        lcalc(id) = .true.
        SwLbl(id) = SewLab
      end do

    case (14)
      !---- NOCO ------------------------------------------------------*
      do i=1,nDisp
        NTPert(i) = ibclr(nTPert(i),2)
      end do

    case (15)
      !---- NOTW ------------------------------------------------------*
      do i=1,nDisp
        NTPert(i) = ibclr(nTPert(i),3)
      end do

    case (16)
      !---- SPIN ------------------------------------------------------*
      SPINPOL = .true.
      ispop = 1
      if (debug) write(u6,*) 'RHF lagrangian, not supported'

    case (17)
      !---- PRIN ------------------------------------------------------*
      read(u5,*) kprint
      if (debug) write(u6,*) 'Print level: ',kprint

    case (18,19)
      !---- PCGD, RESI ------------------------------------------------*
      if (jCom == 18) then
        iBreak = 1
      else if (jCom == 19) then
        iBreak = 2
      end if
      read(u5,*) Eps
      Epsilon_Undef = .false.
      if (debug) write(u6,*) 'Threshold:',Eps

    case (20)
      !---- NOTO ------------------------------------------------------*
      newpre = .false.
      if (debug) write(u6,*) 'New conditioner'

    case (21)
      !---- EXPD ------------------------------------------------------*
      read(u5,*) nexp_max
      if (debug) write(u6,*) 'Maximum explicit preconditioner',nexp_max

    case (22)
      !---- NEGP ------------------------------------------------------*
      NGP = .true.
      if (debug) write(u6,*) 'NGP set to true'

    case (23)
      !---- LOWM ------------------------------------------------------*
      page = .true.
      if (debug) write(u6,*) 'Page memory'

    case (25)
      !---- SAVE ------------------------------------------------------*
      lSAVE = .true.
      if (debug) write(u6,*) 'old integrals, not supported'

    case (26)
      !---- RASS ------------------------------------------------------*
      RASSI = .true.
      if (debug) write(u6,*) 'Output for RASSI'

    case (27)
      !---- DISO ------------------------------------------------------*
      double = .true.    ! Make double isotope substitutions

    case (28)
      !---- CASI ------------------------------------------------------*
      CASINT = .true.
      if (debug) write(u6,*) 'CASPT2 integrals'

    case (29)
      !---- SALA ------------------------------------------------------*
      SA = .true.
      read(u5,*) istate
      override = .true.
      if (debug) write(u6,*) 'Lagrangian for state: ',istate

    case (30)
      !---- NODE ------------------------------------------------------*
      FANCY_PRECONDITIONER = .false.
      if (debug) write(u6,*) 'Turned of the fancy pcg'

    case (31)
      !---- ESTE ------------------------------------------------------*
      esterr = .true.

    case (33)
      !---- MASS ------------------------------------------------------*
      iMass = 0
      call Get_Name_All(Element)

      ! Find out how many different elements are present in the molecule.

      do i=1,nAtoms
        if (Element(i) /= '  ') iMass = iMass+1
        do j=i+1,nAtoms
          if (Element(j) == Element(i)) Element(j) = '  '
        end do
      end do
      call mma_allocate(cmass,iMass,label='cmass')
      call mma_allocate(umass,iMass,label='umass')
      do i=1,iMass
        read(u5,'(A3)') cmass(i)
        read(u5,'(F15.8)') umass(i)
      end do

      ! Put the Info on the run file.

      call Put_iScalar('iMass',iMass)
      call Put_cArray('cmass',cmass(1),3*iMass)
      call Put_dArray('umass',umass,iMass)
      call mma_deallocate(cmass)
      call mma_deallocate(umass)

    case (34)
      !---- NAC  ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if ((Line(1:1) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) NACstates(1),NACstates(2)
      if (istatus /= 0) call Error(istatus)
      isNAC = .true.
      override = .true.
      if (debug) write(u6,*) 'Non-adiabatic couplings for states: ',NACstates(1),NACstates(2)

    case (36)
      !---- THER ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) nsRot
      if (istatus /= 0) call Error(istatus)
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) UserP
      if (istatus /= 0) call Error(istatus)
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        Line = adjustl(Line)
        if (Line(1:1) == '*') cycle
        call UpCase(Line)
        if (Line(1:4) == 'END ') then
          if (nUserPT == 0) then
            nUserPT = 1
            UserT(1) = 298.15_wp
          end if
          exit
        end if
        nUserPT = nUserPT+1
        read(Line,*,iostat=istatus) UserT(nUserPT)
        if (istatus /= 0) call Error(istatus)
      end do

    case (37)
      !---- CHOF ------------------------------------------------------*
      NewCho = .true.

    case (38)
      !---- TWOS ------------------------------------------------------*
      do
        read(u5,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error(istatus)
        call UpCase(Line)
        Line = adjustl(Line)
        if (Line(1:1) /= '*') exit
      end do
      read(Line,*,iostat=istatus) StepType
      if (istatus /= 0) call Error(istatus)
      if (debug) write(u6,*) 'TWOSTEP kind: '//StepType
      if ((StepType(1:4) /= 'FIRS') .and. (StepType(1:4) /= 'SECO') .and. (StepType(1:4) /= 'RUN1') .and. &
          (StepType(1:4) /= 'RUN2')) then
        call WarningMessage(2,'TWOStep: input error!')
        call Quit_OnUserError()
      end if
      if (StepType(1:4) == 'FIRS') StepType(1:4) = 'RUN1'
      if (StepType(1:4) == 'SECO') StepType(1:4) = 'RUN2'
      TwoStep = .true.
      if (debug) write(u6,*) 'TWOSTEP kind: '//StepType

    case default
      jCom = 0
      Skip = .true.

  end select
end do outer

!----------------------------------------------------------------------*
!     "End of input"                                                   *
!----------------------------------------------------------------------*
do i=1,3
  isym = irrfnc(ibset(0,i-1))+1
  ipp = sum(lDisp(isym+1:nsym))
  DspVec(nDisp-ipp+2:nDisp+1) = dspVec(nDisp-ipp+1:nDisp)
  ntpert(nDisp-ipp+2:nDisp+1) = ntpert(nDisp-ipp+1:nDisp)
  !lcalc(nDisp-ipp+2:nDisp+1) = lcalc(nDisp-ipp+1:nDisp)
  Swlbl(nDisp-ipp+2:nDisp+1) = Swlbl(nDisp-ipp+1:nDisp)
  id = ndisp-ipp+1
  DspVec(id) = i
  ldisp(isym) = ldisp(isym)+1
  ndisp = ndisp+1
  ntpert(id) = 2
  !lcalc(id) = .true.
  write(Swlbl(id),'(a,i2)') 'MLTPL ',1
  iRc = -1
  iOpt = ibset(0,sOpSiz)
  call iRdOne(iRc,iOpt,swlbl(id),DspVec(id),idum,iSyLbl)
end do

if (Timedep) then
  do i=1,ndisp
    ntpert(i) = ibset(ntpert(i),5)
  end do
end if

if (Epsilon_Undef) then
  !if (SA) then
  !  Eps = 1.0e-6_wp
  !else
  Eps = 1.0e-4_wp
  ! This I need to change back
  !end if
end if

if (debug) write(u6,*) 'FINITO'
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*

contains

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*
subroutine Error(rc)

  integer(kind=iwp), intent(in) :: rc

  if (rc > 0) then
    write(u6,*) 'RdInp: Error while reading input'
  else
    write(u6,*) 'RdInp: Premature end of input file'
  end if
  write(u6,'(A,A)') 'Last command:',Line
  call Abend()

end subroutine Error

end subroutine RdInp_MCLR
