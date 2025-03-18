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

use Basis_Info, only: Basis_Info_Get
use Center_Info, only: Center_Info_Get
use OneDat, only: sOpSiz
use Exp, only: NewPre, nexp_max
use negpre, only: nGP
use Fock_util_global, only: Deco, dmpk, Estimate, Nscreen, Update
use MCLR_Data, only: ISTATE, OVERRIDE, SA, ESTERR, ISNAC, ISMECIMSPD, FANCY_PRECONDITIONER, NSSA, NACSTATES
use MCLR_Data, only: DspVec, SwLbl, lDisp
use MCLR_Data, only: NoFile
use input_mclr, only: Debug, lRoots, kPrint, mTit, Omega, TimeDep, Page, iBreak, nIter, RASSI, SpinPol, lSave, lCalc, nDisp, &
                      CasInt, NewCho, TwoStep, StepType, double, Eps, IsPop, nSym, nAtoms, ntPert, nsRot, UserP, nUserPT, UserT, &
                      TitleIn
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
#include "rasdim.fh"
character(len=72) Line
character(len=4) Command
character(len=8) Label, SewLab
character(len=2) Element(MxAtom)
logical Epsilon_Undef
integer, parameter :: nCom = 38
character(len=4), parameter :: ComTab(nCom) = ['TITL','DEBU','ROOT','EXTR','PRCI','PROR','ITER','THRE','END ','TIME', &
                                               'CALC','NOFI','SEWA','NOCO','NOTW','SPIN','PRIN','PCGD','RESI','NOTO', &
                                               'EXPD','NEGP','LOWM','ELHE','SAVE','RASS','DISO','CASI','SALA','NODE', &
                                               'ESTE','MOUT','MASS','NAC ','$$$$','THER','CHOF','TWOS']
integer iDum(1), I, JCOM, ICOM, ITIT, ISYM, IP, IRC, IOPT, ICOMP, ISYLBL, IPP, IS, ID, J, IRRFNC, iMass
real*8, allocatable :: umass(:)
character(len=3), allocatable :: cmass(:)

!----------------------------------------------------------------------*
!     Locate "start of input"                                          *
!----------------------------------------------------------------------*
call RdNLst(5,'MCLR')
!----------------------------------------------------------------------*
!     Define default values                                            *
!----------------------------------------------------------------------*
debug = .false.
Epsilon_Undef = .true.
call Basis_Info_Get()
call Center_Info_Get()
call Get_info_Static()
istate = 1     ! State for which the Lagrangian is calc.
override = .false.
if (debug) write(6,*) 'Got Basis_Info and Center_Info'
lRoots = -1
kprint = 0
ngp = .false.
NoFile = .false.
mTit = 0
Omega = 0.0d0
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
do i=1,nDisp
  DspVec(i) = i
end do
CASINT = .true.
NACstates(1) = 0
NACstates(2) = 0
NSSA(1) = 0
NSSA(2) = 0
isNAC = .false.
isMECIMSPD = .false.
NewCho = .false.
! Cholesky. Cannot modify it in the input (yet?)
dmpk = 1.0d-2
Nscreen = 10
Deco = .true.
Update = .true.
Estimate = .false.
TwoStep = .false.
StepType = 'xxxx'
!----------------------------------------------------------------------*
!     Read the input stream line by line and identify key command      *
!----------------------------------------------------------------------*
100 continue
read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 100
call StdFmt(Line,Command)
jCom = 0
do iCom=1,nCom
  if (Command == ComTab(iCom)) jCom = iCom
end do
if (jCom == 0) then
  write(6,'(A,A)') 'RdInp: illegal command:',Command
  call Abend()
end if
!----------------------------------------------------------------------*
!     Branch to the processing of the command sections                 *
!----------------------------------------------------------------------*
110 continue
select case (jCom)
  case (1)
    Go to 10
  case (2)
    Go to 16
  case (3)
    Go to 20
  case (4)
    Go to 30
  case (5)
    Go to 40
  case (6)
    Go to 50
  case (7)
    Go to 60
  case (8)
    Go to 70
  case (9)
    Go to 99
  case (10)
    Go to 80
  case (11)
    Go to 55
  case (12)
    Go to 175
  case (13)
    Go to 185
  case (14)
    Go to 165
  case (15)
    Go to 166
  case (16)
    Go to 177
  case (17)
    Go to 178
  case (18)
    Go to 179
  case (19)
    Go to 180
  case (20)
    Go to 191
  case (21)
    Go to 192
  case (22)
    Go to 193
  case (23)
    Go to 194
  case (24)
    Go to 195
  case (25)
    Go to 196
  case (26)
    Go to 788
  case (27)
    Go to 789
  case (28)
    Go to 198
  case (29)
    Go to 199
  case (30)
    Go to 200
  case (31)
    Go to 201
  case (32)
    Go to 202
  case (33)
    Go to 203
  case (34)
    Go to 204
  case (35)
    Go to 205
  case (36)
    Go to 206
  case (37)
    Go to 210
  case (38)
    Go to 220
end select
!---  TITL ------------------------------------------------------------*
10 continue
15 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if (Line(1:1) == '*') goto 15
call StdFmt(Line,Command)
jCom = 0
do iCom=1,nCom
  if (Command == ComTab(iCom)) jCom = iCom
end do
if (jCom /= 0) goto 110
mTit = mTit+1
if (mTit <= mxTit) then
  read(Line,'(18A4)') (TitleIN(iTit),iTit=(mTit-1)*18+1,mTit*18)
  goto 15
end if
goto 100

!----      ------------------------------------------------------------*
788 RASSI = .true.
if (debug) write(6,*) 'Output for RASSI'
goto 100
!----      ------------------------------------------------------------*
789 double = .true.    ! Make double isotope substitutions
goto 100
!---- DEBU ------------------------------------------------------------*
16 debug = .true.
goto 100
!----      ------------------------------------------------------------*
195 continue
write(6,*) 'ELHE is disabled!'
goto 100
!---- LOWM ------------------------------------------------------------*
194 page = .true.
if (debug) write(6,*) 'Page memory'
goto 100
!----      ------------------------------------------------------------*
191 newpre = .false.
if (debug) write(6,*) 'New conditioner'
goto 100
!----      ------------------------------------------------------------*
196 lSAVE = .true.
if (debug) write(6,*) 'old integrals, not supported'
goto 100
!----      ------------------------------------------------------------*
198 CASINT = .true.
if (debug) write(6,*) 'CASPT2 integrals'
goto 100
!---- EXPD ------------------------------------------------------------*
192 read(5,*) nexp_max
if (debug) write(6,*) 'Maximum explicit preconditioner',nexp_max
goto 100
!----      ------------------------------------------------------------*
179 iBreak = 1
read(5,*) Eps
Epsilon_Undef = .false.
if (debug) write(6,*) 'Threshold:',Eps
goto 100
!----      ------------------------------------------------------------*
180 iBreak = 2
read(5,*) Eps
Epsilon_Undef = .false.
if (debug) write(6,*) 'Threshold:',Eps
goto 100
!----      ------------------------------------------------------------*
178 read(5,*) kprint
if (debug) write(6,*) 'Print level: ',kprint
goto 100
!----      ------------------------------------------------------------*
193 NGP = .true.
if (debug) write(6,*) 'NGP set to true'
goto 100
!---- SALA ------------------------------------------------------------*
199 SA = .true.
read(5,*) istate
override = .true.
if (debug) write(6,*) 'Lagrangian for state: ',istate
goto 100
!----      ------------------------------------------------------------*
201 esterr = .true.
goto 100
!---- NODE ------------------------------------------------------------*
200 FANCY_PRECONDITIONER = .false.
if (debug) write(6,*) 'Turned of the fancy pcg'
goto 100
!----      ------------------------------------------------------------*
177 SPINPOL = .true.
ispop = 1
if (debug) write(6,*) 'RHF lagrangian, not supported'
goto 100
!----      ------------------------------------------------------------*
166 do i=1,nDisp
  NTPert(i) = iand(nTPert(i),247)
end do
goto 100
!----      ------------------------------------------------------------*
165 do i=1,nDisp
  NTPert(i) = iand(nTPert(i),251)
end do
goto 100
!----      ------------------------------------------------------------*

185 continue
186 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
call StdFmt(Line,Command)
if ((Command(1:1) == ' ') .or. (Line(1:1) == '*')) goto 186
if (debug) write(6,*) 'SEWARD INPUT'
if ((Command(1:4) == 'END ') .or. (Command(1:4) == 'ENDS')) goto 100
read(Line,'(A8,I2,I2)',Err=998,end=999) SewLab,isym,ip
iRc = -1
iOpt = ibset(0,sOpSiz)
iComp = ip
iSyLbl = 2**isym
Label = SewLab
call iRdOne(iRc,iOpt,Label,iComp,idum,iSyLbl)
if (iRc /= 0) then
  write(6,*) 'RdInp: Error reading ONEINT'
  write(6,'(A,A)') 'Label=',Label
  call Abend()
end if

!---  read number of symm. species ------------------------------------*

ipp = 0
do is=isym+1,nsym
  ipp = ipp+ldisp(is)
end do
do id=nDisp,ndisp-ipp+1,-1
  DspVec(id+1) = dspVec(id)
  ntpert(id+1) = ntpert(id)
  lcalc(id+1) = lcalc(id)
end do
id = ndisp-ipp+1
DspVec(id) = ip
ldisp(isym) = ldisp(isym)+1
ndisp = ndisp+1
ntpert(id) = 2
lcalc(id) = .true.
SwLbl(id) = SewLab
goto 185

175 Nofile = .true.
if (debug) write(6,*) 'NOFILE'
goto 100

!---  Process the "root" input card -----------------------------------*
20 continue
25 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 25
read(Line,*,Err=998,end=999) lRoots
if (debug) write(6,*) 'LROOT'
goto 100

!---  CALC ------------------------------------------------------------*
55 continue
write(6,*) 'CALC is disabled!'
goto 100
!---  Process the "extract" input card --------------------------------*
30 write(6,*) 'RdInp: EXTRACT option is redundant and is ignored!'
goto 100
!---  Process the "PrCI" input card -----------------------------------*
40 continue
write(6,*) 'PRCI is disabled!'
goto 100
!---  Process the "PrOr" input card -----------------------------------*
50 continue
write(6,*) 'PROR is disabled!'
goto 100
!---  Process the "ITER" input card -----------------------------------*
60 continue
65 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 65
read(Line,*,Err=998,end=999) nIter
goto 100
!---  Process the "THRE" input card -----------------------------------*
70 continue
75 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 75
read(Line,*,Err=998,end=999) Eps
Epsilon_Undef = .false.
goto 100
!---  Process the "TIME" input card -----------------------------------*
80 continue
read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 80
read(Line,*,Err=998,end=999) Omega
TimeDep = .true.
nIter = 100
goto 100
!---  Process the "MOUT" input card -----------------------------------*
202 continue
write(6,*) 'MOUT is disabled!'
goto 100
!---  Process the "MASS" input card -----------------------------------*
203 continue
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
  read(5,'(A3)') cmass(i)
  read(5,'(F15.8)') umass(i)
end do

! Put the Info on the run file.

call Put_iScalar('iMass',iMass)
call Put_cArray('cmass',cmass(1),3*iMass)
call Put_dArray('umass',umass,iMass)
call mma_deallocate(cmass)
call mma_deallocate(umass)

goto 100
!---  Process the "NAC " input card -----------------------------------*
204 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if ((Line(1:1) == ' ') .or. (Line(1:1) == '*')) goto 25
read(Line,*,Err=998,end=999) NACstates(1),NACstates(2)
isNAC = .true.
override = .true.
if (debug) write(6,*) 'Non-adiabatic couplings for states: ',NACstates(1),NACstates(2)
goto 100
!---  Process the "$$$$" input card -----------------------------------*
205 continue
! not used
goto 100
!---  Process the "THERmochemistry input card -------------------------*
206 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if (Line(1:1) == '*') goto 206
read(Line,*,Err=998,end=999) nsRot
2060 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if (Line(1:1) == '*') goto 2060
read(Line,*,Err=998,end=999) UserP
2061 read(5,'(A)',Err=998,end=999) Line
Line = adjustl(Line)
if (Line(1:1) == '*') goto 2061
call UpCase(Line)
if (Line(1:4) == 'END ') then
  if (nUserPT == 0) then
    nUserPT = 1
    UserT(1) = 298.15d0
  end if
  goto 100
end if
nUserPT = nUserPT+1
read(Line,*,Err=998,end=999) UserT(nUserPT)
goto 2061
!---  Process the "NEWCho input card ----------------------------------*
210 NewCho = .true.
goto 100
!---  Process the "TWOStep" input card --------------------------------*
220 read(5,'(A)',Err=998,end=999) Line
call UpCase(Line)
Line = adjustl(Line)
if (Line(1:1) == '*') goto 220
read(Line,*,Err=998,end=999) StepType
if (debug) write(6,*) 'TWOSTEP kind: '//StepType
if ((StepType(1:4) /= 'FIRS') .and. (StepType(1:4) /= 'SECO') .and. (StepType(1:4) /= 'RUN1') .and. (StepType(1:4) /= 'RUN2')) then
  call WarningMessage(2,'TWOStep: input error!')
  call Quit_OnUserError()
end if
if (StepType(1:4) == 'FIRS') StepType(1:4) = 'RUN1'
if (StepType(1:4) == 'SECO') StepType(1:4) = 'RUN2'
TwoStep = .true.
if (debug) write(6,*) 'TWOSTEP kind: '//StepType
goto 100
!----------------------------------------------------------------------*
!     "End of input"                                                   *
!----------------------------------------------------------------------*
99 continue
do i=1,3
  isym = irrfnc(2**(i-1))+1
  ipp = 0
  do is=isym+1,nsym
    ipp = ipp+ldisp(is)
  end do
  do id=nDisp,ndisp-ipp+1,-1
    DspVec(id+1) = dspVec(id)
    ntpert(id+1) = ntpert(id)
    !lcalc(id+1) = lcalc(id)
    Swlbl(id+1) = Swlbl(id)
  end do
  id = ndisp-ipp+1
  DspVec(id) = i
  ldisp(isym) = ldisp(isym)+1
  ndisp = ndisp+1
  ntpert(id) = 2
  !lcalc(id) = .true.
  write(Swlbl(id),'(a,i2)') 'MLTPL ',1
  iRc = -1
  iOpt = ibset(0,sOpSiz)
  call iRdOne(iRc,iOpt,swlbl(id),dspvec(id),idum,iSyLbl)
end do

if (Timedep) then
  do i=1,ndisp
    ntpert(i) = ior(ntpert(i),32)
  end do
end if

if (Epsilon_Undef) then
  !if (SA) then
  !  Eps = 1.0D-6
  !else
  Eps = 1.0D-4
  ! This I need to change back
  !end if
end if

if (debug) write(6,*) 'FINITO'
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*

return
!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*
998 write(6,*) 'RdInp: Error while reading input'
write(6,'(A,A)') 'Last command:',Line
call Abend()
999 write(6,*) 'RdInp: Premature end of input file'
write(6,'(A,A)') 'Last command:',Line
call Abend()

end subroutine RdInp_MCLR
