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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine RdInp_scf()
!***********************************************************************
!                                                                      *
!     purpose: Read input options                                      *
!                                                                      *
!     called from: ReadIn                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!     modified by M.Schuetz @teokem.lu.se, 1995                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: UHF- V.Veryazov, 2003                                   *
!                                                                      *
!***********************************************************************

use KSDFT_Info, only: CoefR, CoefX
use OFembed, only: dfmd, Do_OFemb, KEonly, OFE_KSDFT, ThrFThaw, XSigma
use Functionals, only: Custom_File, Custom_Func
use IOBuf, only: lDaRec, nSect
use Fock_util_global, only: Deco, DensityCheck, Estimate, Update
use SpinAV, only: Do_SpinAV
use InfSCF, only: Addc_KSDFT, AddFragments, ALGO, Aufb, C1DIIS, Cho_Aufb, Damping, dmpk, DDnOff, DelThr, DIIS, DIISTh, DltnTh, &
                  Do_addc, Do_Tw, DoCholesky, DoHLgap, DSCF, DThr, EThr, ExFac, Falcon, FckAuf, FckAuf, FlipThr, FThr, HLgap, &
                  iAu_ab, iCoCo, iDKeep, indxc, InVec, iPrForm, iPrint, iPrOrb, isHDF5, iStatPRN, Iter2run, IterPrlv, jPrint, &
                  jVOut, kIVO, klockan, kOptim_Max, KSDFT, LKon, LstVec, MaxFlip, MiniDn, MSYMON, MxConstr, MxIter, MxOptm, nAufb, &
                  nBas, nConstr, nCore, nD, nDel, nDisc, Neg2_Action, nFro, nIter, nOcc, NoExchange, NoProp, nOrb, nScreen, nSym, &
                  nTit, OccSet_e, OccSet_m, One_Grid, OnlyProp, PmTime, PreSch, QNRTh, QudThr, ReOrd, RFPert, RGEK, RotFac, &
                  RotLev, RotMax, RSRFO, RTemp, SCF_FileOrb, ScrFac, Scrmbl, Teee, TemFac, Thize, ThrEne, Title, Tot_Charge, &
                  Tot_El_Charge, Tot_Nuc_Charge, TStop, WrOutD
use Cholesky, only: ChFracMem, timings
#ifdef _HDF5_
use mh5, only: mh5_is_hdf5, mh5_open_file_r
use InfSCF, only: FileOrb_ID
#endif
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Two, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "hfc_logical.fh"
integer(kind=iwp) :: i, iArray(32), iAuf, iD, iFroz, iOccu, iOrbi, iPri, iStatus, iSym, j, KeywNo, lthSet_a, lthSet_b, LuCF, &
                     LuSpool, nOccSet_e, nOccSet_m, Mode(1), nFunc, nnn, nSqrSum
real(kind=wp) :: Tot_Ml_Charge
logical(kind=iwp) :: CharSet, Chol, FermSet, IfAufChg, lTtl, OccSet, SpinSet, TDen_UsrDef, UHFSet
character(len=180) :: Key, Line
character(len=8) :: Method
integer(kind=iwp), external :: Allocdisk, IsFreeUnit
real(kind=wp), external :: Get_ExFac
character(len=180), external :: Get_Ln

! copy input from standard input to a local scratch file

call SpoolInp(LuSpool)

OccSet = .false.
FermSet = .false.
CharSet = .false.
SpinSet = .false.

Neg2_Action = 'STOP'

! Some initialization

Chol = .false.
ALGO = 4
REORD = .false.
DECO = .true.
DensityCheck = .false.
timings = .false.
UHFSet = .false.
Nscreen = 10    ! default screening interval (# of red sets)
dmpk = 0.045_wp ! default damping of the screening threshold
Estimate = .false.
Update = .true.
#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Half
#endif
SCF_FileOrb = 'INPORB'
isHDF5 = .false.
! Constrained SCF initialization
do i=1,nSym
  nConstr(i) = 0
end do
MxConstr = 0
klockan = 1
Do_Addc = .false.
iTer2run = 2
! Delta_Tw correlation energy calculation
Do_Tw = .false.
! Read Cholesky info from runfile and save in infscf.fh
call DecideOnCholesky(DoCholesky)
if (DoCholesky) then
  Chol = .true.
  call Cho_scf_rdInp(.true.,LuSpool)
  if (ALGO >= 2) then
    DDnOFF = .true. !do not use diffential density
    MiniDn = .false.
  end if
  !tbp, may 2013: no thr modification with Cholesky
  !tbp call Get_dScalar('Cholesky Threshold',ThrCom)
  !tbp EThr = Max(EThr,ThrCom)
else
  DDnOFF = .false. !default for conventional
  MiniDn = .true.
end if
TDen_UsrDef = .false.

! Set up number of orbitals
do iSym=1,nSym
  nOrb(iSym) = nBas(iSym)
end do
! Set up some counters
nSqrSum = 0
do iSym=1,nSym
  nSqrSum = nSqrSum+nBas(iSym)*nBas(iSym)
end do

iOrbi = 0
iFroz = 0
iOccu = 0
nTit = 0
iPrForm = -1
iterprlv = 0
ScrFac = Zero

! Parameters that control how new orbitals are generated by neworb.

RotLev = Zero
RotFac = One
RotMax = Ten
!HLgap = -One
HLgap = 0.2_wp
DoHLgap = .false.
MaxFlip = 10
FlipThr = 0.1_wp

! Skip exchange when building Fock matrix
! (for debugging purposes)

NoExchange = .false.

! Default value for starting orbitals
! Invec=-1 indicate decision in sorb

InVec = -1

! Default to aufbau for neutral species

Aufb = .true.
Teee = .true.
RTemp = Half
TemFac = 0.46_wp
TStop = 0.01_wp
nAufb(1) = -1
nAufb(2) = -1
call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
Tot_El_Charge = Zero
Tot_Ml_Charge = Zero
Tot_Charge = Zero
iAu_ab = 0
iStatPRN = 0

IfAufChg = .false.

Falcon = .false.
MSYMON = .false.

nD = 1

! Locate "start of input"
rewind(LuSpool)
call RdNLst(LuSpool,'SCF')

KeywNo = 0
1000 lTtl = .false.
999 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)
KeywNo = KeywNo+1

if (Line(1:4) == 'TITL') Go To 1100
if (Line(1:4) == 'ITER') Go To 1200
if (Line(1:4) == 'OCCU') Go To 1300
if (Line(1:4) == 'ORBI') Go To 1400
if (Line(1:4) == 'FROZ') Go To 1500
if (Line(1:4) == 'OVLD') Go To 1700
if (Line(1:4) == 'PRLS') Go To 1800
if (Line(1:4) == 'PROR') Go To 1900
if (Line(1:4) == 'KEEP') Go To 2000
if (Line(1:4) == 'STAR') Go To 2100

! Section for Cholesky input
if (Line(1:4) == 'CHOL') Go To 20000
if (Line(1:4) == 'CHOI') Go To 20001

if (Line(1:4) == 'CONS') Go To 20002
if (Line(1:4) == 'CORE') Go To 21000
if (Line(1:4) == 'NDDO') Go To 21001
if (Line(1:4) == 'LUMO') Go To 21002
if (Line(1:4) == 'GSSR') Go To 21003
if (Line(1:4) == 'REST') Go To 21004
if (Line(1:4) == 'THRE') Go To 2200
if (Line(1:4) == 'NODI') Go To 2300
if (Line(1:4) == 'DIIS') Go To 2400
if (Line(1:4) == 'OCCN') Go To 2500
if (Line(1:4) == 'MCCN') Go To 2510
if (Line(1:4) == 'IVO ') Go To 2600
if (Line(1:4) == 'UHF ') Go To 2700
if (Line(1:4) == 'HFC ') Go To 2701
if (Line(1:4) == 'NODA') Go To 2900
if (Line(1:4) == 'CONV') Go To 3000
if (Line(1:4) == 'DISK') Go To 3100
if (Line(1:4) == 'THIZ') Go To 3200
if (Line(1:4) == 'SIMP') Go To 3300
if (Line(1:4) == 'NOMI') Go To 3400
if (Line(1:4) == 'TDEN') Go To 3401
if (Line(1:4) == 'WRDE') Go To 3500
if (Line(1:4) == 'C1DI') Go To 3600
if (Line(1:4) == 'QUAD') Go To 3700
if (Line(1:4) == 'RS-R') Go To 3750
if (Line(1:4) == 'S-GE') Go To 3751
if (Line(1:4) == 'SCRA') Go To 3800
if (Line(1:4) == 'EXTR') Go To 3900
if (Line(1:4) == 'RFPE') Go To 4000
if (Line(1:4) == 'QNRT') Go To 4100
if (Line(1:4) == 'AUFB') Go To 4200
if (Line(1:4) == 'FERM') Go To 4250
if (Line(1:4) == 'TEEE') Go To 4300
if (Line(1:4) == 'CHAR') Go To 4400
if (Line(1:4) == 'NOTE') Go To 4500
if (Line(1:4) == 'KSDF') Go To 4600
if (Line(1:4) == 'DFCF') Go To 4605
if (Line(1:4) == 'OFEM') Go To 4650
if (Line(1:4) == 'FTHA') Go To 4655
if (Line(1:4) == 'DFMD') Go To 4656
if (Line(1:4) == 'KEON') Go To 4660
if (Line(1:4) == 'TWCO') Go To 4661
if (Line(1:4) == 'ADDC') Go To 4662
if (Line(1:4) == 'SAVE') Go To 4663
if (Line(1:4) == 'DEBU') Go To 4700
if (Line(1:4) == 'ZSPI') Go To 4800
if (Line(1:4) == 'SPIN') Go To 4850
if (Line(1:4) == 'EXFA') Go To 4900
if (Line(1:4) == 'ONEG') Go To 4901
if (Line(1:4) == 'ROTP') Go To 5000
if (Line(1:4) == 'HLGA') Go To 5002
if (Line(1:4) == 'FLIP') Go To 5020
if (Line(1:4) == 'PMTI') Go To 6000
if (Line(1:4) == 'STAT') Go To 6010
if (Line(1:4) == 'ADDF') Go To 6020
if (Line(1:4) == 'FILE') Go To 6030
if (Line(1:4) == 'ITPR') Go To 7100
if (Line(1:4) == 'PROP') Go To 7200
if (Line(1:4) == 'NOPR') Go To 7201
if (Line(1:4) == 'NOX ') Go To 8700
if (Line(1:4) == 'PSDC') Go To 8900
if (Line(1:4) == 'USEX') Go To 8901
if (Line(1:4) == 'NEG2') Go To 8903
if (Line(1:4) == 'MSYM') Go To 8904
if (Line(1:4) == 'ITDI') Go To 8905
if (Line(1:4) == 'FCKA') Go To 8906
if (Line(1:4) == 'DEPT') Go To 8907

if (Line(1:4) == 'FALC') Go To 30000

if (Line(1:4) == 'END ') Go To 9000

if (lTtl) Go To 1101
write(u6,*) 'Unidentified key word:',Key
call FindErrorLine()
call Quit_OnUserError()

!>>>>>>>>>>>>> TITL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1100 continue
Line = Get_Ln(LuSpool)
lTtl = .true.
1101 continue
nTit = nTit+1
if (nTit == 1) then
  Title(1) = Line(1:72)
else
  if (nTit == 2) call WarningMessage(1,'More than one title line!')
end if
!if (nTit <= MxTit) Title(nTit) = Line(1:72)
goto 999

!>>>>>>>>>>>>> ITER <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1200 continue
Line = Get_Ln(LuSpool)
Line(179:180) = '-1'
call Put_Ln(Line)
call Get_I(1,nIter,2)
if (nIter(1) == -1) nIter(1) = nIter(0)
goto 1000

!>>>>>>>>>>>>> OCCU <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1300 continue
if (FermSet) then
  call WarningMessage(2,'Options OCCUpied and FERMi are mutually exclusive')
  call Abend()
end if
if (SpinSet) then
  call WarningMessage(2,'Keyword OCCUpied and ZSPIn/SPIN are mutually exclusive')
  call Abend()
end if
if (CharSet) then
  call WarningMessage(2,'Options OCCUpied and CHARGE are mutually exclusive')
  call Abend()
end if
do iD=1,nD
  Line = Get_Ln(LuSpool)
  call Get_I(1,nOcc(1,iD),nSym)
  if (iD == 1) then
    call Put_iArray('nIsh',nOcc(1,iD),nSym)
  else
    call Put_iArray('nIsh beta',nOcc(1,iD),nSym)
  end if
end do

iOccu = 1
if (nD == 1) then
  do i=1,nSym
    Tot_El_Charge = Tot_El_Charge-real(2*nOcc(i,1),kind=wp)
  end do
else
  do i=1,nSym
    Tot_El_Charge = Tot_El_Charge-real(nOcc(i,1)+nOcc(i,2),kind=wp)
  end do
end if
Aufb = .false. ! Disable default action.
Teee = .false.
OccSet = .true.
Cho_Aufb = .false.
UHFSet = .true.
goto 1000

!>>>>>>>>>>>>> ORBI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1400 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nOrb,nSym)
iOrbi = 1
do iSym=1,nSym
  nDel(iSym) = nBas(iSym)-nOrb(iSym)
end do
goto 1000

!>>>>>>>>>>>>> FROZ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1500 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nFro,nSym)
call Put_iArray('nFro',nFro,nSym)
iFroz = 1
goto 1000

!>>>>>>>>>>>>> OVLD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1700 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,DelThr)
goto 1000

!>>>>>>>>>>>>> PRLS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1800 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,iPri)
iPrint = max(iPri,iPrint)
goto 1000

!>>>>>>>>>>>>> PROR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
1900 continue
Line = Get_Ln(LuSpool)
Line(179:180) = '-1'
call Put_Ln(Line)
call Get_I1(1,iPrOrb)
if (iPrOrb >= 2) then
  call Get_F1(2,ThrEne)
  call Get_I1(3,iPrForm)
else
  call Get_I1(2,iPrForm)
end if
goto 1000

!>>>>>>>>>>>>> KEEP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2000 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,iDKeep)
goto 1000

!>>>>>>>>>>>>> STAR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2100 continue
Line = Get_Ln(LuSpool)
Line(157:180) = '-1 -1 -1 -1 -1 -1 -1 -1'
call Put_Ln(Line)
call Get_I(1,LstVec,7)
! temporary hack to use density
if (LstVec(1) == 3) then
  InVec = 3
  RTemp = 0.1_wp
  TemFac = 0.1_wp
  TStop = 0.005_wp
  !Aufb = .false.
  !Teee = .false.
end if
goto 1000

! Explicit Start Options

!>>>>>>>>>>>>> CHOL <<<<<<< Cholesky Default Input <<<<<<<<<<<<<<<<<

20000 continue
Chol = .true.
call Cho_scf_rdInp(.true.,LuSpool)
if (ALGO >= 2) then
  DDnOFF = .true. !do not use diffential density
  MiniDn = .false.
end if
goto 1000

!>>>>>>>>>>>>> CHOI <<<<<<< Cholesky Custom  Input <<<<<<<<<<<<<<<<<

! Cholesky with user-defined settings.

20001 continue
Chol = .true.
DDnOFF = .false. ! reset to default value
MiniDn = .true.  ! reset to default value
call Cho_scf_rdInp(.false.,LuSpool)
if (ALGO >= 2) then
  DDnOFF = .true. !do not use diffential density
  MiniDn = .false.
else
  DECO = .false.
end if
goto 1000

!>>>>>>>>>>>>> CONS <<<<<<<<<<<< Constrained SCF <<<<<<<<<<<
20002 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nConstr,nSym)
call Izero(indxC,16*2*8)
do i=1,nSym
  MxConstr = max(MxConstr,nConstr(i))
  Line = Get_Ln(LuSpool)
  call Get_I(1,indxC(1,1,i),nConstr(i))
  if (nConstr(i) > 16) then
    write(u6,*) ' Max nConstr=16. Increase 1st dim of indxC and recompile'
    call Abend()
  end if
  do j=1,nConstr(i)
    if ((indxC(j,1,i) /= 1) .and. (indxC(j,1,i) /= -1)) then
      write(u6,*) ' Only 1 and -1 are accepted values'
      call Abend()
    else if (indxC(j,1,i) == 1) then
      indxC(j,2,i) = 2
    else
      indxC(j,1,i) = 2
      indxC(j,2,i) = 1
    end if
  end do
end do
goto 1000

!>>>>>>>>>>>>> CORE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21000 continue
InVec = 0
LstVec(1) = 4
LstVec(2) = -1
goto 1000
!>>>>>>>>>>>>> NDDO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21001 continue
if (Chol) then
  call WarningMessage(1,'RdInp: NDDO and Cholesky are incompatible!!; NDDO Option ignored')
else
  InVec = 1
  LstVec(1) = 5
  LstVec(2) = -1
end if
goto 1000
!>>>>>>>>>>>>> LUMO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21002 continue
InVec = 2
One_Grid = .true.
LstVec(1) = 2
LstVec(2) = -1
goto 1000
!>>>>>>>>>>>>> FILE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
6030 continue
InVec = 2
One_Grid = .true.
LstVec(1) = 2
LstVec(2) = -1
Line = Get_Ln(LuSpool)
call fileorb(Line,SCF_FileOrb)
#ifdef _HDF5_
if (mh5_is_hdf5(SCF_FileOrb)) then
  isHDF5 = .true.
  fileorb_id = mh5_open_file_r(SCF_FileOrb)
end if
#endif
goto 1000
!>>>>>>>>>>>>> GSSR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21003 continue
InVec = 9
One_Grid = .true.
LstVec(1) = 1
LstVec(2) = -1
goto 1000
!>>>>>>>>>>>>> REST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
21004 continue
write(u6,*) 'WARNING: Option REST is redundant.'
goto 1000
!>>>>>>>>>>>>> THRE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2200 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,EThr)
call Get_F1(2,DThr)
call Get_F1(3,FThr)
call Get_F1(4,DltNTh)
!tbp, may 2013: no thr modification with Cholesky
!tbp if (DoCholesky) then
!tbp   write(ww,'(a,es20.8)') 'Detected Cholesky or RI/DF calculation BUT user specified EThr will be used. Ethr = ',EThr
!tbp   call WarningMessage(1,ww)
!tbp end if
!if (  DThr*1.0e-2_wp > EThr) then
!  write(u6,*)
!  write(u6,*) '----> WARNING! <----'
!  write(u6,*)
!  write(u6,*) ' The value of DThr is inconsistent with'
!  write(u6,*) ' with the value of EThr. The code will'
!  write(u6,*) ' automatically reset the value to something'
!  write(u6,*) ' more reasonable based on the requested'
!  write(u6,*) ' threshold of the energy.'
!  DThr = 100.0_wp*EThr
!end if
!if (  FThr*1.0e-2_wp > EThr) then
!  write(u6,*)
!  write(u6,*) '----> WARNING! <----'
!  write(u6,*)
!  write(u6,*) ' The value of FThr is inconsistent with'
!  write(u6,*) ' with the value of EThr. The code will'
!  write(u6,*) ' automatically reset the value to something'
!  write(u6,*) ' more reasonable based on the requested'
!  write(u6,*) ' threshold of the energy.'
!  FThr = 100.0_wp*EThr
!end if
!if (DltNTh*1.0e-2_wp > sqrt(EThr)) then
!  write(u6,*)
!  write(u6,*) '----> WARNING! <----'
!  write(u6,*)
!  write(u6,*) ' The value of DltNTh is inconsistent with'
!  write(u6,*) ' with the value of EThr. The code will'
!  write(u6,*) ' automatically reset the value to something'
!  write(u6,*) ' more reasonable based on the requested'
!  write(u6,*) ' threshold of the energy.'
!  DltNTh = 100.0_wp*Sqrt(EThr)
!end if
goto 1000

!>>>>>>>>>>>>> NODI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2300 continue
Diis = .false.
goto 1000

!>>>>>>>>>>>>> DIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2400 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,DiisTh)
goto 1000

!>>>>>>>>>>>>> OCCN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2500 continue
if (iOccu /= 1) then
  call WarningMessage(2,' Input Error!;The OCCNumber option works only if  the OCCUpied option has been specified!')
  call Abend()
end if
lthSet_a = 0
lthSet_b = 0
do iSym=1,nSym
  lthSet_a = lthSet_a+nOcc(iSym,1)
  if (nD == 2) lthSet_b = lthSet_b+nOcc(iSym,2)
end do
nOccSet_e = max(lthSet_a,lthSet_b)
!write(u6,'(a,i5)') 'rdinp: lthset_a ',lthset_a
!write(u6,'(a,i5)') 'rdinp: lthset_b ',lthset_b

! Note, it is dangerous to read Line first. There may be many
! lines with occupation numbers.

call mma_allocate(OccSet_e,nOccSet_e,nD,Label='OccSet_e')
call FZero(OccSet_e,nOccSet_e*nD)
read(LuSpool,*,end=902,Err=903) (OccSet_e(i,1),i=1,lthSet_a)
if (nD == 2) read(LuSpool,*,end=902,Err=903) (OccSet_e(i,2),i=1,lthSet_b)

Tot_El_Charge = Zero
do iD=1,nD
  do i=1,nOccSet_e
    Tot_El_Charge = Tot_El_Charge-OccSet_e(i,iD)
  end do
end do
UHFSet = .true.
iCoCo = 1
goto 1000

!>>>>>>>>>>>>> MCCN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Just as OCCN, but for muons
2510 continue
if (iOccu /= 1) then
  call WarningMessage(2,' Input Error!;The OCCNumber option works only if the OCCUpied option has been specified!')
  call Abend()
end if
lthSet_a = 0
lthSet_b = 0
do iSym=1,nSym
  lthSet_a = lthSet_a+nOcc(iSym,1)
  if (nD == 2) lthSet_b = lthSet_b+nOcc(iSym,2)
end do
nOccset_m = max(lthSet_a,lthSet_b)
!write(u6,'(a,i5)') 'rdinp: lthset_a ',lthset_a
!write(u6,'(a,i5)') 'rdinp: lthset_b ',lthset_b

! Note, it is dangerous to read Line first. There may be many
! lines with occupation numbers.

call mma_allocate(OccSet_m,nOccSet_m,nD,Label='OccSet_m')
call FZero(OccSet_m,nOccSet_m*nD)
read(LuSpool,*,end=902,Err=903) (OccSet_m(i,1),i=1,lthSet_a)
if (nD == 2) read(LuSpool,*,end=902,Err=903) (OccSet_m(i,2),i=1,lthSet_b)

Tot_Ml_Charge = Zero
do iD=1,nD
  do i=1,nOccSet_m
    Tot_Ml_Charge = Tot_Ml_Charge-OccSet_m(i,iD)
  end do
end do
UHFSet = .true.
iCoCo = 1
goto 1000

!>>>>>>>>>>>>> IVO  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2600 continue
kIvo = 1
goto 1000

!>>>>>>>>>>>>> UHF  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2700 continue
if (UHFSet) call sysAbendMsg('rdinp','incorrect input','UHF keyword should be placed before others')
nD = 2
MiniDn = .false.
goto 1000

!>>>>>>>>>>>>> HFC  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2701 continue
UHF_HFC = .true.
goto 1000

!>>>>>>>>>>>>> NODA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
2900 continue
Damping = .false.
goto 1000

!>>>>>>>>>>>>> CONV <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3000 continue
DSCF = .false.
goto 1000

!>>>>>>>>>>>>> DISK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3100 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,nDisc)
call Get_I1(2,nCore)
goto 1000

!>>>>>>>>>>>>> THIZ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3200 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,Thize)
goto 1000

!>>>>>>>>>>>>> SIMP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3300 continue
PreSch = .true.
goto 1000

!>>>>>>>>>>>>> NOMI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3400 continue
MiniDn = .false.  !do not use minimized density diff
goto 1000

!>>>>>>>>>>>>> TDEN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3401 continue
DDnOFF = .true. !do not use diffential density
MiniDn = .false.
TDen_UsrDef = .true.
goto 1000

!>>>>>>>>>>>>> WRDE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3500 continue
WrOutD = .true.
goto 1000

!>>>>>>>>>>>>> C1DI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3600 continue
c1Diis = .true.
goto 1000

!>>>>>>>>>>>>> QUAD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3700 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,QudThr)
goto 1000

!>>>>>>>>>>>>> RS-R <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3750 continue
RSRFO = .true.
RGEK = .false.
goto 1000

!>>>>>>>>>>>>> S-GE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3751 continue
RGEK = .true.
RSRFO = .false.
goto 1000

!>>>>>>>>>>>>> SCRA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3800 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,ScrFac)
Scrmbl = .true.
goto 1000

!>>>>>>>>>>>>> EXTR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
3900 continue
call WarningMessage(1,'EXTRACT option is redundant and is ignored!')
goto 1000

!>>>>>>>>>>>>> RFPE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4000 continue
RFpert = .true.
goto 1000

!>>>>>>>>>>>>> QNRT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4100 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,QNRTh)
goto 1000

!>>>>>>>>>>>>> AUFB <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4200 continue
call WarningMessage(2,' Error: Keyword AUFBau is now obsolete!;Use keyword CHARge')
call Abend()
if (Chol) then
  DDnOFF = .true.
  DECO = .true.
  MiniDn = .false.
end if
Line = Get_Ln(LuSpool)
Line(180:180) = '2'
call Put_Ln(Line)
call Get_I(1,iArray,nD+1)
do i=1,nD
  nAufb(i) = iArray(i)
end do
iAuf = iArray(nD+2)
if (IfAufChg) then
  call WarningMessage(2,'Option AUFBau is mutually exclusive CHARge')
  call Abend()
end if
if (nD == 1) then
  Tot_El_Charge = -Two*real(nAufb(1),kind=wp)
else
  Tot_El_Charge = -real(nAufb(1)+nAufb(2),kind=wp)
end if
IfAufChg = .true.
4210 continue
select case (iAuf)
  case (0)
    Teee = .false.
  case (1)
    RTemp = Half
    TemFac = 0.4_wp
    TStop = 0.005_wp
  case (2)
    RTemp = Half
    TemFac = 0.46_wp
    TStop = 0.01_wp
  case (3)
    RTemp = Half
    TemFac = 0.61_wp
    TStop = 0.025_wp
  case (4)
    RTemp = One
    TemFac = 0.73_wp
    TStop = 0.06_wp
  case default
    if (iAuf /= 5) call WarningMessage(1,' RdInp: Aufbau case must be in the range 0-5;Using case 5!')
    RTemp = One
    TemFac = 0.87_wp
    TStop = 0.15_wp
end select
UHFSet = .true.
goto 1000
!>>>>>>>>>>>>> FERM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4250 continue
if (OccSet) then
  call WarningMessage(2,'Options OCCUpied and FERMi are mutually exclusive')
  call Abend()
end if
if (Chol) then
  DDnOFF = .true.
  DECO = .true.
  MiniDn = .false.
end if
Line = Get_Ln(LuSpool)
call Get_I1(1,iAuf)
FermSet = .true.
goto 4210
!>>>>>>>>>>>>> TEEE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4300 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,RTemp)
call Get_F1(2,TemFac)
call Get_F1(3,TStop)
if (TStop < Zero) then
  call WarningMessage(2,'Input Error!; End temperture < 0.0 ')
  call Abend()
end if
if (TemFac < Zero) then
  call WarningMessage(2,'Input Error!; Temperture factor < 0.0 ')
  call Abend()
end if
if (RTemp < Zero) then
  call WarningMessage(2,'Input Error!; Start temperature < 0.0 ')
  call Abend()
end if
if (RTemp < TStop) then
  call WarningMessage(2,'Input Error!; End temperture > start temperature ')
  call Abend()
end if
goto 1000
!>>>>>>>>>>>>> CHAR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4400 continue
Line = Get_Ln(LuSpool)
call Get_I(1,iArray,1)
Tot_Charge = real(iArray(1),kind=wp)
nAufb(1) = -1
nAufb(2) = -1
if (IfAufChg) then
  call WarningMessage(2,'Input Error!; Option CHARge is mutually exclusive to AUFBau')
  call Abend()
end if
if (OccSet) then
  call WarningMessage(2,'Input Error!; Option CHARge is mutually exclusive to OCCUpied')
  call Abend()
end if
IfAufChg = .true.
CharSet = .true.
goto 1000
!>>>>>>>>>>>>> NOTE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4500 continue
Teee = .false.
goto 1000

!>>>>>>>>>>>>> KSDF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4600 continue
Line = Get_Ln(LuSpool)
call UpCase(Line)
Line = adjustl(Line)
KSDFT = Line(1:80)
nFunc = 0
read(Line,*,iostat=istatus) nFunc
if ((istatus == 0) .and. (nFunc > 0)) then
  KSDFT = Custom_Func
  LuCF = IsFreeUnit(10)
  call molcas_open(LuCF,Custom_File)
  write(LuCF,*) trim(KSDFT),nFunc
  do i=1,nFunc
    Line = Get_Ln(LuSpool)
    write(LuCF,*) trim(Line)
  end do
  close(LuCF)
end if
goto 1000

!>>>>>>>>>>>>> DFCF <<<< Factors to scale exch. and corr. <<
4605 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,CoefX)
call Get_F1(2,CoefR)
!call put_dscalar('DFT exch coeff',CoefX)
!call put_dscalar('DFT corr coeff',CoefR)
goto 1000

!>>>>>>>>>>>>> OFEM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4650 continue
Do_OFemb = .true.
Line = Get_Ln(LuSpool)
call UpCase(Line)
Line = adjustl(Line)
OFE_KSDFT = Line(1:16)
write(u6,*) '  --------------------------------------'
write(u6,*) '   Orbital-Free Embedding Calculation'
write(u6,*) '  --------------------------------------'
if (OFE_KSDFT(1:4) == 'LDTF') then
  write(u6,*) '    T_nad potential   : Thomas-Fermi    '
else
  write(u6,*) '    T_nad potential   : ',OFE_KSDFT(1:4)
end if
if (KEonly) then
  write(u6,*) '    Exc_nad potential :  None           '
else
  write(u6,*) '    Exc_nad potential : ',OFE_KSDFT(6:10)
end if
write(u6,*) '  --------------------------------------'
write(u6,*)
goto 1000

!>>>>>>>>>>>>> FTHA <<<< threshold for Freeze-n-Thaw <<<<<<<
4655 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,ThrFThaw)
goto 1000

!>>>>>>>>>>>>> DFMD <<<< fraction of correlation potential <
4656 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,dFMD)
call Get_F1(2,Xsigma)
if (dFMD+Xsigma < Zero) then
  write(u6,*) ' *** Warning: arguments to DFMD must be nonnegative!'
  write(u6,*) ' ***          I will take their ABS !!! '
  dFMD = abs(dFMD)
  Xsigma = abs(Xsigma)
end if
goto 1000

!>>>>>>>>>>>>> KEON <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4660 continue
KEonly = .true.
if (.not. Do_OFemb) then
  write(u6,*) ' *** Warning:  KEonly works in OFembedding runs!'
else
  write(u6,*) ' *** Warning:  Exc_nad set to NONE at this point ***'
end if
write(u6,*)
goto 1000

!>>>>>>>>>>>>> TWCO <<<<< activate Tw correlation <<<<<<<<<<
4661 continue
Do_Tw = .true.
goto 1000

!>>>>>>>>>>>>> ADDC << add correlation energy (CONStraint) <
4662 continue
Do_Addc = .true.
Line = Get_Ln(LuSpool)
call UpCase(Line)
Line = adjustl(Line)
ADDC_KSDFT = Line(1:80)
goto 1000

!>>>>>>>>>>>>> SAVE << Spin-Averaged wavelets (CONStraint) <
4663 continue
Do_SpinAV = .true.
goto 1000

!>>>>>>>>>>>>> DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4700 continue
Diis = .false.
MiniDn = .false.
Damping = .false.
goto 1000
!>>>>>>>>>>>>> ZSPI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4800 continue
if (SpinSet) then
  call WarningMessage(2,'Multiple definition of SPIN/ZSPIn')
  call Abend()
end if
if (OccSet) then
  call WarningMessage(2,'Input Error!; Keywords OCCUpied and ZSPIn are mutually exclusive')
  call Abend()
end if
Line = Get_Ln(LuSpool)
call Get_I(1,iArray,1)
iAu_ab = iArray(1)
if ((nD /= 2) .and. (iAu_ab /= 0)) then
  call WarningMessage(2,'ZSPIn different from 0 requires UHF before it')
  call Abend()
end if
SpinSet = .true.
goto 1000
!>>>>>>>>>>>>> SPIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4850 continue
if (SpinSet) then
  call WarningMessage(2,'Multiple definition of SPIN/ZSPIn')
  call Abend()
end if
if (OccSet) then
  call WarningMessage(2,'Input Error!; Keywords OCCUpied and SPIN are mutually exclusive')
  call Abend()
end if
Line = Get_Ln(LuSpool)
call Get_I(1,iArray,1)
iAu_ab = iArray(1)-1
if (iAu_ab < 0) then
  call WarningMessage(2,'SPIN must be a positive integer')
  call Abend()
end if
if (iAu_ab /= 0) then
  nD = 2
  MiniDn = .false.
end if
if ((nD /= 2) .and. (iAu_ab /= 0)) then
  call WarningMessage(2,'SPIN greater than 1 requires UHF before it')
  call Abend()
end if
SpinSet = .true.
goto 1000
!>>>>>>>>>>>>> EXFA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4900 continue
ExFac = Zero
goto 1000

!>>>>>>>>>>>>> ONEG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
4901 continue
One_Grid = .true.
goto 1000

!>>>>>>>>>>>>> ROTP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
5000 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,RotLev)
call Get_F1(2,RotFac)
call Get_F1(3,RotMax)
write(u6,'(a,ES15.3)') 'Fock matrix levelshift   ',RotLev
write(u6,'(a,ES15.3)') 'Fock matrix scaling      ',RotFac
write(u6,'(a,ES15.3)') 'Fock matrix max rotation ',RotMax
goto 1000
!>>>>>>>>>>>>> HLGA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
5002 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,HLgap)
write(u6,'(a,ES15.3)') 'Minimum HOMO-LUMO gap    ',HLgap
DoHLgap = .true.
QNRTh = Zero
goto 1000
!>>>>>>>>>>>>> FLIP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
5020 continue
Line = Get_Ln(LuSpool)
Line(178:180) = '0.1'
call Put_Ln(Line)
call Get_I(1,iArray,1)
call Get_F1(2,FlipThr)
MaxFlip = iArray(1)
!write(u6,*) 'MaxFlip:',MaxFlip
!write(u6,*) 'FlipThr:',FlipThr
goto 1000
!>>>>>>>>>>>>> PMTI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Time (CPU *and* wall) subroutine pmat_scf (2-el Fock matrix)
6000 continue
PmTime = .true.
goto 1000
!>>>>>>>>>>>>> STAT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Print Statistic information
6010 continue
iStatPRN = 1
goto 1000
!>>>>>>>>>>>>> ADDF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Add the fragment atoms to the MOLDEN file
6020 continue
AddFragments = .true.
goto 1000
!>>>>>>>>>>>>> ITPR <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
7100 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,iterprlv)
goto 1000
!>>>>>>>>>>>>> PROP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
7200 continue
InVec = 2
One_Grid = .true.
LstVec(1) = 2
LstVec(2) = -1
OnlyProp = .true.
goto 1000
!>>>>>>>>>>>>> SKIP PROP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
7201 continue
NoProp = .true.
goto 1000
!>>>>>>>>>>>>> NOX  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: skip exchange in Fock matrix build
8700 continue
NoExchange = .true.
goto 1000
!>>>>>>>>>>>>> PSDC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: check that full integral matrix is PSD
! input=1: diagonalization
! input=2: Cholesky decomposition
8900 continue
Line = Get_Ln(LuSpool)
call Get_I(1,Mode,1)
goto 1000
!>>>>>>>>>>>>> USEX <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Debug option: use exact diagonal (1) or off-dagonal (2) blocks when
! checking PSD
8901 continue
Line = Get_Ln(LuSpool)
call Get_I(1,Mode,1)
goto 1000
!>>>>>>>>>>>>> NEG2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Specify action when negative two-electron energies are
! encountered (stop, warn, continue).
8903 continue
Line = Get_Ln(LuSpool)
call UpCase(Line)
Neg2_Action = Line(1:4)
if ((Neg2_Action /= 'STOP') .and. (Neg2_Action /= 'WARN') .and. (Neg2_Action /= 'CONT')) then
  write(u6,'(A,A)') 'Illegal input for NEG2 keyword: ',Line(1:4)
  !call FindErrorLine()
  call Quit_OnUserError()
end if
goto 1000
!>>>>>>>>>>>>> MSYM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8904 continue
MSYMON = .true.
goto 1000
!>>>>>>>>>>>>> ITDIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8905 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,iTer2run)
goto 1000
!>>>>>>>>>>>>> ITDIIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8906 continue
Line = Get_Ln(LuSpool)
call UpCase(Line)
if (index(Line,'TRUE') /= 0) then
  FckAuf = .true.
else
  FckAuf = .false.
end if
goto 1000
!>>>>>>>>>>>>> DEPTH  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
8907 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,kOptim_Max)
if (kOptim_Max > MxOptm) then
  write(u6,*) 'kOptim_Max>MxOptm'
  write(u6,*) 'kOptim_Max=',kOptim_Max
  write(u6,*) 'MxOptm=',MxOptm
  write(u6,*) 'Modify infscf.F90 and recompile!'
  call Abend()
end if
goto 1000
!>>>>>>>>>>>>> FALC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
30000 continue
Falcon = .true.
goto 1000

!>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
9000 continue

if (iPrint >= 3) iStatPRN = 1

Tot_El_charge = Tot_El_Charge+Tot_Ml_Charge

! xml tag method

if (nD == 1) then
  if (KSDFT == 'SCF') then
    call xml_cDump('method','','',0,'rhf',1,1)
  else
    call xml_cDump('method','','',0,'rdft',1,1)
  end if
else
  if (KSDFT == 'SCF') then
    call xml_cDump('method','','',0,'uhf',1,1)
  else
    call xml_cDump('method','','',0,'udft',1,1)
  end if
end if

! Even or odd number of electrons

if (.not. SpinSet) then
  nnn = int(tot_nuc_charge-tot_charge+half)
  if ((nnn/2)*2 /= nnn) iAu_ab = 1
end if

! Check start orbital priority list

if ((.not. OccSet) .and. (.not. FermSet)) then
  !write(u6,*) 'rdinp: Checking OCCU/FERM'
  call VecFind(OccSet,FermSet,CharSet,SpinSet)
  if (OccSet .and. (.not. FermSet)) then
    !write(u6,*) 'Using OCCU'
    Aufb = .false.
    Teee = .false.
    Cho_Aufb = .false.
  else if (FermSet .and. (.not. OccSet)) then
    !write(u6,*) 'Using FERM'
    Aufb = .true.
    Teee = .true.
    Cho_Aufb = .true.
  else
    call WarningMessage(2,'Internal (logical) error in rdinp_scf')
    call Abend()
  end if
end if
if (SpinSet) then
  nnn = int(Tot_Nuc_Charge-Tot_Charge-iAu_ab+half)
  if ((nnn/2)*2 /= nnn) then
    call WarningMessage(2,'Input error!;zSpin inconsistent with number of electrons')
    call Abend()
  end if
end if

! Check if certain input parameters are not in conflict

if ((iCoCo == 1) .and. (iOccu == 0)) then
  call WarningMessage(2,'Input error!; The OCCNumber option require that the OCCUpied option is specified!')
  call Abend()
end if

if (max(nIter(0),nIter(1)) > MxIter) then
  call WarningMessage(1,'Input error!;Number of iterations specified in input exceeds allowed maximum!')
  !write(u6,*) 'nIter(0)=',nIter(0)
  !write(u6,*) 'nIter(1)=',nIter(1)
  !write(u6,*) 'MxIter=',MxIter
  !write(u6,*)
  nIter(0) = MxIter
  nIter(1) = MxIter

  !write(u6,*) 'Number of iteration reset to ',MxIter
  !write(u6,*)
end if

if ((iOrbi == 1) .and. ((InVec /= 2) .and. (InVec /= 4))) then
  call WarningMessage(2,'Input error!; The ORBITAL option can only be used with input orbitals!;'// &
                      'Possible exclusive options are: LUMORB RESTART')
  call Abend()
end if

if ((iFroz == 1) .and. (InVec == 3)) then
  call WarningMessage(2,'Input error!; The FROZEN option does not work with an density matrix input')
  call Abend()
end if

if ((InVec == 1) .and. Aufb .and. (.not. Chol)) then
  call WarningMessage(2,'Input error!; Aufbau not compatible with NDDO option')
  call Abend()
end if

if ((Invec < -1) .or. (InVec > 9)) then
  call WarningMessage(2,'Input error!; inappropriate value for Invec')
  call Abend()
end if

if ((nD == 1) .and. UHF_HFC) call sysAbendMsg('rdinp','incorrect input','HFC keyword should be used with UHF')

! Print out warning informations (if any)

if ((iFroz == 1) .and. ((InVec == 2) .or. (InVec == 4))) &
  call WarningMessage(1,'RdInp: Warning!;Option FROZEN is used together with input orbitals as initial guess.;'// &
                      'Freezing of orbitals is performed at MO level.;First orbitals in each symmetry will not be modified.')

if (Aufb .and. (iFroz == 1)) then
  do iSym=1,nSym
    nAufb(1) = nAufb(1)+nFro(iSym)
    if (nD == 2) nAufb(2) = nAufb(2)+nFro(iSym)
  end do
  call ICopy(nSym,[0],0,nFro,1)
  call WarningMessage(2,'Input error!;Aufbau not allowed with frozen orbitals')
  call Abend()
end if

! Check parameters for semi direct SCF

! If semi-direct and I/O buffer value not specified set to default value.
if ((nCore == 0) .and. (nDisc /= 0)) nCore = lDaRec*nSect*2*8/1024
nCore = ((nCore+7)/8)*8

! Adjust disk size to multiple of I/O buffer
if (nDisc /= 0) nDisc = (nDisc*1024+nCore-1)/1024

nDisc = min(10*Allocdisk(),nDisc)

! Set up parameters that follow from others

if ((.not. Diis) .and. (.not. Damping)) iDKeep = 0

if ((InVec == 3) .and. (max(nIter(0),nIter(1)) == 0)) then
  iPrOrb = 0
  jVOut = 0
end if

! Check Cholesky vs. Aufbau
if (Aufb .and. DoCholesky) then
  Cho_Aufb = .true.
  !write(u6,*)
  !write(u6,*) ' ********** WARNING! *********'
  !write(u6,*) ' Cholesky SCF runs much slower with AufBau !'
  !write(u6,*) ' *** Do you really need AufBau in this case? ***'
  !write(u6,*)
end if

! Check CONS vs. UHF+OCCU
if ((MxConstr > 0) .and. (nD-1+iOCCU /= 2)) then
  call WarningMessage(2,'For CONStraints, keywords UHF and OCCUpied are compulsory!')
  call Abend()
end if

! Check CONS vs. ADDC
if ((MxConstr == 0) .and. Do_Addc) then
  call WarningMessage(0,' In absence of CONStraints, ADDCorrelation is de-activated!')
  Do_Addc = .false.
end if

! Check CONS vs. SAVE
if ((MxConstr == 0) .and. Do_SpinAV) then
  call WarningMessage(0,' In absence of CONStraints, SAVErage is de-activated!')
  Do_SpinAV = .false.
end if

if (Do_SpinAV) DThr = 1.0e4_wp ! reset because it is not meaningful
if (MxConstr > 0) InVec = 6

! Check parameters of KS-DFT

if (KSDFT /= 'SCF') then
  if (MiniDn) then
    if (jPrint >= 2) call WarningMessage(0,' Minimized-density-differences option turned off!')
    MiniDn = .false.
  end if
  if (Do_OFemb .and. (dFMD /= Zero)) then
    call WarningMessage(0,' KSDFT/OFE requires DFMD=0 for correlation potential!')
    write(u6,*) ' dFMD = ',dFMD
  end if
else
  if (Do_OFemb .and. (dFMD /= One)) then
    call WarningMessage(0,' HF/OFE may require DFMD=1 for correlation potential!')
    write(u6,*) ' dFMD = ',dFMD
  end if
end if

call Put_iScalar('SCF mode',nD-1)

LKon = (ALGO == 4)

Method = 'RHF-SCF '
if (nD == 2) Method = 'UHF-SCF '
if (kIvo /= 0) Method = 'IVO-SCF '
if (KSDFT /= 'SCF') Method = 'KS-DFT  '
call Put_cArray('Relax Method',Method,8)

ExFac = Get_ExFac(KSDFT)
if ((ExFac == Zero) .and. (.not. TDen_UsrDef) .and. (.not. Do_OFemb)) DDnOFF = .false. ! use always differential density
!DDnOFF = .true.
!MiniDn = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! remove local copy of standard input

call Close_LuSpool(LuSpool)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

! Error exit
902 continue
call WarningMessage(2,'Input error!;Error reading input file for OCCNO option')
call Abend()
903 continue
call WarningMessage(2,'Input error!;End of input file for OCCNO option')
call Abend()

end subroutine RdInp_scf
