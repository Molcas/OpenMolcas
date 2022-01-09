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

subroutine RdCtl_Seward(LuRd,lOPTO,Do_OneEl)

use SW_File
use AMFI_Info
use Basis_Info
use Center_Info
use Her_RW
use Period
use MpmC
use EFP_Module
use Real_Spherical, only: Sphere
use fortran_strings, only: str
use External_Centers
use Symmetry_Info, only: Symmetry_Info_Setup, iSkip, nIrrep, VarR, VarT
use Temporary_Parameters
use Integral_Parameters
use Sizes_of_Seward, only: S
use Real_Info, only: ThrInt, Rtrnc, CutInt, PkAcc, Thrs, E1, E2, RPQMin, SadStep, Shake, kVector, CoM
use DKH_Info
use RICD_Info, only: iRI_Type, LDF, Do_RI, Cholesky, Do_acCD_Basis, Skip_High_AC, DiagCheck, LocalDF, Do_nacCD_Basis, Thrshld_CD
use Logical_Info
use Gateway_Interfaces, only: GetBS
use Gateway_global, only: Run_Mode, G_Mode, S_Mode
#ifdef _FDE_
use Embedding_Global, only: embPot, embPotInBasis, embPotPath, outGridPathGiven, embWriteDens, embWriteEsp, embWriteGrad, &
                            embWriteHess
#endif
#ifndef _HAVE_EXTRA_
use XYZ
#endif
use Para_Info, only: MyRank
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif

implicit real*8(a-h,o-z)
external NucExp
#include "Molcas.fh"
#include "angtp.fh"
#include "constants.fh"
#include "constants2.fh"
#include "SysDef.fh"
#include "notab.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "rmat.fh"
#include "real.fh"
#include "print.fh"
#ifdef _HAVE_EXTRA_
#include "hyper.fh"
#endif
real*8 Lambda
character Key*180, KWord*180, Oper(3)*3, BSLbl*80, Fname*256, DefNm*13, Ref(2)*180, ChSkip*80, AngTyp(0:iTabMx)*1, dbas*(LENIN), &
          filename*180, KeepBasis*256, KeepGroup*180, Previous_Command*12, CtrLDK(10)*(LENIN), Directory*256, BasLib*256, &
          ExtBasDir*256
character(LEN=72) :: Header(2) = ['','']
character(LEN=80) :: Title(10) = ['','','','','','','','','','']
character(LEN=14) :: Vrsn = 'Gateway/Seward'
character(LEN=512) :: Align_Weights = 'MASS'
#include "cgetl.fh"
character*180 Get_Ln
external Get_Ln
logical lTtl, lSkip, lMltpl, DoRys, RF_read, Convert, IfTest, Exist, CutInt_UsrDef, ThrInt_UsrDef, MolWgh_UsrDef, CholeskyWasSet, &
        GWInput, NoAMFI, lOPTO, Do_OneEl
logical do1CCD
logical :: CSPF = .false.
logical APThr_UsrDef, Write_BasLib
integer Cho_MolWgh, BasisTypes(4), BasisTypes_Save(4), iGeoInfo(2), iOpt_XYZ, RC
parameter(Cho_CutInt=1.0D-40,Cho_ThrInt=1.0D-40, Cho_MolWgh=2)
real*8 NucExp, WellCff(3), WellExp(3), WellRad(3), OAMt(3), OMQt(3)
real*8, allocatable :: RTmp(:,:), EFt(:,:), DMSt(:,:), OrigTrans(:,:), OrigRot(:,:,:), mIsot(:)
integer, allocatable :: ITmp(:), nIsot(:,:), iScratch(:)
! Temporary buffer
integer, parameter :: nBuff = 10000
real*8, allocatable :: Buffer(:), Isotopes(:)
character*180, allocatable :: STDINP(:)
character Basis_lib*256, CHAR4*4
character*256 Project, GeoDir, temp1, temp2
integer StrnLn
external StrnLn
logical SymmSet
logical CoordSet, RPSet
logical BasisSet
logical GroupSet
logical DoneCoord
logical NoZMAT
logical ForceZMAT
logical NoDKroll
logical DoTinker
logical DoGromacs
logical OrigInput
logical OriginSet
logical FragSet
logical HyperParSet
logical WriteZMat
#ifdef _HAVE_EXTRA_
logical geoInput, oldZmat, zConstraints
#endif
logical EFgiven
logical Invert
real*8 HypParam(3), RandVect(3)
logical Vlct_, nmwarn, FOUND
logical Basis_test, lECP, lPP
logical :: lDMS = .false., lOAM = .false., lOMQ = .false., lXF = .false., lFAIEMP = .false.
#include "embpcharg.fh"
#ifdef _GROMACS_
integer, dimension(:), allocatable :: CastMM
integer, dimension(:,:), allocatable :: DefLA
real*8, dimension(:), allocatable :: FactLA
#endif
character*256 Message
parameter(MAX_XBAS=20)
character*12 xb_label(MAX_XBAS)
character*128 xb_bas(MAX_XBAS)
data WellCff/.35d0,0.25d0,5.2d0/
data WellExp/4.0d0,3.0d0,2.0d0/
data WellRad/-1.22d0,-3.20d0,-6.20d0/
#include "angstr.fh"
data DefNm/'basis_library'/
data IfTest/.false./
interface
  subroutine datimx(TimeStamp) bind(C,name='datimx_')
    use, intrinsic :: iso_c_binding, only: c_char
    character(kind=c_char) :: TimeStamp(*)
  end subroutine
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 3
iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
IfTest = .true.
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Buffer,nBuff,Label='Buffer')
!                                                                      *
!***********************************************************************
!                                                                      *
do i=0,iTabMx
  AngTyp(i) = Angtp(i)
  call UpCase(AngTyp(i))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call WhichMolcas(Basis_lib)
if (Basis_lib(1:1) /= ' ') then
  ib = index(Basis_lib,' ')-1
  if (ib < 1) call SysAbendMsg('rdCtl','Too long PATH to MOLCAS',' ')
  BasLib = Basis_lib(1:ib)//'/basis_library'
else
  BasLib = DefNm
end if
Write_BasLib = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call Qpg_cArray('Align_Weights',Found,lAW)
if (Found) call Get_cArray('Align_Weights',Align_Weights,512)
!                                                                      *
!***********************************************************************
!                                                                      *
CutInt_UsrDef = .false.
ThrInt_UsrDef = .false.
MolWgh_UsrDef = .false.
APThr_UsrDef = .false.
NoAMFI = .false.

iChk_RI = 0
iChk_CH = 0
iChk_DC = 0
iOpt_XYZ = -1

isXfield = 0
CholeskyThr = -9.99d9

Basis_Test = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
LuWr = 6
LuFS = -1
LuRdSave = -1

nTemp = 0
lMltpl = .false.

CholeskyWasSet = .false.
do1CCD = .false.
spanCD = -9.9d9
lTtl = .false.
RF_read = .false.
lSkip = .false.
NoDKroll = .false.
SymmSet = .false.
BasisSet = .false.
GroupSet = .false.
RPSet = .false.
CoordSet = .false.
DoneCoord = .false.
KeepGroup = 'FULL'
NoZMAT = .false.
ForceZMAT = .false.
DoTinker = .false.
DoGromacs = .false.
OrigInput = .false.
#ifdef _HAVE_EXTRA_
origin_input = .false.
geoInput = .false.
OldZmat = .false.
isHold = -1
nCoord = 0
ZConstraints = .false.
#endif
WriteZMat = .false.
nFragment = 0
iFrag = 0
FragSet = .false.
OriginSet = .false.
HyperParSet = .false.
stepFac1 = 60.0d0
iOptimType = 1
gradLim = 0.0d0
Do_OneEl = .true.
Vlct_ = .false.
#ifdef _FDE_
! Embedding
embPot = .false.
embPotInBasis = .false.
embPotPath = 'EMBPOT'
outGridPathGiven = .false.
embWriteDens = .false.
embWriteEsp = .false.
embWriteGrad = .false.
embWriteHess = .false.
call Put_iScalar('embpot',0)
#endif
DoEmPC = .false.
EFgiven = .false.
Invert = .false.
call Put_iScalar('agrad',0)

ScaleFactor = 1.0d0
lSTDINP = 0
iCoord = 0
iBSSE = -1
SymThr = 0.01d0
nTtl = 0

imix = 0
ifnr = -1
ign = 0
itype = 0
ExtBasDir = ' '
isxbas = 0
nmwarn = .true.

! Selective initialization

if (Run_Mode == S_Mode) then
  iShll = S%Mx_Shll
  mdc = S%Mx_mdc
else
  iShll = 0
  mdc = 0

  ChSkip = ' '
  do i=1,3
    Oper(i) = ' '
  end do
  nOper = 0
end if

iDNG = 0
BasisTypes(:) = 0
KeepBasis = ' '
! period
lthCell = 0
Cell_l = .false.
call izero(ispread,3)
call fzero(VCell,9)
!     Set local DF variables (dummy)
call LDF_SetInc()
rewind(LuRd)
! Count the number of calls to coord to set the number of fragments
! for a geo/hyper-calculation
400 Key = Get_Ln(LuRd)
KWord = Key
call UpCase(KWord)
if (KWord(1:4) == 'COOR') nFragment = nFragment+1
if ((KWord(1:4) == 'HYPE') .or. (KWord(1:4) == 'GEO ')) then
  call getenvf('Project',Project)
  call getenvf('GeoDir',GeoDir)
  temp1 = Project(1:index(Project,' ')-1)//'.Gateway.Input'
  temp2 = GeoDir(1:index(GeoDir,' ')-1)//'/'//Project(1:index(Project,' ')-1)//'.gwcopy.input'
  call fCopy(temp1,temp2,ierr)
  if (ierr /= 0) then
    write(6,*) '*** Detect Hyper input, but no GEO loop'
    call Quit_OnUserError()
  end if
end if
if (KWord(1:4) == 'END ') Go To 401
goto 400
401 call Put_iScalar('nCoordFiles',nFragment)
rewind(LuRd)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Run_Mode == G_Mode) then
  call RdNLst(LuRd,'GATEWAY')
  Previous_Command = '&Gateway'
else
  call RdNLst(LuRd,'SEWARD')
  Previous_Command = '&Seward'
end if

! Default setting of GWInput

GWInput = Run_Mode == G_Mode

! GWInput is a logical flag which tells if a keyword is Gateway or
! Seward specific. This is done down below for each keyword in the
! case that the keyword is not Seward or Gateway specific enter the
! line,
!
!   GWInput = Run_Mode == G_Mode
!
! in the section of the particular keyword
!                                                                      *
!***********************************************************************
!                                                                      *

! KeyWord directed input

call mma_allocate(STDINP,MxAtom*2,label='STDINP')

nDone = 0
998 lTtl = .false.
if (Basis_Test .and. (nDone == 1)) then
  nDone = 0
  Basis_Test = .false.
end if
9988 continue
Key = Get_Ln(LuRd)

9989 if ((Run_Mode == G_Mode) .and. (.not. GWInput)) then
  call WarningMessage(2,'Gateway input error!')
  write(LuWr,*) 'The keyword : "',Previous_Command,'" is not allowed in the Gateway input!'
  write(LuWr,*) 'This keyword most likely belongs in the Seward input section!.'
  write(LuWr,*)
  call Quit_OnUserError()
else if ((Run_Mode == S_Mode) .and. GWInput) then
  call WarningMessage(2,'Seward input error!')
  write(LuWr,*) 'The keyword : "',Previous_Command,'" is not allowed in the Seward input when the Gateway is used!'
  write(LuWr,*) ' Try putting the keyword in the Gateway input section!'
  write(LuWr,*)
  call Quit_OnUserError()
end if
GWInput = .false.

if (IfTest) write(LuWr,*) ' RdCtl: Processing:',Key
KWord = Key
call UpCase(KWord)
Previous_Command = KWord(1:4)
if (KWord(1:1) == '*') Go To 998
if (KWord == '') Go To 998
if (Basis_Test) nDone = 1

! KEYWORDs in ALPHABETIC ORDER!

if (KWord(1:4) == '1C-C') Go To 9022
if (KWord(1:4) == '1CCD') Go To 9022
if (KWord(1:4) == 'ACCD') Go To 8005
if (KWord(1:4) == 'ACD ') Go To 8004
if (KWord(1:4) == 'ALIG') Go To 8013
if (KWord(1:4) == 'AMFI') Go To 9761
if (KWord(1:4) == 'AMF1') Go To 8761
if (KWord(1:4) == 'AMF2') Go To 8762
if (KWord(1:4) == 'AMF3') Go To 8763
if (KWord(1:4) == 'AMPR') Go To 9951
if (KWord(1:4) == 'ANGM') Go To 995
if (KWord(1:4) == 'APTH') Go To 38
if (KWord(1:4) == 'AUXS') Go To 912
if (KWord(1:4) == 'BASD') Go To 7700
if (KWord(1:4) == 'BASI') Go To 920
if (KWord(1:4) == 'BASL') Go To 8029
if (KWord(1:4) == 'BSSE') Go To 6020
if (KWord(1:4) == 'BSSH') Go To 9121
if (KWord(1:4) == 'BSSM') Go To 9760
if (KWord(1:4) == 'CDTH') Go To 8001
if (KWord(1:4) == 'CELL') Go To 887
if (KWord(1:4) == 'CENT') Go To 973
if (KWord(1:4) == 'CHEC') Go To 39
if (KWord(1:4) == 'CLDF') Go To 42
if (KWord(1:4) == 'CHOL') Go To 9091
if (KWord(1:4) == 'CHOI') Go To 9092
if (KWord(1:4) == 'CLIG') Go To 9000
if (KWord(1:4) == 'CONS') Go To 8010
if (KWord(1:4) == 'COOR') Go To 6000
if (KWord(1:4) == 'CSPF') Go To 9110
if (KWord(1:4) == 'CUTO') Go To 942
if (KWord(1:4) == 'DCRN') Go To 958
if (KWord(1:4) == 'DIAG') Go To 9087
if (KWord(1:4) == 'DIRE') Go To 9770
if (KWord(1:4) == 'DK1H') Go To 9001
if (KWord(1:4) == 'DK2H') Go To 9002
if (KWord(1:4) == 'DK3F') Go To 9004
if (KWord(1:4) == 'DK3H') Go To 9003
if (KWord(1:4) == 'DOAN') Go To 8019
if (KWord(1:4) == 'DOFM') Go To 8006
if (KWord(1:4) == 'DOUG') Go To 976
if (KWord(1:4) == 'DOWN') Go To 1002
if (KWord(1:4) == 'DSHD') Go To 996
if (KWord(1:4) == 'ECPS') Go To 912
if (KWord(1:4) == 'EFLD') Go To 993
if (KWord(1:4) == 'EFP ') Go To 9088
#ifdef _FDE_
if (KWord(1:4) == 'EMBE') Go To 666
#endif
if (KWord(1:4) == 'EMFR') Go To 8035
if (KWord(1:4) == 'EMPC') Go To 974
if (KWord(1:4) == 'EPOT') Go To 9932
if (KWord(1:4) == 'EXPE') Go To 9771
if (KWord(1:4) == 'EXTR') Go To 9960
if (KWord(1:4) == 'FAKE') Go To 9759
if (KWord(1:4) == 'FAT-') Go To 8004
if (KWord(1:4) == 'FCD ') Go To 9091
if (KWord(1:4) == 'FILE') Go To 904
if (KWord(1:4) == 'FINI') Go To 9762
if (KWord(1:4) == 'FLDG') Go To 994
if (KWord(1:4) == 'FNMC') Go To 9086
if (KWord(1:4) == 'FOOC') Go To 8000
if (KWord(1:4) == 'FPCO') Go To 9764
if (KWord(1:4) == 'FPPR') Go To 9765
if (KWord(1:4) == 'FRGM') Go To 8025
if (KWord(1:4) == 'GEO ') Go To 8024
if (KWord(1:4) == 'GEOE') Go To 8020
if (KWord(1:4) == 'GIAO') Go To 9020
if (KWord(1:4) == 'GRID') Go To 9773
if (KWord(1:4) == 'GROM') Go To 8034
if (KWord(1:4) == 'GROU') Go To 6010
if (KWord(1:4) == 'HIGH') Go To 9096
if (KWord(1:4) == 'HYPE') Go To 8016
if (KWord(1:4) == 'ISOT') Go To 7654
if (KWord(1:4) == 'JMAX') Go To 971
if (KWord(1:4) == 'KHAC') Go To 8003
if (KWord(1:4) == 'LDF ') Go To 35
if (KWord(1:4) == 'LDF1') Go To 35
if (KWord(1:4) == 'LDF2') Go To 36
if (KWord(1:4) == 'RLOC') Go To 658
if (KWord(1:4) == 'LINK') Go To 8036
if (KWord(1:4) == 'LOCA') Go To 35
if (KWord(1:4) == 'LOW ') Go To 9094
if (KWord(1:4) == 'MEDI') Go To 9095
if (KWord(1:4) == 'MGAU') Go To 8009
if (KWord(1:4) == 'MOLC') Go To 957
if (KWord(1:4) == 'MOLE') Go To 960
if (KWord(1:4) == 'MOLP') Go To 959
if (KWord(1:4) == 'MOVE') Go To 6040
if (KWord(1:4) == 'MULT') Go To 972
if (KWord(1:4) == 'NACC') Go To 8030
if (KWord(1:4) == 'NEMO') Go To 800
if (KWord(1:4) == 'NGEX') Go To 501
if (KWord(1:4) == 'NOAL') Go To 7013
if (KWord(1:4) == 'NOAM') Go To 8007
if (KWord(1:4) == 'NOCD') Go To 9084
if (KWord(1:4) == 'NOCH') Go To 9089
if (KWord(1:4) == 'NODE') Go To 7070
if (KWord(1:4) == 'NODK') Go To 989
if (KWord(1:4) == 'NOGU') Go To 9100
if (KWord(1:4) == 'NOMO') Go To 6050
if (KWord(1:4) == 'NOON') Go To 8023
if (KWord(1:4) == 'NOPA') Go To 9910
if (KWord(1:4) == 'NOTA') Go To 980
if (KWord(1:4) == 'NOUN') Go To 46
if (KWord(1:4) == 'NUME') Go To 8031
if (KWord(1:4) == 'OLDD') Go To 8121
if (KWord(1:4) == 'OLDZ') Go To 8021
if (KWord(1:4) == 'OMQI') Go To 999
if (KWord(1:4) == 'ONEO') Go To 990
if (KWord(1:4) == 'OPTH') Go To 8022
if (KWord(1:4) == 'OPTO') Go To 940
if (KWord(1:4) == 'ORBA') Go To 913
if (KWord(1:4) == 'ORBC') Go To 906
if (KWord(1:4) == 'ORIG') Go To 8015
if (KWord(1:4) == 'OVER') Go To 41
if (KWord(1:4) == 'PAMF') Go To 8060
if (KWord(1:4) == 'PART') Go To 9763
if (KWord(1:4) == 'PKTH') Go To 9940
if (KWord(1:4) == 'MXTC') Go To 9023
if (KWord(1:4) == 'PRIN') Go To 930
!if ((KWord(1:1) == 'R') .and. &
!    ((KWord(2:2) >= '0') .and. (KWord(2:2) <= '9')) .and. ((KWord(3:3) >= '0') .and. (KWord(3:3) <= '9')) .and. &
!    ((KWord(4:4) == 'O') .or. (KWord(4:4) == 'E') .or. (KWord(4:4) == 'S') .or. (KWord(4:4) == 'M') .or. &
!     (KWord(4:4) == 'C'))) Go To 657
if (KWord(1:4) == 'RA0F') Go To 9012
if (KWord(1:4) == 'RA0H') Go To 9011
if (KWord(1:4) == 'RADI') Go To 909
if (KWord(1:4) == 'RAIH') Go To 9013
if (KWord(1:4) == 'RBSS') Go To 9015
if (KWord(1:4) == 'RELA') Go To 657
if (KWord(1:4) == 'RELI') Go To 962
if (KWord(1:4) == 'RESC') Go To 978
if (KWord(1:4) == 'RF-I') Go To 9970
if (KWord(1:4) == 'RLDF') Go To 8012
if (KWord(1:4) == 'RIC ') Go To 9097
if (KWord(1:4) == 'RICD') Go To 9080
if (KWord(1:4) == 'RIJ ') Go To 9098
if (KWord(1:4) == 'RIJK') Go To 9099
if (KWord(1:4) == 'RMAT') Go To 880
if (KWord(1:4) == 'RMBP') Go To 886
if (KWord(1:4) == 'RMDI') Go To 884
if (KWord(1:4) == 'RMEA') Go To 881
if (KWord(1:4) == 'RMER') Go To 882
if (KWord(1:4) == 'RMEQ') Go To 885
if (KWord(1:4) == 'RMQC') Go To 883
if (KWord(1:4) == 'ROT ') Go To 8027
if (KWord(1:4) == 'RP-C') Go To 9093
if (KWord(1:4) == 'RPQM') Go To 8008
if (KWord(1:4) == 'RTRN') Go To 950
if (KWord(1:4) == 'RX2C') Go To 9014
if (KWord(1:4) == 'SADD') Go To 9081
if (KWord(1:4) == 'SCAL') Go To 8018
if (KWord(1:4) == 'SHAK') Go To 8050
if (KWord(1:4) == 'SDEL') Go To 7071
if (KWord(1:4) == 'SDIP') Go To 992
if (KWord(1:4) == 'SHAC') Go To 8002
if (KWord(1:4) == 'SKIP') Go To 9950
if (KWord(1:4) == 'SLIM') Go To 8005
if (KWord(1:4) == 'SPAN') Go To 890
if (KWord(1:4) == 'SPRE') Go To 889
if (KWord(1:4) == 'STDO') Go To 9930
if (KWord(1:4) == 'SYMM') Go To 900
if (KWord(1:4) == 'SYMT') Go To 6060
if (KWord(1:4) == 'TARG') Go To 37
if (KWord(1:4) == 'TDEL') Go To 7072
if (KWord(1:4) == 'TEST') Go To 991
if (KWord(1:4) == 'THRC') Go To 9021
if (KWord(1:4) == 'THRE') Go To 941
if (KWord(1:4) == 'THRL') Go To 37
if (KWord(1:4) == 'THRS') Go To 907
if (KWord(1:4) == 'TINK') Go To 8014
if (KWord(1:4) == 'TITL') Go To 910
if (KWord(1:4) == 'TRAN') Go To 8026
if (KWord(1:4) == 'UNCO') Go To 43
if (KWord(1:4) == 'UNIQ') Go To 45
if (KWord(1:4) == 'UNNO') Go To 908
if (KWord(1:4) == 'UPON') Go To 1001
if (KWord(1:4) == 'VECT') Go To 905
if (KWord(1:4) == 'VART') Go To 8032
if (KWord(1:4) == 'VARR') Go To 8033
if (KWord(1:4) == 'VERB') Go To 9122
if (KWord(1:4) == 'VERI') Go To 40
if (KWord(1:4) == 'WEIG') Go To 7014
if (KWord(1:4) == 'WELL') Go To 986
if (KWord(1:4) == 'WRUC') Go To 44
if (KWord(1:4) == 'XBAS') Go To 1924
if (KWord(1:4) == 'XFIE') Go To 975
if (KWord(1:4) == 'XRIC') Go To 9085
if (KWord(1:4) == 'XYZ ') Go To 1917
if (KWord(1:4) == 'ZCON') Go To 8017
if (KWord(1:4) == 'ZMAT') Go To 1920
if (KWord(1:4) == 'ZONL') Go To 8028

if (KWord(1:4) == 'END ') Go To 997
if (lTtl) Go To 911

if (Basis_test) then

  ! So the Basis keyword was in the native format.
  ! We have to back step until we find the command line!

  backspace(LuRd)
  backspace(LuRd)
  read(LuRd,'(A)') Key
  call UpCase(Key)
  do while (index(Key(1:4),'BASI') == 0)
    backspace(LuRd)
    backspace(LuRd)
    read(LuRd,'(A)') Key
    call UpCase(Key)
  end do
  Basis_test = .false.
  nDone = 0
  Go To 9201
end if
iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(LuWr,*)
call WarningMessage(2,KWord(1:Last)//' is not a keyword!, Error in keyword.')
call Quit_OnUserError()

!977 call WarningMessage(2,' Premature end of input file.')
!call Quit_OnUserError()

#ifdef _FDE_
!                                                                      *
!***** EMBE ************************************************************
!                                                                      *
! Read in the name of a file for an embedding potential on a grid

666 embPot = .true.
call Put_iScalar('embpot',1)
667 KWord = Get_Ln(LuRd)
if (KWord(1:4) == 'ENDE') then
  ! Sanity check
  if (embPotInBasis .and. (.not. outGridPathGiven) .and. (embWriteDens .or. embWriteESP .or. embWriteGrad .or. embWriteHess)) then
    call WarningMessage(2,' No grid to write output in embedding.')
    call Quit_OnUserError()
  end if
  ! Write the embpot runfile
  call EmbPotWrRun()
  Go To 998
end if
! Switch whether the potential is given in basis set
! representation instead of grid representation
if (KWord(1:4) == 'BASI') then
  embPotInBasis = .true.
  write(LuWr,*) 'Set embPotInBasis to ',embPotInBasis
  Go To 667
end if
! Get the EMBInfile path containing an embedding pot. on a grid
if (KWord(1:4) == 'EMBI') then
  KWord = Get_Ln(LuRd)
  call Get_S(1,embPotPath,1)
  Go To 667
end if
! If the output grid is different from the input grid a path to
! a grid file is specified here
if (KWord(1:4) == 'OUTG') then
  KWord = Get_Ln(LuRd)
  call Get_S(1,outGridPath,1)
  outGridPathGiven = .true.
  Go To 667
end if
! If given, the final density is written out on a grid, the path
! to the file where to write it to must be given as well.
if (KWord(1:4) == 'WRDE') then
  embWriteDens = .true.
  KWord = Get_Ln(LuRd)
  call Get_S(1,embOutDensPath,1)
  Go To 667
end if
! If given, the final electrostatic potential is written out on
! a grid, the path to the file where to write it to must be given
! as well.
if (KWord(1:4) == 'WREP') then
  KWord = Get_Ln(LuRd)
  call Get_S(1,embOutEspPath,1)
  embWriteESP = .true.
  Go To 667
end if
! If given, the density gradient is written out on a grid, the
! path to the file where to write it to must be given as well.
if (KWord(1:4) == 'WRGR') then
  KWord = Get_Ln(LuRd)
  call Get_S(1,embOutGradPath,1)
  embWriteGrad = .true.
  Go To 667
end if
! If given, the Hessian of the density is written out on a grid,
! the path to the file where to write it to must be given as
! well.
if (KWord(1:4) == 'WRHE') then
  embWriteHess = .true.
  KWord = Get_Ln(LuRd)
  call Get_S(1,embOutHessPath,1)
  Go To 667
end if
! Should not be reached.
call WarningMessage(2,'Error in input of EMBEdding block!')
call Quit_OnUserError()
#endif
!                                                                      *
!***** SYMM ************************************************************
!                                                                      *
! Read distinct symmetry operators apart for the unity operator

900 SymmSet = .true.
KWord = Get_Ln(LuRd)
GWInput = .true.
if ((.not. DoneCoord) .and. (iCoord /= 0)) then
  call WarningMessage(2,'SYMMETRY keyword is not compatible with COORD')
  call Quit_OnUserError()
end if
if (Run_Mode == S_Mode) then
  call WarningMessage(2,'Seward input error!')
  write(LuWr,'(A,A,A)') 'The command : "',Previous_Command,'" is not allowed in Seward input when Gateway is used!'
  write(LuWr,*)
  call Quit_OnUserError()
end if
call UpCase(KWord)
iChrct = len(KWord)
901 Last = iCLast(KWord,iChrct)
iFrst = iCFrst(KWord,iChrct)
if (iFrst <= Last) then
  nOper = nOper+1
  Oper(nOper)(1:1) = KWord(iFrst:iFrst)
  KWord(iFrst:iFrst) = ' '
  iFrst = iFrst+1
  if (KWord(iFrst:iFrst) == ' ') Go To 901
  Oper(nOper)(2:2) = KWord(iFrst:iFrst)
  KWord(iFrst:iFrst) = ' '
  iFrst = iFrst+1
  if (KWord(iFrst:iFrst) == ' ') Go To 901
  Oper(nOper)(3:3) = KWord(iFrst:iFrst)
  KWord(iFrst:iFrst) = ' '
  Go to 901
else
  Go To 998
end if
!                                                                      *
!***** FILE ************************************************************
!                                                                      *
! Specify filename for input orbitals

904 Line = Get_Ln(LuRd)
call FileOrb(Line,SW_FileOrb)
Go To 998
!                                                                      *
!***** VECT ************************************************************
!                                                                      *
! Change mode of Seward to property calculations.

905 Prprt = .true.
Go To 998
!                                                                      *
!***** ORBC ************************************************************
!                                                                      *
! Request property output with explicit listing of orbital contributions.

906 Short = .false.
Go To 998
!                                                                      *
!***** THRS ************************************************************
!                                                                      *
! Change default for non-zero occupation numbers

907 KWord = Get_Ln(LuRd)
call Get_F1(1,Thrs)
Thrs = abs(Thrs)
Go To 998
!                                                                      *
!***** UNNO ************************************************************
!                                                                      *
! Change default to unnormalized integrals

908 UnNorm = .true.
Go To 998
!                                                                      *
!***** RADI ************************************************************
!                                                                      *
! Integral cutoff to be based on radial overlap integrals

909 lSchw = .false.
Go To 998
!                                                                      *
!***** TITL ************************************************************
!                                                                      *
! Read the Title card

910 Key = Get_Ln(LuRd)
911 continue
GWInput = Run_Mode == G_Mode
lTtl = .true.
nTtl = nTtl+1
if (nTtl > 10) then
  call WarningMessage(2,' Too many title cards')
  call Quit_OnUserError()
end if
i1 = iCFrst(Key,80)
i2 = iCLast(Key,80)
nc = 80-(i2-i1+1)
nc2 = nc/2
Title(nTtl) = ''
Title(nTtl)(nc2+1:nc2+i2-i1+1) = Key(i1:i2)
Go To 9988
!                                                                      *
!***** ECPS **** or ****** AUXS ****************************************
!                                                                      *
! Allow printing of ECP data

912 nPrint(2) = max(10,nPrint(2))
GWInput = .true.
Go To 998
!                                                                      *
!***** BSSH ************************************************************
!                                                                      *
! Allow printing of basis set data

9121 nPrint(2) = max(6,nPrint(2))
GWInput = .true.
Go To 998
!                                                                      *
!***** VERB ************************************************************
!                                                                      *
! Verbose printing

9122 nPrint(2) = max(10,nPrint(2))
nPrint(117) = 6
nPrint(80) = 6
nPrint(1) = 6
GWInput = Run_Mode == G_Mode
Go To 998
!                                                                      *
!***** ORBA ************************************************************
!                                                                      *
! Request property output with explicit listing of properties of
! all orbitals, including all unoccupied (ignoring THRS), and not
! weighted by occupation numbers. (S.S.Dong, 2018)

913 ifallorb = .true.
Go To 998
!                                                                      *
!***** ZMAT ************************************************************
!                                                                      *
! Read Basis Sets & Coordinates in Z-Matrix format

1920 continue
if (isxbas == 0) call Quit_OnUserError()
call ZMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
if (iErr /= 0) call Quit_OnUserError()
GWInput = .true.
call StdSewInput(LuRd,ifnr,mdc,iShll,BasisTypes,STDINP,lSTDINP,iErr)
if (iErr /= 0) call Quit_OnUserError()
Go To 998
!                                                                      *
!***** XBAS ************************************************************
!                                                                      *
1924 continue
call read_xbas(LuRd,iglobal,nxbas,xb_label,xb_bas,ierr)
GWInput = .true.
isxbas = 1
if (ierr == 1) call Quit_OnUserError()
goto 998
!                                                                      *
!***** XYZ  ************************************************************
!                                                                      *
1917 continue
if (isxbas == 0) call Quit_OnUserError()
call XMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
if (iErr /= 0) call Quit_OnUserError()
GWInput = .true.
call StdSewInput(LuRd,ifnr,mdc,iShll,BasisTypes,STDINP,lSTDINP,iErr)
if (iErr /= 0) call Quit_OnUserError()
!if (SymmSet) then
!  call WarningMessage(2,'SYMMETRY keyword is not compatible with XYZ')
!  call Quit_OnUserError()
!end if
if (GroupSet) then
  call WarningMessage(2,'GROUP keyword is not compatible with XYZ')
  call Quit_OnUserError()
end if
Go To 998

!                                                                      *
!***** COOR ************************************************************
!                                                                      *
! Read Basis Sets & Coordinates in xyz format

6000 continue
if (SymmSet) then
  call WarningMessage(2,'SYMMETRY keyword is not compatible with COORD')
  call Quit_OnUserError()
end if
if (iOpt_XYZ == 0) then
  call WarningMessage(2,'COORD keyword is not compatible with non-XYZ format')
  call Quit_OnUserError()
end if
!if (CoordSet) then
!  call WarningMessage(2,'There is more than one COORD keyword!')
!  ! Should this be an error and stop? Only if no HYPER or GEO?
!end if
if (isxbas == 1) call Quit_OnUserError()
iCoord = iCoord+1

CoordSet = .true.
GWInput = .true.
#ifdef _HAVE_EXTRA_
call XYZread(LuRd,ForceZMAT,nCoord,iErr)
if (iErr /= 0) call Quit_OnUserError()
call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
call Read_XYZ(LuRd,OrigRot,OrigTrans)
#endif
Go To 998
!                                                                      *
!***** GROUP ***********************************************************
!                                                                      *
! Read information for a group

6010 continue
if (SymmSet) then
  call WarningMessage(2,'SYMMETRY keyword is not compatible with GROUP')
  call Quit_OnUserError()
end if
if (.not. CoordSet) then
  call WarningMessage(2,'COORD keyword is not found')
  call Quit_OnUserError()
end if
KeepGroup = Get_Ln(LuRd)
! Simplistic validity check for value
temp1 = KeepGroup
call UpCase(temp1)
if (temp1(1:4) == 'FULL') goto 6015
if (temp1(1:1) == 'E') goto 6015
if (temp1(1:2) == 'C1') goto 6015
if (temp1(1:5) == 'NOSYM') goto 6015
do i=1,StrnLn(temp1)
  if ((temp1(i:i) /= 'X') .and. (temp1(i:i) /= 'Y') .and. (temp1(i:i) /= 'Z') .and. (temp1(i:i) /= ' ')) then
    call WarningMessage(2,'Illegal symmetry group or operator: '//temp1(:StrnLn(temp1)))
    call Quit_OnUserError()
  end if
end do
6015 continue
GroupSet = .true.
GWInput = .true.
DoneCoord = .true.
goto 998
!                                                                      *
!***** BSSE ************************************************************
!                                                                      *
! Read information for BSSE

6020 continue
if (.not. CoordSet) then
  call WarningMessage(2,'COORD keyword is not found')
  call Quit_OnUserError()
end if
KWord = Get_Ln(LuRd)
read(Kword,*,end=6666,err=6666) iBSSE
GWInput = .true.
goto 998
!                                                                      *
!***** MOVE ************************************************************
!                                                                      *
! allow to MOVE coordinates

6040 continue
if (.not. CoordSet) then
  call WarningMessage(2,'COORD keyword is not found')
  call Quit_OnUserError()
end if
#ifdef _HAVE_EXTRA_
isHold = 0
#endif
GWInput = .true.
goto 998
!                                                                      *
!***** NOMOVE **********************************************************
!                                                                      *
! Do NOT allow to MOVE coordinates

6050 continue
if (.not. CoordSet) then
  call WarningMessage(2,'COORD keyword is not found')
  call Quit_OnUserError()
end if
#ifdef _HAVE_EXTRA_
isHold = 1
#endif
GWInput = .true.
goto 998
!                                                                      *
!***** SYMT ************************************************************
!                                                                      *
! Threshold for findsym

6060 continue
if (.not. CoordSet) then
  call WarningMessage(2,'COORD keyword is not found')
  call Quit_OnUserError()
end if
KWord = Get_Ln(LuRd)
read(Kword,*,end=6666,err=6666) SymThr
GWInput = .true.
goto 998
!                                                                      *
!***** NODE ************************************************************
!                                                                      *
! Set global delete parameters to disable orbital deleting
! Keywords: NODElete at 7070
!           SDELete  at 7071
!           TDELete  at 7072

7070 continue
call Put_dScalar('S delete thr',0.0d0)
call Put_dScalar('T delete thr',1.0d15)
goto 998
7071 continue
KWord = Get_Ln(LuRd)
call Get_F1(1,sDel)
call Put_dScalar('S delete thr',sDel)
goto 998
7072 continue
KWord = Get_Ln(LuRd)
call Get_F1(1,tDel)
call Put_dScalar('T delete thr',tDel)
goto 998
7700 continue
KWord = Get_Ln(LuRd)
ExtBasDir = KWord
GWInput = .true.
goto 998
!                                                                      *
!***** BASI ************************************************************
!                                                                      *
! Read information for a basis set

920 continue

! Check if the format is old or new style. Damn the person who used
! the same keword for two different styles of input and making the
! input require a specific order of the keyword. Comrade 55?

Basis_Test = .true.

GWInput = .true.
Key = Get_Ln(LuRd)
BSLbl = Key(1:80)
if (BasisSet) then
  KeepBasis = KeepBasis(1:index(KeepBasis,' '))//','//BSLbl
else
  KeepBasis = BSLbl
  BasisSet = .true.
end if
temp1 = KeepBasis
call UpCase(temp1)
!if (INDEX(temp1,'INLINE') /= 0) then
!  write(LuWr,*) 'XYZ input and Inline basis set are not compatible'
!  write(LuWr,*) 'Consult the manual how to change inline basis set'
!  write(LuWr,*) ' into basis set library'
!  call Quit_OnUserError()
!end if
iOpt_XYZ = 1
goto 998

9201 continue
iOpt_XYZ = 0
GWInput = .true.
nCnttp = nCnttp+1
if (Run_Mode == S_Mode) then
  call WarningMessage(2,'Seward input error!')
  write(LuWr,*) 'The command : "',Previous_Command,'" is not allowed in the Seward input when the Gateway is used!'
  write(LuWr,*)
  call Quit_OnUserError()
end if
if (nCnttp > Mxdbsc) then
  call WarningMessage(2,' Increase Mxdbsc')
  call Quit_OnUserError()
end if
NoZMAT = .true.

! Read the basis set label

Key = Get_Ln(LuRd)
BSLbl = Key(1:80)

! If dummy atom point at the ANO-RCC set where the specification is.

call UpCase(BSLbl)
iDummy_basis = 0
call ICopy(4,BasisTypes,1,BasisTypes_save,1)
if ((BSLbl(1:2) == 'X.') .and. (index(BSLbl,'INLINE') == 0) .and. (index(BSLbl,'RYDBERG') == 0)) then
  BSLbl = 'X.ANO-RCC.'
  do i=11,80
    BSLbl(i:i) = '.'
  end do
  iDummy_basis = 1
end if
!call UpCase(BSLbl)
LenBSL = len(BSLbl)
Last = iCLast(BSLbl,LenBSL)
GWInput = .true.
Indx = index(BSLbl,'/')
if (Indx == 0) then
  Fname = BasLib
  Indx = Last+1
  dbsc(nCnttp)%Bsl = BSLbl
else
  Fname = BSLbl(Indx+2:Last)
  if (Fname == ' ') then
    call WarningMessage(2,' No basis set library specified for BSLbl='//BSLbl(1:Indx-1)//' Fname='//Fname)
    call Quit_OnUserError()
  end if
1919 if (Fname(1:1) == ' ') then
    Fname(1:79) = Fname(2:80)
    Fname(80:80) = ' '
    Go To 1919
  end if
  dbsc(nCnttp)%Bsl = BSLbl(1:Indx-1)
end if

n = index(dbsc(nCnttp)%Bsl,' ')
if (n == 0) n = 81
do i=n,80
  dbsc(nCnttp)%Bsl(i:i) = '.'
end do

if ((Show .and. (nPrint(2) >= 6)) .or. Write_BasLib) then
  write(LuWr,*)
  write(LuWr,*)
  write(LuWr,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLbl(1:Indx-1)
  write(LuWr,'(1X,A,A)') 'Basis set is read from library:',Fname(1:index(Fname,' '))
end if

jShll = iShll
dbsc(nCnttp)%Bsl_old = dbsc(nCnttp)%Bsl
dbsc(nCnttp)%mdci = mdc
call GetBS(Fname,dbsc(nCnttp)%Bsl,iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,lSTDINP,.false.,Expert,ExtBasDir)

Do_FckInt = Do_FckInt .and. dbsc(nCnttp)%FOp .and. (dbsc(nCnttp)%AtmNr <= 96)
#ifdef _DEMO_
Do_GuessOrb = .false.
#else
Do_GuessOrb = Do_GuessOrb .and. (dbsc(nCnttp)%AtmNr <= 96)
#endif

if (iDummy_Basis == 1) call ICopy(4,BasisTypes_Save,1,BasisTypes,1)
if (itype == 0) then
  if ((BasisTypes(3) == 1) .or. (BasisTypes(3) == 2) .or. (BasisTypes(3) == 14)) iType = BasisTypes(3)
else
  if ((BasisTypes(3) == 1) .or. (BasisTypes(3) == 2) .or. (BasisTypes(3) == 14)) then
    if (BasisTypes(3) /= iType) then
      imix = 1
      BasisTypes(3) = -1
    end if
    iType = BasisTypes(3)
  end if
end if
if (itype == 1) ifnr = 1
if ((itype == 2) .or. (itype == 14)) ifnr = 0

if (ign == 0) then
  ign = BasisTypes(4)
else if (abs(BasisTypes(4)) /= abs(ign)) then
  if (nmwarn) then
    call WarningMessage(1,'SEWARD found basis sets of mixed nuclear charge model. The most advanced one will be used.')
  end if
  nmwarn = .false.
  ign = max(ign,BasisTypes(4))
  BasisTypes(4) = ign
end if

if (Show .and. (nPrint(2) >= 6) .and. (Ref(1) /= '') .and. (Ref(2) /= '')) then
  write(LuWr,'(1x,a)') 'Basis Set Reference(s):'
  if (Ref(1) /= '') write(LuWr,'(5x,a)') trim(Ref(1))
  if (Ref(2) /= '') write(LuWr,'(5x,a)') trim(Ref(2))
  write(LuWr,*)
  write(LuWr,*)
end if
dbsc(nCnttp)%ECP = (dbsc(nCnttp)%nPP+dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nM1+dbsc(nCnttp)%nM2) /= 0

lAng = max(dbsc(nCnttp)%nVal,dbsc(nCnttp)%nSRO,dbsc(nCnttp)%nPrj)-1
S%iAngMx = max(S%iAngMx,lAng)
! No transformation needed for s and p shells
Shells(jShll+1)%Transf = .false.
Shells(jShll+1)%Prjct = .false.
Shells(jShll+2)%Transf = .false.
Shells(jShll+2)%Prjct = .false.
dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal+dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nPP
nCnt = 0
if (dbsc(nCnttp)%Aux) then
  do iSh=jShll+1,iShll
    Shells(iSh)%Aux = .true.
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Set Cartesian functions if specified by the basis type

if (BasisTypes(1) == 9) then
  do iSh=jShll+3,iShll
    Shells(iSh)%Transf = .false.
    Shells(iSh)%Prjct = .false.
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Automatic onset of muonic charge if the basis type is muonic.
! This will also automatically activate finite nuclear mass correction.

KWord = ''
KWord(1:Indx-1) = BSLbl(1:Indx-1)
call UpCase(KWord)
if (index(KWord,'MUONIC') /= 0) then
  dbsc(nCnttp)%fMass = CONST_MUON_MASS_IN_SI_/CONST_ELECTRON_MASS_IN_SI_
  FNMC = .true.
  tDel = 1.0d50
  call Put_dScalar('T delete thr',tDel)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Update BasisTypes

do i=1,4
  if (BasisTypes_save(i) == 0) cycle
  if (BasisTypes(i) /= BasisTypes_save(i)) BasisTypes(i) = -1
end do
!                                                                      *
!***********************************************************************
!                                                                      *

777 KWord = Get_Ln(LuRd)
call UpCase(KWord)
call LeftAd(KWord)
if (KWord(1:4) == 'PSEU') then
  dbsc(nCnttp)%pChrg = .true.
  dbsc(nCnttp)%Fixed = .true.
  Go To 777
end if
if (KWord(1:4) == 'ACDT') then
  KWord = Get_Ln(LuRd)
  call Get_F1(1,dbsc(nCnttp)%aCD_Thr)
  Go To 777
end if
if (KWord(1:4) == 'MUON') then
  dbsc(nCnttp)%fMass = CONST_MUON_MASS_IN_SI_/CONST_ELECTRON_MASS_IN_SI_
  Go To 777
end if
if (KWord(1:4) == 'NUCL') then
  KWord = Get_Ln(LuRd)
  call Get_F1(1,dbsc(nCnttp)%ExpNuc)
  Go To 777
end if
if (KWord(1:4) == 'FIXE') then
  dbsc(nCnttp)%Fixed = .true.
  Go To 777
end if
if (KWord(1:4) == 'SPHE') then
  if (index(KWord,'ALL') /= 0) then
    do iSh=jShll+3,iShll
      Shells(iSh)%Transf = .true.
      Shells(iSh)%Prjct = .true.
    end do
    Go To 777
  end if
  ist = index(KWord,' ')
  iAng = 2
  do iSh=jShll+3,iShll
    if (index(KWord(ist:80),AngTyp(iAng)) /= 0) then
      Shells(iSh)%Transf = .true.
      Shells(iSh)%Prjct = .true.
    end if
    iAng = iAng+1
  end do
  Go To 777
end if
if (KWord(1:4) == 'CART') then
  if (index(KWord,'ALL') /= 0) then
    do iSh=jShll+1,iShll
      Shells(iSh)%Transf = .false.
      Shells(iSh)%Prjct = .false.
    end do
    Go To 777
  end if
  ist = index(KWord,' ')
  iAng = 0
  do iSh=jShll+1,iShll
    if (index(KWord(ist:80),AngTyp(iAng)) /= 0) then
      Shells(iSh)%Transf = .false.
      Shells(iSh)%Prjct = .false.
    end if
    iAng = iAng+1
  end do
  Go To 777
end if
if (KWord(1:4) == 'CONT') then
  if (index(KWord,'ALL') /= 0) then
    do iSh=jShll+1,iShll
      Shells(iSh)%Prjct = .false.
    end do
    Go To 777
  end if
  ist = index(KWord,' ')
  iAng = 0
  do iSh=jShll+1,iShll
    if (index(KWord(ist:80),AngTyp(iAng)) /= 0) Shells(iSh)%Prjct = .false.
    iAng = iAng+1
  end do
  Go To 777
end if
if (KWord(1:4) == 'CHAR') then
  KWord = Get_Ln(LuRd)
  call UpCase(KWord)
  call Get_F1(1,dbsc(nCnttp)%Charge)
  ist = index(KWord,' ')
  if (dbsc(nCnttp)%IsMM /= 0) then
    call WarningMessage(1,' Found a charge associated with a MM atom. Ignore it')
    dbsc(nCnttp)%Charge = Zero
  end if
  Go To 777
end if
if (KWord(1:4) == 'FRAG') then
  dbsc(nCnttp)%pChrg = .true.
  dbsc(nCnttp)%Fixed = .true.
  lFAIEMP = .true.
  Go To 777
end if
if (KWord(1:4) == 'END ') then
  if (nCnt == 0) then
    call WarningMessage(2,' Input error, no center specified!')
    call Quit_OnUserError()
  end if
  dbsc(nCnttp)%nCntr = nCnt
  mdc = mdc+nCnt
  ! Now allocate the array for the coordinates and copy them over.
  ! Call Allocate(dbsc(nCnttp)%Coor(1:3,1:nCnt)
  call mma_Allocate(dbsc(nCnttp)%Coor_Hidden,3,nCnt,Label='dbsc:C')
  dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
  call DCopy_(3*nCnt,Buffer,1,dbsc(nCnttp)%Coor,1)

  Go To 998
end if

! Read Coordinates

nCnt = nCnt+1
n_dc = max(mdc+nCnt,n_dc)
if (mdc+nCnt > MxAtom) then
  call WarningMessage(2,' RdCtl: Increase MxAtom')
  write(LuWr,*) '        MxAtom=',MxAtom
  call Quit_OnUserError()
end if
jend = index(KWord,' ')
if (jEnd > LENIN+1) then
  write(6,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LENIN,' characters!'
end if
dc(mdc+nCnt)%LblCnt = KWord(1:min(LENIN,jend-1))
dbas = dc(mdc+nCnt)%LblCnt(1:LENIN)
call Upcase(dbas)
if (dbas == 'DBAS') then
  RMat_On = .true.
end if
if (mdc+nCnt > 1) then
  call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)
end if
iOff = 1+(nCnt-1)*3
call Get_F(2,Buffer(iOff),3)
if (index(KWord,'ANGSTROM') /= 0) then
  do i=0,2
    Buffer(iOff+i) = Buffer(iOff+i)/angstr
  end do
end if

if (Cell_l) then
  nCnt0 = nCnt
  iOff0 = iOff
  lthCell = lthCell+1
  AdCell(lthCell) = mdc+nCnt0  ! the sequence atom No
  ii = 0
  do n1=-ispread(1),ispread(1)
    do n2=-ispread(2),ispread(2)
      do n3=-ispread(3),ispread(3)
        if ((n1 == 0) .and. (n2 == 0) .and. (n3 == 0)) goto 110

        ii = ii+1

        if (ii >= 10000) then
          call WarningMessage(2,' Too many atoms in Seward')
          call Quit_OnUserError()
        else
          if (ii < 1000) then
            CHAR4 = '_'//str(ii)
          else
            CHAR4 = str(ii)
          end if
        end if

        nCnt = nCnt+1

        n_dc = max(mdc+nCnt,n_dc)
        if (mdc+nCnt > MxAtom) then
          call WarningMessage(2,' RdCtl: Increase MxAtom')
          write(LuWr,*) '        MxAtom=',MxAtom
          call Quit_OnUserError()
        end if

        jend = index(KWord,' ')
        if (jEnd > 5) then
          write(6,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LENIN,' characters!'
        end if
        dc(mdc+nCnt)%LblCnt = KWord(1:min(LENIN,jend-1))//CHAR4

        call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)

        iOff = 1+(nCnt-1)*3

        ! Copy old coordinate  first
        call DCOPY_(3,Buffer(iOff0),1,Buffer(iOff),1)
        call DAXPY_(3,dble(n1),VCell(1,1),1,Buffer(iOff),1)
        call DAXPY_(3,dble(n2),VCell(1,2),1,Buffer(iOff),1)
        call DAXPY_(3,dble(n3),VCell(1,3),1,Buffer(iOff),1)

110     continue

      end do
    end do
  end do
end if
Go To 777
!                                                                      *
!***** PRIN ************************************************************
!                                                                      *
! Print level

930 KWord = Get_Ln(LuRd)
GWInput = Run_Mode == G_Mode
call Get_I1(1,n)
do i=1,n
  KWord = Get_Ln(LuRd)
  call Get_I1(1,jRout)
  call Get_I1(2,iPrint)
  nPrint(jRout) = iPrint
end do
Go To 998
!                                                                      *
!***** OPTO ************************************************************
!                                                                      *
! Reduce the output for optimizations

940 lOPTO = .true.
Go To 998
!                                                                      *
!***** THRE ************************************************************
!                                                                      *
! Threshold for writing integrals to disk

941 KWord = Get_Ln(LuRd)
call Get_F1(1,ThrInt)
ThrInt = abs(ThrInt)
ThrInt_UsrDef = .true.
Go To 998
!                                                                      *
!***** CUTO ************************************************************
!                                                                      *
! Cutoff for computing primitive integrals [a0|c0]

942 KWord = Get_Ln(LuRd)
call Get_F1(1,CutInt)
CutInt = abs(CutInt)
CutInt_UsrDef = .true.
Go To 998
!                                                                      *
!***** RTRN ************************************************************
!                                                                      *
! Defining max bond distance for bonds, angles and dihedrals
! Define max number of atoms to list for

950 KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_I1(1,S%Max_Center)
call Get_F1(2,rtrnc)
if (index(KWord,'ANGSTROM') /= 0) Rtrnc = Rtrnc/angstr
GWInput = .true.
Go To 998
!                                                                      *
!***** DIRE ************************************************************
!                                                                      *
! Force direct calculations & disable two-electron integrals

9770 DirInt = .true.
Onenly = .true.
iChk_DC = 1
if ((iChk_RI+iChk_CH) > 0) then
  call WarningMessage(2,'Direct is incompatible with RI and Cholesky keywords')
  call Quit_OnUserError()
end if
Do_RI = .false.
Go To 998
!                                                                      *
!***** CSPF ************************************************************
!                                                                      *
! Turn on the use of Condon-Shortley phase factors

9110 CSPF = .true.
GWInput = Run_Mode == G_Mode
Go To 998
!                                                                      *
!***** EXPE ************************************************************
!                                                                      *
! Expert mode

9771 Expert = .true.
GWInput = Run_Mode == G_Mode
call WarningMessage(1,' EXPERT option is ON!')
Go To 998
!                                                                      *
!***** MOLC or DCRN ****************************************************
!                                                                      *
! Weight for DCR summation (this is the default for conventional
! calculations)

957 continue
958 MolWgh = 0
MolWgh_UsrDef = .true.
Go To 998
!                                                                      *
!***** MOLP ************************************************************
!                                                                      *
! Weight for DCR summation modified to MOLPRO format

959 MolWgh = 2
MolWgh_UsrDef = .true.
Go To 998
!                                                                      *
!***** MOLE ************************************************************
!                                                                      *
! Weight for DCR summation modified to MOLECULE format

960 MolWgh = 1
MolWgh_UsrDef = .true.
Go To 998
!                                                                      *
!***** RELI ************************************************************
!                                                                      *
! Compute integrals for first order relativistic corrections
! of the energy, i.e. the mass-velocity integrals and the
! one-electron Darwin contract term integrals.

962 lRel = .true.
Go To 998
!                                                                      *
!***** JMAX ************************************************************
!                                                                      *
! Change max j quantum number for the rigid rotor analysis

971 KWord = Get_Ln(LuRd)
call Get_I1(1,S%jMax)
Go To 998
!                                                                      *
!***** MULT ************************************************************
!                                                                      *
! Read order of highest multipole to be computed

972 KWord = Get_Ln(LuRd)
call Get_I1(1,S%nMltpl)
Go To 998
!                                                                      *
!***** CENT ************************************************************
!                                                                      *
! User specified centers of multipole moment operators.

973 KWord = Get_Ln(LuRd)
call Get_I1(1,nTemp)
if (lMltpl) then
  call WarningMessage(2,' Abend: User specified centers already defined; Abend: Correct input error!')
  call Quit_OnUserError()
else
  lMltpl = .true.
end if
! Allocate temporary memory for user defined centers
call mma_allocate(RTmp,3,nTemp,label='RTmp')
call mma_allocate(ITmp,nTemp,label='ITmp')
do i=1,nTemp
  KWord = Get_Ln(LuRd)
  call Upcase(KWord)
  call Get_I1(1,iMltpl)
  call Get_F(2,RTmp(1,i),3)
  if (index(KWord,'ANGSTROM') /= 0) call DScal_(3,One/angstr,RTmp(1,i),1)
  ITmp(i) = iMltpl
end do
Go To 998
!                                                                      *
!***** EMPC ************************************************************
!                                                                      *
! Compute Orbital-Free Embedding integrals from Point Charges
! specified in XFIEld

974 DoEmPC = .true.
#ifdef _HAVE_EXTRA_
isHold = 1 ! avoid coordinate moving
#endif
GWInput = .true.
Go To 998
!                                                                      *
!***** XFIE ************************************************************
!                                                                      *
! User specified external field

975 lXF = .true.
GWInput = .true.
KWord = Get_Ln(LuRd)
! Open external file if the line does not start with an integer
! Note that the "Err" signal cannot be completely trusted, since
! a slash (i.e., an absolute path) marks end of input and gives
! no error
LuRd_saved = LuRd
ibla = -1
read(KWord,*,Err=9751) ibla
if (ibla < 0) goto 9751
goto 9752
9751 LuRd = 1
call Get_S(1,filename,1)
9753 call OpnFl(filename(1:(index(filename,' ')-1)),LuRd,Exist)
if (.not. Exist) then
  call WarningMessage(2,'Error! File not found: '//filename(1:(index(filename,' ')-1)))
  call Quit_OnUserError()
end if
write(LuWr,*) 'Reading external field from file: ',filename(1:(index(filename,' ')-1))
KWord = Get_Ln(LuRd)
9752 call Get_I1(1,nXF)
Convert = .false.
call Upcase(KWord)
if (index(KWord,'ANGSTROM') /= 0) then
  Convert = .true.
  ix = index(KWord,'ANGSTROM')
  KWord(ix:ix+7) = '        '
end if

KWord(170:180) = '-2 -2 -2 -2'
call Put_Ln(KWord)
call Get_I1(2,nOrd_XF)
call Get_I1(3,iXPolType)
call Get_I1(4,nXMolnr)
call Get_I1(5,nReadEle)

! Set defaults: ch+dip, no polarisabilities,
!               exclude only its own multipole,
!               no element read

if (nOrd_XF == -2) nOrd_XF = 1
if (iXPolType == -2) iXPolType = 0
if (nXMolnr == -2) nXMolnr = 0
if (nReadEle == -2) nReadEle = 0

if ((nOrd_XF > 2) .or. (nOrd_XF < -1)) then
  call WarningMessage(2,'Error! Illegal value of nOrd_XF')
  write(LuWr,*) 'nOrd_XF= ',nOrd_XF
  call Quit_OnUserError()
end if
if ((iXPolType > 2) .or. (iXPolType < 0)) then
  call WarningMessage(2,'Error! Illegal value of iXPolType')
  write(LuWr,*) 'iXPolType= ',iXPolType
  call Quit_OnUserError()
end if
if ((nXMolnr > 100) .or. (nXMolnr < 0)) then
  call WarningMessage(2,'Error! Illegal value of nXMolnr')
  write(LuWr,*) 'nXMolnr= ',nXMolnr
  call Quit_OnUserError()
end if
if ((nReadEle > 1) .or. (nReadEle < 0)) then
  call WarningMessage(2,'Error! Illegal value of nReadEle')
  write(LuWr,*) 'nReadEle= ',nReadEle
  call Quit_OnUserError()
end if

nData_XF = 3
do iOrd_XF=0,nOrd_XF
  nData_XF = nData_XF+(iOrd_XF+1)*(iOrd_XF+2)/2
end do

if (iXPolType > 0) then
  nData_XF = nData_XF+6
  lRF = .true.  ! Polarisabilities treated as Langevin
end if

if (iXPolType == 1) then
  nDataRead = nData_XF-5 ! Read only one pol value
else
  nDataRead = nData_XF
end if

call mma_allocate(XF,nData_XF,nXF,Label='XF')
call mma_allocate(XMolnr,nXMolnr,nXF,Label='XMolnr')
call mma_allocate(XEle,nXF,Label='XEle')

call Upcase(KWord)

do iXF=1,nXF
  XEle(iXF) = 0  ! default: no element spec.

  ! If reading from external file, use free format to allow
  ! long lines of input. On the other hand, comments are
  ! not allowed in external files.

  if (LuRd /= LuRd_saved) then
    call mma_Allocate(iScratch,nXMolnr+nReadEle,Label='iScratch')
    read(LuRd,*) (iScratch(k),k=1,nXMolnr),(iScratch(nXMolnr+k),k=1,nReadEle),(XF(k,iXF),k=1,nDataRead)
    do i=1,nXMolnr
      XMolnr(i,iXF) = iScratch(i)
    end do
    do i=1,nReadEle
      XEle(iXF+(i-1)) = iScratch(nXMolnr+i)
    end do
    call mma_deallocate(iScratch)
  else
    KWord = Get_Ln(LuRd)
    KWord(170:180) = ' 0.0 0.0 0.0'
    call Put_Ln(KWord)

    do i=1,nXMolnr
      call Get_I1(i,iTemp)
      XMolnr(i,iXF) = iTemp
    end do
    do i=1,nReadEle
      call Get_I1(nXMolnr+i,iTemp)
      XEle(iXF+(i-1)) = iTemp
    end do
    call Get_F(nXMolnr+nReadEle+1,XF(1,iXF),nDataRead)
  end if

  XF(1:3,iXF) = XF(1:3,iXF)*ScaleFactor
  if (Convert) XF(1:3,iXF) = XF(1:3,iXF)/angstr

end do

! Close file and reset LuRd if external file was used
if (LuRd /= LuRd_saved) then
  close(LuRd)
  LuRd = LuRd_saved
end if
if (isXfield == 1) then
  goto 9755
end if
Go To 998
!                                                                      *
!***** DOUG ************************************************************
!                                                                      *
! Full Douglas-Kroll operator

976 continue
IRELAE = 0
DKroll = .true.
GWInput = .true.
go to 998

! Full Douglas-Kroll (DK1) operator

9001 continue
IRELAE = 1
DKroll = .true.
GWInput = .true.
go to 998

! Full Douglas-Kroll (DK2) operator

9002 continue
IRELAE = 2
DKroll = .true.
GWInput = .true.
go to 998

! Douglas-Kroll (DK3) operator

9003 continue
IRELAE = 3
DKroll = .true.
GWInput = .true.
go to 998

! Full Douglas-Kroll (DK3) operator

9004 continue
IRELAE = 4
DKroll = .true.
GWInput = .true.
go to 998

! Full RESC operator

978 continue
IRELAE = 11
DKroll = .true.
GWInput = .true.
go to 998

! Full ZORA

9011 continue
IRELAE = 21
DKroll = .true.
GWInput = .true.
go to 998

! Full ZORA-FP

9012 continue
IRELAE = 22
DKroll = .true.
GWInput = .true.
go to 998

! Full IORA

9013 continue
IRELAE = 23
DKroll = .true.
GWInput = .true.
go to 998

! Exact decoupling X2C method

9014 continue
IRELAE = 101
DKroll = .true.
GWInput = .true.
go to 998

! Exact decoupling BSS method

9015 continue
IRELAE = 102
DKroll = .true.
GWInput = .true.
go to 998

! Switch on the old DKH routine

8121 continue
IRFLAG1 = 1
GWInput = .true.
goto 998
!                                                                      *
!***** BSSM ************************************************************
!                                                                      *
! BSS method

9760 continue
BSS = .true.
GWInput = .true.
go to 976
!                                                                      *
!***** AMFI ************************************************************
!                                                                      *
! AMFI integrals

9761 continue
lAMFI = .true.
GWInput = .true.
Go To 998
!                                                                      *
!***** AMF1 ************************************************************
!                                                                      *
! AMFI integrals

8761 continue
lAMFI = .true.
GWInput = .true.
Go To 998
!                                                                      *
!***** AMF2 ************************************************************
!                                                                      *
! AMFI integrals (including 2nd-order)

8762 continue
lAMFI = .true.
GWInput = .true.
Go To 998
!                                                                      *
!***** AMF3 ************************************************************
!                                                                      *
! AMFI integrals (including 3rd-order)

8763 continue
lAMFI = .true.
GWInput = .true.
Go To 998
!                                                                      *
!***** FAKE ************************************************************
!                                                                      *
! Fake run : conventional, RI or CD ERIs not computed
!            but some info is set to runfile (e.g. CD thrs)
!            Not the same as ONEOnly !!

9759 continue
Fake_ERIs = .true.
Go To 998
!                                                                      *
!***** FINI ************************************************************
!                                                                      *
! Finite nuclei - Gaussian type

9762 continue
Nuclear_Model = Gaussian_Type
GWInput = .true.
Go To 998
!                                                                      *
!***** MGAU ************************************************************
!                                                                      *
! Finite nuclei - modified Gaussian type

8009 continue
Nuclear_Model = mGaussian_Type
GWInput = .true.
Go To 998
!                                                                      *
!***** PART ************************************************************
!                                                                      *
! Show partitioning statistics

9763 continue
nPrint(10) = 6
Go To 998
!                                                                      *
!***** FPCO ************************************************************
!                                                                      *
! Force partitioning for contracted functions

9764 continue
force_part_c = .true.
Go To 998
!                                                                      *
!***** FPPR ************************************************************
!                                                                      *
! Force partitioning for primitive functions

9765 continue
force_part_p = .true.
Go To 998
!                                                                      *
!***** NOTA ************************************************************
!                                                                      *
! Do not use tables for the roots and weights of the
! Rys polynomials.

980 NoTab = .true.
Go To 998
!                                                                      *
!***** WELL ************************************************************
!                                                                      *
! Read radius and exponents for spherical well integrals and coefficient

986 KWord = Get_Ln(LuRd)
GWInput = .true.
call Get_I1(1,nWel)
if (nWel <= 0) then
  ! Use automatic set up for well integrals
  nWel = 3
  call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
  do iWel=1,nWel
    Wel_Info(3,iWel) = WellCff(iWel)
    Wel_Info(2,iWel) = WellExp(iWel)
    Wel_Info(1,iWel) = WellRad(iWel)
  end do
else
  call mma_allocate(Wel_Info,3,nWel,Label='Wel_Info')
  do iWel=1,nWel
    ! Read the Coefficient, Exponent, and Radius
    KWord = Get_Ln(LuRd)
    call Get_F1(1,Wel_Info(3,iWel))
    call Get_F1(2,Wel_Info(2,iWel))
    call Get_F1(3,Wel_Info(1,iWel))
    call Upcase(KWord)
    if (index(KWord,'ANGSTROM') /= 0) then
      Wel_Info(1,iWel) = Wel_Info(1,iWel)/angstr
      Wel_Info(2,iWel) = Wel_Info(2,iWel)*angstr
    end if
  end do
end if
Go To 998
!                                                                      *
!***** NODK ************************************************************
!                                                                      *
! Do not compute Douglas-Kroll integrals.

989 NoDKroll = .true.
Go To 998
!                                                                      *
!***** ONEO ************************************************************
!                                                                      *
! Do not compute two electron integrals.

990 Onenly = .true.
Go To 998
!                                                                      *
!***** TEST ************************************************************
!                                                                      *
! Process only the input.

991 Test = .true.
GWInput = Run_Mode == G_Mode
Go To 998
!                                                                      *
!***** SDIP ************************************************************
!                                                                      *
! Compute integrals for transition dipole moment

992 Vlct_ = .true.
GWInput = .true.
Go To 998
!                                                                      *
!***** EPOT ************************************************************
!                                                                      *
! Compute the electric potential for a number of points.  If nEF is
! set to 0 this will cause the points to coincide with the unique
! centers.

9932 nOrdEF = max(nOrdEF,0)
GWInput = .true.
if (EFgiven) then
  call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
  call Quit_OnUserError()
end if
EFgiven = .true.
Go To 9931
!                                                                      *
!***** EFLD ************************************************************
!                                                                      *
! Compute the electric potential and electric field for a number of
! points. If nEF is set to 0 this will cause the points to coincide
! with the unique centers.

993 nOrdEF = max(nOrdEF,1)
GWInput = .true.
if (EFgiven) then
  call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
  call Quit_OnUserError()
end if
EFgiven = .true.
Go To 9931
!                                                                      *
!***** FLDG ************************************************************
!                                                                      *
! Compute the electric potential, electric field, and electric field
! gradient for a number of points. If nEF is set to 0 this will
! cause the points to coincide with the Unique centers.

994 nOrdEF = max(nOrdEF,2)
GWInput = .true.
if (EFgiven) then
  call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
  call Quit_OnUserError()
end if
EFgiven = .true.
Go To 9931

9931 KWord = Get_Ln(LuRd)
call Get_I1(1,nEF)
if (nEF < 0) nEF = 0
if (nEF == 0) Go To 998
call mma_allocate(EFt,3,nEF,label='nEF')
do iEF=1,nEF
  KWord = Get_Ln(LuRd)
  ! Check whether a label is specified instead of a coordinate
  call Get_S(1,Key,1)
  call LeftAd(Key)
  call Upcase(Key)
  jTmp = ichar(Key(1:1))
  if ((jTmp >= 65) .and. (jTmp <= 90)) then
    jEnd = index(Key,' ')-1
    iOff = 0
    iFound_Label = 0
    do iCnttp=1,nCnttp
      do iCnt=iOff+1,iOff+dbsc(iCnttp)%nCntr
        if (Key(1:jEnd) == dc(iCnt)%LblCnt(1:jEnd)) then
          iFound_Label = 1
          EFt(1:3,iEF) = dbsc(iCnttp)%Coor(1:3,iCnt-iOff)
        end if
      end do
      iOff = iOff+dbsc(iCnttp)%nCntr
    end do
    if (iFound_Label == 0) then
      call WarningMessage(2,'; Error in processing the keyword FLDG.; The label "'//Key(1:jEnd)// &
                          '" could not be found among the centers.; Remember to specify the atom center before specifying the '// &
                          'FLDG keyword.')
      call Quit_OnUserError()
    end if
  else
    call Get_F(1,EFt(1,iEF),3)
    call Upcase(KWord)
    if (index(KWord,'ANGSTROM') /= 0) call DScal_(3,One/angstr,EFt(1,iEF),1)
  end if
end do
Go To 998
!                                                                      *
!***** ANGM ************************************************************
!                                                                      *
! Orbital angular momentum

995 lOAM = .true.
GWInput = .true.
KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_F(1,OAMt,3)
if (index(KWord,'ANGSTROM') /= 0) call DScal_(3,One/angstr,OAMt,1)
Go To 998
!                                                                      *
!***** ANGM derivative restriction *************************************
!                                                                      *
! Orbital angular momentum restriction

1001 lUPONLY = .true.
Go To 998

1002 lDOWNONLY = .true.
Go To 998
!                                                                      *
!***** OMQI ************************************************************
!                                                                      *
! Orbital magnetic quadrupole

999 lOMQ = .true.
GWInput = .true.
KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_F(1,OMQt,3)
if (index(KWord,'ANGSTROM') /= 0) call DScal_(3,One/angstr,OMQt,1)
Go To 998
!                                                                      *
!***** AMPR ************************************************************
!                                                                      *
! Angular momentum products

9951 GWInput = .true.
if ((Run_Mode == S_Mode) .and. GWInput) Go To 9989
call mma_allocate(AMP_Center,3,Label='AMP_Center')
KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_F(1,AMP_Center,3)
if (index(KWord,'ANGSTROM') /= 0) AMP_Center(:) = (One/angstr)*AMP_Center(:)
Go To 998
!                                                                      *
!***** DSHD ************************************************************
!                                                                      *
! Compute the diamagnetic shielding for a number of points. If nDMS
! is set to 0 this will cause the points to coincide with the
! unique centers.

996 lDMS = .true.
GWInput = .true.
KWord = Get_Ln(LuRd)
call Get_F(1,Dxyz,3)
KWord = Get_Ln(LuRd)
call Get_I1(1,nDMS)
if (nDMS < 0) nDMS = 0
if (nDMS == 0) Go To 998
call mma_allocate(DMSt,3,nDMS,label='DMSt')
do iDMS=1,nDMS
  KWord = Get_Ln(LuRd)
  call Upcase(KWord)
  call Get_F(1,DMSt(1,iDMS),3)
  if (index(KWord,'ANGSTROM') /= 0) call DScal_(3,One/angstr,DMSt(1,iDMS),1)
end do
Go To 998
!                                                                      *
!***** NOPA ************************************************************
!                                                                      *
! Set integral packing flag
! Note      : this flag is only active if iWRopt=0
! iPack=0   : pack 2el integrals (= Default)
! iPack=1   : do not pack 2el integrals

9910 iPack = 1
Go To 998
!                                                                      *
!***** STDO ************************************************************
!                                                                      *
! Set integral write option for 2 el integrals
! iWRopt=0  : 2 el integrals are written in the MOLCAS2 format,
!             i.e., in canonical order, no labels and packed format
!             (= Default)
! iWRopt=1  : 2 el integrals are written in a format identical
!             to MOLECULE, i.e., values and labels

9930 iWRopt = 1
Go To 998
!                                                                      *
!***** PKTH ************************************************************
!                                                                      *
! Read desired packing accuracy ( Default = 1.0D-14 )
! Note      : this flag is only active if iWRopt=0

9940 KWord = Get_Ln(LuRd)
call Get_F1(1,PkAcc)
PkAcc = abs(PkAcc)
Go To 998
!                                                                      *
!***** SKIP ************************************************************
!                                                                      *
! Read skip parameters,i.e.,
! if 2el integral symmetry blocks containing a given symmetry
! will not be needed in subsequent calculations their computation
! ans storage can be omitted.
! ( Default = 0,0,0,0,0,0,0,0 )
! Note      : this flag is only activ if iWRopt=0

9950 KWord = Get_Ln(LuRd)
lSkip = .true.
ChSkip = KWord(1:80)
Go To 998
!                                                                      *
!***** EXTR ************************************************************
!                                                                      *
! Put the program name and the time stamp onto the extract file

9960 write(LuWr,*) 'RdCtl: keyword EXTRACT is obsolete and is ignored!'
Go To 998
!                                                                      *
!***** REAC ************************************************************
!                                                                      *
! Read reaction field input.

9970 continue
GWInput = .true.
if (.not. RF_read) then
  call InpRct(LuRd)
  if (lLangevin .or. PCM) Go To 9971

  ! Add a center corresponding to the center of the RF cavity.

  RF_read = .true.
9971 continue

else
  call WarningMessage(2,'RdCtl: A second RF-input block discovered!')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** GRID *************************************************************
!                                                                      *
9773 continue
call Funi_input(LuRd)
Go To 998
!                                                                      *
!***** CLIG ************************************************************
!                                                                      *
! Speed of light (in au)

9000 KWord = Get_Ln(LuRd)
call Get_F1(1,CLightAU)
CLightAU = abs(CLightAU)
write(LuWr,*) 'The speed of light in this calculation =',CLightAU
Go To 998
!                                                                      *
!**** NEMO *************************************************************
!                                                                      *
800 continue
NEMO = .true.
Go To 998
!                                                                      *
!**** RMAT *************************************************************
!                                                                      *
! RmatR    : radius of the R-matrix sphere (bohr)

880 KWord = Get_Ln(LuRd)
call Get_F1(1,RMatR)
Go To 998
!                                                                      *
!**** RMEA *************************************************************
!                                                                      *
! Epsabs   : absolute precision of numerical radial integration

881 KWord = Get_Ln(LuRd)
call Get_F1(1,Epsabs)
Go To 998
!                                                                      *
!**** RMER *************************************************************
!                                                                      *
! Epsrel   : relative precision of numerical radial integration

882 KWord = Get_Ln(LuRd)
call Get_F1(1,Epsrel)
Go To 998
!                                                                      *
!**** RMQC *************************************************************
!                                                                      *
! qCoul    : effective charge of the target molecule

883 KWord = Get_Ln(LuRd)
call Get_F1(1,qCoul)
Go To 998
!                                                                      *
!**** RMDI *************************************************************
!                                                                      *
! dipol(3) : effective dipole moment of the target molecule

884 KWord = Get_Ln(LuRd)
call Get_F(1,dipol,3)
dipol1 = abs(dipol(1))+abs(dipol(2))+abs(dipol(3))
Go To 998
!                                                                      *
!**** RMEQ *************************************************************
!                                                                      *
! epsq     : minimal value of qCoul and/or dipol1 to be considered

885 KWord = Get_Ln(LuRd)
call Get_F1(1,epsq)
Go To 998
!                                                                      *
!**** RMBP *************************************************************
!                                                                      *
! bParm    : Bloch term parameter

886 KWord = Get_Ln(LuRd)
call Get_F1(1,bParm)
Go To 998
!                                                                      *
!**** GIAO *************************************************************
!                                                                      *
! Enable GIAO integrals.

9020 GIAO = .true.
Go To 998
!                                                                      *
!**** NOCH *************************************************************
!                                                                      *
! Deactivate Cholesky decomposition.

9089 continue
Cholesky = .false.
CholeskyWasSet = .true.
Do_RI = .false.
Go To 998
!                                                                      *
!**** CHOL *************************************************************
!                                                                      *
! Activate Cholesky decomposition with default settings.
! This section can only be executed once.

9091 continue
Do_RI = .false.
if (.not. CholeskyWasSet) then
  CholeskyWasSet = .true.
  Cholesky = .true.
  Do_RI = .false.
  DirInt = .true.
  call Cho_Inp(.true.,-1,6)
  iChk_CH = 1
end if
if ((iChk_RI+iChk_DC) > 0) then
  call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
  call Quit_OnUserError()
end if
GWInput = .false.
Go To 998
!                                                                      *
!**** THRC *************************************************************
!                                                                      *
! Set Cholesky decomposition threshold to specified value.

9021 continue
KWord = Get_Ln(LuRd)
call Get_F1(1,CholeskyThr)
Go To 998
!                                                                      *
!**** 1CCD *************************************************************
!                                                                      *
! Use one-center Cholesky.

9022 continue
do1CCD = .true.
Do_RI = .false.
iChk_Ch = 1
if ((iChk_RI+iChk_DC) > 0) then
  call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** CHOI *************************************************************
!                                                                      *
! Activate Cholesky decomposition with user-defined settings.
! This section can be executed any number of times.
! User-defined settings will be preferred to defaults also if the
! keywords appear in "wrong" order,
!
! ChoInput
! ...
! End ChoInput
! Cholesky

9092 continue
Do_RI = .false.
CholeskyWasSet = .true.
Cholesky = .true.
DirInt = .true.
call Cho_Inp(.false.,LuRd,6)
iChk_CH = 1
if ((iChk_RI+iChk_DC) > 0) then
  call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** RP-C *************************************************************
!                                                                      *
9093 lRP = .true.
KWord = Get_Ln(LuRd)
jTmp = 0
nRP_prev = -1
read(KWord,*,err=9082) nRP
isnumber = 1
ifile = index(KWord,' ')
do i=1,ifile
  ii = index(' 0123456789',KWord(i:i))
  if (ii == 0) then
    isnumber = 0
  end if
end do
!
if (isnumber == 0) goto 9082
nRP = 3*nRP
call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')

! Inline input

call UpCase(KWord)
if (index(KWord,'ANGSTROM') /= 0) then
  Fact = One/Angstr
else
  Fact = One
end if

KWord = Get_Ln(LuRd)
call Get_F1(1,E1)
call Read_v(LuRd,RP_Centers(1,1,1),1,nRP,1,iErr)
KWord = Get_Ln(LuRd)
call Get_F1(1,E2)
call Read_v(LuRd,RP_Centers(1,1,2),1,nRP,1,iErr)
RP_Centers(:,:,:) = Fact*RP_Centers(:,:,:)
GWInput = .true.
Go To 998

! Files

9082 continue
RPSet = .true.
jTmp = jTmp+1
ifile = index(KWord,' ')
if (KWord(1:1) == '/') then
  call f_inquire(KWord(1:ifile-1),Exist)
  Key = KWord
else
  call getenvf('MOLCAS_SUBMIT_DIR',Directory)
  if (Directory(1:1) /= ' ') then
    i = index(Directory,' ')
    Key = Directory(1:i-1)//'/'//KWord(1:ifile-1)
    ifile = i+ifile
    call f_inquire(Key(1:iFile-1),Exist)
  else
    Exist = .false.
  end if
  if (.not. Exist) then
    Key = Key(i+1:iFile-1)
    ifile = ifile-i
    call f_inquire(Key(1:iFile-1),Exist)
  end if
end if
if (.not. Exist) then
  call WarningMessage(2,'File '//Key(1:ifile)//' is not found')
  call Quit_OnUserError()
end if
LuIn = 8
LuIn = isFreeUnit(LuIn)
call molcas_open(LuIn,Key(1:iFile-1))

KWord = Get_Ln(LuIn)
read(KWord,*,err=9083) nRP
if ((nRP_prev >= 0) .and. (nRP /= nRP_prev)) then
  call WarningMessage(2,'The numbers of atoms in the two RP structures do not match.')
  call Quit_OnUserError()
end if
nRP_prev = nRP
nRP = 3*nRP
if (.not. allocated(RP_Centers)) call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')
KWord = Get_Ln(LuIn)
call UpCase(KWord)
if (index(KWord,'BOHR') /= 0) then
  Fact = One
else
  Fact = One/Angstr
end if
if (jTmp == 1) then
  LuRP = 10
  LuRP = isFreeUnit(LuRP)
  call molcas_open(LuRP,'findsym.RP1')
  read(KWord,*,err=9083) E1

  ! write a separate file for findsym

  write(LuRP,*) nRP/3
# ifdef _HAVE_EXTRA_
  write(LuRP,'(a)')
# else
  write(LuRP,'(a)') 'bohr'
# endif
  do i=1,nRP/3
    KWord = Get_Ln(LuIn)
    read(KWord,*,err=9083) Key,(RP_Centers(j,i,1),j=1,3)
    write(LuRP,'(A,3F20.12)') Key(1:LENIN),(RP_Centers(j,i,1)*Fact,j=1,3)
  end do
  KWord = Get_Ln(LuRd)
  close(LuIn)
  close(LuRP)
  Go To 9082
else
  LuRP = 10
  LuRP = isFreeUnit(LuRP)
  call molcas_open(LuRP,'findsym.RP2')
  write(LuRP,*) nRP/3
# ifdef _HAVE_EXTRA_
  write(LuRP,'(a)')
# else
  write(LuRP,'(a)') 'bohr'
# endif
  read(KWord,*,err=9083) E2
  do i=1,nRP/3
    KWord = Get_Ln(LuIn)
    read(KWord,*,err=9083) Key,(RP_Centers(j,i,2),j=1,3)
    write(LuRP,'(A,3F20.12)') Key(1:LENIN),(RP_Centers(j,i,2)*Fact,j=1,3)
  end do
  close(LuRP)
end if
RP_Centers(:,:,:) = Fact*RP_Centers(:,:,:)

close(LuIn)
GWInput = Run_Mode == G_Mode
Go To 998

! Error

9083 continue
write(6,'(a,a)') 'Error reading from file ',Key(1:iFile-1)
write(6,'(a,a)') 'unable to process line: ',KWord
call Quit_OnUserError()
!                                                                      *
!**** SADD *************************************************************
!                                                                      *
! Saddle options
9081 Key = Get_Ln(LuRd)
call Get_F1(1,SadStep)
GWInput = .true.
Go To 998
!                                                                      *
!**** CELL *************************************************************
!                                                                      *
! VCell(3,3)    : the vectors of the cell

887 Key = Get_Ln(LuRd)
call Upcase(Key)
if (index(Key,'ANGSTROM') /= 0) KWord = Get_Ln(LuRd)
call Get_F(1,VCell(1,1),3)
KWord = Get_Ln(LuRd)
call Get_F(1,VCell(1,2),3)
KWord = Get_Ln(LuRd)
call Get_F(1,VCell(1,3),3)
if (index(Key,'ANGSTROM') /= 0) call DScal_(9,One/angstr,VCell,1)
Cell_l = .true.
call mma_allocate(AdCell,MxAtom)
Go To 998
!                                                                      *
!**** SPAN *************************************************************
!                                                                      *
! Set span factor in Cholesky decomposition (0 < span < 1).
! The span decides the smallest diagonal element that can be
! treated as span*max(Diag). Span=1 thus implies full pivoting.

890 continue
KWord = Get_Ln(LuRd)
call Get_F1(1,spanCD)
spanCD = abs(spanCD)
Go To 998
!                                                                      *
!**** SPREAD ***********************************************************
!                                                                      *
! ispread(3)    : the number of cells to spread in different directions

889 KWord = Get_Ln(LuRd)
call Get_I(1,ispread,3)
Go To 998
!                                                                      *
!**** LOW  *************************************************************
!                                                                      *
! Activate low-accuracy Cholesky decomposition.

9094 continue
Do_RI = .false.
if (.not. CholeskyWasSet) then
  CholeskyWasSet = .true.
  Cholesky = .true.
  DirInt = .true.
  call Cho_Inp(.true.,-1,6)
  call Cho_InpMod('LOW ')
  Thrshld_CD = 1.0D-4
end if
Go To 998
!                                                                      *
!**** MEDI *************************************************************
!                                                                      *
! Activate medium-accuracy Cholesky decomposition.

9095 continue
Do_RI = .false.
if (.not. CholeskyWasSet) then
  CholeskyWasSet = .true.
  Cholesky = .true.
  DirInt = .true.
  call Cho_Inp(.true.,-1,6)
  call Cho_InpMod('MEDI')
  Thrshld_CD = 1.0D-6
end if
Go To 998
!                                                                      *
!**** HIGH *************************************************************
!                                                                      *
! Activate high-accuracy Cholesky decomposition.

9096 continue
Do_RI = .false.
if (.not. CholeskyWasSet) then
  CholeskyWasSet = .true.
  Cholesky = .true.
  DirInt = .true.
  call Cho_Inp(.true.,-1,6)
  call Cho_InpMod('HIGH')
  Thrshld_CD = 1.0D-8
end if
Go To 998
!                                                                      *
!**** DIAG *************************************************************
!                                                                      *
9087 continue
DiagCheck = .true.
Go To 998
!                                                                      *
!**** RI   *************************************************************
!                                                                      *
! Activate RI approach

9097 continue
Do_RI = .true.
GWInput = .true.
iRI_Type = 3
iChk_RI = 1
if ((iChk_DC+iChk_CH) > 0) then
  call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
  call Quit_OnUserError()
end if
Go To 998
9098 continue
Do_RI = .true.
GWInput = .true.
iRI_Type = 1
iChk_RI = 1
if ((iChk_DC+iChk_CH) > 0) then
  call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
  call Quit_OnUserError()
end if
Go To 998
9099 continue
Do_RI = .true.
GWInput = .true.
iRI_Type = 2
iChk_RI = 1
if ((iChk_DC+iChk_CH) > 0) then
  call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
  call Quit_OnUserError()
end if
Go To 998
9080 continue
Do_RI = .true.
GWInput = .true.
iRI_Type = 4
iChk_RI = 1
if ((iChk_DC+iChk_CH) > 0) then
  call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
  call Quit_OnUserError()
end if
Go To 998
9085 continue
Do_RI = .true.
GWInput = .true.
iRI_Type = 5
iChk_RI = 1
if ((iChk_DC+iChk_CH) > 0) then
  call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** NOGU *************************************************************
!                                                                      *
! Disable atomatic execution of GuessOrb

9100 Do_GuessOrb = .false.
Go To 998
!                                                                      *
!**** RELA *************************************************************
!                                                                      *
! DKH option: order and parameterization.
! xx: order of Hamiltonian
!  y: parameterization
! zz: order of properties

657 continue
kWord = Get_Ln(LuRd)
if ((KWord(1:1) == 'R') .and. &
    ((KWord(2:2) >= '0') .and. (KWord(2:2) <= '9')) .and. ((KWord(3:3) >= '0') .and. (KWord(3:3) <= '9')) .and. &
    ((KWord(4:4) == 'O') .or. (KWord(4:4) == 'E') .or. (KWord(4:4) == 'S') .or. (KWord(4:4) == 'M') .or. (KWord(4:4) == 'C'))) then
else
  call WarningMessage(2,'Error in RELA keyword')
  call Quit_OnUserError()
end if
DKroll = .true.

! DKH order in the Hamiltonian

read(KWord(2:3),*) idk_ord
IRELAE = 1000+idk_ord*10

! Method of parametrization

if (kWord(4:4) == 'O') IRELAE = IRELAE+1
if (kWord(4:4) == 'E') IRELAE = IRELAE+2
if (kWord(4:4) == 'S') IRELAE = IRELAE+3
if (kWord(4:4) == 'M') IRELAE = IRELAE+4
if (kWord(4:4) == 'C') IRELAE = IRELAE+5

! DKH order in the property integrals

if ((KWord(5:5) >= '0') .and. (KWord(5:5) <= '9') .and. (KWord(6:6) >= '0') .and. (KWord(6:6) <= '9')) then
  read(KWord(5:6),*) iprop_ord
  IRELAE = IRELAE+iprop_ord*10000
else
  IRELAE = IRELAE+idk_ord*10000
end if

Go To 998
!                                                                      *
!**** LDKH *************************************************************
!                                                                      *
! Local Douglas-Kroll-Hess/X2C/BSS

658 LDKroll = .true.
!GWInput = .True.
nCtrLD = 0
radiLD = 5.5d0

KWord = Get_Ln(LuRd)
call Upcase(KWord)
if (KWord(1:3) == 'DLU') Go To 998
if (KWord(1:3) == 'DLH') then
  radiLD = 0.0d0
  Go To 998
end if
read(KWord,*,end=6582,err=6582) nCtrLD,radiLD
if (nCtrLD > 10) then
  call WarningMessage(2,'The number of centers for LDKH is limited to 10')
  call Quit_OnUserError()
end if
if (index(KWord,'ANGSTROM') /= 0) then
  radiLD = radiLD/angstr
end if

KWord = Get_Ln(LuRd)
call Upcase(KWord)
read(KWord,*,end=6666,err=6581) (iCtrLD(i),i=1,nCtrLD)
Go To 998
6581 read(Kword,*,end=6666,err=6666) (CtrLDK(i),i=1,nCtrLD)
call Get_nAtoms_all(nAtom)
k = 0
do i=1,nAtom
  do j=1,nCtrLD
    if (CtrLDK(j) == dc(i)%LblCnt(1:LENIN)) then
      iCtrLD(j) = i
      k = k+1
    end if
  end do
end do
if (k /= nCtrLD) then
  call WarningMessage(2,'Error in LDKH Centers definitions')
  call Quit_OnUserError()
end if
Go To 998

! Automatic choice: all heavy elements (from K 19)

6582 continue
!DP if (nCtrLD == 0) radiLD = 0.0d0
Key = KWord
Go To 9989
!                                                                      *
!**** FOOC *************************************************************
!                                                                      *
! Force the use of the out-of-core RI algorithm.

8000 Force_Out_of_Core = .true.
Go To 998
!                                                                      *
!**** CDTH *************************************************************
!                                                                      *
! Threshold for CD to generate RICD auxiliary basis sets

8001 Key = Get_Ln(LuRd)
call Get_F1(1,Thrshld_CD)
GWInput = .true.
Go To 998
!                                                                      *
!**** SHAC *************************************************************
!                                                                      *
! Skip high angular combinations when constructing RICD aux basis.

8002 Skip_High_AC = .true.
GWInput = .true.
Go To 998
!                                                                      *
!**** KHAC *************************************************************
!                                                                      *
! Keep high angular combinations when constructing RICD aux basis.

8003 Skip_High_AC = .false.
GWInput = .true.
Go To 998
!                                                                      *
!**** ACD  *************************************************************
!                                                                      *
! Generate a aCD basis.

8004 Do_acCD_Basis = .false.
GWInput = .true.
Go To 998
!                                                                      *
!**** ACCD *************************************************************
!                                                                      *
! Generate a acCD basis.

8005 Do_acCD_Basis = .true.
GWInput = .true.
Go To 998
!                                                                      *
!**** NACC *************************************************************
!                                                                      *
! Generate a nacCD basis.

8030 Do_acCD_Basis = .false.
Do_nacCD_Basis = .true.
GWInput = .true.
Go To 998
!                                                                      *
!**** DOFM *************************************************************
!                                                                      *
! DoFMM: activate FMM option

8006 DoFMM = .true.
Go To 998
!                                                                      *
!**** NOAM *************************************************************
!                                                                      *
! No computation of AMFI integrals

8007 NoAMFI = .true.
GWInput = .true.
Go To 998
!                                                                      *
!**** RPQM *************************************************************
!                                                                      *
! Set RPQMin for FMM option

8008 Key = Get_Ln(LuRd)
call Get_F1(1,RPQMin)
Go To 998
!                                                                      *
!**** CONS *************************************************************
!                                                                      *
! Have the Gateway read the constraints for Slapaf

8010 continue
GWInput = .true.
Lu_UDC = 97
Lu_UDC = IsFreeUnit(Lu_UDC)
call Molcas_Open(Lu_UDC,'UDC.Gateway')
8011 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
write(Lu_UDC,'(A)') trim(Key)
if (Key(1:4) /= 'END ') Go To 8011
! This rather obscure feature seems to be needed to to make Intel
! compilers behave like the others when detecting EOF
endfile(Lu_UDC)
close(Lu_UDC)
Go To 998
!                                                                      *
!**** NGEX *************************************************************
!                                                                      *
! Have the Gateway read the constraints for Numerical_gradient

501 continue
GWInput = .true.
Lu_UDC = 97
Lu_UDC = IsFreeUnit(Lu_UDC)
call Molcas_Open(Lu_UDC,'UDC.NG')
502 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
if (adjustl(Key) == 'INVERT') then
  Invert = .true.
  Go To 502
end if
write(Lu_UDC,'(A)') trim(Key)
if (Key(1:4) /= 'END ') Go To 502
! This rather obscure feature seems to be needed to to make Intel
! compilers behave like the others when detecting EOF
endfile(Lu_UDC)
close(Lu_UDC)
Go To 998
!                                                                      *
!**** LOCA or LDF1 or LDF **********************************************
!                                                                      *
! Activate Local Density Fitting.

35 continue
LocalDF = .true.
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** LDF2 *************************************************************
!                                                                      *
! Activate Local Density Fitting with 2-center functions included
! when needed to achieve target accuracy.

36 continue
LocalDF = .true.
call LDF_SetLDF2(.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** TARG or THRL *****************************************************
!                                                                      *
! Set target accuracy for Local Density Fitting.
! This implies inclusion of 2-center functions (the only way we can
! affect accuracy).

37 continue
Key = Get_Ln(LuRd)
call Get_F1(1,Target_Accuracy)
call LDF_SetThrs(Target_Accuracy)
LocalDF = .true.
call LDF_SetLDF2(.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** APTH *************************************************************
!                                                                      *
! Set screening threshold for LDF - i.e. threshold for defining
! significant atom pairs.

38 continue
Key = Get_Ln(LuRd)
call Get_F1(1,APThr)
call LDF_SetPrescreen(APThr)
LocalDF = .true.
APThr_UsrDef = .true.
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** CHEC *************************************************************
!                                                                      *
! LDF debug option: check pair integrals.

39 continue
call LDF_SetOptionFlag('CHEC',.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** VERI *************************************************************
!                                                                      *
! LDF debug option: verify fit for each atom pair.

40 continue
call LDF_SetOptionFlag('VERI',.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** OVER *************************************************************
!                                                                      *
! LDF debug option: check overlap integrals (i.e. charge)

41 continue
call LDF_SetOptionFlag('OVER',.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** CLDF *************************************************************
!                                                                      *
! Constrained LDF - read constraint order
! order=-1 --- unconstrained
! order=0  --- charge (i.e. overlap)

42 continue
Key = Get_Ln(LuRd)
call Get_I1(1,iCLDF)
call LDF_AddConstraint(iCLDF)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** UNCO *************************************************************
!                                                                      *
! Unconstrained LDF (same as CLDF=-1)

43 continue
call LDF_AddConstraint(-1)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** WRUC *************************************************************
!                                                                      *
! Write unconstrained coefficients to disk.
! Only meaningful along with constrained fitting.
! For debugging purposes: enables constrained fit verification in
! modules other than Seward.

44 continue
call LDF_SetOptionFlag('WRUC',.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** UNIQ *************************************************************
!                                                                      *
! LDF: use unique atom pairs.

45 continue
call LDF_SetOptionFlag('UNIQ',.true.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** NOUN *************************************************************
!                                                                      *
! LDF: do not use unique atom pairs.

46 continue
call LDF_SetOptionFlag('UNIQ',.false.)
GWInput = .false. ! Only in Seward
Go To 998
!                                                                      *
!**** RLDF *************************************************************
!                                                                      *
! Activate local DF/RI, Roland's original LDF test implementation

8012 LDF = .true.
GWInput = .true.
Go To 998
!                                                                      *
!**** NOAL *************************************************************
!                                                                      *
! Do not align reactants and products

7013 Do_Align = .false.
if (Align_Only) then
  call WarningMessage(2,'Keywords ALIG and NOAL are not compatible')
  call Quit_OnUserError()
end if
GWInput = .true.
Go To 998
!                                                                      *
!**** WEIG *************************************************************
!                                                                      *
! Weights for alignment of reactants and products

7014 Align_Weights = Get_Ln(LuRd)
call UpCase(Align_Weights)
GWInput = .true.
Go To 998
!                                                                      *
!**** ALIG *************************************************************
!                                                                      *
! Align reactants and products

8013 Align_Only = .true.
if (.not. Do_Align) then
  call WarningMessage(2,'Keywords ALIG and NOAL are not compatible')
  call Quit_OnUserError()
end if
GWInput = .true.
Go To 998
!                                                                      *
!***** TINK ************************************************************
!                                                                      *
! Read Coordinates in Tinker's xyz format

8014 if (SymmSet) then
  call WarningMessage(2,'SYMMETRY keyword is not compatible with TINKER')
  call Quit_OnUserError()
end if
DoTinker = .true.
ITkQMMM = 1
if (MyRank == 0) then
  ITkQMMM = IsFreeUnit(ITkQMMM)
  call Molcas_Open(ITkQMMM,'QMMM')
  write(ITkQMMM,'(A)') 'Molcas -1 0'
  close(ITkQMMM)

  call Getenvf('TINKER ',Key)
  mLine = len(Key)
  iLast = iCLast(Key,mLine)
  if (iLast == 0) then
    call Getenvf('MOLCAS',Key)
    mLine = len(Key)
    iLast = iCLast(Key,mLine)
    Key = Key(1:iLast)//'/tinker/bin'
  end if
  iLast = iCLast(Key,mLine)
  call Getenvf('Project',Project)
  mLine = len(Project)
  jLast = iCLast(Project,mLine)
  Key = Key(1:iLast)//'/tkr2qm_s '//Project(1:jLast)//'.xyz>'//Project(1:jLast)//'.Tinker.log'
  mLine = len(Key)
  iLast = iCLast(Key,mLine)
  write(6,*) 'TINKER keyword found, run ',Key(1:iLast)
  call StatusLine(' Gateway:',' Read input from Tinker')
  RC = 0
  call Systemf(Key(1:iLast),RC)
  if (RC /= 0) then
    Key = 'RdCtl_Seward: Tinker call terminated abnormally'
    call WarningMessage(2,Key)
    call Abend()
  end if
end if
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  call PFGet_ASCII('QMMM')
  call GA_Sync()
end if
#endif
iCoord = iCoord+1
CoordSet = .true.
ITkQMMM = IsFreeUnit(ITkQMMM)
call Molcas_Open(ITkQMMM,'QMMM')
#ifdef _HAVE_EXTRA_
if (Expert) then
  if (iCoord > 1) then
    call WarningMessage(1,'TINKER and COORD keywords cannot be combined with molcas_extra')
  end if
end if
call XYZread(ITkQMMM,ForceZMAT,nCoord,iErr)
if (iErr /= 0) then
  Key = 'RdCtl_Seward: Tinker+XYZread failed: check Tinker input files'
  call WarningMessage(2,Key)
  call Abend()
end if
call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
if (Expert) then
  if (iCoord > 1) then
    call WarningMessage(1,'TINKER coordinates replacing COORD')
  end if
  call Read_XYZ(ITkQMMM,OrigRot,OrigTrans,Replace=(iCoord > 1))
else
  call Read_XYZ(ITkQMMM,OrigRot,OrigTrans)
end if
#endif
close(ITkQMMM)
GWInput = .true.
Go To 998
!                                                                      *
!**** ORIG *************************************************************
!                                                                      *
! Defines translation and rotation for each xyz-file

8015 OrigInput = .true.
#ifdef _HAVE_EXTRA_
Origin_input = .true.
#endif
if (FragSet) then
  write(6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
  call Quit_OnUserError()
end if
if (.not. OriginSet) then
  call mma_allocate(OrigTrans,3,nFragment,label='OrigTrans')
  call mma_allocate(OrigRot,3,3,nFragment,label='OrigRot')
  OriginSet = .true.
end if
do iFrag=1,nFragment
  KWord = Get_Ln(LuRd)
  call Get_F(1,OrigTrans(1,iFrag),3)
  KWord = Get_Ln(LuRd)
  call Get_F(1,OrigRot(1,1,iFrag),9)
end do
GWinput = .true.
Go To 998
!                                                                      *
!**** HYPE *************************************************************
!                                                                      *
8016 KWord = Get_Ln(LuRd)
#ifdef _HAVE_EXTRA_
geoInput = .true.
#endif
writeZmat = .true.
call Get_F(1,HypParam,3)
GWinput = .true.
HyperParSet = .true.
Go To 998
!                                                                      *
!**** ZCON *************************************************************
!                                                                      *
8017 writeZMat = .true.
#ifdef _HAVE_EXTRA_
ZConstraints = .true.
#endif
GWinput = .true.
Go To 998
!                                                                      *
!**** SCAL *************************************************************
!                                                                      *
8018 KWord = Get_Ln(LuRd)
call Get_F1(1,ScaleFactor)
GWinput = .true.
if (.not. CoordSet) then
  call WarningMessage(2,'Scale can be used only with xyz input')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** DOAN *************************************************************
!                                                                      *
8019 call Put_iScalar('agrad',1)
Go To 998
!                                                                      *
!**** GEOE *************************************************************
!                                                                      *
8020 Kword = Get_Ln(LuRd)
call Get_I1(1,iGeoInfo(2))
GWinput = .true.
iGeoInfo(1) = 1
call Put_iArray('GeoInfo',iGeoInfo,2)
Go To 998
!                                                                      *
!**** OLDZ *************************************************************
!                                                                      *
8021 GWinput = .true.
#ifdef _HAVE_EXTRA_
oldZmat = .true.
#endif
Go To 998
!                                                                      *
!**** OPTH *************************************************************
!                                                                      *
8022 GWinput = .true.
Kword = Get_Ln(LuRd)
call Get_I1(1,iOptimType)
Kword = Get_Ln(LuRd)
call Get_F1(1,StepFac1)
if (iOptimType == 2) then
  KWord = Get_Ln(LuRd)
  call Get_F1(1,gradLim)
end if
Go To 998
!                                                                      *
!**** NOON *************************************************************
!                                                                      *
8023 Do_OneEl = .false.
Go To 998
!                                                                      *
!**** GEO  *************************************************************
!                                                                      *
8024 continue
#ifdef _HAVE_EXTRA_
geoInput = .true.
#endif
writeZMat = .true.
! Parameters for the gridsize is set to default-values if geo is
! used instead of hyper
if (.not. HyperParSet) then
  HypParam(1) = 0.15d0
  HypParam(2) = 2.5d0
  HypParam(3) = 2.5d0
end if
GWinput = .true.
Go To 998
!                                                                      *
!**** GEN1INT **********************************************************
!                                                                      *
! GEN1INT integrals

9023 if (IRELAE == 101) then
  lMXTC = .true.
else
  write(6,*) 'Keyword MXTC must be preceded by keyword RX2C!'
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!**** FRGM *************************************************************
!                                                                      *
8025 OrigInput = .true.
#ifdef _HAVE_EXTRA_
Origin_input = .true.
#endif
GWinput = .true.
if (OriginSet) then
  write(6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
  call Quit_OnUserError()
end if
if (.not. FragSet) then
  call mma_allocate(OrigTrans,3,nFragment,label='OrigTrans')
  call mma_allocate(OrigRot,3,3,nFragment,label='OrogRot')
  ! Set up no translation and no rotation as default
  call FZero(OrigTrans,3*nFragment)
  call FZero(OrigRot,9*nFragment)
  do i=1,nFragment
    OrigRot(1,1,i) = 1.0d0
    OrigRot(2,2,i) = 1.0d0
    OrigRot(3,3,i) = 1.0d0
  end do
  FragSet = .true.
end if
Kword = Get_Ln(LuRd)
call Get_I1(1,iFrag)
Go To 998
!                                                                      *
!**** TRAN *************************************************************
!                                                                      *
8026 if (.not. FragSet) then
  write(6,*) 'Keyword TRANS must be preceded by keyword FRAG!'
  call Quit_OnUserError()
end if
GWinput = .true.
Kword = Get_Ln(LuRd)
call Get_F(1,OrigTrans(1,iFrag),3)
Go To 998
!                                                                      *
!***** ROT  ************************************************************
!                                                                      *
8027 if (.not. FragSet) then
  write(6,*) 'Keyword ROT must be preceded by keyword FRAG!'

end if
GWinput = .true.
Kword = Get_Ln(LuRd)
call Get_F(1,OrigRot(1,1,iFrag),9)
Go To 998
!                                                                      *
!****** ZONL ***********************************************************
!                                                                      *
8028 GWinput = .true.
WriteZMat = .true.
Go To 998
!                                                                      *
!****** BASL ***********************************************************
!                                                                      *
8029 GWinput = .true.
BasLib = Get_Ln(LuRd)
Write_BasLib = .true.
Go To 998
!                                                                      *
!****** NUME ***********************************************************
!                                                                      *
8031 GWinput = .true.
iDNG = 1
Go To 998
!                                                                      *
!****** VART ***********************************************************
!                                                                      *
8032 GWinput = .true.
VarT = .true.
Go To 998
!                                                                      *
!****** VARR ***********************************************************
!                                                                      *
8033 GWinput = .true.
VarR = .true.
Go To 998
!                                                                      *
!****** SHAK ***********************************************************
!                                                                      *
8050 continue
GWinput = .true.
KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_F1(1,Shake)
if (index(KWord,'ANGSTROM') /= 0) Shake = Shake/angstr
Go To 998
!                                                                      *
!***** PAMF ************************************************************
!                                                                      *
! Disable AMFI for an atom type

8060 KWord = Get_Ln(LuRd)
call Get_I1(1,iAtom_Number)
No_AMFI(iAtom_Number) = .true.
Go To 998
!                                                                      *
!****** GROM ***********************************************************
!                                                                      *
! Import definition of QMMM system from Gromacs

8034 continue
#ifdef _GROMACS_
if (SymmSet) then
  Message = 'SYMMETRY keyword is not compatible with GROMACS'
  call WarningMessage(2,Message)
  call Quit_OnUserError()
end if
DoGromacs = .true.
GWInput = .true.
VarT = .true.
VarR = .true.

! Check for options
KWord = Get_Ln(LuRd)
call UpCase(KWord)
if (KWord(1:4) == 'SIMP') then
  nCastMM = 0
  call mma_allocate(CastMM,nCastMM)
else if (KWord(1:4) == 'CAST') then
  KWord = Get_Ln(LuRd)
  call Get_I(1,nCastMM,1)
  if (nCastMM <= 0) then
    Message = 'nCastMM is zero or negative'
    call WarningMessage(2,Message)
    call Quit_OnUserError()
  end if
  call mma_allocate(CastMM,nCastMM)
  KWord = Get_Ln(LuRd)
  call Get_I(1,CastMM,nCastMM)
  do iCastMM=1,nCastMM
    if (CastMM(iCastMM) <= 0) then
      Message = 'Impossible, MM index < 1'
      call WarningMessage(2,Message)
      call Quit_OnUserError()
    end if
  end do
else
  Message = 'GROMACS keyword found, but no valid option'
  call WarningMessage(2,Message)
  call Quit_OnUserError()
end if

! After the call to Fetch_QMMM, the inner subsystem is in a temporary
! xyz file and the outer subsystem is on the runfile
call Fetch_QMMM(CastMM,nCastMM)

call mma_deallocate(CastMM)

! Let Molcas read the xyz file
iCoord = iCoord+1
CoordSet = .true.
LuXYZ = 1
LuXYZ = isFreeUnit(LuXYZ)
call molcas_open(LuXYZ,'GMX.XYZ')
#ifdef _HAVE_EXTRA_
if (Expert) then
  if (iCoord > 1) then
    call WarningMessage(1,'GROMACS and COORD keywords cannot be combined with molcas_extra')
  end if
end if
call XYZread(LuXYZ,ForceZMAT,nCoord,iErr)
if (iErr /= 0) then
  Message = 'RdCtl_Seward: XYZread returned non-zero error code'
  call WarningMessage(2,Message)
  call Abend()
end if
call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#else
if (Expert) then
  if (iCoord > 1) then
    call WarningMessage(1,'TINKER coordinates replacing COORD')
  end if
  call Read_XYZ(LuXYZ,OrigRot,OrigTrans,Replace=(iCoord > 1))
else
  call Read_XYZ(LuXYZ,OrigRot,OrigTrans)
end if
#endif
close(LuXYZ)
#else
Message = 'Interface to Gromacs not installed'
call WarningMessage(2,Message)
call Quit_OnUserError()
#endif
Go To 998
!                                                                      *
!****** LINK ***********************************************************
!                                                                      *
! Define link atoms for a Molcas/Gromacs run

8036 continue
#ifdef _GROMACS_
GWInput = .true.
KWord = Get_Ln(LuRd)
call Get_I(1,nLA,1)
if (nLA <= 0) then
  Message = 'LA definition: nLA is zero or negative'
  call WarningMessage(2,Message)
  call Quit_OnUserError()
end if
#ifdef _DEBUGPRINT_
write(LuWr,'(/,a)') ' Link atoms (Gromacs numbering):'
write(LuWr,'(/,a)') '      LA     QM     MM     Scaling factor'
#endif
call mma_allocate(DefLA,3,nLA)
call mma_allocate(FactLA,nLA)
do iLA=1,nLA
  KWord = Get_Ln(LuRd)
  call Get_I(1,DefLA(1,iLA),3)
  call Get_F(4,FactLA(iLA),1)
# ifdef _DEBUGPRINT_
  write(LuWr,'(i8,2i7,F19.8)') (DefLA(i,iLA),i=1,3),FactLA(iLA)
# endif
  if (DefLA(1,iLA) <= 0) then
    call WarningMessage(2,'LA definition: index of LA atom < 1')
    call Quit_OnUserError()
  else if (DefLA(2,iLA) <= 0) then
    call WarningMessage(2,'LA definition: index of QM atom < 1')
    call Quit_OnUserError()
  else if (DefLA(3,iLA) <= 0) then
    call WarningMessage(2,'LA definition: index of MM atom < 1')
    call Quit_OnUserError()
  else if ((FactLA(iLA) <= Zero) .or. (FactLA(iLA) >= One)) then
    call WarningMessage(2,'LA definition: bad scaling factor')
    call Quit_OnUserError()
  end if
end do
call Put_iArray('LA Def',DefLA,3*nLA)
call Put_dArray('LA Fact',FactLA,nLA)
call mma_deallocate(DefLA)
call mma_deallocate(FactLA)
#else
Message = 'Interface to Gromacs not installed'
call WarningMessage(2,Message)
call Quit_OnUserError()
#endif
Go To 998
!                                                                      *
!***** EMFR ************************************************************
!                                                                      *
8035 GWinput = .true.
Kword = Get_Ln(LuRd)
call Upcase(KWord)
EMFR = .true.
call Get_F(1,KVector,3)
Temp = sqrt(KVector(1)**2+KVector(2)**2+KVector(3)**2)
KVector(1) = KVector(1)/Temp
KVector(2) = KVector(2)/Temp
KVector(3) = KVector(3)/Temp
! Get the wavelength in atomic units.
call Get_F1(4,Lambda)
if (index(KWord,'ANGSTROM') /= 0) Lambda = Lambda/angstr
if (index(KWord,'NANOMETER') /= 0) then
  Lambda = Ten*Lambda/angstr
end if
KVector(1) = ((Two*Pi)/Lambda)*KVector(1)
KVector(2) = ((Two*Pi)/Lambda)*KVector(2)
KVector(3) = ((Two*Pi)/Lambda)*KVector(3)
Go To 998
!                                                                      *
!***** NOCD ************************************************************
!                                                                      *
9084 GWinput = .true.
if (.not. CholeskyWasSet) then
  Do_RI = .false.
  iRI_Type = 0
  Cholesky = .false.
  CholeskyWasSet = .true.
end if
Go To 998
!                                                                      *
!***** FNMC ************************************************************
!                                                                      *
9086 GWinput = .true.
FNMC = .true.
Go To 998
!                                                                      *
!***** ISOT ************************************************************
!                                                                      *
7654 GWinput = .true.
KWord = Get_Ln(LuRd)
call Upcase(KWord)
call Get_I1(1,nIsotopes)
call mma_allocate(nIsot,nIsotopes,2)
call mma_allocate(mIsot,nIsotopes)
do i=1,nIsotopes
  KWord = Get_Ln(LuRd)
  call Upcase(KWord)
  call Get_I1(1,iAt)
  nIsot(i,1) = iAt
  if (index(KWord,'DALTON') /= 0) then
    call Get_F1(2,dMass)
    nIsot(i,2) = -1
    mIsot(i) = dMass*UToAU
  else
    call Get_I1(2,iIso)
    nIsot(i,2) = iIso
    mIsot(i) = -One
  end if
end do
Go To 998
!                                                                      *
!***** EFP  ************************************************************
!                                                                      *
9088 GWinput = .true.
Kword = Get_Ln(LuRd)
call Get_I1(1,nEFP_fragments)
allocate(FRAG_TYPE(nEFP_fragments))
allocate(ABC(3,nEFP_fragments))
Kword = Get_Ln(LuRd)
call Upcase(kWord)
if (KWord == 'XYZABC') then
  Coor_Type = XYZABC_type
  nEFP_Coor = 6
  allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
  write(LuWr,*) 'XYZABC option to be implemented'
  call Abend()
else if (KWord == 'POINTS') then
  Coor_Type = POINTS_type
  nEFP_Coor = 9
  allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
  do iFrag=1,nEFP_fragments
    KWord = Get_Ln(LuRd)
    FRAG_Type(iFrag) = KWord
    do i=1,3
      KWord = Get_Ln(LuRd)
      jend = index(KWord,' ')
      if (jEnd > LENIN+1) then
        write(LuWr,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LENIN,' characters!'
      end if
      ABC(i,iFrag) = KWord(1:min(LENIN,jend-1))
      call Get_F(2,EFP_COORS((i-1)*3+1,iFrag),3)
    end do
  end do
else if (KWord == 'ROTMAT') then
  Coor_Type = ROTMAT_type
  nEFP_Coor = 12
  allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
  write(LuWr,*) 'ROTMAT option to be implemented'
  call Abend()
else
  write(LuWr,*) 'Illegal EFP format :',KWord
  write(LuWr,*)
  write(LuWr,*) 'Allowed format: XYZABC,'
  write(LuWr,*) '                POINTS, and'
  write(LuWr,*) '                ROTMAT'
end if
lEFP = .true.
Go To 998
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     P O S T   P R O C E S S I N G                                    *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *

997 continue
! Postprocessing for COORD
!  ik=index(KeepBasis,'....')
!  if (ik /= 0) then
!    KeepBasis = KeepBasis(1:ik-1)
!  end if
if (CoordSet) then
  do ik=len(KeepBasis),1,-1
    if ((KeepBasis(ik:ik) /= ' ') .and. (KeepBasis(ik:ik) /= '.')) Go To 1997
  end do
1997 continue
  KeepBasis = KeepBasis(1:ik)
# ifdef _HAVE_EXTRA_
  call ProcessXYZ(BasisSet,KeepBasis,KeepGroup,iBSSE,SymThr,isHold,ScaleFactor,HyperParSet,isXfield)
# else
  call Parse_Basis(KeepBasis)
  call Parse_Group(KeepGroup,SymThr)
  call Write_SewInp('COORD',[iBSSE])
# endif
  if (writeZmat) then
#   ifdef _HAVE_EXTRA_
    stepFactor = stepFac1/(hypParam(1)*hypParam(1))
    call Geo_Setup_Drv(ishold,oldZMat,zConstraints,geoInput,hypParam,nFragment,iOptimType,stepFactor,gradLim)
#   else
    call WarningMessage(2,'molcas-extra not installed')
    call Quit_OnUserError()
#   endif
  end if
  DoneCoord = .true.
  if (isXfield == 1) then
    LuRd_saved = LuRd
    filename = 'findsym.xfield'
    lXF = .true.
    goto 9753
  end if
end if

9755 continue
if (CoordSet) then
  CoordSet = .false.
  LuRdSave = LuRd
  LuFS = IsFreeUnit(1)
# ifdef _HAVE_EXTRA_
  call Molcas_Open(LuFS,'FS.std')
# else
  call Molcas_Open(LuFS,'COORD')
# endif
  LuRd = LuFS
  GWInput = .true.
  Go To 998
else
  if (DoneCoord) then
    close(LuFS)
    LuRd = LuRdSave
#   ifndef _HAVE_EXTRA_
    call Clear_XYZ()
#   endif
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Put_lScalar('CSPF',CSPF)
call Put_cArray('Align_Weights',Align_Weights,512)
!                                                                      *
!***********************************************************************
!                                                                      *
! Isotopic specifications

if (.not. allocated(nIsot)) call mma_allocate(nIsot,0,2)

if (Run_Mode /= S_Mode) then
  ! Loop over unique centers
  iUnique = 0
  call mma_allocate(Isotopes,nCnttp,Label='Isotopes')
  do iCnttp=1,nCnttp
    nCnt = dbsc(iCnttp)%nCntr
    do iCnt=1,nCnt
      iUnique = iUnique+1
      ! Get the mass for this center
      dm = rMass(dbsc(iCnttp)%AtmNr)
      do j=1,size(nIsot,1)
        if (nIsot(j,1) == iUnique) then
          if (nIsot(j,2) >= 0) then
            dm = rMassx(dbsc(iCnttp)%AtmNr,nIsot(j,2))
          else
            dm = mIsot(j)
          end if
          exit
        end if
      end do
      if (iCnt == 1) then
        dbsc(iCnttp)%CntMass = dm
      else
        if (dm /= dbsc(iCnttp)%CntMass) then
          call WarningMessage(2,'Error: All centers of the same type must have the same mass')
          call Quit_OnUserError()
        end if
      end if
    end do
    Isotopes(iCnttp) = dbsc(iCnttp)%CntMass
  end do ! iCnttp
  call Put_dArray('Isotopes',Isotopes,nCnttp)
  call mma_deallocate(Isotopes)

  ! Find errors
  do j=1,size(nIsot,1)
    if (nIsot(j,1) > iUnique) then
      call WarningMessage(2,'Error: Isotope specification index larger than the number of unique centers')
      call Quit_OnUserError()
    end if
  end do
end if

! Deallocate
if (allocated(nIsot)) call mma_deallocate(nIsot)
if (allocated(mIsot)) call mma_deallocate(mIsot)
!                                                                      *
!***********************************************************************
!                                                                      *
! post-processing for RP-Coord

if (lRP .and. RPset) call processRP(KeepGroup,SymThr)

lAMFI = lAMFI .and. (.not. NoAMFI)

! Disable the RI flag if only one-electron integrals are requested

Do_RI = (.not. Onenly) .and. Do_RI
if (Do_RI) then
  if (LDF .and. LocalDF) then
    call WarningMessage(2,'LDF and LocalDF are incompatible')
    call Quit_OnUserError()
  end if
end if

iPrint = nPrint(iRout)

S%Mx_Shll = iShll+1
Max_Shells = S%Mx_Shll

if (nCnttp == 0) then
  call WarningMessage(2,'Input does not contain any basis sets')
  call Quit_OnUserError()
end if
if (mdc == 0) then
  call WarningMessage(2,'Input does not contain coordinates')
  call Quit_OnUserError()
end if
if (S%iAngMx < 0) then
  call WarningMessage(2,' There is an error somewhere in the input!;S%iAngMx < 0')
  call Quit_OnUserError()
end if
if (S%iAngMx > iTabMx) then
  call WarningMessage(2,' Too High angular momentum !!!')
  call Quit_OnUserError()
end if
if (DoTinker .and. DoGromacs) then
  call WarningMessage(2,'TINKER and GROMACS keywords cannot be used together')
  call Quit_OnUserError()
end if
if (DoTinker .and. (iCoord > 1)) then
  if (.not. Expert) then
    call WarningMessage(2,'TINKER and COORD keywords cannot be used together')
    call Quit_OnUserError()
  end if
end if
if (DoGromacs .and. (iCoord > 1)) then
  if (.not. Expert) then
    call WarningMessage(2,'GROMACS and COORD keywords cannot be used together')
    call Quit_OnUserError()
  end if
end if

if (Test) then
  Do_GuessOrb = .false.
  Do_FckInt = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Automatic DK and AMFI option for relativistic basis sets
!                                                                      *
!***********************************************************************
!                                                                      *
if (BSS .and. (.not. DKroll)) then
  call WarningMessage(2,';BSSM GOES ALWAYS WITH DOUGLAS. THE OPPOSITE IS NOT TRUE')
  call Abend()
end if

lECP = .false.
lPP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
end do
if ((lECP .or. lPP) .and. DKroll .and. (.not. Expert)) then
  call WarningMessage(2,' ECP option not compatible with Douglas-Kroll!')
  call Quit_OnUserError()
end if

if (imix == 1) then
  if (Expert) then
    call WarningMessage(1,' input is inconsistent!;SEWARD found basis sets of mixed relativistic (or non-relativistic) types!;'// &
                        'No relativistic option will be automatically enabled')
  else
    call WarningMessage(2,' input is inconsistent!;SEWARD found basis sets of mixed relativistic (or non-relativistic) types!')
    call Quit_OnUserError()
  end if
end if
if ((ifnr == 1) .and. (.not. Expert)) then
  if (DKroll) then
    call WarningMessage(1,';you requested the DK-option for;a non-relativistic basis.;This request will be ignored')
  end if
  DKroll = .false.
else if (ifnr == 0) then
  lAMFI = .not. NoAMFI
  if (.not. DKroll) then
    DKroll = .true.
    !if (iRELAE == -1) IRELAE = 201022
    if (iRELAE == -1) then
      if (itype == 2) then
        IRELAE = 1022
      else if (itype == 14) then
        IRELAE = 101
      end if
    end if
  end if
  if ((MolWgh /= 0) .and. (MolWgh /= 2)) MolWgh = 2
end if

if (NoDKroll) DKroll = .false.
#ifdef _HAVE_EXTRA_
if (DoEMPC) isHold = 1
#endif

if ((lECP .or. lPP) .and. lAMFI .and. (.not. Expert)) then
  call WarningMessage(2,' ECP option not compatible with AMFI!')
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Activate Finite Nucleus parameters

if (Nuclear_Model == Point_Charge) then
  if (ign == 2) Nuclear_Model = Gaussian_Type
  if (ign == 3) Nuclear_Model = mGaussian_Type
end if

do iCnttp=1,nCnttp
  if (Nuclear_Model == Gaussian_Type) then

    ! If ExpNuc not explicitly defined use default value.

    nMass = nint(dbsc(iCnttp)%CntMass/UToAU)
    if (dbsc(iCnttp)%ExpNuc < Zero) dbsc(iCnttp)%ExpNuc = NucExp(nMass)
  else if (Nuclear_Model == mGaussian_Type) then

    ! Get parameters for the Modified Gaussian Nuclear
    ! charge distribution.

    jAtmNr = dbsc(iCnttp)%AtmNr
    nMass = nint(dbsc(iCnttp)%CntMass/UToAU)
    call ModGauss(dble(jAtmNr),nMass,dbsc(iCnttp)%ExpNuc,dbsc(iCnttp)%w_mGauss)

  else

    ! Nothing to do for point charges!

  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Cholesky-specific postprocessing:
! 0) if 1C-CD is requested, do it.
! 1) reset integral prescreening thresholds (if not user-defined).
! 2) use default Cholesky normalization (if not user-defined);
!    for AMFI or Douglas-Kroll, use Molcas normalization (again,
!    if not user-defined).
! 3) Turn off Onenly flag, which might be set through the Direct
!    keyword. Thus, specifying Cholesky will force 2-el. int.
!    processing even with Direct specifed as well!
! 4) Turn off Dist flag (makes no sense with Cholesky).
! 5) Integral format on disk is irrelevant, except that Aces II is
!    not allowed. So, reset iWrOpt or quit (for Aces II).
! 6) if Cholesky threshold is specified, use it.
! 7) if span factor is specified, use it.

if (do1CCD) then
  if (.not. Cholesky) then
    DirInt = .true.
    call Cho_Inp(.true.,-1,6)
  end if
  Cholesky = .true.
  call Cho_InpMod('1CCD')
end if
if (Cholesky) then
  if (Onenly) then
    Cholesky = .false. ! we gotta be lazy in such cases
  else
    if (.not. CutInt_UsrDef) CutInt = Cho_CutInt
    if (.not. ThrInt_UsrDef) ThrInt = Cho_ThrInt
    if (.not. MolWgh_UsrDef) then
      if (lAMFI .or. DKroll) then
        MolWgh = 2
      else
        MolWgh = Cho_MolWgh
      end if
    end if
    if (iWrOpt == 2) then
      write(LuWr,*) 'Acess II format not allowed with Cholesky!!'
      call Quit_OnUserError()
    else if ((iWrOpt /= 0) .and. (iWrOpt /= 3)) then
      iWrOpt = 0
    end if
    if (CholeskyThr >= 0.0d0) then
      Thrshld_CD = CholeskyThr
      call Cho_SetDecompositionThreshold(Thrshld_CD)
      call Put_Thr_Cho(Thrshld_CD)
    end if
    if (spanCD >= 0.0d0) then
      v = min(spanCD,1.0d0)
      call Cho_SetSpan(v)
    end if
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Prprt) then
  Onenly = .true.
  Vlct = .false.
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Post processing for FAIEMP fragment data

if (lFAIEMP .and. (Run_Mode /= S_Mode)) call FragExpand(LuRd)
!                                                                      *
!***********************************************************************
!                                                                      *
! Post processing for RI and RI/CD option

if (Do_RI .and. (Run_Mode /= S_Mode)) then
  if (iRI_Type == 4) then

    ! Generate on-the-fly aCD or aTrue.cCD auxiliary basis sets.

    call Mk_RICD_Shells()

  else

    ! Pick up an externally defined auxiliary basis set.

    call Mk_RI_Shells(LuRd)

  end if
end if
if (Do_RI .and. LocalDF .and. (Run_Mode == S_Mode)) then
  call SetTargetAccuracy_LDF()
  if (CutInt_UsrDef .and. (.not. APThr_UsrDef)) then
    call LDF_SetPrescreen(CutInt)
    call LDF_CheckThrs()
  end if
  call LDF_CheckConfig()
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Post processing for Well integrals

do iWel=1,nWel
  if (Wel_Info(1,iWel) < Zero) then
    if (.not. lRF) then
      call WarningMessage(2,'; Input inconsistency!; ;Relative positions of well integrals can only be used if the cavity '// &
                          'radius has been specified!')
      call Quit_OnUserError()
    end if
    Wel_Info(1,iWel) = rds+abs(Wel_Info(1,iWel))
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Put up list for point at which the electric field will be
! evaluated. If nEF=0 the default points will be the unique
! centers.

if ((nOrdEF >= 0) .and. (Run_Mode /= S_Mode)) then
  if (nEF /= 0) then
    call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')
    EF_Centers(:,:) = EFt(:,:)
    call mma_deallocate(EFt)
  else
    nEF = 0
    do iCnttp=1,nCnttp
      if (.not. (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) nEF = nEF+dbsc(iCnttp)%nCntr
    end do
    call mma_allocate(EF_Centers,3,nEF,Label='EF_Centers')

    iEF = 1
    do iCnttp=1,nCnttp
      if (.not. (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) then
        call dcopy_(3*dbsc(iCnttp)%nCntr,dbsc(iCnttp)%Coor,1,EF_Centers(1,iEF),1)
        iEF = iEF+dbsc(iCnttp)%nCntr
      end if
    end do
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put up list for point at which the diamagnetic shielding will
! be evaluated. If nDMS=0 the default points will be the unique
! centers.

if (lDMS .and. (Run_Mode /= S_Mode)) then
  if (nDMS /= 0) then
    call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
    DMS_Centers(:,:) = DMSt(:,:)
    call mma_deallocate(DMSt)
  else
    nDMS = 0
    do iCnttp=1,nCnttp
      nDMS = nDMS+dbsc(iCnttp)%nCntr
    end do
    call mma_allocate(DMS_Centers,3,nDMS,Label='DMS_Centers')
    iDMS = 1
    do iCnttp=1,nCnttp
      call dcopy_(3*dbsc(iCnttp)%nCntr,dbsc(iCnttp)%Coor,1,DMS_Centers(1,iDMS),1)
      iDMS = iDMS+dbsc(iCnttp)%nCntr
    end do
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! If no multipole moment integrals are requested turn also of the
! computation of the velocity integrals.

if (S%nMltpl == 0) Vlct = .false.

! But turn it on again if explicitly requested

if (Vlct_) Vlct = .true.

! This is the highest order of any property operator.
! The default value of 4 is due to the mass-velocity operator
! which is computed by default.

nPrp = max(4,S%nMltpl)

! Setup of tables for coefficients of the Rys roots and weights.

nDiff = 0
if (S%iAngMx == 0) nDiff = 2
DoRys = .true.
if (DKroll .and. (nOrdEF > 0)) nDiff = nDiff+nOrdEF
if ((.not. Test) .and. (Run_Mode /= S_Mode)) call SetUp_RW(DoRys,nDiff)
!                                                                      *
!***********************************************************************
!                                                                      *
! Store information for the Douglas-Kroll code.

if (DKroll .or. NEMO) call Fill_rInfo1()
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute kOffAO and lOffAO

call Setup_OffAO()
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate labels for Cartesian and spherical basis sets.
! Generate the transformation matrix for cartesian to sphericals
! and contaminants. This has to be done after adding auxiliary or
! fragment basis sets.

call Sphere(S%iAngMx)
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
! Set up Symmetry_Info

call Symmetry_Info_Setup(nOper,Oper,max(S%iAngMx,3))

if (lSkip) then
  call Put_Ln(ChSkip)
  call Get_I(1,iSkip,nIrrep)
  Do_GuessOrb = .false.
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
! Generate list of Stabilizers , Stabilizer Index and distinct cosets

S%mCentr = 0
S%mCentr_Aux = 0
S%mCentr_Frag = 0
do iCnttp=1,nCnttp
  nCnt = dbsc(iCnttp)%nCntr
  do iCnt=1,nCnt
    mdc = iCnt+dbsc(iCnttp)%mdci
    S%Mx_mdc = max(S%Mx_mdc,mdc)
    n_dc = max(mdc,n_dc)
    if (mdc > MxAtom) then
      call WarningMessage(2,' mdc > MxAtom!; Increase MxAtom in Molcas.fh.')
      write(LuWr,*) ' MxAtom=',MxAtom
      call Abend()
    end if

    ! The symmetry operators of the fragment's atoms should
    ! always be identical to that of the fragment's
    ! pseudocenter/placeholder

    if (dbsc(iCnttp)%Frag) then
      !  Check the FragExpand routine!
      if (abs(dbsc(iCnttp)%nFragCoor) > mdc) then
        write(6,*) 'rdctl_seward: incorrect mdc index'
        call Abend()
      end if
      iChxyz = dc(abs(dbsc(iCnttp)%nFragCoor))%iChCnt
    else

      ! To assign the character of a center we need to find
      ! the cartesian components that are permutable. We
      ! will only need to loop over the generators of the
      ! group. We will use the three first bits to indicate if
      ! the cartesian component is affected by any symmetry
      ! operation.

      iChxyz = iChAtm(dbsc(iCnttp)%Coor(:,iCnt))
    end if
    dc(mdc)%iChCnt = iChxyz
    call Stblz(iChxyz,dc(mdc)%nStab,dc(mdc)%iStab,nIrrep,dc(mdc)%iCoSet)

    ! Perturb the initial geometry if the SHAKE keyword was given,
    ! but maintain the symmetry

    if ((Shake > Zero) .and. (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux))) then
      jTmp = 0
      do j=1,dc(mdc)%nStab-1
        jTmp = ior(jTmp,dc(mdc)%iStab(j))
      end do
      S%nDim = 0
      do j=0,2
        if (iand(jTmp,2**j) == 0) S%nDim = S%nDim+1
      end do
      if (S%nDim > 0) then
        call Random_Vector(S%nDim,RandVect(1:S%nDim),.false.)
        jDim = 0
        do j=0,2
          if (iand(jTmp,2**j) == 0) then
            jDim = jDim+1
            dbsc(iCnttp)%Coor(j+1,iCnt) = dbsc(iCnttp)%Coor(j+1,iCnt)+Shake*RandVect(jDim)
          end if
        end do
      end if
    end if
    if (dbsc(iCnttp)%Frag) then
      S%mCentr_Frag = S%mCentr_Frag+nIrrep/dc(mdc)%nStab
    else if (dbsc(iCnttp)%Aux) then
      S%mCentr_Aux = S%mCentr_Aux+nIrrep/dc(mdc)%nStab
    else
      S%mCentr = S%mCentr+nIrrep/dc(mdc)%nStab
    end if
  end do
end do
if (S%mCentr > MxAtom) then
  call WarningMessage(2,'RdCtl: S%mCentr > MxAtom')
  write(6,*) 'S%mCentr=',S%mCentr
  write(6,*) 'Edit src/Include/Molcas.fh'
  write(6,*) 'Set MxAtom to the value of S%mCentr.'
  write(6,*) 'Recompile MOLCAS and try again!'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((SymmSet .or. (nIrrep > 1)) .and. LDKroll) then
  call WarningMessage(2,'Local DKH approach is not yet implemented with SYMMETRY')
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_GuessOrb .and. (Run_Mode /= S_Mode)) call Fix_FockOp(LuRd)
!                                                                      *
!***********************************************************************
!                                                                      *
if ((iXPolType /= 0) .and. (nIrrep /= 1)) then
  call WarningMessage(2,'Polarizabilities are not compatible with symmetry.')
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((nTtl /= 0) .and. (Run_Mode == G_Mode)) then
  if (iPrint >= 6) then
    write(LuWr,*)
    write(LuWr,'(15X,88A)') ('*',i=1,88)
    write(LuWr,'(15X,88A)') '*',(' ',i=1,86),'*'
    do iTtl=1,nTtl
      write(LuWr,'(15X,A,A,A)') '*   ',Title(iTtl),'   *'
    end do
    write(LuWr,'(15X,88A)') '*',(' ',i=1,86),'*'
    write(LuWr,'(15X,88A)') ('*',i=1,88)
  else
    write(LuWr,*)
    write(LuWr,'(A)') ' Title:'
    do iTtl=1,nTtl
      write(LuWr,'(8X,A)') Title(iTtl)
    end do
    write(LuWr,*)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the weights used for alignment and distance measurement

call Process_Weights(iPrint)
!
!***********************************************************************
!                                                                      *
! Set structures for TS optimization according to the Saddle method.

if (Run_Mode /= G_Mode) then
  call Saddle()
!                                                                      *
!***********************************************************************
!                                                                      *
! Read coordinates from run file (if any), ditto for external field.
! Do not do this in the Gateway!

  call GeoNew(Show)
  if (lXF) call GeoNew_PC()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Gen_GeoList()
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of centers for multipole operators

call SetMltplCenters()

! Put in user specified centers if any

if (lMltpl) then
  do i=1,nTemp
    iMltpl = ITmp(i)
    if (iMltpl <= S%nMltpl) call dcopy_(3,RTmp(1,i),1,Coor_MPM(1,iMltpl+1),1)
  end do
  call mma_deallocate(RTmp)
  call mma_deallocate(ITmp)
end if
#ifdef _DEBUGPRINT_
call RecPrt(' Multipole centers',' ',Coor_MPM,3,S%nMltpl+1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Put up list for point at which the orbital angular momentum
! will be computed.

if (lOAM .and. (Run_Mode /= S_Mode)) then
  call mma_allocate(OAM_Center,3,Label='OAM_Center')
  call dcopy_(3,OAMt,1,OAM_Center,1)
else if (.not. allocated(OAM_Center)) then
  call mma_allocate(OAM_Center,3,Label='OAM_Center')
  call dcopy_(3,CoM,1,OAM_Center,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put up list for point at which the orbital magnetic quadrupole
! will be computed.

if (lOMQ .and. (Run_Mode /= S_Mode)) then
  call mma_allocate(OMQ_Center,3,Label='OMQ_Center')
  call DCopy_(3,OMQt,1,OMQ_Center,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate fields from keyword ORIGIN

if (OrigInput) then
  call mma_deallocate(OrigRot)
  call mma_deallocate(OrigTrans)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (NoZMAT .and. (.not. ForceZMAT)) call Put_iScalar('N ZMAT',0)
if (Run_Mode /= S_Mode) call Put_iArray('BasType',BasisTypes,4)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Run_Mode == G_Mode) call Put_lScalar('Invert constraints',Invert)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Buffer)
!                                                                      *
!***********************************************************************
!                                                                      *
call Datimx(KWord)
Header(1) = Title(1)(5:76)
write(Header(2),'(4A)') ' Integrals generated by ',Vrsn,', ',KWord(1:24)
call Put_cArray('Seward Title',Header(1),144)
if (nTtl > 0) call Put_cArray('SewardXTitle',Title(1),nTtl*80)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Run_Mode == G_Mode) call Put_iScalar('DNG',iDNG)
!                                                                      *
!***********************************************************************
!                                                                      *

call mma_deallocate(STDINP)

return

6666 call WarningMessage(2,'Unable to read data from '//KWord)
call Quit_OnUserError()

end subroutine RdCtl_Seward
