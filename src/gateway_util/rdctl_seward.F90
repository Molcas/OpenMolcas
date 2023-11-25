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

subroutine RdCtl_Seward(LuRd_,lOPTO,Do_OneEl)

use AMFI_Info, only: No_AMFI
use Basis_Info, only: dbsc, Gaussian_Type, Max_Shells, mGaussian_Type, MolWgh, nCnttp, Nuclear_Model, Point_Charge, Shells
use Center_Info, only: dc, n_dc
use Her_RW, only: nPrp
use Period, only: AdCell, Cell_l, lthCell, ispread, VCell
use MpmC, only: Coor_MPM
use EFP_Module, only: ABC, Coor_Type, EFP_COORS, FRAG_TYPE, lEFP, nEFP_Coor, nEFP_fragments, POINTS_TYPE, ROTMAT_type, XYZABC_type
use Real_Spherical, only: Sphere
use fortran_strings, only: str
use External_Centers, only: AMP_Center, DMS_Centers, Dxyz, EF_Centers, iXPolType, nData_XF, nDMS, nEF, nOrd_XF, nOrdEF, nRP, nWel, &
                            nXF, nXMolnr, OAM_Center, OMQ_Center, RP_Centers, Wel_Info, XEle, XF, XMolnr
use Symmetry_Info, only: iSkip, nIrrep, Symmetry_Info_Setup, VarR, VarT
use Sizes_of_Seward, only: S
use Gateway_Info, only: Align_Only, CoM, CutInt, Do_Align, Do_FckInt, Do_GuessOrb, DoFMM, E1, E2, EMFR, FNMC, GIAO, kVector, &
                        lAMFI, lDOWNONLY, lMXTC, lRel, lRP, lSchw, lUPONLY, NEMO, PkAcc, RPQMin, Rtrnc, SadStep, Shake, ThrInt, &
                        Thrs, UnNorm, Vlct
use DKH_Info, only: iCtrLD, BSS, cLightAU, DKroll, IRELAE, LDKRoll, nCtrlD, radiLD
use Cholesky, only: Span, ThrCom
use RICD_Info, only: Chol => Cholesky, DiagCheck, Do_acCD_Basis, Do_DCCD, Do_RI, iRI_Type, Skip_High_AC, Thrshld_CD
use Gateway_global, only: DirInt, Expert, Fake_ERIs, Force_Out_of_Core, force_part_c, force_part_p, G_Mode, ifallorb, iPack, &
                          NoTab, Onenly, Prprt, Run_Mode, S_Mode, Short, SW_FileOrb, Test
use rctfld_module, only: lLangevin, lRF, PCM, RDS
use rmat, only: bParm, Dipol, Dipol1, EpsAbs, EpsQ, EpsRel, QCoul, RMat_On, RMatR
use define_af, only: AngTp, iTabMx
use getline_mod, only: Line
#ifdef _FDE_
use Embedding_Global, only: embOutDensPath, embOutEspPath, embOutGradPath, embOutHessPath, embPot, embPotInBasis, embPotPath, &
                            embWriteDens, embWriteEsp, embWriteGrad, embWriteHess, outGridPath, outGridPathGiven
#endif
#ifndef _HAVE_EXTRA_
use XYZ, only: Clear_XYZ, Parse_Basis, Parse_Group, Read_XYZ, Write_SewInp
#endif
use Para_Info, only: MyRank
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Ten, Pi, Angstrom, mu2elmass, UtoAU
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuRd_
logical(kind=iwp), intent(inout) :: lOPTO
logical(kind=iwp), intent(out) :: Do_OneEl
#include "Molcas.fh"
#include "print.fh"
#include "embpcharg.fh"
#ifdef _HAVE_EXTRA_
#include "hyper.fh"
#endif
integer(kind=iwp), parameter :: MAX_XBAS = 20
integer(kind=iwp) :: BasisTypes(4), BasisTypes_Save(4), i, i1, i2, iAng, iAt, iAtom_Number, ib, ibla, iBSSE, iChk_CH, iChk_DC, &
                     iChk_RI, iChrct, iChxyz, iCnt, iCnttp, iCoord, idk_ord, iDMS, iDNG, iDummy_basis, iEF, ierr, ifile, ifnr, &
                     iFound_Label, iFrag, iFrst, iGeoInfo(2), iglobal, ign, ii, iIso, ik, imix, iMltpl, Indx, iOff, iOff0, &
                     iOpt_XYZ, iOptimType, iPrint, iprop_ord, iRout, iShll, ist, istatus, isxbas, isXfield, iTemp, ITkQMMM, iTtl, &
                     itype, iUnique, iWel, ix, j, jAtmNr, jDim, jend, jRout, jShll, jTmp, k, lAng, Last, lAW, lSTDINP, Lu_UDC, &
                     LuFS, LuIn, LuRd, LuRd_saved, LuRdSave, LuRP, mdc, n, nAtom, nc, nc2, nCnt, nCnt0, nDataRead, nDiff, nDone, &
                     nFragment, nIsotopes, nMass, nOper, nReadEle, nRP_prev, nTemp, nTtl, nxbas, RC
real(kind=wp) :: CholeskyThr, dm, dMass, Fact, gradLim, HypParam(3), Lambda, OAMt(3), OMQt(3), RandVect(3), ScaleFactor, sDel, &
                 spanCD, stepFac1, SymThr, tDel, Temp

logical(kind=iwp) :: AnyMode, Basis_test, BasisSet, CholeskyWasSet, Convert, CoordSet, CSPF = .false., &
                     CutInt_UsrDef, DoGromacs, DoneCoord, DoRys, DoTinker, EFgiven, Exists, FinishBasis, ForceZMAT, FOUND, &
                     FragSet, GroupSet, GWInput, HyperParSet, Invert, isnumber, lDMS = .false., lECP, lFAIEMP = .false., lMltpl, &
                     lOAM = .false., lOMQ = .false., lPP, lSkip, lTtl, lXF = .false., MolWgh_UsrDef, nmwarn, NoAMFI, NoDKroll, &
                     NoZMAT, OrigInput, OriginSet, RF_read, RPSet, Skip1, Skip2, SymmSet, ThrInt_UsrDef, Vlct_, Write_BasLib, &
                     WriteZMat
character(len=LenIn) :: CtrLDK(10), dbas
character(len=512) :: Align_Weights = 'MASS'
character(len=256) :: BasLib, Basis_lib, Directory, ExtBasDir, Fname, GeoDir, KeepBasis, Message, Project, temp1, temp2
character(len=180) :: filename, KeepGroup, Key, KWord, Ref(2)
character(len=80) :: BSLbl, ChSkip, Title(10) = ''
character(len=72) :: Header(2) = ''
character(len=14) :: Vrsn = 'Gateway/Seward'
character(len=12) :: Previous_Command
character(len=4) :: CHAR4
character(len=3) :: Oper(3)
character :: AngTyp(0:iTabMx)
integer(kind=iwp), allocatable :: ITmp(:), iScratch(:), nIsot(:,:)
real(kind=wp), allocatable :: Buffer(:), DMSt(:,:), EFt(:,:), Isotopes(:), mIsot(:), OrigRot(:,:,:), OrigTrans(:,:), RTmp(:,:)
character(len=180), allocatable :: STDINP(:)
character(len=128), allocatable :: xb_bas(:)
character(len=12), allocatable :: xb_label(:)
#ifdef _HAVE_EXTRA_
logical(kind=iwp) :: geoInput, oldZmat, zConstraints
#endif
#ifdef _GROMACS_
integer(kind=iwp) :: iCastMM, iLA, LuXYZ, nCastMM, nLA
integer(kind=iwp), allocatable :: CastMM(:), DefLA(:,:)
real(kind=wp), allocatable :: FactLA(:)
#endif
#ifdef _DEBUGPRINT_
#define _TEST_ .true.
#else
#define _TEST_ .false.
#endif
integer(kind=iwp), parameter :: Cho_MolWgh = 2, nBuff = 10000
real(kind=wp), parameter :: Cho_CutInt = 1.0e-40_wp, Cho_ThrInt = 1.0e-40_wp, &
                            WellCff(3) = [0.35_wp,0.25_wp,5.2_wp], &
                            WellExp(3) = [Four,Three,Two], &
                            WellRad(3) = [-1.22_wp,-3.20_wp,-6.20_wp]
logical(kind=iwp), parameter :: IfTest = _TEST_
character(len=*), parameter :: DefNm = 'basis_library'
! Note: blank keywords have been removed and can be reused
character(len=*), parameter :: KeyW(188) = ['END ','EMBE','SYMM','FILE','VECT','ORBC','THRS','UNNO','RADI','TITL','ECPS','AUXS', &
                                            'BSSH','VERB','ORBA','ZMAT','XBAS','XYZ ','COOR','GROU','BSSE','MOVE','NOMO','SYMT', &
                                            'NODE','SDEL','TDEL','BASD','BASI','PRIN','OPTO','THRE','CUTO','RTRN','DIRE','CSPF', &
                                            'EXPE','MOLC','DCRN','MOLP','MOLE','RELI','JMAX','MULT','CENT','EMPC','XFIE','DOUG', &
                                            'DK1H','DK2H','DK3H','DK3F','RESC','RA0H','RA0F','RAIH','RX2C','RBSS','DCCD','BSSM', &
                                            'AMFI','AMF1','AMF2','AMF3','FAKE','FINI','MGAU','PART','FPCO','FPPR','NOTA','WELL', &
                                            'NODK','ONEO','TEST','SDIP','EPOT','EFLD','FLDG','ANGM','UPON','DOWN','OMQI','AMPR', &
                                            'DSHD','NOPA','    ','PKTH','SKIP','EXTR','RF-I','GRID','CLIG','NEMO','RMAT','RMEA', &
                                            'RMER','RMQC','RMDI','RMEQ','RMBP','GIAO','NOCH','CHOL','FCD ','THRC','1CCD','1C-C', &
                                            'CHOI','RP-C','SADD','CELL','SPAN','SPRE','LOW ','MEDI','HIGH','DIAG','RIC ','RIJ ', &
                                            'RIJK','RICD','XRIC','NOGU','RELA','RLOC','FOOC','CDTH','SHAC','KHAC','ACD ','FAT-', &
                                            'ACCD','SLIM','    ','DOFM','NOAM','RPQM','CONS','NGEX','LOCA','    ','    ','    ', &
                                            '    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ','    ', &
                                            'NOAL','WEIG','ALIG','TINK','ORIG','HYPE','ZCON','SCAL','DOAN','GEOE','OLDZ','OPTH', &
                                            'NOON','GEO ','MXTC','FRGM','TRAN','ROT ','ZONL','BASL','NUME','VART','VARR','SHAK', &
                                            'PAMF','GROM','LINK','EMFR','NOCD','FNMC','ISOT','EFP ']
integer(kind=iwp), external :: iCFrst, iChAtm, IsFreeUnit
real(kind=wp), external :: NucExp, rMass, rMassx
character(len=180), external :: Get_Ln

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 3
iPrint = nPrint(iRout)
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
NoAMFI = .false.

iChk_RI = 0
iChk_CH = 0
iChk_DC = 0
iOpt_XYZ = -1

isXfield = 0
CholeskyThr = -huge(CholeskyThr)

Basis_Test = .false.
FinishBasis = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
LuRd = LuRd_
LuFS = -1
LuRdSave = -1

nTemp = 0
lMltpl = .false.

CholeskyWasSet = .false.
spanCD = -huge(spanCD)
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
stepFac1 = 60.0_wp
iOptimType = 1
gradLim = Zero
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

ScaleFactor = One
lSTDINP = 0
iCoord = 0
iBSSE = -1
SymThr = 0.01_wp
nTtl = 0

imix = 0
ifnr = -1
ign = 0
itype = 0
ExtBasDir = ' '
isxbas = 0
nmwarn = .true.

call mma_allocate(xb_bas,MAX_XBAS,label='xb_bas')
call mma_allocate(xb_label,MAX_XBAS,label='xb_label')

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
ispread(:) = 0
VCell(:,:) = Zero
rewind(LuRd)
! Count the number of calls to coord to set the number of fragments
! for a geo/hyper-calculation
do
  Key = Get_Ln(LuRd)
  KWord = Key
  call UpCase(KWord)
  select case (KWord(1:4))
    case (KeyW(19)) ! COOR
      nFragment = nFragment+1
    case (KeyW(162),KeyW(170)) ! HYPE, GEO
      call getenvf('Project',Project)
      call getenvf('GeoDir',GeoDir)
      temp1 = Project(1:index(Project,' ')-1)//'.Gateway.Input'
      temp2 = GeoDir(1:index(GeoDir,' ')-1)//'/'//Project(1:index(Project,' ')-1)//'.gwcopy.input'
      call fCopy(temp1,temp2,ierr)
      if (ierr /= 0) then
        write(u6,*) '*** Detect Hyper input, but no GEO loop'
        call Quit_OnUserError()
      end if
    case (KeyW(1)) ! END
      exit
  end select
end do
call Put_iScalar('nCoordFiles',nFragment)
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

AnyMode = Run_Mode == G_Mode

! GWInput is a logical flag which tells if a keyword is Gateway or
! Seward specific. This is done down below for each keyword:
!
!   GWInput = .true.    ! keyword is Gateway specific
!   GWInput = .false.   ! keyword is Seward specific
!   GWInput = AnyMode   ! keyword is not specific
!
! in the section of the particular keyword
!                                                                      *
!***********************************************************************
!                                                                      *

! KeyWord directed input

call mma_allocate(STDINP,MxAtom*2,label='STDINP')

GWInput = AnyMode

nDone = 0
Skip1 = .false.
Skip2 = .false.
! This loop is here because if the COORD keyword is given, we will need
! to read the generated native BASIS input, and rather than rewrite
! a reduced version of the input parsing code, we just process the file
! as a full input
do
  do
    if (Skip2) then
      Skip2 = .false.
    else
      if (Skip1) then
        Skip1 = .false.
      else
        lTtl = .false.
        if (Basis_Test .and. (nDone == 1)) then
          nDone = 0
          Basis_Test = .false.
        end if
      end if
      Key = Get_Ln(LuRd)
    end if

    if ((Run_Mode == G_Mode) .and. (.not. GWInput)) then
      call WarningMessage(2,'Gateway input error!')
      write(u6,*) 'The keyword : "',Previous_Command,'" is not allowed in the Gateway input!'
      write(u6,*) 'This keyword most likely belongs in the Seward input section!.'
      write(u6,*)
      call Quit_OnUserError()
    else if ((Run_Mode == S_Mode) .and. GWInput) then
      call WarningMessage(2,'Seward input error!')
      write(u6,*) 'The keyword : "',Previous_Command,'" is not allowed in the Seward input when the Gateway is used!'
      write(u6,*) ' Try putting the keyword in the Gateway input section!'
      write(u6,*)
      call Quit_OnUserError()
    end if
    GWInput = .false.

    if (IfTest) write(u6,*) ' RdCtl: Processing:',Key
    KWord = Key
    call UpCase(KWord)
    Previous_Command = KWord(1:4)
    if (KWord(1:1) == '*') cycle
    if (KWord == '') cycle
    if (Basis_Test) nDone = 1

    select case (KWord(1:4))
      case (KeyW(1))
        !                                                              *
        !***** END  ****************************************************
        !                                                              *
        FinishBasis = .true.
        exit

#     ifdef _FDE_
      case (KeyW(2))
        !                                                              *
        !***** EMBE ****************************************************
        !                                                              *
        ! Read in the name of a file for an embedding potential on a grid

        embPot = .true.
        call Put_iScalar('embpot',1)
        do
          KWord = Get_Ln(LuRd)
          select case (KWord(1:4))
            case ('ENDE')
              ! Sanity check
              if (embPotInBasis .and. (.not. outGridPathGiven) .and. &
                  (embWriteDens .or. embWriteESP .or. embWriteGrad .or. embWriteHess)) then
                call WarningMessage(2,' No grid to write output in embedding.')
                call Quit_OnUserError()
              end if
              ! Write the embpot runfile
              call EmbPotWrRun()
              exit
            case ('BASI')
              ! Switch whether the potential is given in basis set
              ! representation instead of grid representation
              embPotInBasis = .true.
              write(u6,*) 'Set embPotInBasis to ',embPotInBasis
            case ('EMBI')
              ! Get the EMBInfile path containing an embedding pot. on a grid
              KWord = Get_Ln(LuRd)
              call Get_S(1,embPotPath,1)
            case ('OUTG')
              ! If the output grid is different from the input grid a path to
              ! a grid file is specified here
              KWord = Get_Ln(LuRd)
              call Get_S(1,outGridPath,1)
              outGridPathGiven = .true.
            case ('WRDE')
              ! If given, the final density is written out on a grid, the path
              ! to the file where to write it to must be given as well.
              embWriteDens = .true.
              KWord = Get_Ln(LuRd)
              call Get_S(1,embOutDensPath,1)
            case ('WREP')
              ! If given, the final electrostatic potential is written out on a
              ! grid, the path to the file where to write it to must be given as well.
              KWord = Get_Ln(LuRd)
              call Get_S(1,embOutEspPath,1)
              embWriteESP = .true.
            case ('WRGR')
              ! If given, the density gradient is written out on a grid, the
              ! path to the file where to write it to must be given as well.
              KWord = Get_Ln(LuRd)
              call Get_S(1,embOutGradPath,1)
              embWriteGrad = .true.
            case ('WRHE')
              ! If given, the Hessian of the density is written out on a grid,
              ! the path to the file where to write it to must be given as well.
              embWriteHess = .true.
              KWord = Get_Ln(LuRd)
              call Get_S(1,embOutHessPath,1)
            case default
              ! Should not be reached.
              call WarningMessage(2,'Error in input of EMBEdding block!')
              call Quit_OnUserError()
          end select
        end do
#     endif

      case (KeyW(3))
        !                                                              *
        !***** SYMM ****************************************************
        !                                                              *
        ! Read distinct symmetry operators apart for the unity operator

        SymmSet = .true.
        KWord = Get_Ln(LuRd)
        GWInput = .true.
        if ((.not. DoneCoord) .and. (iCoord /= 0)) then
          call WarningMessage(2,'SYMMETRY keyword is not compatible with COORD')
          call Quit_OnUserError()
        end if
        if (Run_Mode == S_Mode) then
          call WarningMessage(2,'Seward input error!')
          write(u6,'(A,A,A)') 'The command : "',Previous_Command,'" is not allowed in Seward input when Gateway is used!'
          write(u6,*)
          call Quit_OnUserError()
        end if
        call UpCase(KWord)
        iChrct = len(KWord)
        do
          Last = len_trim(KWord)
          iFrst = iCFrst(KWord,iChrct)
          if (iFrst > Last) exit
          nOper = nOper+1
          Oper(nOper)(1:1) = KWord(iFrst:iFrst)
          KWord(iFrst:iFrst) = ' '
          iFrst = iFrst+1
          if (KWord(iFrst:iFrst) == ' ') cycle
          Oper(nOper)(2:2) = KWord(iFrst:iFrst)
          KWord(iFrst:iFrst) = ' '
          iFrst = iFrst+1
          if (KWord(iFrst:iFrst) == ' ') cycle
          Oper(nOper)(3:3) = KWord(iFrst:iFrst)
          KWord(iFrst:iFrst) = ' '
        end do

      case (KeyW(4))
        !                                                              *
        !***** FILE ****************************************************
        !                                                              *
        ! Specify filename for input orbitals

        Line = Get_Ln(LuRd)
        call FileOrb(Line,SW_FileOrb)

      case (KeyW(5))
        !                                                              *
        !***** VECT ****************************************************
        !                                                              *
        ! Change mode of Seward to property calculations.

        Prprt = .true.

      case (KeyW(6))
        !                                                              *
        !***** ORBC ****************************************************
        !                                                              *
        ! Request property output with explicit listing of orbital contributions.

        Short = .false.

      case (KeyW(7))
        !                                                              *
        !***** THRS ****************************************************
        !                                                              *
        ! Change default for non-zero occupation numbers

        KWord = Get_Ln(LuRd)
        call Get_F1(1,Thrs)
        Thrs = abs(Thrs)

      case (KeyW(8))
        !                                                              *
        !***** UNNO ****************************************************
        !                                                              *
        ! Change default to unnormalized integrals

        UnNorm = .true.

      case (KeyW(9))
        !                                                              *
        !***** RADI ****************************************************
        !                                                              *
        ! Integral cutoff to be based on radial overlap integrals

        lSchw = .false.

      case (KeyW(10))
        !                                                              *
        !***** TITL ****************************************************
        !                                                              *
        ! Read the Title card

        Key = Get_Ln(LuRd)
        call ProcessTitle()

      case (KeyW(11),KeyW(12))
        !                                                              *
        !***** ECPS or AUXS ********************************************
        !                                                              *
        ! Allow printing of ECP data

        nPrint(2) = max(10,nPrint(2))
        GWInput = .true.

      case (KeyW(13))
        !                                                              *
        !***** BSSH ****************************************************
        !                                                              *
        ! Allow printing of basis set data

        nPrint(2) = max(6,nPrint(2))
        GWInput = .true.

      case (KeyW(14))
        !                                                              *
        !***** VERB ****************************************************
        !                                                              *
        ! Verbose printing

        nPrint(2) = max(10,nPrint(2))
        nPrint(117) = 6
        nPrint(80) = 6
        nPrint(1) = 6
        GWInput = AnyMode

      case (KeyW(15))
        !                                                              *
        !***** ORBA ****************************************************
        !                                                              *
        ! Request property output with explicit listing of properties of
        ! all orbitals, including all unoccupied (ignoring THRS), and not
        ! weighted by occupation numbers. (S.S.Dong, 2018)

        ifallorb = .true.

      case (KeyW(16))
        !                                                              *
        !***** ZMAT ****************************************************
        !                                                              *
        ! Read Basis Sets & Coordinates in Z-Matrix format

        if (isxbas == 0) call Quit_OnUserError()
        call ZMatrixConverter(LuRd,u6,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
        if (iErr /= 0) call Quit_OnUserError()
        GWInput = .true.
        call StdSewInput(LuRd,ifnr,mdc,iShll,BasisTypes,STDINP,lSTDINP,iErr)
        if (iErr /= 0) call Quit_OnUserError()

      case (KeyW(17))
        !                                                              *
        !***** XBAS ****************************************************
        !                                                              *

        call read_xbas(LuRd,iglobal,nxbas,xb_label,xb_bas,ierr)
        GWInput = .true.
        isxbas = 1
        if (ierr == 1) call Quit_OnUserError()

      case (KeyW(18))
        !                                                              *
        !***** XYZ  ****************************************************
        !                                                              *

        if (isxbas == 0) call Quit_OnUserError()
        call XMatrixConverter(LuRd,u6,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
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

      case (KeyW(19))
        !                                                              *
        !***** COOR ****************************************************
        !                                                              *
        ! Read Basis Sets & Coordinates in xyz format

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
#       ifdef _HAVE_EXTRA_
        call XYZread(LuRd,ForceZMAT,nCoord,iErr)
        if (iErr /= 0) call Quit_OnUserError()
        call XYZcollect(iCoord,nCoord,OrigTrans,OrigRot,nFragment)
#       else
        call Read_XYZ(LuRd,OrigRot,OrigTrans)
#       endif

      case (KeyW(20))
        !                                                              *
        !***** GROUP ***************************************************
        !                                                              *
        ! Read information for a group

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
        if ((temp1(1:4) /= 'FULL') .and. (temp1(1:1) /= 'E') .and. (temp1(1:2) /= 'C1') .and. (temp1(1:5) /= 'NOSYM')) then
          do i=1,len_trim(temp1)
            if ((temp1(i:i) /= 'X') .and. (temp1(i:i) /= 'Y') .and. (temp1(i:i) /= 'Z') .and. (temp1(i:i) /= ' ')) then
              call WarningMessage(2,'Illegal symmetry group or operator: '//trim(temp1))
              call Quit_OnUserError()
            end if
          end do
        end if
        GroupSet = .true.
        GWInput = .true.
        DoneCoord = .true.

      case (KeyW(21))
        !                                                              *
        !***** BSSE ****************************************************
        !                                                              *
        ! Read information for BSSE

        if (.not. CoordSet) then
          call WarningMessage(2,'COORD keyword is not found')
          call Quit_OnUserError()
        end if
        KWord = Get_Ln(LuRd)
        read(Kword,*,iostat=istatus) iBSSE
        if (istatus /= 0) call Error(1)
        GWInput = .true.

      case (KeyW(22))
        !                                                              *
        !***** MOVE ****************************************************
        !                                                              *
        ! allow to MOVE coordinates

        if (.not. CoordSet) then
          call WarningMessage(2,'COORD keyword is not found')
          call Quit_OnUserError()
        end if
#       ifdef _HAVE_EXTRA_
        isHold = 0
#       endif
        GWInput = .true.

      case (KeyW(23))
        !                                                              *
        !***** NOMOVE **************************************************
        !                                                              *
        ! Do NOT allow to MOVE coordinates

        if (.not. CoordSet) then
          call WarningMessage(2,'COORD keyword is not found')
          call Quit_OnUserError()
        end if
#       ifdef _HAVE_EXTRA_
        isHold = 1
#       endif
        GWInput = .true.

      case (KeyW(24))
        !                                                              *
        !***** SYMT ****************************************************
        !                                                              *
        ! Threshold for findsym

        if (.not. CoordSet) then
          call WarningMessage(2,'COORD keyword is not found')
          call Quit_OnUserError()
        end if
        KWord = Get_Ln(LuRd)
        read(Kword,*,iostat=istatus) SymThr
        if (istatus /= 0) call Error(1)
        GWInput = .true.

      case (KeyW(25))
        !                                                              *
        !***** NODE ****************************************************
        !                                                              *
        ! Set global delete parameters to disable orbital deleting
        ! See also: SDELete, TDELete

        call Put_dScalar('S delete thr',Zero)
        call Put_dScalar('T delete thr',1.0e15_wp)

      case (KeyW(26))
        !                                                              *
        !***** SDEL ****************************************************
        !                                                              *

        KWord = Get_Ln(LuRd)
        call Get_F1(1,sDel)
        call Put_dScalar('S delete thr',sDel)

      case (KeyW(27))
        !                                                              *
        !***** TDEL ****************************************************
        !                                                              *

        KWord = Get_Ln(LuRd)
        call Get_F1(1,tDel)
        call Put_dScalar('T delete thr',tDel)

      case (KeyW(28))
        !                                                              *
        !***** BASD ****************************************************
        !                                                              *

        KWord = Get_Ln(LuRd)
        ExtBasDir = KWord
        GWInput = .true.

      case (KeyW(29))
        !                                                              *
        !***** BASI ****************************************************
        !                                                              *
        ! Read information for a basis set

        ! Check if the format is old or new style. Damn the person who used
        ! the same keword for two different styles of input and making the
        ! input require a specific order of the keyword. Comrade 55?

        Basis_Test = .true.

        GWInput = .true.
        Key = Get_Ln(LuRd)
        if (len_trim(Key) > len(BSLbl)) then
          call WarningMessage(2,'BASIS keyword: line too long')
          call Quit_OnUserError()
        end if
        BSLbl = Key(1:len(BSLbl))
        if (BasisSet) then
          ! iOpt_XYZ=1 means we're dealing with "new style" and we need to accumulate all the labels
          if ((.not. FinishBasis) .and. (iOpt_XYZ == 1)) then
            if (len_trim(KeepBasis)+len_trim(BsLbl) > len(KeepBasis)-1) then
              call WarningMessage(2,'BASIS keyword: total basis too long')
              call Quit_OnUserError()
            end if
            KeepBasis = trim(KeepBasis)//','//BsLbl
          end if
        else
          KeepBasis = BSLbl
          BasisSet = .true.
        end if
        !temp1 = KeepBasis
        !call UpCase(temp1)
        !if (INDEX(temp1,'INLINE') /= 0) then
        !  write(u6,*) 'XYZ input and Inline basis set are not compatible'
        !  write(u6,*) 'Consult the manual how to change inline basis set'
        !  write(u6,*) ' into basis set library'
        !  call Quit_OnUserError()
        !end if
        iOpt_XYZ = 1

      case (KeyW(30))
        !                                                              *
        !***** PRIN ****************************************************
        !                                                              *
        ! Print level

        KWord = Get_Ln(LuRd)
        GWInput = AnyMode
        call Get_I1(1,n)
        do i=1,n
          KWord = Get_Ln(LuRd)
          call Get_I1(1,jRout)
          call Get_I1(2,iPrint)
          nPrint(jRout) = iPrint
        end do

      case (KeyW(31))
        !                                                              *
        !***** OPTO ****************************************************
        !                                                              *
        ! Reduce the output for optimizations

        lOPTO = .true.

      case (KeyW(32))
        !                                                              *
        !***** THRE ****************************************************
        !                                                              *
        ! Threshold for writing integrals to disk

        KWord = Get_Ln(LuRd)
        call Get_F1(1,ThrInt)
        ThrInt = abs(ThrInt)
        ThrInt_UsrDef = .true.

      case (KeyW(33))
        !                                                              *
        !***** CUTO ****************************************************
        !                                                              *
        ! Cutoff for computing primitive integrals [a0|c0]

        KWord = Get_Ln(LuRd)
        call Get_F1(1,CutInt)
        CutInt = abs(CutInt)
        CutInt_UsrDef = .true.

      case (KeyW(34))
        !                                                              *
        !***** RTRN ****************************************************
        !                                                              *
        ! Defining max bond distance for bonds, angles and dihedrals
        ! Define max number of atoms to list for

        KWord = Get_Ln(LuRd)
        call Upcase(KWord)
        call Get_I1(1,S%Max_Center)
        call Get_F1(2,rtrnc)
        if (index(KWord,'ANGSTROM') /= 0) Rtrnc = Rtrnc/Angstrom
        GWInput = .true.

      case (KeyW(35))
        !                                                              *
        !***** DIRE ****************************************************
        !                                                              *
        ! Force direct calculations & disable two-electron integrals

        DirInt = .true.
        Onenly = .true.
        iChk_DC = 1
        if ((iChk_RI+iChk_CH) > 0) then
          call WarningMessage(2,'Direct is incompatible with RI and Cholesky keywords')
          call Quit_OnUserError()
        end if
        Do_RI = .false.

      case (KeyW(36))
        !                                                              *
        !***** CSPF ****************************************************
        !                                                              *
        ! Turn on the use of Condon-Shortley phase factors

        CSPF = .true.
        GWInput = AnyMode

      case (KeyW(37))
        !                                                              *
        !***** EXPE ****************************************************
        !                                                              *
        ! Expert mode

        Expert = .true.
        GWInput = AnyMode
        call WarningMessage(1,' EXPERT option is ON!')

      case (KeyW(38),KeyW(39))
        !                                                              *
        !***** MOLC or DCRN ********************************************
        !                                                              *
        ! Weight for DCR summation (this is the default for conventional
        ! calculations)

        MolWgh = 0
        MolWgh_UsrDef = .true.

      case (KeyW(40))
        !                                                              *
        !***** MOLP ****************************************************
        !                                                              *
        ! Weight for DCR summation modified to MOLPRO format

        MolWgh = 2
        MolWgh_UsrDef = .true.

      case (KeyW(41))
        !                                                              *
        !***** MOLE ****************************************************
        !                                                              *
        ! Weight for DCR summation modified to MOLECULE format

        MolWgh = 1
        MolWgh_UsrDef = .true.

      case (KeyW(42))
        !                                                              *
        !***** RELI ****************************************************
        !                                                              *
        ! Compute integrals for first order relativistic corrections
        ! of the energy, i.e. the mass-velocity integrals and the
        ! one-electron Darwin contract term integrals.

        lRel = .true.

      case (KeyW(43))
        !                                                              *
        !***** JMAX ****************************************************
        !                                                              *
        ! Change max j quantum number for the rigid rotor analysis

        KWord = Get_Ln(LuRd)
        call Get_I1(1,S%jMax)

      case (KeyW(44))
        !                                                              *
        !***** MULT ****************************************************
        !                                                              *
        ! Read order of highest multipole to be computed

        KWord = Get_Ln(LuRd)
        call Get_I1(1,S%nMltpl)

      case (KeyW(45))
        !                                                              *
        !***** CENT ****************************************************
        !                                                              *
        ! User specified centers of multipole moment operators.

        KWord = Get_Ln(LuRd)
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
          call Get_F(2,RTmp(:,i),3)
          if (index(KWord,'ANGSTROM') /= 0) RTmp(:,i) = RTmp(:,i)/Angstrom
          ITmp(i) = iMltpl
        end do

      case (KeyW(46))
        !                                                              *
        !***** EMPC ****************************************************
        !                                                              *
        ! Compute Orbital-Free Embedding integrals from Point Charges
        ! specified in XFIEld

        DoEmPC = .true.
#       ifdef _HAVE_EXTRA_
        isHold = 1 ! avoid coordinate moving
#       endif
        GWInput = .true.

      case (KeyW(47))
        !                                                              *
        !***** XFIE ****************************************************
        !                                                              *
        ! User specified external field

        lXF = .true.
        GWInput = .true.
        KWord = Get_Ln(LuRd)
        ! Open external file if the line does not start with an integer
        ! Note that the "iostat" value cannot be completely trusted, since
        ! a slash (i.e., an absolute path) marks end of input and gives
        ! no error
        LuRd_saved = LuRd
        ibla = -1
        read(KWord,*,iostat=istatus) ibla
        if ((ibla < 0) .or. (istatus > 0)) then
          LuRd = 1
          call Get_S(1,filename,1)
          call OpnFl(filename(1:(index(filename,' ')-1)),LuRd,Exists)
          if (.not. Exists) then
            call WarningMessage(2,'Error! File not found: '//filename(1:(index(filename,' ')-1)))
            call Quit_OnUserError()
          end if
          write(u6,*) 'Reading external field from file: ',filename(1:(index(filename,' ')-1))
          KWord = Get_Ln(LuRd)
        end if
        call ProcessXF()

      case (KeyW(48))
        !                                                              *
        !***** DOUG ****************************************************
        !                                                              *
        ! Full Douglas-Kroll operator

        IRELAE = 0
        DKroll = .true.
        GWInput = .true.

      case (KeyW(49))
        !                                                              *
        !***** DK1H ****************************************************
        !                                                              *
        ! Full Douglas-Kroll (DK1) operator

        IRELAE = 1
        DKroll = .true.
        GWInput = .true.

      case (KeyW(50))
        !                                                              *
        !***** DK2H ****************************************************
        !                                                              *
        ! Full Douglas-Kroll (DK2) operator

        IRELAE = 2
        DKroll = .true.
        GWInput = .true.

      case (KeyW(51))
        !                                                              *
        !***** DK3H ****************************************************
        !                                                              *
        ! Douglas-Kroll (DK3) operator

        IRELAE = 3
        DKroll = .true.
        GWInput = .true.

      case (KeyW(52))
        !                                                              *
        !***** DK3F ****************************************************
        !                                                              *
        ! Full Douglas-Kroll (DK3) operator

        IRELAE = 4
        DKroll = .true.
        GWInput = .true.

      case (KeyW(53))
        !                                                              *
        !***** RESC ****************************************************
        !                                                              *
        ! Full RESC operator

        IRELAE = 11
        DKroll = .true.
        GWInput = .true.

      case (KeyW(54))
        !                                                              *
        !***** RA0H ****************************************************
        !                                                              *
        ! Full ZORA

        IRELAE = 21
        DKroll = .true.
        GWInput = .true.

      case (KeyW(55))
        !                                                              *
        !***** RA0F ****************************************************
        !                                                              *
        ! Full ZORA-FP

        IRELAE = 22
        DKroll = .true.
        GWInput = .true.

      case (KeyW(56))
        !                                                              *
        !***** RAIH ****************************************************
        !                                                              *
        ! Full IORA

        IRELAE = 23
        DKroll = .true.
        GWInput = .true.

      case (KeyW(57))
        !                                                              *
        !***** RX2C ****************************************************
        !                                                              *
        ! Exact decoupling X2C method

        IRELAE = 101
        DKroll = .true.
        GWInput = .true.

      case (KeyW(58))
        !                                                              *
        !***** RBSS ****************************************************
        !                                                              *
        ! Exact decoupling BSS method

        IRELAE = 102
        DKroll = .true.
        GWInput = .true.

      case (KeyW(59))
        !                                                              *
        !**** DCCD *****************************************************
        !                                                              *

        Do_DCCD = .true.

        ! RICD
        Do_RI = .true. !ORDINT ERROR
        GWInput = .true.
        if (iChk_RI == 0) then
          call WarningMessage(2,'DCCD option set without RI type defined.')
          call Quit_OnUserError()
        end if

        ! DIRE
        DirInt = .true.

      case (KeyW(60))
        !                                                              *
        !***** BSSM ****************************************************
        !                                                              *
        ! BSS method

        IRELAE = 0
        DKroll = .true.
        BSS = .true.
        GWInput = .true.

      case (KeyW(61))
        !                                                              *
        !***** AMFI ****************************************************
        !                                                              *
        ! AMFI integrals

        lAMFI = .true.
        GWInput = .true.

      case (KeyW(62))
        !                                                              *
        !***** AMF1 ****************************************************
        !                                                              *
        ! AMFI integrals

        lAMFI = .true.
        GWInput = .true.

      case (KeyW(63))
        !                                                              *
        !***** AMF2 ****************************************************
        !                                                              *
        ! AMFI integrals (including 2nd-order)

        lAMFI = .true.
        GWInput = .true.

      case (KeyW(64))
        !                                                              *
        !***** AMF3 ****************************************************
        !                                                              *
        ! AMFI integrals (including 3rd-order)

        lAMFI = .true.
        GWInput = .true.

      case (KeyW(65))
        !                                                              *
        !***** FAKE ****************************************************
        !                                                              *
        ! Fake run : conventional, RI or CD ERIs not computed
        !            but some info is set to runfile (e.g. CD thrs)
        !            Not the same as ONEOnly !!

        Fake_ERIs = .true.

      case (KeyW(66))
        !                                                              *
        !***** FINI ****************************************************
        !                                                              *
        ! Finite nuclei - Gaussian type

        Nuclear_Model = Gaussian_Type
        GWInput = .true.

      case (KeyW(67))
        !                                                              *
        !***** MGAU ****************************************************
        !                                                              *
        ! Finite nuclei - modified Gaussian type

        Nuclear_Model = mGaussian_Type
        GWInput = .true.

      case (KeyW(68))
        !                                                              *
        !***** PART ****************************************************
        !                                                              *
        ! Show partitioning statistics

        nPrint(10) = 6

      case (KeyW(69))
        !                                                              *
        !***** FPCO ****************************************************
        !                                                              *
        ! Force partitioning for contracted functions

        force_part_c = .true.

      case (KeyW(70))
        !                                                              *
        !***** FPPR ****************************************************
        !                                                              *
        ! Force partitioning for primitive functions

        force_part_p = .true.

      case (KeyW(71))
        !                                                              *
        !***** NOTA ****************************************************
        !                                                              *
        ! Do not use tables for the roots and weights of the Rys polynomials.

        NoTab = .true.

      case (KeyW(72))
        !                                                              *
        !***** WELL ****************************************************
        !                                                              *
        ! Read radius and exponents for spherical well integrals and coefficient

        KWord = Get_Ln(LuRd)
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
              Wel_Info(1,iWel) = Wel_Info(1,iWel)/Angstrom
              Wel_Info(2,iWel) = Wel_Info(2,iWel)*Angstrom
            end if
          end do
        end if

      case (KeyW(73))
        !                                                              *
        !***** NODK ****************************************************
        !                                                              *
        ! Do not compute Douglas-Kroll integrals.

        NoDKroll = .true.

      case (KeyW(74))
        !                                                              *
        !***** ONEO ****************************************************
        !                                                              *
        ! Do not compute two electron integrals.

        Onenly = .true.

      case (KeyW(75))
        !                                                              *
        !***** TEST ****************************************************
        !                                                              *
        ! Process only the input.

        Test = .true.
        GWInput = AnyMode

      case (KeyW(76))
        !                                                              *
        !***** SDIP ****************************************************
        !                                                              *
        ! Compute integrals for transition dipole moment

        Vlct_ = .true.
        GWInput = .true.

      case (KeyW(77))
        !                                                              *
        !***** EPOT ****************************************************
        !                                                              *
        ! Compute the electric potential for a number of points.  If nEF is
        ! set to 0 this will cause the points to coincide with the unique
        ! centers.

        nOrdEF = max(nOrdEF,0)
        GWInput = .true.
        if (EFgiven) then
          call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
          call Quit_OnUserError()
        end if
        EFgiven = .true.
        call ProcessEF(KWord(1:4))

      case (KeyW(78))
        !                                                              *
        !***** EFLD ****************************************************
        !                                                              *
        ! Compute the electric potential and electric field for a number of
        ! points. If nEF is set to 0 this will cause the points to coincide
        ! with the unique centers.

        nOrdEF = max(nOrdEF,1)
        GWInput = .true.
        if (EFgiven) then
          call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
          call Quit_OnUserError()
        end if
        EFgiven = .true.
        call ProcessEF(KWord(1:4))

      case (KeyW(79))
        !                                                              *
        !***** FLDG ****************************************************
        !                                                              *
        ! Compute the electric potential, electric field, and electric field
        ! gradient for a number of points. If nEF is set to 0 this will
        ! cause the points to coincide with the Unique centers.

        nOrdEF = max(nOrdEF,2)
        GWInput = .true.
        if (EFgiven) then
          call WarningMessage(2,'Only one of EPOT,EFLD,FLDG may be given')
          call Quit_OnUserError()
        end if
        EFgiven = .true.
        call ProcessEF(KWord(1:4))

      case (KeyW(80))
        !                                                              *
        !***** ANGM ****************************************************
        !                                                              *
        ! Orbital angular momentum

        lOAM = .true.
        GWInput = .true.
        KWord = Get_Ln(LuRd)
        call Upcase(KWord)
        call Get_F(1,OAMt,3)
        if (index(KWord,'ANGSTROM') /= 0) OAMt(:) = OAMt/Angstrom

      case (KeyW(81))
        !                                                              *
        !***** UPON ****************************************************
        !                                                              *
        ! Orbital angular momentum restriction

        lUPONLY = .true.

      case (KeyW(82))
        !                                                              *
        !***** DOWN ****************************************************
        !                                                              *
        ! Orbital angular momentum restriction

        lDOWNONLY = .true.

      case (KeyW(83))
        !                                                              *
        !***** OMQI ****************************************************
        !                                                              *
        ! Orbital magnetic quadrupole

        lOMQ = .true.
        GWInput = .true.
        KWord = Get_Ln(LuRd)
        call Upcase(KWord)
        call Get_F(1,OMQt,3)
        if (index(KWord,'ANGSTROM') /= 0) OMQt(:) = OMQt/Angstrom

      case (KeyW(84))
        !                                                              *
        !***** AMPR ****************************************************
        !                                                              *
        ! Angular momentum products

        GWInput = .true.
        if (Run_Mode == S_Mode) then
          Skip2 = .true.
        else
          call mma_allocate(AMP_Center,3,Label='AMP_Center')
          KWord = Get_Ln(LuRd)
          call Upcase(KWord)
          call Get_F(1,AMP_Center,3)
          if (index(KWord,'ANGSTROM') /= 0) AMP_Center(:) = (One/Angstrom)*AMP_Center(:)
        end if

      case (KeyW(85))
        !                                                              *
        !***** DSHD ****************************************************
        !                                                              *
        ! Compute the diamagnetic shielding for a number of points. If nDMS
        ! is set to 0 this will cause the points to coincide with the
        ! unique centers.

        lDMS = .true.
        GWInput = .true.
        KWord = Get_Ln(LuRd)
        call Get_F(1,Dxyz,3)
        KWord = Get_Ln(LuRd)
        call Get_I1(1,nDMS)
        if (nDMS < 0) nDMS = 0
        if (nDMS /= 0) then
          call mma_allocate(DMSt,3,nDMS,label='DMSt')
          do iDMS=1,nDMS
            KWord = Get_Ln(LuRd)
            call Upcase(KWord)
            call Get_F(1,DMSt(1,iDMS),3)
            if (index(KWord,'ANGSTROM') /= 0) DMSt(:,iDMS) = DMSt(:,iDMS)/Angstrom
          end do
        end if

      case (KeyW(86))
        !                                                              *
        !***** NOPA ****************************************************
        !                                                              *
        ! Set integral packing flag
        ! iPack=0   : pack 2el integrals (= Default)
        ! iPack=1   : do not pack 2el integrals

        iPack = 1

      case (KeyW(88))
        !                                                              *
        !***** PKTH ****************************************************
        !                                                              *
        ! Read desired packing accuracy ( Default = 1.0e-14 )

        KWord = Get_Ln(LuRd)
        call Get_F1(1,PkAcc)
        PkAcc = abs(PkAcc)

      case (KeyW(89))
        !                                                              *
        !***** SKIP ****************************************************
        !                                                              *
        ! Read skip parameters,i.e.,
        ! if 2el integral symmetry blocks containing a given symmetry
        ! will not be needed in subsequent calculations their computation
        ! ans storage can be omitted.
        ! ( Default = 0,0,0,0,0,0,0,0 )

        KWord = Get_Ln(LuRd)
        lSkip = .true.
        ChSkip = KWord(1:80)

      case (KeyW(90))
        !                                                              *
        !***** EXTR ****************************************************
        !                                                              *
        ! Put the program name and the time stamp onto the extract file

        write(u6,*) 'RdCtl: keyword EXTRACT is obsolete and is ignored!'

      case (KeyW(91))
        !                                                              *
        !***** RF-I ****************************************************
        !                                                              *
        ! Read reaction field input.

        GWInput = .true.
        if (.not. RF_read) then
          call InpRct(LuRd)
          if (.not. (lLangevin .or. PCM)) then

            ! Add a center corresponding to the center of the RF cavity.

            RF_read = .true.
          end if
        else
          call WarningMessage(2,'RdCtl: A second RF-input block discovered!')
          call Quit_OnUserError()
        end if

      case (KeyW(92))
        !                                                              *
        !**** GRID *****************************************************
        !                                                              *

        call Funi_input(LuRd)

      case (KeyW(93))
        !                                                              *
        !***** CLIG ****************************************************
        !                                                              *
        ! Speed of light (in au)

        KWord = Get_Ln(LuRd)
        call Get_F1(1,cLightAU)
        cLightAU = abs(cLightAU)
        write(u6,*) 'The speed of light in this calculation =',cLightAU

      case (KeyW(94))
        !                                                              *
        !**** NEMO *****************************************************
        !                                                              *

        NEMO = .true.

      case (KeyW(95))
        !                                                              *
        !**** RMAT *****************************************************
        !                                                              *
        ! RmatR    : radius of the R-matrix sphere (bohr)

        KWord = Get_Ln(LuRd)
        call Get_F1(1,RMatR)

      case (KeyW(96))
        !                                                              *
        !**** RMEA *****************************************************
        !                                                              *
        ! Epsabs   : absolute precision of numerical radial integration

        KWord = Get_Ln(LuRd)
        call Get_F1(1,Epsabs)

      case (KeyW(97))
        !                                                              *
        !**** RMER *****************************************************
        !                                                              *
        ! Epsrel   : relative precision of numerical radial integration

        KWord = Get_Ln(LuRd)
        call Get_F1(1,Epsrel)

      case (KeyW(98))
        !                                                              *
        !**** RMQC *****************************************************
        !                                                              *
        ! qCoul    : effective charge of the target molecule

        KWord = Get_Ln(LuRd)
        call Get_F1(1,qCoul)

      case (KeyW(99))
        !                                                              *
        !**** RMDI *****************************************************
        !                                                              *
        ! dipol(3) : effective dipole moment of the target molecule

        KWord = Get_Ln(LuRd)
        call Get_F(1,dipol,3)
        dipol1 = abs(dipol(1))+abs(dipol(2))+abs(dipol(3))

      case (KeyW(100))
        !                                                              *
        !**** RMEQ *****************************************************
        !                                                              *
        ! epsq     : minimal value of qCoul and/or dipol1 to be considered

        KWord = Get_Ln(LuRd)
        call Get_F1(1,epsq)

      case (KeyW(101))
        !                                                              *
        !**** RMBP *****************************************************
        !                                                              *
        ! bParm    : Bloch term parameter

        KWord = Get_Ln(LuRd)
        call Get_F1(1,bParm)

      case (KeyW(102))
        !                                                              *
        !**** GIAO *****************************************************
        !                                                              *
        ! Enable GIAO integrals.

        GIAO = .true.

      case (KeyW(103))
        !                                                              *
        !**** NOCH *****************************************************
        !                                                              *
        ! Deactivate Cholesky decomposition.

        Chol = .false.
        CholeskyWasSet = .true.
        Do_RI = .false.

      case (KeyW(104),KeyW(105))
        !                                                              *
        !**** CHOL *****************************************************
        !                                                              *
        ! Activate Cholesky decomposition with default settings.
        ! This section can only be executed once.

        Do_RI = .false.
        if (.not. CholeskyWasSet) then
          CholeskyWasSet = .true.
          Chol = .true.
          Do_RI = .false.
          DirInt = .true.
          call Cho_Inp(.true.,-1,u6)
          iChk_CH = 1
        end if
        if ((iChk_RI+iChk_DC) > 0) then
          call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
          call Quit_OnUserError()
        end if
        GWInput = .false.

      case (KeyW(106))
        !                                                              *
        !**** THRC *****************************************************
        !                                                              *
        ! Set Cholesky decomposition threshold to specified value.

        KWord = Get_Ln(LuRd)
        call Get_F1(1,CholeskyThr)

      case (KeyW(107),KeyW(108))
        !                                                              *
        !**** 1CCD or 1C-C *********************************************
        !                                                              *
        ! Use one-center Cholesky.

        Do_RI = .false.
        iChk_Ch = 1
        DirInt = .true.
        call Cho_Inp(.true.,-1,u6)
        Chol = .true.
        call Cho_InpMod('1CCD')
        if ((iChk_RI+iChk_DC) > 0) then
          call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(109))
        !                                                              *
        !**** CHOI *****************************************************
        !                                                              *
        ! Activate Cholesky decomposition with user-defined settings.
        ! This section can be executed any number of times.
        ! User-defined settings will be preferred to defaults also if the
        ! keywords appear in "wrong" order,
        !
        ! ChoInput
        ! ...
        ! End ChoInput
        ! Cholesky

        Do_RI = .false.
        CholeskyWasSet = .true.
        Chol = .true.
        DirInt = .true.
        call Cho_Inp(.false.,LuRd,u6)
        iChk_CH = 1
        if ((iChk_DC) > 0) then
          call WarningMessage(2,'Cholesky is incompatible with RI and Direct keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(110))
        !                                                              *
        !**** RP-C *****************************************************
        !                                                              *

        lRP = .true.
        KWord = Get_Ln(LuRd)
        nRP_prev = -1
        read(KWord,*,iostat=istatus) nRP
        if (istatus > 0) then
          isnumber = .false.
        else
          isnumber = .true.
          ifile = index(KWord,' ')
          do i=1,ifile
            ii = index(' 0123456789',KWord(i:i))
            if (ii == 0) then
              isnumber = .false.
            end if
          end do
        end if

        if (isnumber) then
          nRP = 3*nRP
          call mma_allocate(RP_Centers,3,nRP/3,2,Label='RP_Centers')

          ! Inline input

          call UpCase(KWord)
          if (index(KWord,'ANGSTROM') /= 0) then
            Fact = One/Angstrom
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
        else

          ! Files

          do jTmp=1,2
            RPSet = .true.
            ifile = index(KWord,' ')
            if (KWord(1:1) == '/') then
              call f_inquire(KWord(1:ifile-1),Exists)
              Key = KWord
            else
              call getenvf('MOLCAS_SUBMIT_DIR',Directory)
              if (Directory(1:1) /= ' ') then
                i = index(Directory,' ')
                Key = Directory(1:i-1)//'/'//KWord(1:ifile-1)
                ifile = i+ifile
                call f_inquire(Key(1:iFile-1),Exists)
              else
                Exists = .false.
              end if
              if (.not. Exists) then
                Key = Key(i+1:iFile-1)
                ifile = ifile-i
                call f_inquire(Key(1:iFile-1),Exists)
              end if
            end if
            if (.not. Exists) then
              call WarningMessage(2,'File '//Key(1:ifile)//' is not found')
              call Quit_OnUserError()
            end if
            LuIn = isFreeUnit(8)
            call molcas_open(LuIn,Key(1:iFile-1))

            KWord = Get_Ln(LuIn)
            read(KWord,*,iostat=istatus) nRP
            if (istatus > 0) call Error(2)
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
              Fact = One/Angstrom
            end if
            LuRP = isFreeUnit(10)
            if (jTmp == 1) then
              call molcas_open(LuRP,'findsym.RP1')
              read(KWord,*,iostat=istatus) E1
            else
              call molcas_open(LuRP,'findsym.RP2')
              read(KWord,*,iostat=istatus) E2
            end if
            if (istatus > 0) call Error(2)

            ! write a separate file for findsym

            write(LuRP,*) nRP/3
#          ifdef _HAVE_EXTRA_
            write(LuRP,'(a)')
#          else
            write(LuRP,'(a)') 'bohr'
#          endif
            do i=1,nRP/3
              KWord = Get_Ln(LuIn)
              read(KWord,*,iostat=istatus) Key,(RP_Centers(j,i,jTmp),j=1,3)
              if (istatus > 0) call Error(2)
              write(LuRP,'(A,3F20.12)') Key(1:LenIn),(RP_Centers(j,i,jTmp)*Fact,j=1,3)
            end do
            close(LuRP)
            if (jTmp == 1) then
              KWord = Get_Ln(LuRd)
              close(LuIn)
            end if
          end do
          RP_Centers(:,:,:) = Fact*RP_Centers(:,:,:)

          close(LuIn)
          GWInput = AnyMode
        end if

      case (KeyW(111))
        !                                                              *
        !**** SADD *****************************************************
        !                                                              *
        ! Saddle options

        Key = Get_Ln(LuRd)
        call Get_F1(1,SadStep)
        GWInput = .true.

      case (KeyW(112))
        !                                                              *
        !**** CELL *****************************************************
        !                                                              *
        ! VCell(3,3)    : the vectors of the cell

        Key = Get_Ln(LuRd)
        call Upcase(Key)
        if (index(Key,'ANGSTROM') /= 0) KWord = Get_Ln(LuRd)
        call Get_F(1,VCell(1,1),3)
        KWord = Get_Ln(LuRd)
        call Get_F(1,VCell(1,2),3)
        KWord = Get_Ln(LuRd)
        call Get_F(1,VCell(1,3),3)
        if (index(Key,'ANGSTROM') /= 0) VCell(:,:) = VCell/Angstrom
        Cell_l = .true.
        call mma_allocate(AdCell,MxAtom)

      case (KeyW(113))
        !                                                              *
        !**** SPAN *****************************************************
        !                                                              *
        ! Set span factor in Cholesky decomposition (0 < span < 1).
        ! The span decides the smallest diagonal element that can be
        ! treated as span*max(Diag). Span=1 thus implies full pivoting.

        KWord = Get_Ln(LuRd)
        call Get_F1(1,spanCD)
        spanCD = abs(spanCD)

      case (KeyW(114))
        !                                                              *
        !**** SPREAD ***************************************************
        !                                                              *
        ! ispread(3)    : the number of cells to spread in different directions

        KWord = Get_Ln(LuRd)
        call Get_I(1,ispread,3)

      case (KeyW(115))
        !                                                              *
        !**** LOW  *****************************************************
        !                                                              *
        ! Activate low-accuracy Cholesky decomposition.

        Do_RI = .false.
        if (.not. CholeskyWasSet) then
          CholeskyWasSet = .true.
          Chol = .true.
          DirInt = .true.
          call Cho_Inp(.true.,-1,u6)
          call Cho_InpMod('LOW ')
          Thrshld_CD = 1.0e-4_wp
        end if

      case (KeyW(116))
        !                                                              *
        !**** MEDI *****************************************************
        !                                                              *
        ! Activate medium-accuracy Cholesky decomposition.

        Do_RI = .false.
        if (.not. CholeskyWasSet) then
          CholeskyWasSet = .true.
          Chol = .true.
          DirInt = .true.
          call Cho_Inp(.true.,-1,u6)
          call Cho_InpMod('MEDI')
          Thrshld_CD = 1.0e-6_wp
        end if

      case (KeyW(117))
        !                                                              *
        !**** HIGH *****************************************************
        !                                                              *
        ! Activate high-accuracy Cholesky decomposition.

        Do_RI = .false.
        if (.not. CholeskyWasSet) then
          CholeskyWasSet = .true.
          Chol = .true.
          DirInt = .true.
          call Cho_Inp(.true.,-1,u6)
          call Cho_InpMod('HIGH')
          Thrshld_CD = 1.0e-8_wp
        end if

      case (KeyW(118))
        !                                                              *
        !**** DIAG *****************************************************
        !                                                              *

        DiagCheck = .true.

      case (KeyW(119))
        !                                                              *
        !**** RIC  *****************************************************
        !                                                              *
        ! Activate RI approach

        Do_RI = .true.
        GWInput = .true.
        iRI_Type = 3
        if (iChk_RI == 1) then
          call WarningMessage(2,'RI basis already defined.')
          call Quit_OnUserError()
        end if
        iChk_RI = 1
        if ((iChk_DC+iChk_CH) > 0) then
          call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(120))
        !                                                              *
        !**** RIJ  *****************************************************
        !                                                              *

        Do_RI = .true.
        GWInput = .true.
        iRI_Type = 1
        if (iChk_RI == 1) then
          call WarningMessage(2,'RI basis already defined.')
          call Quit_OnUserError()
        end if
        iChk_RI = 1
        if ((iChk_DC+iChk_CH) > 0) then
          call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(121))
        !                                                              *
        !**** RIJK *****************************************************
        !                                                              *
        Do_RI = .true.
        GWInput = .true.
        iRI_Type = 2
        if (iChk_RI == 1) then
          call WarningMessage(2,'RI basis already defined.')
          call Quit_OnUserError()
        end if
        iChk_RI = 1
        if ((iChk_DC+iChk_CH) > 0) then
          call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(122))
        !                                                              *
        !**** RICD *****************************************************
        !                                                              *
        Do_RI = .true.
        GWInput = .true.
        iRI_Type = 4
        if (iChk_RI == 1) then
          call WarningMessage(2,'RI basis already defined.')
          call Quit_OnUserError()
        end if
        iChk_RI = 1
        if ((iChk_DC+iChk_CH) > 0) then
          call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(123))
        !                                                              *
        !**** XRIC *****************************************************
        !                                                              *
        Do_RI = .true.
        GWInput = .true.
        iRI_Type = 5
        if (iChk_RI == 1) then
          call WarningMessage(2,'RI basis already defined.')
          call Quit_OnUserError()
        end if
        iChk_RI = 1
        if ((iChk_DC+iChk_CH) > 0) then
          call WarningMessage(2,'RI is incompatible with Direct and Cholesky keywords')
          call Quit_OnUserError()
        end if

      case (KeyW(124))
        !                                                              *
        !**** NOGU *****************************************************
        !                                                              *
        ! Disable automatic execution of GuessOrb

        Do_GuessOrb = .false.

      case (KeyW(125))
        !                                                              *
        !**** RELA *****************************************************
        !                                                              *
        ! DKH option: order and parameterization.
        ! xx: order of Hamiltonian
        !  y: parameterization
        ! zz: order of properties

        kWord = Get_Ln(LuRd)
        if ((KWord(1:1) == 'R') .and. &
            ((KWord(2:2) >= '0') .and. (KWord(2:2) <= '9')) .and. ((KWord(3:3) >= '0') .and. (KWord(3:3) <= '9')) .and. &
            ((KWord(4:4) == 'O') .or. (KWord(4:4) == 'E') .or. (KWord(4:4) == 'S') .or. (KWord(4:4) == 'M') .or. &
             (KWord(4:4) == 'C'))) then
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

      case (KeyW(126))
        !                                                              *
        !**** RLOC *****************************************************
        !                                                              *
        ! Local Douglas-Kroll-Hess/X2C/BSS

        LDKroll = .true.
        !GWInput = .true.
        nCtrLD = 0
        radiLD = 5.5_wp

        KWord = Get_Ln(LuRd)
        call Upcase(KWord)
        if (KWord(1:3) == 'DLU') then
          !continue
        else if (KWord(1:3) == 'DLH') then
          radiLD = Zero
        else
          read(KWord,*,iostat=istatus) nCtrLD,radiLD
          if (istatus == 0) then
            if (nCtrLD > 10) then
              call WarningMessage(2,'The number of centers for LDKH is limited to 10')
              call Quit_OnUserError()
            end if
            if (index(KWord,'ANGSTROM') /= 0) then
              radiLD = radiLD/Angstrom
            end if

            KWord = Get_Ln(LuRd)
            call Upcase(KWord)
            read(KWord,*,iostat=istatus) (iCtrLD(i),i=1,nCtrLD)
            if (istatus < 0) then
              call Error(1)
            else if (istatus > 0) then
              read(Kword,*,iostat=istatus) (CtrLDK(i),i=1,nCtrLD)
              if (istatus /= 0) call Error(1)
              call Get_nAtoms_all(nAtom)
              k = 0
              do i=1,nAtom
                do j=1,nCtrLD
                  if (CtrLDK(j) == dc(i)%LblCnt(1:LenIn)) then
                    iCtrLD(j) = i
                    k = k+1
                  end if
                end do
              end do
              if (k /= nCtrLD) then
                call WarningMessage(2,'Error in LDKH Centers definitions')
                call Quit_OnUserError()
              end if
            end if
          else

            ! Automatic choice: all heavy elements (from K 19)

            !DP if (nCtrLD == 0) radiLD = Zero
            Key = KWord
            Skip2 = .true.
          end if
        end if

      case (KeyW(127))
        !                                                              *
        !**** FOOC *****************************************************
        !                                                              *
        ! Force the use of the out-of-core RI algorithm.

        Force_Out_of_Core = .true.

      case (KeyW(128))
        !                                                              *
        !**** CDTH *****************************************************
        !                                                              *
        ! Threshold for CD to generate RICD auxiliary basis sets

        Key = Get_Ln(LuRd)
        call Get_F1(1,Thrshld_CD)
        GWInput = .true.

      case (KeyW(129))
        !                                                              *
        !**** SHAC *****************************************************
        !                                                              *
        ! Skip high angular combinations when constructing RICD aux basis.

        Skip_High_AC = .true.
        GWInput = .true.

      case (KeyW(130))
        !                                                              *
        !**** KHAC *****************************************************
        !                                                              *
        ! Keep high angular combinations when constructing RICD aux basis.

        Skip_High_AC = .false.
        GWInput = .true.

      case (KeyW(131),KeyW(132))
        !                                                              *
        !**** ACD  or FAT- *********************************************
        !                                                              *
        ! Generate a aCD basis.

        Do_acCD_Basis = .false.
        GWInput = .true.

      case (KeyW(133),KeyW(134))
        !                                                              *
        !**** ACCD or SLIM *********************************************
        !                                                              *
        ! Generate a acCD basis.

        Do_acCD_Basis = .true.
        GWInput = .true.

      case (KeyW(136))
        !                                                              *
        !**** DOFM *****************************************************
        !                                                              *
        ! DoFMM: activate FMM option

        DoFMM = .true.

      case (KeyW(137))
        !                                                              *
        !**** NOAM *****************************************************
        !                                                              *
        ! No computation of AMFI integrals

        NoAMFI = .true.
        GWInput = .true.

      case (KeyW(138))
        !                                                              *
        !**** RPQM *****************************************************
        !                                                              *
        ! Set RPQMin for FMM option

        Key = Get_Ln(LuRd)
        call Get_F1(1,RPQMin)

      case (KeyW(139))
        !                                                              *
        !**** CONS *****************************************************
        !                                                              *
        ! Have the Gateway read the constraints for Slapaf

        GWInput = .true.
        Lu_UDC = 97
        Lu_UDC = IsFreeUnit(Lu_UDC)
        call Molcas_Open(Lu_UDC,'UDC.Gateway')
        do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          write(Lu_UDC,'(A)') trim(Key)
          if (Key(1:4) == 'END ') exit
        end do
        ! This rather obscure feature seems to be needed to to make Intel
        ! compilers behave like the others when detecting EOF
        endfile(Lu_UDC)
        close(Lu_UDC)

      case (KeyW(140))
        !                                                              *
        !**** NGEX *****************************************************
        !                                                              *
        ! Have the Gateway read the constraints for Numerical_gradient

        GWInput = .true.
        Lu_UDC = 97
        Lu_UDC = IsFreeUnit(Lu_UDC)
        call Molcas_Open(Lu_UDC,'UDC.NG')
        do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          if (adjustl(Key) == 'INVERT') then
            Invert = .true.
          else
            write(Lu_UDC,'(A)') trim(Key)
            if (Key(1:4) == 'END ') exit
          end if
        end do
        ! This rather obscure feature seems to be needed to to make Intel
        ! compilers behave like the others when detecting EOF
        endfile(Lu_UDC)
        close(Lu_UDC)

      case (KeyW(157))
        !                                                              *
        !**** NOAL *****************************************************
        !                                                              *
        ! Do not align reactants and products

        Do_Align = .false.
        if (Align_Only) then
          call WarningMessage(2,'Keywords ALIG and NOAL are not compatible')
          call Quit_OnUserError()
        end if
        GWInput = .true.

      case (KeyW(158))
        !                                                              *
        !**** WEIG *****************************************************
        !                                                              *
        ! Weights for alignment of reactants and products

        Align_Weights = Get_Ln(LuRd)
        call UpCase(Align_Weights)
        GWInput = .true.

      case (KeyW(159))
        !                                                              *
        !**** ALIG *****************************************************
        !                                                              *
        ! Align reactants and products

        Align_Only = .true.
        if (.not. Do_Align) then
          call WarningMessage(2,'Keywords ALIG and NOAL are not compatible')
          call Quit_OnUserError()
        end if
        GWInput = .true.

      case (KeyW(160))
        !                                                              *
        !***** TINK ****************************************************
        !                                                              *
        ! Read Coordinates in Tinker's xyz format

        if (SymmSet) then
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
          if (Key == '') then
            call Getenvf('MOLCAS',Key)
            Key = trim(Key)//'/tinker/bin'
          end if
          call Getenvf('Project',Project)
          Project = Project(1:index(Project,' ')-1)
          Key = trim(Key)//'/tkr2qm_s '//trim(Project)//'.xyz>'//trim(Project)//'.Tinker.log'
          write(u6,*) 'TINKER keyword found, run ',trim(Key)
          call StatusLine(' Gateway:',' Read input from Tinker')
          RC = 0
          call Systemf(trim(Key),RC)
          if (RC /= 0) then
            Key = 'RdCtl_Seward: Tinker call terminated abnormally'
            call WarningMessage(2,Key)
            call Abend()
          end if
        end if
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          call GA_Sync()
          call PFGet_ASCII('QMMM')
          call GA_Sync()
        end if
#       endif
        iCoord = iCoord+1
        CoordSet = .true.
        ITkQMMM = IsFreeUnit(ITkQMMM)
        call Molcas_Open(ITkQMMM,'QMMM')
#       ifdef _HAVE_EXTRA_
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
#       else
        if (Expert) then
          if (iCoord > 1) then
            call WarningMessage(1,'TINKER coordinates replacing COORD')
          end if
          call Read_XYZ(ITkQMMM,OrigRot,OrigTrans,Replace=(iCoord > 1))
        else
          call Read_XYZ(ITkQMMM,OrigRot,OrigTrans)
        end if
#       endif
        close(ITkQMMM)
        GWInput = .true.

      case (KeyW(161))
        !                                                              *
        !**** ORIG *****************************************************
        !                                                              *
        ! Defines translation and rotation for each xyz-file

        OrigInput = .true.
#       ifdef _HAVE_EXTRA_
        Origin_input = .true.
#       endif
        if (FragSet) then
          write(u6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
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

      case (KeyW(162))
        !                                                              *
        !**** HYPE *****************************************************
        !                                                              *

        KWord = Get_Ln(LuRd)
#       ifdef _HAVE_EXTRA_
        geoInput = .true.
#       endif
        writeZmat = .true.
        call Get_F(1,HypParam,3)
        GWinput = .true.
        HyperParSet = .true.

      case (KeyW(163))
        !                                                              *
        !**** ZCON *****************************************************
        !                                                              *

        writeZMat = .true.
#       ifdef _HAVE_EXTRA_
        ZConstraints = .true.
#       endif
        GWinput = .true.

      case (KeyW(164))
        !                                                              *
        !**** SCAL *****************************************************
        !                                                              *

        KWord = Get_Ln(LuRd)
        call Get_F1(1,ScaleFactor)
        GWinput = .true.
        if (.not. CoordSet) then
          call WarningMessage(2,'Scale can be used only with xyz input')
          call Quit_OnUserError()
        end if

      case (KeyW(165))
        !                                                              *
        !**** DOAN *****************************************************
        !                                                              *

        write(u6,*) 'RdCtl: keyword DOAN is obsolete and is ignored!'

      case (KeyW(166))
        !                                                              *
        !**** GEOE *****************************************************
        !                                                              *

        Kword = Get_Ln(LuRd)
        call Get_I1(1,iGeoInfo(2))
        GWinput = .true.
        iGeoInfo(1) = 1
        call Put_iArray('GeoInfo',iGeoInfo,2)

      case (KeyW(167))
        !                                                              *
        !**** OLDZ *****************************************************
        !                                                              *

        GWinput = .true.
#       ifdef _HAVE_EXTRA_
        oldZmat = .true.
#       endif

      case (KeyW(168))
        !                                                              *
        !**** OPTH *****************************************************
        !                                                              *

        GWinput = .true.
        Kword = Get_Ln(LuRd)
        call Get_I1(1,iOptimType)
        Kword = Get_Ln(LuRd)
        call Get_F1(1,StepFac1)
        if (iOptimType == 2) then
          KWord = Get_Ln(LuRd)
          call Get_F1(1,gradLim)
        end if

      case (KeyW(169))
        !                                                              *
        !**** NOON *****************************************************
        !                                                              *

        Do_OneEl = .false.

      case (KeyW(170))
        !                                                              *
        !**** GEO  *****************************************************
        !                                                              *

#       ifdef _HAVE_EXTRA_
        geoInput = .true.
#       endif
        writeZMat = .true.
        ! Parameters for the gridsize is set to default-values if geo is
        ! used instead of hyper
        if (.not. HyperParSet) then
          HypParam(1) = 0.15_wp
          HypParam(2) = 2.5_wp
          HypParam(3) = 2.5_wp
        end if
        GWinput = .true.

      case (KeyW(171))
        !                                                              *
        !**** MXTC *****************************************************
        !                                                              *
        ! GEN1INT integrals

        if (IRELAE == 101) then
          lMXTC = .true.
        else
          write(u6,*) 'Keyword MXTC must be preceded by keyword RX2C!'
          call Quit_OnUserError()
        end if

      case (KeyW(172))
        !                                                              *
        !**** FRGM *****************************************************
        !                                                              *

        OrigInput = .true.
#       ifdef _HAVE_EXTRA_
        Origin_input = .true.
#       endif
        GWinput = .true.
        if (OriginSet) then
          write(u6,*) 'Keywords FRGM and ORIG are mutually exclusive!'
          call Quit_OnUserError()
        end if
        if (.not. FragSet) then
          call mma_allocate(OrigTrans,3,nFragment,label='OrigTrans')
          call mma_allocate(OrigRot,3,3,nFragment,label='OrogRot')
          ! Set up no translation and no rotation as default
          OrigTrans(:,:) = Zero
          OrigRot(:,:,:) = Zero
          OrigRot(1,1,:) = One
          OrigRot(2,2,:) = One
          OrigRot(3,3,:) = One
          FragSet = .true.
        end if
        Kword = Get_Ln(LuRd)
        call Get_I1(1,iFrag)

      case (KeyW(173))
        !                                                              *
        !**** TRAN *****************************************************
        !                                                              *

        if (.not. FragSet) then
          write(u6,*) 'Keyword TRANS must be preceded by keyword FRAG!'
          call Quit_OnUserError()
        end if
        GWinput = .true.
        Kword = Get_Ln(LuRd)
        call Get_F(1,OrigTrans(1,iFrag),3)

      case (KeyW(174))
        !                                                              *
        !***** ROT  ****************************************************
        !                                                              *

        if (.not. FragSet) then
          write(u6,*) 'Keyword ROT must be preceded by keyword FRAG!'
          call Quit_OnUserError()
        end if
        GWinput = .true.
        Kword = Get_Ln(LuRd)
        call Get_F(1,OrigRot(1,1,iFrag),9)

      case (KeyW(175))
        !                                                              *
        !****** ZONL ***************************************************
        !                                                              *

        GWinput = .true.
        WriteZMat = .true.

      case (KeyW(176))
        !                                                              *
        !****** BASL ***************************************************
        !                                                              *

        GWinput = .true.
        BasLib = Get_Ln(LuRd)
        Write_BasLib = .true.

      case (KeyW(177))
        !                                                              *
        !****** NUME ***************************************************
        !                                                              *

        GWinput = .true.
        iDNG = 1

      case (KeyW(178))
        !                                                              *
        !****** VART ***************************************************
        !                                                              *

        GWinput = .true.
        VarT = .true.

      case (KeyW(179))
        !                                                              *
        !****** VARR ***************************************************
        !                                                              *

        GWinput = .true.
        VarR = .true.

      case (KeyW(180))
        !                                                              *
        !****** SHAK ***************************************************
        !                                                              *

        GWinput = .true.
        KWord = Get_Ln(LuRd)
        call Upcase(KWord)
        call Get_F1(1,Shake)
        if (index(KWord,'ANGSTROM') /= 0) Shake = Shake/Angstrom

      case (KeyW(181))
        !                                                              *
        !***** PAMF ****************************************************
        !                                                              *
        ! Disable AMFI for an atom type

        KWord = Get_Ln(LuRd)
        call Get_I1(1,iAtom_Number)
        No_AMFI(iAtom_Number) = .true.

      case (KeyW(182))
        !                                                              *
        !****** GROM ***************************************************
        !                                                              *
        ! Import definition of QMMM system from Gromacs

#       ifdef _GROMACS_
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
        select case (KWord(1:4))
          case ('SIMP')
            nCastMM = 0
            call mma_allocate(CastMM,nCastMM)
          case ('CAST')
            KWord = Get_Ln(LuRd)
            call Get_I1(1,nCastMM)
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
          case default
            Message = 'GROMACS keyword found, but no valid option'
            call WarningMessage(2,Message)
            call Quit_OnUserError()
        end select

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
#       ifdef _HAVE_EXTRA_
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
#       else
        if (Expert) then
          if (iCoord > 1) then
            call WarningMessage(1,'TINKER coordinates replacing COORD')
          end if
          call Read_XYZ(LuXYZ,OrigRot,OrigTrans,Replace=(iCoord > 1))
        else
          call Read_XYZ(LuXYZ,OrigRot,OrigTrans)
        end if
#       endif
        close(LuXYZ)
#       else
        Message = 'Interface to Gromacs not installed'
        call WarningMessage(2,Message)
        call Quit_OnUserError()
#       endif

      case (KeyW(183))
        !                                                              *
        !****** LINK ***************************************************
        !                                                              *
        ! Define link atoms for a Molcas/Gromacs run

#       ifdef _GROMACS_
        GWInput = .true.
        KWord = Get_Ln(LuRd)
        call Get_I1(1,nLA)
        if (nLA <= 0) then
          Message = 'LA definition: nLA is zero or negative'
          call WarningMessage(2,Message)
          call Quit_OnUserError()
        end if
#       ifdef _DEBUGPRINT_
        write(u6,'(/,a)') ' Link atoms (Gromacs numbering):'
        write(u6,'(/,a)') '      LA     QM     MM     Scaling factor'
#       endif
        call mma_allocate(DefLA,3,nLA)
        call mma_allocate(FactLA,nLA)
        do iLA=1,nLA
          KWord = Get_Ln(LuRd)
          call Get_I(1,DefLA(1,iLA),3)
          call Get_F(4,FactLA(iLA),1)
#        ifdef _DEBUGPRINT_
          write(u6,'(i8,2i7,F19.8)') (DefLA(i,iLA),i=1,3),FactLA(iLA)
#        endif
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
#       else
        Message = 'Interface to Gromacs not installed'
        call WarningMessage(2,Message)
        call Quit_OnUserError()
#       endif

      case (KeyW(184))
        !                                                              *
        !***** EMFR ****************************************************
        !                                                              *

        GWinput = .true.
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
        if (index(KWord,'ANGSTROM') /= 0) Lambda = Lambda/Angstrom
        if (index(KWord,'NANOMETER') /= 0) then
          Lambda = Ten*Lambda/Angstrom
        end if
        KVector(1) = ((Two*Pi)/Lambda)*KVector(1)
        KVector(2) = ((Two*Pi)/Lambda)*KVector(2)
        KVector(3) = ((Two*Pi)/Lambda)*KVector(3)

      case (KeyW(185))
        !                                                              *
        !***** NOCD ****************************************************
        !                                                              *

        GWinput = .true.
        if (.not. CholeskyWasSet) then
          Do_RI = .false.
          iRI_Type = 0
          Chol = .false.
          CholeskyWasSet = .true.
        end if

      case (KeyW(186))
        !                                                              *
        !***** FNMC ****************************************************
        !                                                              *

        GWinput = .true.
        FNMC = .true.

      case (KeyW(187))
        !                                                              *
        !***** ISOT ****************************************************
        !                                                              *

        GWinput = .true.
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

      case (KeyW(188))
        !                                                              *
        !***** EFP  ****************************************************
        !                                                              *

        GWinput = .true.
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
          write(u6,*) 'XYZABC option to be implemented'
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
              if (jEnd > LenIn+1) then
                write(u6,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LenIn,' characters!'
              end if
              ABC(i,iFrag) = KWord(1:min(LenIn,jend-1))
              call Get_F(2,EFP_COORS((i-1)*3+1,iFrag),3)
            end do
          end do
        else if (KWord == 'ROTMAT') then
          Coor_Type = ROTMAT_type
          nEFP_Coor = 12
          allocate(EFP_COORS(nEFP_Coor,nEFP_fragments))
          write(u6,*) 'ROTMAT option to be implemented'
          call Abend()
        else
          write(u6,*) 'Illegal EFP format :',KWord
          write(u6,*)
          write(u6,*) 'Allowed format: XYZABC,'
          write(u6,*) '                POINTS, and'
          write(u6,*) '                ROTMAT'
        end if
        lEFP = .true.

      case default
        if (lTtl) then
          call ProcessTitle()
        else if (Basis_test) then
          ! So the Basis keyword was in the native format.
          ! We have to back step until we find the command line!

          backspace(LuRd)
          backspace(LuRd)
          read(LuRd,'(A)') Key
          call UpCase(Key)
          do while (index(Key(1:4),KeyW(29)) == 0) ! BASI
            backspace(LuRd)
            backspace(LuRd)
            read(LuRd,'(A)') Key
            call UpCase(Key)
          end do
          Basis_test = .false.
          nDone = 0
          call ProcessBasis()
        else
          Last = len_trim(KWord)
          write(u6,*)
          call WarningMessage(2,KWord(1:Last)//' is not a keyword!, Error in keyword.')
          call Quit_OnUserError()

          !call WarningMessage(2,' Premature end of input file.')
          !call Quit_OnUserError()
          exit
        end if
    end select

  end do

  ! Postprocessing for COORD
  !ik = index(KeepBasis,'....')
  !if (ik /= 0) then
  !    KeepBasis = KeepBasis(1:ik-1)
  !end if
  if (CoordSet) then
    do ik=len(KeepBasis),1,-1
      if ((KeepBasis(ik:ik) /= ' ') .and. (KeepBasis(ik:ik) /= '.')) exit
    end do
    KeepBasis = KeepBasis(1:ik)
#   ifdef _HAVE_EXTRA_
    call ProcessXYZ(BasisSet,KeepBasis,KeepGroup,iBSSE,SymThr,isHold,ScaleFactor,HyperParSet,isXfield)
#   else
    call Parse_Basis(KeepBasis)
    call Parse_Group(KeepGroup,SymThr)
    call Write_SewInp('COORD',[iBSSE])
#   endif
    if (writeZmat) then
#     ifdef _HAVE_EXTRA_
      stepFactor = stepFac1/(hypParam(1)*hypParam(1))
      call Geo_Setup_Drv(ishold,oldZMat,zConstraints,geoInput,hypParam,nFragment,iOptimType,stepFactor,gradLim)
#     else
      call WarningMessage(2,'molcas-extra not installed')
      call Quit_OnUserError()
#     endif
    end if
    DoneCoord = .true.
    if (isXfield == 1) then
      LuRd_saved = LuRd
      filename = 'findsym.xfield'
      lXF = .true.
      call OpnFl(filename(1:(index(filename,' ')-1)),LuRd,Exists)
      if (.not. Exists) then
        call WarningMessage(2,'Error! File not found: '//filename(1:(index(filename,' ')-1)))
        call Quit_OnUserError()
      end if
      write(u6,*) 'Reading external field from file: ',filename(1:(index(filename,' ')-1))
      KWord = Get_Ln(LuRd)
      call ProcessXF()
    end if

    CoordSet = .false.
    LuRdSave = LuRd
    LuFS = IsFreeUnit(1)
#   ifdef _HAVE_EXTRA_
    call Molcas_Open(LuFS,'FS.std')
#   else
    call Molcas_Open(LuFS,'COORD')
#   endif
    LuRd = LuFS
    GWInput = .true.
  else
    if (DoneCoord) then
      close(LuFS)
      LuRd = LuRdSave
#     ifndef _HAVE_EXTRA_
      call Clear_XYZ()
#     endif
    end if
    exit
  end if
end do

!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     P O S T   P R O C E S S I N G                                    *
!                                                                      *
!***********************************************************************
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
    call ModGauss(real(jAtmNr,kind=wp),nMass,dbsc(iCnttp)%ExpNuc,dbsc(iCnttp)%w_mGauss)

  else

    ! Nothing to do for point charges!

  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Cholesky-specific postprocessing:
! 1) reset integral prescreening thresholds (if not user-defined).
! 2) use default Cholesky normalization (if not user-defined);
!    for AMFI or Douglas-Kroll, use Molcas normalization (again,
!    if not user-defined).
! 3) Turn off Onenly flag, which might be set through the Direct
!    keyword. Thus, specifying Cholesky will force 2-el. int.
!    processing even with Direct specifed as well!
! 4) Turn off Dist flag (makes no sense with Cholesky).
! 5) Integral format on disk is irrelevant.
! 6) if Cholesky threshold is specified, use it.
! 7) if span factor is specified, use it.

if (Chol) then
  if (Onenly) then
    Chol = .false. ! we gotta be lazy in such cases
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
    if (CholeskyThr >= Zero) then
      Thrshld_CD = CholeskyThr
      ThrCom = Thrshld_CD
    end if
    if (spanCD >= Zero) Span = min(spanCD,One)
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
      write(u6,*) ' MxAtom=',MxAtom
      call Abend()
    end if

    ! The symmetry operators of the fragment's atoms should
    ! always be identical to that of the fragment's
    ! pseudocenter/placeholder

    if (dbsc(iCnttp)%Frag) then
      !  Check the FragExpand routine!
      if (abs(dbsc(iCnttp)%nFragCoor) > mdc) then
        write(u6,*) 'rdctl_seward: incorrect mdc index'
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
        if (.not. btest(jTmp,j)) S%nDim = S%nDim+1
      end do
      if (S%nDim > 0) then
        call Random_Vector(S%nDim,RandVect(1:S%nDim),.false.)
        jDim = 0
        do j=0,2
          if (.not. btest(jTmp,j)) then
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
  write(u6,*) 'S%mCentr=',S%mCentr
  write(u6,*) 'Edit src/Include/Molcas.fh'
  write(u6,*) 'Set MxAtom to the value of S%mCentr.'
  write(u6,*) 'Recompile MOLCAS and try again!'
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
    write(u6,*)
    write(u6,'(15X,A)') repeat('*',88)
    write(u6,'(15X,A,A,A)') '*',repeat(' ',86),'*'
    do iTtl=1,nTtl
      write(u6,'(15X,A,A,A)') '*   ',Title(iTtl),'   *'
    end do
    write(u6,'(15X,A,A,A)') '*',repeat(' ',86),'*'
    write(u6,'(15X,A)') repeat('*',88)
  else
    write(u6,*)
    write(u6,'(A)') ' Title:'
    do iTtl=1,nTtl
      write(u6,'(8X,A)') Title(iTtl)
    end do
    write(u6,*)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the weights used for alignment and distance measurement

call Process_Weights(iPrint)
!                                                                      *
!***********************************************************************
!                                                                      *
! Set structures for TS optimization according to the Saddle method.

if (Run_Mode /= G_Mode) then
  call Saddle()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Read coordinates from run file (if any), ditto for external field.
  ! Do not do this in the Gateway!

  call GeoNew(Show)
  if (lXF) call GeoNew_PC()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Gen_GeoList()
call mma_deallocate(xb_bas)
call mma_deallocate(xb_label)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of centers for multipole operators

call SetMltplCenters()

! Put in user specified centers if any

if (lMltpl) then
  do i=1,nTemp
    iMltpl = ITmp(i)
    if (iMltpl <= S%nMltpl) Coor_MPM(:,iMltpl+1) = RTmp(:,i)
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
  OAM_Center(:) = OAMt
else if (.not. allocated(OAM_Center)) then
  call mma_allocate(OAM_Center,3,Label='OAM_Center')
  OAM_Center(:) = CoM
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put up list for point at which the orbital magnetic quadrupole
! will be computed.

if (lOMQ .and. (Run_Mode /= S_Mode)) then
  call mma_allocate(OMQ_Center,3,Label='OMQ_Center')
  OMQ_Center(:) = OMQt
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

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      call WarningMessage(2,'Unable to read data from '//KWord)
    case (2)
      write(u6,'(a,a)') 'Error reading from file ',Key(1:iFile-1)
      write(u6,'(a,a)') 'unable to process line: ',KWord
  end select
  call Quit_OnUserError()

end subroutine Error

subroutine ProcessTitle()

  GWInput = AnyMode
  lTtl = .true.
  nTtl = nTtl+1
  if (nTtl > 10) then
    call WarningMessage(2,' Too many title cards')
    call Quit_OnUserError()
  end if
  i1 = iCFrst(Key,80)
  i2 = len_trim(Key(1:len(Title)))
  nc = 80-(i2-i1+1)
  nc2 = nc/2
  Title(nTtl) = ''
  Title(nTtl)(nc2+1:nc2+i2-i1+1) = Key(i1:i2)
  Skip1 = .true.

end subroutine ProcessTitle

subroutine ProcessBasis()

  integer(kind=iwp) :: i, iSh, n1, n2, n3

  iOpt_XYZ = 0
  GWInput = .true.
  nCnttp = nCnttp+1
  if (Run_Mode == S_Mode) then
    call WarningMessage(2,'Seward input error!')
    write(u6,*) 'The command : "',Previous_Command,'" is not allowed in the Seward input when the Gateway is used!'
    write(u6,*)
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
  BasisTypes_save(:) = BasisTypes
  if ((BSLbl(1:2) == 'X.') .and. (index(BSLbl,'INLINE') == 0) .and. (index(BSLbl,'RYDBERG') == 0)) then
    BSLbl = 'X.ANO-RCC.'
    do i=11,80
      BSLbl(i:i) = '.'
    end do
    iDummy_basis = 1
  end if
  !call UpCase(BSLbl)
  Last = len_trim(BSLbl)
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
    Fname = adjustl(Fname)
    dbsc(nCnttp)%Bsl = BSLbl(1:Indx-1)
  end if

  n = index(dbsc(nCnttp)%Bsl,' ')
  if (n == 0) n = 81
  do i=n,80
    dbsc(nCnttp)%Bsl(i:i) = '.'
  end do

  if ((Show .and. (nPrint(2) >= 6)) .or. Write_BasLib) then
    write(u6,*)
    write(u6,*)
    write(u6,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLbl(1:Indx-1)
    write(u6,'(1X,A,A)') 'Basis set is read from library:',Fname(1:index(Fname,' '))
  end if

  jShll = iShll
  dbsc(nCnttp)%Bsl_old = dbsc(nCnttp)%Bsl
  dbsc(nCnttp)%mdci = mdc
  call GetBS(Fname,dbsc(nCnttp)%Bsl,iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,lSTDINP,.false.,Expert,ExtBasDir)

  Do_FckInt = Do_FckInt .and. dbsc(nCnttp)%FOp .and. (dbsc(nCnttp)%AtmNr <= 96)
# ifdef _DEMO_
  Do_GuessOrb = .false.
# else
  Do_GuessOrb = Do_GuessOrb .and. (dbsc(nCnttp)%AtmNr <= 96)
# endif

  if (iDummy_Basis == 1) BasisTypes(:) = BasisTypes_save
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
    write(u6,'(1x,a)') 'Basis Set Reference(s):'
    if (Ref(1) /= '') write(u6,'(5x,a)') trim(Ref(1))
    if (Ref(2) /= '') write(u6,'(5x,a)') trim(Ref(2))
    write(u6,*)
    write(u6,*)
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
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Set Cartesian functions if specified by the basis type

  if (BasisTypes(1) == 9) then
    do iSh=jShll+3,iShll
      Shells(iSh)%Transf = .false.
      Shells(iSh)%Prjct = .false.
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Automatic onset of muonic charge if the basis type is muonic.
  ! This will also automatically activate finite nuclear mass correction.

  KWord = ''
  KWord(1:Indx-1) = BSLbl(1:Indx-1)
  call UpCase(KWord)
  if (index(KWord,'MUONIC') /= 0) then
    dbsc(nCnttp)%fMass = mu2elmass
    FNMC = .true.
    tDel = 1.0e50_wp
    call Put_dScalar('T delete thr',tDel)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Update BasisTypes

  do i=1,4
    if (BasisTypes_save(i) == 0) cycle
    if (BasisTypes(i) /= BasisTypes_save(i)) BasisTypes(i) = -1
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  do
    KWord = Get_Ln(LuRd)
    call UpCase(KWord)
    KWord = adjustl(KWord)
    select case (KWord(1:4))
      case ('PSEU')
        dbsc(nCnttp)%pChrg = .true.
        dbsc(nCnttp)%Fixed = .true.
      case ('ACDT')
        KWord = Get_Ln(LuRd)
        call Get_F1(1,dbsc(nCnttp)%aCD_Thr)
      case ('MUON')
        dbsc(nCnttp)%fMass = mu2elmass
      case ('NUCL')
        KWord = Get_Ln(LuRd)
        call Get_F1(1,dbsc(nCnttp)%ExpNuc)
      case ('FIXE')
        dbsc(nCnttp)%Fixed = .true.
      case ('SPHE')
        if (index(KWord,'ALL') /= 0) then
          do iSh=jShll+3,iShll
            Shells(iSh)%Transf = .true.
            Shells(iSh)%Prjct = .true.
          end do
        else
          ist = index(KWord,' ')
          iAng = 2
          do iSh=jShll+3,iShll
            if (index(KWord(ist:80),AngTyp(iAng)) /= 0) then
              Shells(iSh)%Transf = .true.
              Shells(iSh)%Prjct = .true.
            end if
            iAng = iAng+1
          end do
        end if
      case ('CART')
        if (index(KWord,'ALL') /= 0) then
          do iSh=jShll+1,iShll
            Shells(iSh)%Transf = .false.
            Shells(iSh)%Prjct = .false.
          end do
        else
          ist = index(KWord,' ')
          iAng = 0
          do iSh=jShll+1,iShll
            if (index(KWord(ist:80),AngTyp(iAng)) /= 0) then
              Shells(iSh)%Transf = .false.
              Shells(iSh)%Prjct = .false.
            end if
            iAng = iAng+1
          end do
        end if
      case ('CONT')
        if (index(KWord,'ALL') /= 0) then
          do iSh=jShll+1,iShll
            Shells(iSh)%Prjct = .false.
          end do
        else
          ist = index(KWord,' ')
          iAng = 0
          do iSh=jShll+1,iShll
            if (index(KWord(ist:80),AngTyp(iAng)) /= 0) Shells(iSh)%Prjct = .false.
            iAng = iAng+1
          end do
        end if
      case ('CHAR')
        KWord = Get_Ln(LuRd)
        call UpCase(KWord)
        call Get_F1(1,dbsc(nCnttp)%Charge)
        ist = index(KWord,' ')
        if (dbsc(nCnttp)%IsMM /= 0) then
          call WarningMessage(1,' Found a charge associated with a MM atom. Ignore it')
          dbsc(nCnttp)%Charge = Zero
        end if
      case ('FRAG')
        dbsc(nCnttp)%pChrg = .true.
        dbsc(nCnttp)%Fixed = .true.
        lFAIEMP = .true.
      case ('END ')
        if (nCnt == 0) then
          call WarningMessage(2,' Input error, no center specified!')
          call Quit_OnUserError()
        end if
        dbsc(nCnttp)%nCntr = nCnt
        mdc = mdc+nCnt
        ! Now allocate the array for the coordinates and copy them over.
        ! allocate(dbsc(nCnttp)%Coor(1:3,1:nCnt)
        call mma_Allocate(dbsc(nCnttp)%Coor_Hidden,3,nCnt,Label='dbsc:C')
        dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
        call DCopy_(3*nCnt,Buffer,1,dbsc(nCnttp)%Coor,1)
        exit
      case default
        ! Read Coordinates

        nCnt = nCnt+1
        n_dc = max(mdc+nCnt,n_dc)
        if (mdc+nCnt > MxAtom) then
          call WarningMessage(2,' RdCtl: Increase MxAtom')
          write(u6,*) '        MxAtom=',MxAtom
          call Quit_OnUserError()
        end if
        jend = index(KWord,' ')
        if (jEnd > LenIn+1) then
          write(u6,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LenIn,' characters!'
        end if
        dc(mdc+nCnt)%LblCnt = KWord(1:min(LenIn,jend-1))
        dbas = dc(mdc+nCnt)%LblCnt(1:LenIn)
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
            Buffer(iOff+i) = Buffer(iOff+i)/Angstrom
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
                if ((n1 == 0) .and. (n2 == 0) .and. (n3 == 0)) cycle

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
                  write(u6,*) '        MxAtom=',MxAtom
                  call Quit_OnUserError()
                end if

                jend = index(KWord,' ')
                if (jEnd > 5) then
                  write(u6,*) 'Warning: the label ',KWord(1:jEnd),' will be truncated to ',LenIn,' characters!'
                end if
                dc(mdc+nCnt)%LblCnt = KWord(1:min(LenIn,jend-1))//CHAR4

                call Chk_LblCnt(dc(mdc+nCnt)%LblCnt,mdc+nCnt-1)

                iOff = 1+(nCnt-1)*3

                ! Copy old coordinate  first
                Buffer(iOff:iOff+2) = Buffer(iOff0:iOff0+2)+n1*VCell(:,1)+n2*VCell(:,2)+n3*VCell(:,3)

              end do
            end do
          end do
        end if
    end select
  end do

end subroutine ProcessBasis

subroutine ProcessEF(KWName)

  character(len=4), intent(in) :: KWName
  integer(kind=iwp) :: iCnt, iCnttp, iEF

  KWord = Get_Ln(LuRd)
  call Get_I1(1,nEF)
  if (nEF < 0) nEF = 0
  if (nEF /= 0) then
    call mma_allocate(EFt,3,nEF,label='nEF')
    do iEF=1,nEF
      KWord = Get_Ln(LuRd)
      ! Check whether a label is specified instead of a coordinate
      call Get_S(1,Key,1)
      call Upcase(Key)
      Key = adjustl(Key)
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
          call WarningMessage(2,'; Error in processing the keyword '//KWName//'.;'// &
                              'The label "'//Key(1:jEnd)//'" could not be found among the centers.;'// &
                              'Remember to specify the atom center before specifying the '//KWName//' keyword.')
          call Quit_OnUserError()
        end if
      else
        call Get_F(1,EFt(:,iEF),3)
        call Upcase(KWord)
        if (index(KWord,'ANGSTROM') /= 0) EFt(:,i) = EFt(:,i)/Angstrom
      end if
    end do
  end if

end subroutine ProcessEF

subroutine ProcessXF()

  integer(kind=iwp) :: i, iOrd_XF, iXF, k

  call Get_I1(1,nXF)
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
    write(u6,*) 'nOrd_XF= ',nOrd_XF
    call Quit_OnUserError()
  end if
  if ((iXPolType > 2) .or. (iXPolType < 0)) then
    call WarningMessage(2,'Error! Illegal value of iXPolType')
    write(u6,*) 'iXPolType= ',iXPolType
    call Quit_OnUserError()
  end if
  if ((nXMolnr > 100) .or. (nXMolnr < 0)) then
    call WarningMessage(2,'Error! Illegal value of nXMolnr')
    write(u6,*) 'nXMolnr= ',nXMolnr
    call Quit_OnUserError()
  end if
  if ((nReadEle > 1) .or. (nReadEle < 0)) then
    call WarningMessage(2,'Error! Illegal value of nReadEle')
    write(u6,*) 'nReadEle= ',nReadEle
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
    if (Convert) XF(1:3,iXF) = XF(1:3,iXF)/Angstrom

  end do

  ! Close file and reset LuRd if external file was used
  if (LuRd /= LuRd_saved) then
    close(LuRd)
    LuRd = LuRd_saved
  end if

end subroutine ProcessXF

end subroutine RdCtl_Seward
