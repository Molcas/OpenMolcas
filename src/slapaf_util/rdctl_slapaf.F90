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

subroutine RdCtl_Slapaf(LuSpool,Dummy_Call)

use kriging_mod
use ThermoChem
use Symmetry_Info, only: Symmetry_Info_Get
use Slapaf_Info, only: Cx, Gx, Weights, MF, Atom, nSup, RefGeo, GradRef, nStab, Lbl, mRowH, Coor
use Slapaf_Parameters, only: iRow, iRow_c, ddV_Schlegel, HWRS, iOptH, HrmFrq_Show, IRC, Curvilinear, Redundant, FindTS, nBVec, &
                             User_Def, MaxItr, iOptC, rHidden, CnstWght, lOld, Beta, Beta_Disp, Line_Search, TSConstraints, &
                             GNrm_Threshold, Mode, ThrEne, ThrGrd, nLambda, ThrCons, ThrMEP, Baker, eMEPTest, rMEP, MEP, nMEP, &
                             MEPNum, MEPCons, dMEPStep, MEP_Type, MEP_Algo, Max_Center, Delta, RtRnc, rFuzz, lNmHss, Cubic, &
                             Request_Alaska, CallLast, lCtoF, Track, isFalcon, MxItr, nWndw, Iter, WeightedConstraints, NADC, &
                             Fallback
use UnixInfo, only: SuperName

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
integer iDum(1)
logical Found, Dummy_Call
character(len=180) Get_Ln
character*16 FilNam
character*3 MEPLab
character(len=180), parameter :: BLine = ''
character(len=180) :: Key = '', Char = ''
real*8, allocatable :: DIR(:,:), Tmp(:), TmpRx(:)
#include "cgetl.fh"
external Get_Ln
logical External_UDC, Explicit_IRC, Expert, ThrInp, FirstNum, Manual_Beta
#include "angstr.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
Expert = .false.
Lu = 6
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate some parameters

call Symmetry_Info_Get()
call Init_Slapaf()
nsAtom = size(Coor,2)
iPrint = nPrint(iRout)
iSetAll = 2**30-1

call f_Inquire('UDC.Gateway',External_UDC)

iMEP = 0
Explicit_IRC = .false.
WeightedConstraints = .false.
ThrInp = .false.
call Qpg_iScalar('nMEP',Found)
if (Found) call Get_iScalar('nMEP',iMEP)
if (iMEP == 0) then
  iOff_Iter = 0
  call Put_iScalar('iOff_Iter',iOff_Iter)
end if
Manual_Beta = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! When called from outside Slapaf or as a dummy call, process no
! input but proceed with default values only.

if ((SuperName == 'slapaf') .and. (.not. Dummy_Call)) then
  !                                                                    *
  !*********************************************************************
  !************************   Input section   **************************
  !*********************************************************************
  !                                                                    *
  LuRd = LuSpool
  call RdNlst(LuRd,'SLAPAF')
  do
    Char = Get_Ln(LuRd)
    call UpCase(Char)
    !write(Lu,'(A)') Char
    !write(Lu,*) iOptC
    if (char == BLine) cycle
    if (char(1:1) == '*') cycle
    select case (char(1:4))
      !case ('AIL ')
      !  !                                                              *
      !  !***** AIL  ****************************************************
      !  !                                                              *
      !  ! Width limits of the Matern function
      !
      !  Char = Get_Ln(LuRd)
      !  call Get_F(1,lb,3)

      !case ('AIP ')
      !  !                                                              *
      !  !***** AIP  ****************************************************
      !  !                                                              *
      !  ! Parameter of differentiability for Matern function
      !
      !  Char = Get_Ln(LuRd)
      !  call Get_F1(1,pAI)
      !  if ((pAI > 3) .or. (pAI < 1)) anMd = .false.

      !case ('AISP')
      !  !                                                              *
      !  !***** AISP ****************************************************
      !  !                                                              *
      !  ! Defining the number of source points for the AI method
      !
      !  Char = Get_Ln(LuRd)
      !  call Get_I1(1,nspAI)

      !case ('AIMD')
      !  !                                                              *
      !  !***** AIMD ****************************************************
      !  !                                                              *
      !  ! Analytical or numerical Matern derivatives
      !
      !  Char = Get_Ln(LuRd)
      !  if ((Char == 'False') .or. (Char == 'false')) then
      !    anMd = .false.
      !  else
      !    anMd = .true.
      !  end if

      !case ('AIME')
      !  !                                                              *
      !  !***** AIME ****************************************************
      !  !                                                              *
      !  ! Minimum energy differences of the last two iterations
      !  ! (loop exit condition)
      !
      !  Char = Get_Ln(LuRd)
      !  call Get_F1(1,Thr_microiterations)

      !case ('AIBL')
      !  !                                                              *
      !  !***** AIBL ****************************************************
      !  !                                                              *
      !  ! Base line modification value to not ordinary
      !  ! (Trend Function on GEK)
      !
      !  Char = Get_Ln(LuRd)
      !  call Get_F1(1,blvAI)
      !  blAI = .true.

      !case ('AIMB')
      !  !                                                              *
      !  !***** AIMB ****************************************************
      !  !                                                              *
      !  ! Base line modification value to maximum value of the Energy
      !  ! This option supersedes any value assigned to blAI
      !
      !  Char = Get_Ln(LuRd)
      !  mblAI = .true.

      !case ('L-VA')
      !  !                                                              *
      !  !***** L-VA ****************************************************
      !  !                                                              *
      !  ! Change the l value of the GEK.
      !
      !  Char = Get_Ln(LuRd)
      !  Set_l = .true.
      !  call Get_F1(1,Value_l)
      !  call Qpg_dScalar('Value_l',Found)
      !  if (.not. Found) call Put_dScalar('Value_l',Value_l)

      case ('NDEL')
        !                                                              *
        !***** NDELta **************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_I1(1,nD_In)

      case ('BAKE')
        !                                                              *
        !***** BAKE ****************************************************
        !                                                              *
        Baker = .true.

      case ('C1-D')
        !                                                              *
        !***** C1-D ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**0-2**1-2**2-2**3
        iOptC = ior(2**1,iand(iOptC,Mask))

      case ('C2-D')
        !                                                              *
        !***** C2-D ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**0-2**1-2**2-2**3
        iOptC = ior(2**2,iand(iOptC,Mask))

      case ('CART')
        !                                                              *
        !***** CART ****************************************************
        !                                                              *
        CurviLinear = .false.

      case ('CNWE')
        !                                                              *
        !***** CNWE ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_F1(1,CnstWght)

      case ('CONS')
        !                                                              *
        !***** CONS ****************************************************
        !                                                              *
        ! Copy constraints definition into the UDC file, to be read
        ! (after fixing and merging, if necessary) by DefInt2/Cllct2.

        if (.not. Expert) then
          write(Lu,*)
          write(Lu,*) ' ************ ERROR ***************'
          write(Lu,*) ' Obsolete input standard!'
          write(Lu,*) ' The CONSTRAINT section should'
          write(Lu,*) ' be define in the &Gateway input.'
          write(Lu,*)
          write(Lu,*) ' To override add the EXPERT keyword'
          write(Lu,*) ' to the top of the SLAPAF input.'
          write(Lu,*) ' **********************************'
          call Quit_OnUserError()
        end if
        if (External_UDC) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '****************** ERROR *********************'
          write(Lu,*) 'Multiple definitions of constraints.'
          write(Lu,*) 'Check Gateway and Slapaf inputs for conflicts!'
          write(Lu,*) '**********************************************'
          call Quit_OnUserError()
        end if
        Lu_UDC = 20
        FilNam = 'UDC'
        Lu_UDC = IsFreeUnit(Lu_UDC)
        call Molcas_Open(Lu_UDC,FilNam)
        do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          Key = adjustl(Key)
          write(Lu_UDC,'(A)') trim(Key)
          if (Key(1:4) == 'END') exit
        end do
        close(Lu_UDC)

      case ('CTOF')
        !                                                              *
        !***** CTOF ****************************************************
        !                                                              *
        ! Read the internal (C)oordinate specification (TO) be (F)ollowed.

        if (iRow > 0) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '*********** ERROR ***************'
          write(Lu,*) 'CtoF and User-defined Coordinates'
          write(Lu,*) 'are mutually exclusive.          '
          write(Lu,*) '*********************************'
          call Quit_OnUserError()
        end if
        iNull = 0
        New_Line = 1
        lCtoF = .true.
        Lu_UDIC = 91
        FilNam = 'UDIC'
        call molcas_open(Lu_UDIC,FilNam)
        rewind(Lu_UDIC)
        Key = Get_Ln(LuRd)
        call UpCase(Key)
        call FixEqualSign2(Key,LuRd,Lu_UDIC,iNull,New_Line)
        write(Lu_UDIC,'(A)') Key
        close(Lu_UDIC)

      case ('CUBI')
        !                                                              *
        !***** CUBI ****************************************************
        !                                                              *
        Cubic = .true.

      case ('DDVS')
        !                                                              *
        !***** DDVS ****************************************************
        !                                                              *
        DDV_Schlegel = .true.

      case ('DELT')
        !                                                              *
        !***** DELT ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_F1(1,Delta)

      case ('DISO')
        !                                                              *
        !***** DISO ****************************************************
        !                                                              *
        lDoubleIso = .true.

      case ('DXDX')
        !                                                              *
        !***** DXDX ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**4-2**5-2**6
        iOptC = ior(2**4,iand(iOptC,Mask))

      case ('DXG ','GDX ')
        !                                                              *
        !***** DXG  ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**4-2**5-2**6
        iOptC = ior(2**5,iand(iOptC,Mask))

      case ('END ')
        exit

      case ('EXPE')
        !                                                              *
        !***** EXPE ****************************************************
        !                                                              *
        Expert = .true.

      case ('EXTR')
        !                                                              *
        !***** EXTR ****************************************************
        !                                                              *
        ! Put the program name and the time stamp onto the extract file

        write(Lu,*) 'RdCtl_Slapaf: EXTRACT option is redundant and is ignored!'

      case ('FALC')
        !                                                              *
        !***** FALC ****************************************************
        !                                                              *
        isFalcon = .true.

      case ('FIND')
        !                                                              *
        !***** FIND ****************************************************
        !                                                              *
        FindTS = .true.

      case ('FUZZ')
        !                                                              *
        !***** FUZZ ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call UpCase(Char)
        call Get_F1(1,rFuzz)
        if (index(Char,'ANGSTROM') /= 0) rFuzz = rFuzz/angstr
        rFuzz = max(rFuzz,1.0D-3)

      case ('GG  ')
        !                                                              *
        !***** GG   ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**4-2**5-2**6
        iOptC = ior(2**6,iand(iOptC,Mask))

      case ('GNRM')
        !                                                              *
        !***** GNRM ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_F1(1,GNrm_Threshold)

      case ('GRAD')
        !                                                              *
        !***** GRAD ****************************************************
        !                                                              *
        call mma_allocate(GradRef,3,nsAtom,Label='GradRef')

        call Read_v(LuRd,GradRef,1,3*nsAtom,1,iErr)
        if (iErr /= 0) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '*************** ERROR ******************'
          write(Lu,*) 'Error reading reference gradient vector.'
          write(Lu,*) '****************************************'
          call Quit_OnUserError()
        end if

        ! If there is a transverse vector stored, we are not using this one

        call qpg_dArray('Transverse',Found,nRP)
        if (Found) call mma_deallocate(GradRef)

      case ('HRMS')
        !                                                              *
        !***** HRMS ****************************************************
        !                                                              *
        HrmFrq_Show = .true.

      case ('HUPD')
        !                                                              *
        !***** HUPD ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        read(Char,*) Char
        call UpCase(Char)
        if (trim(Char) == 'BFGS') then
          iOptH = 4
        !else if (trim(Char) == 'MEYER') then
        !  iOptH = ior(1,iAnd(iOptH,32))
        !else if (trim(Char) == 'BP') then
        !  iOptH = ior(2,iAnd(iOptH,32))
        else if (trim(Char) == 'NONE') then
          iOptH = ior(8,iand(iOptH,32))
        else if (trim(Char) == 'MSP') then
          iOptH = ior(16,iand(iOptH,32))
        else if (trim(Char) == 'EU') then
          iOptH = ior(64,iand(iOptH,32))
        else if (trim(Char) == 'TS-BFGS') then
          iOptH = ior(128,iand(iOptH,32))
        else
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '************ ERROR ****************'
          write(Lu,*) 'Unsupported Hessian update method: ',trim(Char)
          write(Lu,*) '***********************************'
          call Quit_OnUserError()
        end if

      case ('HWRS')
        !                                                              *
        !***** HWRS ****************************************************
        !                                                              *
        HWRS = .true.

      case ('INTE')
        !                                                              *
        !***** INTE ****************************************************
        !                                                              *
        ! Read the internal coordinate specification.

        New_Line = 1
        Lu_UDIC = 91
        FilNam = 'UDIC'
        call molcas_open(Lu_UDIC,FilNam)
        rewind(Lu_UDIC)

        ! mInt is the number of internal coordinates you will define.
        ! mInt = nDimBC - mTROld
        ! Subroutine DefInt defines the B matrix.
        ! The matrix B relates a shift in an internal coordinate to
        ! shifts in cartesian coordinates,
        !
        !           |dq> = B |dx>
        !                      =
        ! and has the dimension (3*nsAtom x mInt).
        inte: do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          if (Key(1:4) == 'END ') then
            close(Lu_UDIC)
            exit inte
          end if

          ! Here is a fix because auto will break up the lines if there is an
          ! equal sign in the input.

          ! Lines with VARY or FIX doesn't have equal signs

          if (Key(1:4) == 'VARY') nBVec = iRow
          if ((Key(1:4) == 'VARY') .or. (Key(1:3) == 'FIX') .or. (Key(1:4) == 'ROWH')) then
            New_Line = 0
          end if

          do
            if (New_Line /= 1) exit
            if (index(Key,'=') == 0) call FixEqualSign2(Key,LuRd,Lu_UDIC,iRow,New_Line)
            if (New_Line == 2) then
              close(Lu_UDIC)
              exit inte
            end if
          end do

          iRow = iRow+1

          write(Lu_UDIC,'(A)') Key

          ! If this line does not have a continuation the next line should
          ! have a equal sign!
          if (index(Key,'&') == 0) New_Line = 1
        end do inte

      case ('IRC ')
        !                                                              *
        !***** IRC *****************************************************
        !                                                              *
        call Qpg_iScalar('IRC',Found)
        if (Found) then
          call Get_iScalar('IRC',IRC)
        else
          IRC = 1
          call Put_iScalar('IRC',IRC)
        end if
        MEP = .true.
        rMEP = .false.

      case ('ITER')
        !                                                              *
        !***** ITER ****************************************************
        !                                                              *
        ! read max iterations

        Char = Get_Ln(LuRd)
        call Get_I1(1,iTmp)
        MxItr = min(iTmp,MxItr)

      case ('KRIG')
        !                                                              *
        !***** KRIG ****************************************************
        !                                                              *
        ! Activate Kriging

        Kriging = .true.
        Line_Search = .false.

      case ('LAST')
        !                                                              *
        !***** LAST ****************************************************
        !                                                              *
        do
          Char = Get_Ln(LuRd)
          Char = adjustl(Char)
          if ((Char /= BLine) .and. (char(1:1) == '*')) exit
        end do
        call UpCase(Char)
        call Put_cArray('LastEnergyMethod',Char,8)

      case ('LINE')
        !                                                              *
        !***** LINE ****************************************************
        !                                                              *
        Line_Search = .true.

      case ('MAXS')
        !                                                              *
        !***** MAXS ****************************************************
        !                                                              *
        do
          Char = Get_Ln(LuRd)
          if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
        end do
        call Get_F1(1,Beta)

      case ('MAXD')
        !                                                              *
        !***** MAXD ****************************************************
        !                                                              *
        do
          Char = Get_Ln(LuRd)
          if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
        end do
        call Get_F1(1,Beta_Disp)
        Manual_Beta = .true.

      case ('MEP-','MEP ')
        !                                                              *
        !***** MEP-/MEP  ***********************************************
        !                                                              *
        MEP = .true.
        rMEP = .false.

      case ('MEPA','IRCA')
        !                                                              *
        !***** MEPA/IRCA ***********************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call UpCase(Char)
        if (char(1:2) == 'GS') then
          MEP_Algo = 'GS'
        else if (char(1:2) == 'MB') then
          MEP_Algo = 'MB'
        else
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '************* ERROR ************'
          write(Lu,*) ' Unrecognized MEP/IRC algorithm.'
          write(Lu,*) '********************************'
          call Quit_OnUserError()
        end if

      case ('MEPC','IRCC')
        !                                                              *
        !***** MEPC/IRCC ***********************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_F1(1,ThrMEP)
        ThrMEP = max(Zero,ThrMEP)

      case ('MEPS','IRCS')
        !                                                              *
        !***** MEPStep/IRCStep *****************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call UpCase(Char)
        call Get_F1(1,dMEPStep)

        ! Note that according to the Gonzalez-Schlegel method, only half
        ! this step is used in the constraint

        if (index(Char,'ANGSTROM') /= 0) dMEPStep = dMEPStep/angstr

      case ('MEPT','IRCT')
        !                                                              *
        !***** MEPT/IRCT ***********************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call UpCase(Char)
        if (char(1:6) == 'SPHERE') then
          MEP_Type = 'SPHERE'
        else if (char(1:5) == 'PLANE') then
          MEP_Type = 'TRANSVERSE'
        else
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '********** ERROR **********'
          write(Lu,*) ' Unrecognized MEP/IRC type.'
          write(Lu,*) '***************************'
          call Quit_OnUserError()
        end if

      case ('MODE')
        !                                                              *
        !***** MODE ****************************************************
        !                                                              *
        ! Mode following algorithm

        do
          Char = Get_Ln(LuRd)
          if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
        end do
        call Get_I1(1,mode)

      case ('MXMI')
        !                                                              *
        !***** MXMI ****************************************************
        !                                                              *
        ! Maximum number of Iterations for the Kriging method

        Char = Get_Ln(LuRd)
        call Get_I1(1,Max_Microiterations)

      case ('NMEP','NIRC')
        !                                                              *
        !***** NMEP/NIRC ***********************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call Get_I1(1,nMEP)
        nMEP = min(max(nMEP,1),MaxItr)

      case ('NEWT')
        !                                                              *
        !***** NEWT ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**0-2**1-2**2-2**3
        iOptC = ior(2**0,iand(iOptC,Mask))

      case ('NOEM')
        !                                                              *
        !***** NOEM ****************************************************
        !                                                              *
        eMEPTest = .false.

      case ('NOFA')
        !                                                              *
        !***** NOFA ****************************************************
        !                                                              *
        ! Deactivate fallback to conventional

        Fallback = .false.

      case ('NOHW')
        !                                                              *
        !***** NOHW ****************************************************
        !                                                              *
        HWRS = .false.

      case ('NOLA')
        !                                                              *
        !***** NOLA ****************************************************
        !                                                              *
        CallLast = .false.

      case ('NOLI')
        !                                                              *
        !***** NOLI ****************************************************
        !                                                              *
        Line_Search = .false.

      case ('NOWB')
        !                                                              *
        !***** NO VDWB VdW correction both coordinate and Hessian ******
        !                                                              *
        Mask = iSetAll-2**10-2**11
        iOptC = iand(Mask,iOptC)

      case ('NOWC')
        !                                                              *
        !***** NO VDWB VdW correction for coordinate only **************
        !                                                              *
        Mask = iSetAll-2**11
        iOptC = iand(Mask,iOptC)

      case ('NOWH')
        !                                                              *
        !***** NO VDWB VdW correction for Hessian only *****************
        !                                                              *
        Mask = iSetAll-2**10
        iOptC = iand(Mask,iOptC)

      case ('NUME')
        !                                                              *
        !***** NUME ****************************************************
        !                                                              *
        lNmHss = .true.

      case ('OLDF')
        !                                                              *
        !***** OLDF ****************************************************
        !                                                              *
        lOld = .true.

      case ('PRFC')
        !                                                              *
        !***** PRFC ****************************************************
        !                                                              *
        ! set nPrint to print internal coordinates and hessian

        nPrint(21) = 6  ! Eigen-/Values/Vectors of the Hessian (diagmtrx)
        nPrint(116) = 6 ! Internal Forces (rlxctl)
        if (.not. Request_Alaska) nPrint(30) = 6 ! Coords & Forces (defint)
        nPrint(122) = 6 ! Auto-Defined Internal coordinates (printq_sl)

      case ('PRIN')
        !                                                              *
        !***** PRIN ****************************************************
        !                                                              *
        do
          Char = Get_Ln(LuRd)
          call UpCase(Char)
          if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
        end do
        call Get_I1(1,mPrint)
        do i=1,mPrint
          do
            Char = Get_Ln(LuRd)
            call UpCase(Char)
            if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
          end do
          call Get_I1(1,iRout)
          call Get_I1(2,kPrint)
          nPrint(iRout) = kPrint
        end do

      case ('RATI')
        !                                                              *
        !***** RATI ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**0-2**1-2**2-2**3
        iOptC = ior(2**3,iand(iOptC,Mask))

      case ('REDU')
        !                                                              *
        !***** REDU ****************************************************
        !                                                              *
        Redundant = .true.

      case ('REAC')
        !                                                              *
        !***** REAC ****************************************************
        !                                                              *
        Explicit_IRC = .true.
        call mma_allocate(TmpRx,3*nsAtom,Label='TmpRx')
        call Read_v(LuRd,TmpRx,1,3*nsAtom,1,iErr)
        if (IErr /= 0) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '********** ERROR ***********************'
          write(Lu,*) ' Error while reading the Reaction vector'
          write(Lu,*) '****************************************'
          call Quit_OnUserError()
        end if

      case ('REFE')
        !                                                              *
        !***** REFE ****************************************************
        !                                                              *
        call mma_allocate(RefGeo,3,nsAtom,Label='RefGeo')
        call Read_v(LuRd,RefGeo,1,3*nsAtom,1,iErr)
        if (iErr /= 0) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '************ ERROR ***************'
          write(Lu,*) 'Error reading reference structure.'
          write(Lu,*) '**********************************'
          call Quit_OnUserError()
        end if

      case ('RHID')
        !                                                              *
        !***** rHidden *************************************************
        !                                                              *
        Key = Get_Ln(LuRd)
        call UpCase(Key)
        call Get_F1(1,rHidden)
        if (rHidden < Zero) then
          call WarningMessage(2,'Error in RdCtl_Slapaf')
          write(Lu,*)
          write(Lu,*) '************ ERROR *****************'
          write(Lu,*) 'Error reading rHidden. Should be >0.'
          write(Lu,*) '************************************'
          call Quit_OnUserError()
        end if
        if (index(Key,'ANGSTROM') /= 0) rHidden = rHidden/angstr

      case ('RMEP')
        !                                                              *
        !***** rMEP ****************************************************
        !                                                              *
        rMEP = .true.
        MEP = .false.

      case ('RS-P')
        !                                                              *
        !***** RS-P ****************************************************
        !                                                              *
        Mask = iSetAll
        Mask = Mask-2**9
        iOptC = iand(iOptC,Mask)

      case ('RTRN')
        !                                                              *
        !***** RTRN ****************************************************
        !                                                              *
        Char = Get_Ln(LuRd)
        call UpCase(Char)
        call Get_I1(1,Max_Center)
        call Get_F1(2,rtrnc)
        if (index(Char,'ANGSTROM') /= 0) Rtrnc = Rtrnc/angstr

      case ('SUPS')
        !                                                              *
        !***** SUPS ****************************************************
        !                                                              *
        ! Introduce supersymmetry
        ! Input format
        ! nsg                (number of super groups)
        ! Repeat nsg times
        ! nmem, (ind.., i = 1, nmem)

        Char = Get_Ln(LuRd)
        call Get_I1(1,nSupSy)
        call mma_allocate(nSup,NSUPSY,Label='nSup')
        call mma_allocate(Atom,nsAtom,Label='Atom')
        jStrt = 1
        do i=1,nSupSy
          read(LuRd,*,iostat=istatus) nSup(i),(Atom(j),j=jStrt,jStrt+nSup(i)-1)
          if (istatus > 0) then
            call WarningMessage(2,'Error in RdCtl_Slapaf')
            write(Lu,*)
            write(Lu,*) '************ ERROR ****************'
            write(Lu,*) 'Error while reading supersymmetry.'
            write(Lu,*) '***********************************'
            call Quit_OnUserError()
          end if
          jStrt = jStrt+nSup(i)
        end do

      case ('TFOF')
        !                                                              *
        !***** TFOF ****************************************************
        !                                                              *
        ! adding energy to the last energy value of the base line
        ! This option supersedes any value assigned to blAI and mblAI

        Char = Get_Ln(LuRd)
        call Get_F1(1,blavAI)

      case ('THER')
        !                                                              *
        !***** THER ****************************************************
        !                                                              *
        lNmHss = .true.
        lTherm = .true.
        Char = Get_Ln(LuRd)
        call Get_I1(1,nsRot)
        Char = Get_Ln(LuRd)
        call Get_F1(1,UserP)
        do
          Char = Get_Ln(LuRd)
          call UpCase(Char)
          if (char(1:4) == 'END ') then
            if (nUserPT == 0) then
              nUserPT = 1
              UserT(1) = 298.15d0
            end if
            exit
          end if
          nUserPT = nUserPT+1
          call Get_F1(1,UserT(nUserPT))
        end do

      case ('THRS')
        !                                                              *
        !***** THRS ****************************************************
        !                                                              *
        ! read the gradient threshold

        Char = Get_Ln(LuRd)
        call Get_F1(1,ThrEne)
        call Get_F1(2,ThrGrd)
        ThrInp = .true.

      case ('TOLE')
        !                                                              *
        !***** TOLE ****************************************************
        !                                                              *
        ! read the constraints threshold

        Char = Get_Ln(LuRd)
        call Get_F1(1,ThrCons)
        ThrCons = abs(ThrCons)

      case ('TRAC')
        !                                                              *
        !***** TRAC ****************************************************
        !                                                              *
        Track = .true.

      case ('TS  ')
        !                                                              *
        !***** TS   ****************************************************
        !                                                              *
        Mask = iSetAll-2**7
        iOptC = iand(Mask,iOptC)

      case ('TSCO')
        !                                                              *
        !***** TSCO ****************************************************
        !                                                              *
        LuTS = 20
        FilNam = 'TSC'
        LuTS = IsFreeUnit(LuTS)
        call Molcas_Open(LuTS,FilNam)
        do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          Key = adjustl(Key)
          write(LuTS,'(A)') trim(Key)
          if (Key(1:4) == 'END') exit
        end do
        close(LuTS)
        TSConstraints = .true.

      case ('VDWB')
        !                                                              *
        !***** VDWB VdW correction both coordinate and Hessian *********
        !                                                              *
        iOptC = ior(1024,iOptC)
        iOptC = ior(2048,iOptC)

      case ('VDWC')
        !                                                              *
        !***** VDWB VdW correction for coordinate only *****************
        !                                                              *
        iOptC = ior(2048,iOptC)

      case ('VDWH')
        !                                                              *
        !***** VDWB VdW correction for Hessian only ********************
        !                                                              *
        iOptC = ior(1024,iOptC)

      case ('WIND')
        !                                                              *
        !***** WIND ****************************************************
        !                                                              *
        do
          Char = Get_Ln(LuRd)
          call UpCase(Char)
          if ((Char /= BLine) .and. (char(1:1) /= '*')) exit
        end do
        call Get_I1(1,nWndw)

      case default
        call WarningMessage(2,'Error in RdCtl_Slapaf')
        if (char(1:1) == ' ') then
          write(Lu,*) ' RdCtl_Slapaf: Command line starts with a blank.'
        else
          write(Lu,*)
          write(Lu,*) ' *********** ERROR ***********'
          write(Lu,*) ' The program has been supplied'
          write(Lu,*) ' with an unknown command.     '
          write(Lu,*) ' *****************************'
        end if
        write(Lu,'(A)') Char
        call Quit_OnUserError()

    end select
  end do
  !                                                                    *
  !*********************************************************************
  !**********************   End of input section   *********************
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Now start fixing constraints.
! First get the external constraints.

if (External_UDC) then
  call Merge_Constraints('UDC.Gateway','','UDC',nLambda,iRow_c)
else
  call Merge_Constraints('','UDC','UDC',nLambda,iRow_c)
end if

! Initial preprocessing

if (iRow_c > 1) then
  Lu_UDC = IsFreeUnit(20)
  call Molcas_Open(Lu_UDC,'UDC')
  call Preprocess_UDC(Lu_UDC,iPrint)
  close(Lu_UDC)
else
  NADC = .false.
end if

! Add NAC if needed

if (NADC) then
  Lu_UDCTMP = IsFreeUnit(20)
  call Molcas_Open(Lu_UDCTMP,'UDCTMP')
  write(Lu_UDCTMP,*) 'NADC = NAC'
  write(Lu_UDCTMP,*) 'VALUE'
  write(Lu_UDCTMP,*) 'NADC = 0.0'
  write(Lu_UDCTMP,*) 'END'
  close(Lu_UDCTMP)
  call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
end if

! Add MEP/IRC if needed

if (MEP .or. rMEP .or. (abs(IRC) == 1)) then
  if (abs(IRC) == 1) then
    MEPLab = 'IRC'
  else
    MEPLab = 'MEP'
  end if
  if (MEPCons .and. (.not. Expert)) then
    call WarningMessage(2,'Error in RdCtl_Slapaf')
    write(Lu,*)
    write(Lu,*) '***************** ERROR ********************'
    write(Lu,*) ' There is a '//trim(Mep_Type)//' constraint that may'
    write(Lu,*) ' conflict with '//MEPLab//' calculations.'
    write(Lu,*) ' You should not explictly specify this constraint,'
    write(Lu,*) ' but just rely on '//MEPLab//'Step/'//MEPLab//'Type keywords.'
    write(Lu,*) ' If you really know what you are doing, you'
    write(Lu,*) ' can use the EXPERT keyword.'
    write(Lu,*) '********************************************'
    call Quit_OnUserError()
  end if
  WeightedConstraints = .true.
  Valu = dMEPStep
  if (MEP_Type == 'SPHERE') Valu = abs(Valu)
  if (MEP .and. (MEP_Algo == 'GS')) Valu = Half*Valu
  if (rMEP) Valu = max(dble(iMEP+1),One)*Valu
  Lu_UDCTMP = IsFreeUnit(20)
  call Molcas_Open(Lu_UDCTMP,'UDCTMP')
  write(Lu_UDCTMP,*) MEPLab//' = '//MEP_Type
  write(Lu_UDCTMP,*) 'VALUE'
  write(Lu_UDCTMP,*) MEPLab//' = ',Valu
  write(Lu_UDCTMP,*) 'END'
  close(Lu_UDCTMP)
  call Merge_Constraints('UDC','UDCTMP','UDC',nLambda,iRow_c)
  Beta = min(Beta,abs(Valu))
  MEPnum = nLambda
end if

! Final fixes

call Fix_UDC(iRow_c,nLambda,nsAtom,nStab,.true.)
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate some variables which can only be set after the input has
! been read.

if ((.not. ThrInp) .and. (.not. Baker)) ThrEne = Zero

call Init2()

! Gradients are not needed at the first iteration of a numerical
! Hessian procedure (and only that, i.e. MxItr=0)
FirstNum = (allocated(mRowH) .or. lNmHss .or. Cubic) .and. (Iter == 1) .and. (MxItr == 0)

if ((SuperName == 'slapaf') .and. (.not. FirstNum)) then

  if (Track) then
    call Process_Track()
  else
    call Put_iArray('Root Mapping',iDum,0)
  end if

  call Process_Gradients()

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Put in the "Reaction vector" in Cartesians.
! Priority order:
! 1) Explicit by user input (REAC keyword)
! 2) Found on RunOld
! 3) Found on RunFile

if (abs(IRC) == 1) then

  ! If this is the first macro iteration in the IRC search then
  ! pick up the reaction vector.

  if (Explicit_IRC .and. (iMEP == 0)) then
    ! Case 1)
    call dcopy_(3*nsAtom,TmpRx,1,MF,1)
  else if (iMEP == 0) then
    call NameRun('RUNOLD')
    call qpg_dArray('Reaction Vector',Found,nRx)
    !write(6,*) 'RUNOLD: Found=',Found
    if (Found) then
      ! Case 2)
      call Get_dArray('Reaction Vector',MF,3*nsAtom)
      call NameRun('#Pop')
    else
      call NameRun('#Pop')
      call qpg_dArray('Reaction Vector',Found,nRx)
      !write(6,*) 'RUNFILE: Found=',Found
      if (Found) then
        ! Case 3)
        call Get_dArray('Reaction Vector',MF,3*nsAtom)
      else
        call WarningMessage(2,'Error in RdCtl_Slapaf')
        write(6,*)
        write(6,*) '************ ERROR **************'
        write(6,*) 'IRC calculation but no IRC vector'
        write(6,*) '*********************************'
        call Quit_OnUserError()
      end if
    end if
  end if

  ! Fix the direction forward/backwards

  if ((iMEP == 0) .and. (iRC == -1)) call DScal_(3*nsAtom,-1.0d0,MF,1)
  if ((iMEP == 0) .and. (MEP_Type == 'TRANSVERSE')) call Put_dArray('Transverse',MF,3*nsAtom)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (FindTS .and. (.not. TSConstraints)) then
  call SysWarnMsg('RdCtl_Slapaf','WARNING:','FindTS specified, but no TSConstraints. '// &
                  'It is highly recommended to use TSConstraints in SLAPAF instead of (or in addition to) global constraints '// &
                  'when using FindTS. TSConstraints will be lifted in the final TS search.')
end if
TSConstraints = TSConstraints .and. FindTS

if ((MEP .or. rMEP) .and. (.not. Request_Alaska)) then

  ! If no initial direction given, use the gradient (force)

  call qpg_dArray('Transverse',Found,nRP)
  if (.not. Found) then
    ! Assume the initial reaction vector is already massaged
    if (Explicit_IRC) then
      call Put_dArray('Transverse',TmpRx,3*nsAtom)
    else
      ! The direction is given by the gradient, but in weighted coordinates
      call mma_allocate(Dir,3,nsAtom,Label='Dir')
      do iAtom=1,nsAtom
        xWeight = Weights(iAtom)
        do ixyz=1,3
          Dir(ixyz,iAtom) = Gx(ixyz,iAtom,iter)/xWeight
        end do
      end do
      call Put_dArray('Transverse',Dir,3*nsAtom)
      call mma_deallocate(Dir)
    end if
  end if
end if

if (Explicit_IRC) call mma_deallocate(TmpRx)

! Activate MPS update of Hessian if FindTS

if (FindTS) then

  if (iand(iOptH,64) == 64) then
    iOptH = ior(64,iand(iOptH,32)) ! EU
  else if (iand(iOptH,128) == 128) then
    iOptH = ior(128,iand(iOptH,32)) ! TS-BFGS
  else
    iOptH = ior(16,iand(iOptH,32)) ! MSP
  end if
  iOptC = ior(iOptC,4096)

  ! Increase the update window so that we will not lose the update
  ! which generated the negative curvature.

  nWndw = 4*nWndw
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify some options if TS search

if (iand(iOptC,128) /= 128) then
  if (iand(iOptH,8) /= 8) iOptH = ior(16,iand(iOptH,32)) ! MSP
  Line_search = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For TS optimization with the Saddle method set update to MSP

call qpg_dArray('Saddle',Found,nSaddle)
if (Found .and. (nSaddle /= 0)) then
  call mma_allocate(Tmp,nSaddle,Label='Tmp')
  call Get_dArray('Saddle',Tmp,nSaddle)
  HSR0 = Tmp(nSaddle-2)
  HSR = Tmp(nSaddle-1)
  Update = Tmp(nSaddle)
  if (Update == 2.0d0) then

    ! Enable FindTS procedure

    !write(6,*) 'Enable FindTS procedure'
    if (iand(iOptH,8) /= 8) iOptH = ior(16,iand(iOptH,32)) ! MSP
    nWndw = 4*nWndw
    ! make it look as if this were FindTS with constraints
    FindTS = .true.
    TSConstraints = .true.
    iOptC = ior(iOptC,4096)
    iOptC = ior(iOptC,8192)
    Beta = 0.1d0

  else

    ! Normal constrained optimization with a reduced threshold.
    ! Let the threshold be somewhat tighter as we are close to the TS.

    if ((HSR/HSR0 < 0.20d0) .or. (HSR < 0.20d0)) then
      !ThrGrd = 0.0003D0
      Beta = 0.1d0
    else
      !ThrGrd = 0.003D0
      ThrGrd = Ten*ThrGrd
    end if

    ! Add the constraints from the Saddle method

    call Merge_Constraints('UDC','UDC.Saddle','UDC',nLambda,iRow_c)

  end if
  call mma_deallocate(Tmp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
nLbl = max(3*nsAtom+nLambda,iRow+iRow_c)
call mma_allocate(Lbl,nLbl,Label='Lbl')
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify some options if constraints are part of the calculation.

if ((nLambda > 0) .or. TSConstraints) then
  iOptC = ior(iOptC,256) ! Constraints
  Line_search = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! No iterations set iOptC=0

if (MxItr == 0) iOptC = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Activate some additional printing for numerical Hessian

!GGd: Coherency with patch 7.1.615 !      If (lNmHss) nPrint(122) = 10
!                                                                      *
!***********************************************************************
!                                                                      *
! Do some preprocessing due to input choice

if (Request_Alaska) nPrint(51) = 0
call PrePro(nsAtom,Cx(1,1,iter))
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of Kriging we use a sorting step in update_sl. For this
! to work we need the values of the internal coordinates for more
! points than the window size. Here we increase it with a factor of
! 2 temporarily. The sorted list will still be of the original size.
! However, the default window for kriging is twice as large as
! for conventional calculations.

if (Kriging) then
  nWndw = 4*nWndw  ! 2*2=4

  ! No micro iterations the first MEP iteration

  if ((MEP .or. rMEP) .and. (iter == 1)) Max_Microiterations = 0

  ! Reduce default maximum dispersion during the initial
  ! stage of a FindTS calculation: we don't want to fulfil the
  ! constraints too early

  call Qpg_iScalar('TS Search',Found)
  if (Found) call Get_lScalar('TS Search',Found)
  if (FindTS .and. (.not. (Found .or. Manual_Beta))) then
    Beta_Disp = 0.1d0
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out input parameters, No output if we didn't find the proper
! gradients on the runfile. We will come back!!!

if (SuperName == 'slapaf') then
  if (.not. Request_Alaska) call WrInp_sl()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
User_Def = iRow /= 0
if (.not. User_Def) nBVec = 1
!                                                                      *
!***********************************************************************
!                                                                      *
if (lNmHss .and. allocated(mRowH)) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '**************************************************'
  write(Lu,*) ' ERROR: NUMErical and ROWH are mutually exclusive '
  write(Lu,*) '**************************************************'
  call Quit_OnUserError()
end if
if (lCtoF .and. User_Def) then
  call WarningMessage(2,'Error in RdCtl_Slapaf')
  write(Lu,*)
  write(Lu,*) '******************************************'
  write(Lu,*) ' ERROR: CtoF and User-defined Coordinates '
  write(Lu,*) '        are mutually exclusive.           '
  write(Lu,*) '******************************************'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine RdCtl_Slapaf
