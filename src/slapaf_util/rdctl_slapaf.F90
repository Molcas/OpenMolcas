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

use Symmetry_Info, only: Symmetry_Info_Get
use Slapaf_Info, only: Atom, Baker, Beta, Beta_Disp, CallLast, CnstWght, Coor, Cubic, Curvilinear, Cx, ddV_Schlegel, Delta, &
                       dMEPStep, eMEPTest, Fallback, FindTS, GNrm_Threshold, GradRef, Gx, HrmFrq_Show, HWRS, iOptC, iOptH, IRC, &
                       iRow, iRow_c, isFalcon, Iter, Lbl, lCtoF, lDoubleIso, Line_Search, lNmHss, lOld, lTherm, Max_Center, &
                       MaxItr, MEP, MEP_Algo, MEP_Type, MEPCons, MEPNum, MF, Mode, mRowH, MxItr, NADC, nBVec, nLambda, nMEP, &
                       nsRot, nStab, nSup, nUserPT, nWndw, Redundant, RefGeo, Request_Alaska, rFuzz, rHidden, rMEP, RtRnc, &
                       ThrCons, ThrEne, ThrGrd, ThrMEP, Track, TSConstraints, User_Def, UserP, UserT, WeightedConstraints, Weights
use kriging_mod, only: blavAI, Kriging, Max_Microiterations, nD_In
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Ten, Half, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
logical(kind=iwp), intent(in) :: Dummy_Call
#include "print.fh"
integer(kind=iwp) :: i, iAtom, iDum(1), iErr, iMEP, iNull, iOff_Iter, iPrint, iRout, istatus, iTmp, j, jStrt, kPrint, Lu, Lu_UDC, &
                     Lu_UDCTMP, Lu_UDIC, LuRd, LuTS, mPrint, NewLine, nLbl, nRP, nRx, nSaddle, nsAtom, nSupSy
real(kind=wp) :: HSR, HSR0, Update, Valu, xWeight
logical(kind=iwp) :: Expert, Explicit_IRC, External_UDC, FirstNum, Found, Manual_Beta, ThrInp
character(len=180) :: Chr, Key
character(len=16) :: FilNam
character(len=3) :: MEPLab
real(kind=wp), allocatable :: DIR(:,:), Tmp(:), TmpRx(:)
integer(kind=iwp), external :: IsFreeUnit
character(len=180), external :: Get_Ln

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
Expert = .false.
Lu = u6
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate some parameters

call Symmetry_Info_Get()
call Init_Slapaf()
nsAtom = size(Coor,2)
iPrint = nPrint(iRout)

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
    Chr = Get_Ln(LuRd)
    call UpCase(Chr)
    !write(Lu,'(A)') Chr
    !write(Lu,*) iOptC
    if (Chr == '') cycle
    if (Chr(1:1) == '*') cycle
    select case (Chr(1:4))
      !case ('AIL ')
      !  !                                                              *
      !  !***** AIL  ****************************************************
      !  !                                                              *
      !  ! Width limits of the Matern function
      !
      !  Chr = Get_Ln(LuRd)
      !  call Get_F(1,lb,3)

      !case ('AIP ')
      !  !                                                              *
      !  !***** AIP  ****************************************************
      !  !                                                              *
      !  ! Parameter of differentiability for Matern function
      !
      !  Chr = Get_Ln(LuRd)
      !  call Get_F1(1,pAI)
      !  if ((pAI > 3) .or. (pAI < 1)) anMd = .false.

      !case ('AISP')
      !  !                                                              *
      !  !***** AISP ****************************************************
      !  !                                                              *
      !  ! Defining the number of source points for the AI method
      !
      !  Chr = Get_Ln(LuRd)
      !  call Get_I1(1,nspAI)

      !case ('AIMD')
      !  !                                                              *
      !  !***** AIMD ****************************************************
      !  !                                                              *
      !  ! Analytical or numerical Matern derivatives
      !
      !  Chr = Get_Ln(LuRd)
      !  if ((Chr == 'False') .or. (Chr == 'false')) then
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
      !  Chr = Get_Ln(LuRd)
      !  call Get_F1(1,Thr_microiterations)

      !case ('AIBL')
      !  !                                                              *
      !  !***** AIBL ****************************************************
      !  !                                                              *
      !  ! Base line modification value to not ordinary
      !  ! (Trend Function on GEK)
      !
      !  Chr = Get_Ln(LuRd)
      !  call Get_F1(1,blvAI)
      !  blAI = .true.

      !case ('AIMB')
      !  !                                                              *
      !  !***** AIMB ****************************************************
      !  !                                                              *
      !  ! Base line modification value to maximum value of the Energy
      !  ! This option supersedes any value assigned to blAI
      !
      !  Chr = Get_Ln(LuRd)
      !  mblAI = .true.

      !case ('L-VA')
      !  !                                                              *
      !  !***** L-VA ****************************************************
      !  !                                                              *
      !  ! Change the l value of the GEK.
      !
      !  Chr = Get_Ln(LuRd)
      !  Set_l = .true.
      !  call Get_F1(1,Value_l)
      !  call Qpg_dScalar('Value_l',Found)
      !  if (.not. Found) call Put_dScalar('Value_l',Value_l)

      case ('NDEL')
        !                                                              *
        !***** NDELta **************************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
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
        iOptC = ibclr(iOptC,0)
        iOptC = ibset(iOptC,1)
        iOptC = ibclr(iOptC,2)
        iOptC = ibclr(iOptC,3)

      case ('C2-D')
        !                                                              *
        !***** C2-D ****************************************************
        !                                                              *
        iOptC = ibclr(iOptC,0)
        iOptC = ibclr(iOptC,1)
        iOptC = ibset(iOptC,2)
        iOptC = ibclr(iOptC,3)

      case ('CART')
        !                                                              *
        !***** CART ****************************************************
        !                                                              *
        CurviLinear = .false.

      case ('CNWE')
        !                                                              *
        !***** CNWE ****************************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
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
        NewLine = 1
        lCtoF = .true.
        Lu_UDIC = 91
        FilNam = 'UDIC'
        call molcas_open(Lu_UDIC,FilNam)
        rewind(Lu_UDIC)
        Key = Get_Ln(LuRd)
        call UpCase(Key)
        call FixEqualSign2(Key,LuRd,Lu_UDIC,iNull,NewLine)
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
        Chr = Get_Ln(LuRd)
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
        iOptC = ibset(iOptC,4)
        iOptC = ibclr(iOptC,5)
        iOptC = ibclr(iOptC,6)

      case ('DXG ','GDX ')
        !                                                              *
        !***** DXG  ****************************************************
        !                                                              *
        iOptC = ibclr(iOptC,4)
        iOptC = ibset(iOptC,5)
        iOptC = ibclr(iOptC,6)

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
        Chr = Get_Ln(LuRd)
        call UpCase(Chr)
        call Get_F1(1,rFuzz)
        if (index(Chr,'ANGSTROM') /= 0) rFuzz = rFuzz/Angstrom
        rFuzz = max(rFuzz,1.0e-3_wp)

      case ('GG  ')
        !                                                              *
        !***** GG   ****************************************************
        !                                                              *
        iOptC = ibclr(iOptC,4)
        iOptC = ibclr(iOptC,4)
        iOptC = ibset(iOptC,4)

      case ('GNRM')
        !                                                              *
        !***** GNRM ****************************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
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
        Chr = Get_Ln(LuRd)
        read(Chr,*) Chr
        call UpCase(Chr)
        select case (Chr)
          !case ('MEYER') then
          !  iOptH = ibset(0,0)
          !case ('BP') then
          !  iOptH = ibset(0,1)
          case ('BFGS')
            iOptH = ibset(0,2)
          case ('NONE')
            iOptH = ibset(0,3)
          case ('MSP')
            iOptH = ibset(0,4)
          case ('EU')
            iOptH = ibset(0,5)
          case ('TS-BFGS')
            iOptH = ibset(0,6)
          case default
            call WarningMessage(2,'Error in RdCtl_Slapaf')
            write(Lu,*)
            write(Lu,*) '************ ERROR ****************'
            write(Lu,*) 'Unsupported Hessian update method: ',trim(Chr)
            write(Lu,*) '***********************************'
            call Quit_OnUserError()
        end select

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

        NewLine = 1
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
          if ((Key(1:4) == 'VARY') .or. (Key(1:3) == 'FIX') .or. (Key(1:4) == 'ROWH')) NewLine = 0

          do
            if (NewLine /= 1) exit
            if (index(Key,'=') == 0) call FixEqualSign2(Key,LuRd,Lu_UDIC,iRow,NewLine)
            if (NewLine == 2) then
              close(Lu_UDIC)
              exit inte
            end if
          end do

          iRow = iRow+1

          write(Lu_UDIC,'(A)') Key

          ! If this line does not have a continuation the next line should
          ! have a equal sign!
          if (index(Key,'&') == 0) NewLine = 1
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

        Chr = Get_Ln(LuRd)
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
          Chr = Get_Ln(LuRd)
          Chr = adjustl(Chr)
          if ((Chr /= '') .and. (Chr(1:1) == '*')) exit
        end do
        call UpCase(Chr)
        call Put_cArray('LastEnergyMethod',Chr,8)

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
          Chr = Get_Ln(LuRd)
          if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
        end do
        call Get_F1(1,Beta)

      case ('MAXD')
        !                                                              *
        !***** MAXD ****************************************************
        !                                                              *
        do
          Chr = Get_Ln(LuRd)
          if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
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
        Chr = Get_Ln(LuRd)
        call UpCase(Chr)
        if (Chr(1:2) == 'GS') then
          MEP_Algo = 'GS'
        else if (Chr(1:2) == 'MB') then
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
        Chr = Get_Ln(LuRd)
        call Get_F1(1,ThrMEP)
        ThrMEP = max(Zero,ThrMEP)

      case ('MEPS','IRCS')
        !                                                              *
        !***** MEPStep/IRCStep *****************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
        call UpCase(Chr)
        call Get_F1(1,dMEPStep)

        ! Note that according to the Gonzalez-Schlegel method, only half
        ! this step is used in the constraint

        if (index(Chr,'ANGSTROM') /= 0) dMEPStep = dMEPStep/Angstrom

      case ('MEPT','IRCT')
        !                                                              *
        !***** MEPT/IRCT ***********************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
        call UpCase(Chr)
        if (Chr(1:6) == 'SPHERE') then
          MEP_Type = 'SPHERE'
        else if (Chr(1:5) == 'PLANE') then
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
          Chr = Get_Ln(LuRd)
          if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
        end do
        call Get_I1(1,mode)

      case ('MXMI')
        !                                                              *
        !***** MXMI ****************************************************
        !                                                              *
        ! Maximum number of Iterations for the Kriging method

        Chr = Get_Ln(LuRd)
        call Get_I1(1,Max_Microiterations)

      case ('NMEP','NIRC')
        !                                                              *
        !***** NMEP/NIRC ***********************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
        call Get_I1(1,nMEP)
        nMEP = min(max(nMEP,1),MaxItr)

      case ('NEWT')
        !                                                              *
        !***** NEWT ****************************************************
        !                                                              *
        iOptC = ibset(iOptC,0)
        iOptC = ibclr(iOptC,1)
        iOptC = ibclr(iOptC,2)
        iOptC = ibclr(iOptC,3)

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
        iOptC = ibclr(iOptC,10)
        iOptC = ibclr(iOptC,11)

      case ('NOWC')
        !                                                              *
        !***** NO VDWB VdW correction for coordinate only **************
        !                                                              *
        iOptC = ibclr(iOptC,11)

      case ('NOWH')
        !                                                              *
        !***** NO VDWB VdW correction for Hessian only *****************
        !                                                              *
        iOptC = ibclr(iOptC,10)

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
          Chr = Get_Ln(LuRd)
          call UpCase(Chr)
          if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
        end do
        call Get_I1(1,mPrint)
        do i=1,mPrint
          do
            Chr = Get_Ln(LuRd)
            call UpCase(Chr)
            if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
          end do
          call Get_I1(1,iRout)
          call Get_I1(2,kPrint)
          nPrint(iRout) = kPrint
        end do

      case ('RATI')
        !                                                              *
        !***** RATI ****************************************************
        !                                                              *
        iOptC = ibclr(iOptC,0)
        iOptC = ibclr(iOptC,1)
        iOptC = ibclr(iOptC,2)
        iOptC = ibset(iOptC,3)

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
        if (index(Key,'ANGSTROM') /= 0) rHidden = rHidden/Angstrom

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
        iOptC = ibclr(iOptC,9)

      case ('RTRN')
        !                                                              *
        !***** RTRN ****************************************************
        !                                                              *
        Chr = Get_Ln(LuRd)
        call UpCase(Chr)
        call Get_I1(1,Max_Center)
        call Get_F1(2,rtrnc)
        if (index(Chr,'ANGSTROM') /= 0) Rtrnc = Rtrnc/Angstrom

      case ('SUPS')
        !                                                              *
        !***** SUPS ****************************************************
        !                                                              *
        ! Introduce supersymmetry
        ! Input format
        ! nsg                (number of super groups)
        ! Repeat nsg times
        ! nmem, (ind.., i = 1, nmem)

        Chr = Get_Ln(LuRd)
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

        Chr = Get_Ln(LuRd)
        call Get_F1(1,blavAI)

      case ('THER')
        !                                                              *
        !***** THER ****************************************************
        !                                                              *
        lNmHss = .true.
        lTherm = .true.
        Chr = Get_Ln(LuRd)
        call Get_I1(1,nsRot)
        Chr = Get_Ln(LuRd)
        call Get_F1(1,UserP)
        do
          Chr = Get_Ln(LuRd)
          call UpCase(Chr)
          if (Chr(1:4) == 'END ') then
            if (nUserPT == 0) then
              nUserPT = 1
              UserT(1) = 298.15_wp
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

        Chr = Get_Ln(LuRd)
        call Get_F1(1,ThrEne)
        call Get_F1(2,ThrGrd)
        ThrInp = .true.

      case ('TOLE')
        !                                                              *
        !***** TOLE ****************************************************
        !                                                              *
        ! read the constraints threshold

        Chr = Get_Ln(LuRd)
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
        iOptC = ibclr(iOptC,7)

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
        iOptC = ibset(iOptC,10)
        iOptC = ibset(iOptC,11)

      case ('VDWC')
        !                                                              *
        !***** VDWB VdW correction for coordinate only *****************
        !                                                              *
        iOptC = ibset(iOptC,11)

      case ('VDWH')
        !                                                              *
        !***** VDWB VdW correction for Hessian only ********************
        !                                                              *
        iOptC = ibset(iOptC,10)

      case ('WIND')
        !                                                              *
        !***** WIND ****************************************************
        !                                                              *
        do
          Chr = Get_Ln(LuRd)
          call UpCase(Chr)
          if ((Chr /= '') .and. (Chr(1:1) /= '*')) exit
        end do
        call Get_I1(1,nWndw)

      case default
        call WarningMessage(2,'Error in RdCtl_Slapaf')
        if (Chr(1:1) == ' ') then
          write(Lu,*) ' RdCtl_Slapaf: Command line starts with a blank.'
        else
          write(Lu,*)
          write(Lu,*) ' *********** ERROR ***********'
          write(Lu,*) ' The program has been supplied'
          write(Lu,*) ' with an unknown command.     '
          write(Lu,*) ' *****************************'
        end if
        write(Lu,'(A)') Chr
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
  if (rMEP) Valu = max(real(iMEP+1,kind=wp),One)*Valu
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
    MF(:,1:nsAtom) = reshape(TmpRx,[3,nsAtom])
  else if (iMEP == 0) then
    call NameRun('RUNOLD')
    call qpg_dArray('Reaction Vector',Found,nRx)
    !write(u6,*) 'RUNOLD: Found=',Found
    if (Found) then
      ! Case 2)
      call Get_dArray('Reaction Vector',MF,3*nsAtom)
      call NameRun('#Pop')
    else
      call NameRun('#Pop')
      call qpg_dArray('Reaction Vector',Found,nRx)
      !write(u6,*) 'RUNFILE: Found=',Found
      if (Found) then
        ! Case 3)
        call Get_dArray('Reaction Vector',MF,3*nsAtom)
      else
        call WarningMessage(2,'Error in RdCtl_Slapaf')
        write(u6,*)
        write(u6,*) '************ ERROR **************'
        write(u6,*) 'IRC calculation but no IRC vector'
        write(u6,*) '*********************************'
        call Quit_OnUserError()
      end if
    end if
  end if

  ! Fix the direction forward/backwards

  if ((iMEP == 0) .and. (iRC == -1)) MF(:,1:nsAtom) = -MF(:,1:nsAtom)
  if ((iMEP == 0) .and. (MEP_Type == 'TRANSVERSE')) call Put_dArray('Transverse',MF,3*nsAtom)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (FindTS .and. (.not. TSConstraints)) call SysWarnMsg('RdCtl_Slapaf','WARNING:','FindTS specified, but no TSConstraints. '// &
                                                        'It is highly recommended to use TSConstraints in SLAPAF instead of '// &
                                                        '(or in addition to) global constraints when using FindTS. '// &
                                                        'TSConstraints will be lifted in the final TS search.')
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
        Dir(:,iAtom) = Gx(:,iAtom,iter)/xWeight
      end do
      call Put_dArray('Transverse',Dir,3*nsAtom)
      call mma_deallocate(Dir)
    end if
  end if
end if

if (Explicit_IRC) call mma_deallocate(TmpRx)

! Activate MPS update of Hessian if FindTS

if (FindTS) then

  if (btest(iOptH,5)) then
    iOptH = ibset(0,5) ! EU
  else if (btest(iOptH,6)) then
    iOptH = ibset(0,6) ! TS-BFGS
  else
    iOptH = ibset(0,4) ! MSP
  end if
  iOptC = ibset(iOptC,12)

  ! Increase the update window so that we will not lose the update
  ! which generated the negative curvature.

  nWndw = 4*nWndw
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify some options if TS search

if (.not. btest(iOptC,7)) then
  if (.not. btest(iOptH,3)) iOptH = ibset(0,4) ! MSP
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
  if (Update == Two) then

    ! Enable FindTS procedure

    !write(u6,*) 'Enable FindTS procedure'
    if (.not. btest(iOptH,3)) iOptH = ibset(0,4) ! MSP
    nWndw = 4*nWndw
    ! make it look as if this were FindTS with constraints
    FindTS = .true.
    TSConstraints = .true.
    iOptC = ibset(iOptC,12)
    iOptC = ibset(iOptC,13)
    Beta = 0.1_wp

  else

    ! Normal constrained optimization with a reduced threshold.
    ! Let the threshold be somewhat tighter as we are close to the TS.

    if ((HSR/HSR0 < 0.2_wp) .or. (HSR < 0.2_wp)) then
      !ThrGrd = 0.0003_wp
      Beta = 0.1_wp
    else
      !ThrGrd = 0.003_wp
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
  iOptC = ibset(iOptC,8) ! Constraints
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
  if (FindTS .and. (.not. (Found .or. Manual_Beta))) Beta_Disp = 0.1_wp
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out input parameters, No output if we didn't find the proper
! gradients on the runfile. We will come back!!!

if ((SuperName == 'slapaf') .and. (.not. Request_Alaska)) call WrInp_sl()
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
