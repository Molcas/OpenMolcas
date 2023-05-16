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

subroutine Numerical_Gradient(ireturn)

#ifndef _HAVE_EXTRA_
use Prgm, only: PrgmFree
#endif
use Para_Info, only: MyRank, nProcs, Set_Do_Parallel
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
use Para_Info, only: King
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, OneHalf, Angstrom, auTokcalmol
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "LenIn.fh"
#include "standard_iounits.fh"
#include "warnings.h"
real(kind=wp) :: Energy_Ref, FX(3), rDum(1), Dsp, EMinus, EPlus, Grada, Gradb, rDeg, rDelta, rMax, rTest, Sgn, TempX, TempY, &
                 TempZ, x, x0, y, y0, z, z0
integer(kind=iwp) :: iOper(0:7), jStab(0:7), iCoSet(0:7,0:7), iDispXYZ(3), rc, error, i, iAt, iAtom, ibla, iBlabla, iChxyz, iCoor, &
                     id_Tsk, iData, iDisp, iEm, iEp, ii, iMlt, iPL, iPL_Save, IPotFl, iQMChg, iR, iRank, iRC, irlxroot1, &
                     irlxroot2, iRoot, iSave, ITkQMMM, ixyz, j, LuWr_save, MaxDCR, mDisp, mInt, MltOrd, nAll, nAtMM, nAtoms, &
                     nCList, nDisp, nDisp2, nGNew, nGrad, nIrrep, nLambda, nMult, nRoots, nStab, nSym
character(len=8) :: Method
character(len=LenIn) :: Namei
character(len=10) :: ESPFKey
character(len=180) :: Line
logical(kind=iwp) :: DispX, DispY, DispZ, Do_ESPF, DoDirect, DoFirst, DoTinker, DynExtPot, Exists, External_Coor_List, Found, &
                     Is_Roots_Set, KeepOld, NMCart, StandAlone
integer(kind=iwp), allocatable :: IsMM(:)
real(kind=wp), allocatable :: EnergyArray(:,:), GradArray(:,:), OldGrads(:,:), Grad(:), GNew(:), MMGrd(:,:), BMtrx(:,:), &
                              TMtrx(:,:), Coor(:,:), Energies_Ref(:), XYZ(:,:), AllC(:,:), Disp(:), Deg(:,:), Mltp(:), C(:,:), &
                              Tmp2(:), Tmp(:,:)
character(len=LenIn), allocatable :: AtomLbl(:)
real(kind=wp), parameter :: ToHartree = One/auTokcalmol
integer(kind=iwp), external :: Read_Grad, IsFreeUnit, iPrintLevel, iChAtm, iDeg
character(len=180), external :: Get_Ln
logical(kind=iwp), external :: Rsv_Tsk, Reduce_Prt
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
character(len=80) :: SSTMNGR
integer(kind=iwp) :: SSTMODE
logical(kind=iwp), external :: Rsv_Tsk_Even
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the print level.

iPL_Save = iPrintLevel(-1)
iPL = iPL_Save

if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
ireturn = _RC_ALL_IS_WELL_

! Get information regarding the last method used

call Get_cArray('Relax Method',Method,8)
call DecideOnESPF(Do_ESPF)
Is_Roots_Set = .false.
call Qpg_iScalar('Number of roots',Is_Roots_Set)
if (Is_Roots_Set) then
  call Get_iScalar('Number of roots',nRoots)
  if (nRoots == 1) then
    iRoot = 1
  else
    call qpg_iScalar('NumGradRoot',Found)
    if (Found) then
      call Get_iScalar('NumGradRoot',iRoot)
    else
      iRoot = 1
    end if
  end if
else
  nRoots = 1
  iRoot = 1
end if
call mma_allocate(Energies_Ref,nRoots,Label='Energies_Ref')
!write(u6,*) 'Is_Roots_Set, nRoots, iRoot = ',Is_Roots_Set,nRoots,iRoot
if (iRoot > nRoots) then
  write(LuWr,*)
  write(LuWr,*) '****************** ERROR ******************'
  write(LuWr,*) 'It was selected to run numerical gradient'
  write(LuWr,*) 'using energies for root/state number ',iRoot
  write(LuWr,*) 'but only ',nRoots,' exists.'
  write(LuWr,*) '*******************************************'
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nRoots > 1) then
  call Get_dArray('Last energies',Energies_Ref,nRoots)
else
  call Get_dScalar('Last energy',Energies_Ref(iRoot))
end if
call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Coor,3,nAtoms,Label='Coor')
call mma_allocate(AtomLbl,nAtoms,Label='AtomLbl')
call Get_cArray('Unique Atom Names',AtomLbl,LenIn*nAtoms)
call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
if (iPL_Save >= 3) call RecPrt('Original coordinates',' ',Coor,3,nAtoms)
!                                                                      *
!***********************************************************************
!                                                                      *
! If this is a QM/MM calculation, the gradient on the MM atoms
! is analytical (using the ESPF method for the QM/MM electrostatics)
! Accordingly, no need to loop over the MM atoms

DoTinker = .false.
DoDirect = .false.
call F_Inquire('ESPF.DATA',Exists)
if (Exists) then
  IPotFl = IsFreeUnit(15)
  call Molcas_Open(IPotFl,'ESPF.DATA')
  Line = ' '
  do while (index(Line,'ENDOFESPF ') == 0)
    Line = Get_Ln(IPotFl)
    if (index(Line,'TINKER ') /= 0) then
      DoTinker = .true.
    else if (index(Line,'DIRECT ') /= 0) then
      DoDirect = .true.
    else if (index(Line,'MLTORD ') /= 0) then
      call Get_I1(2,MltOrd)
      ibla = 0
      do ii=0,MltOrd
        ibla = ibla+(ii+2)*(ii+1)/2
      end do
      MltOrd = ibla
    else if (index(Line,'MULTIPOLE ') /= 0) then
      call Get_I1(2,nMult)
      call mma_allocate(Mltp,nMult,Label='Mltp')
      do iMlt=1,nMult,MltOrd
        Line = Get_Ln(IPotFl)
        call Get_I1(1,iAt)
        call Get_F(2,Mltp(iMlt),MltOrd)
      end do
    end if
  end do
  close(IPotFl)
end if

nAtMM = 0
if (DoTinker) then
  call mma_allocate(IsMM,nAtoms,Label='IsMM')
  call MMCount(nAtoms,nAtMM,IsMM)
  if (nAtMM > 0) then
    iQMChg = 0
    StandAlone = .false.
    DoFirst = .not. allocated(Mltp)
    call RunTinker(nAtoms,Coor,Mltp,DoFirst,IsMM,MltOrd,DynExtPot,iQMchg,iBlabla,StandAlone,DoDirect)
    call mma_allocate(MMGrd,3,nAtoms,Label='MMGrd')
    MMGrd(:,:) = Zero
    ITkQMMM = IsFreeUnit(15)
    call Molcas_Open(ITkQMMM,'QMMM')
    Line = ' '
    do while (index(Line,'TheEnd') == 0)
      Line = Get_Ln(ITkQMMM)
      if (index(Line,'MMGradient') /= 0) then
        call Get_I1(2,iAtom)
        call Get_F(3,FX,3)
        if (IsMM(iAtom) == 1) MMGrd(:,iAtom) = FX(:)
      end if
    end do
    call DScal_(3*nAtoms,Angstrom*ToHartree,MMGrd,1)
    close(ITkQMMM)
    if (iPL_Save >= 3) call RecPrt('MM Grad:',' ',MMGrd,3,nAtoms)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up rDelta from the runfile
call Get_dScalar('Numerical Gradient rDelta',rDelta)
!                                                                      *
!***********************************************************************
!                                                                      *
! Check if there is a coordinate list from Slapaf on the run file.
! Read the B-matrix and T-matrix. To be used in case of
! differentiation in internal coordinates according to Slapaf.

if (nAtMM == 0) call GenCxCTL(iRC,NMCart,rDelta)

call qpg_dArray('CList',Found,nCList)
if (Found) then
  External_Coor_List = .true.
  call Get_iScalar('No of Internal Coordinates',mInt)
  call Get_iScalar('nLambda',nLambda)
  call mma_allocate(BMtrx,3*nAtoms,mInt,Label='BMtrx')
  call mma_allocate(TMtrx,mInt,mInt,Label='TMtrx')
  call Get_dArray('BMtrx',BMtrx,size(BMtrx))
  call Get_dArray('T-Matrix',TMtrx,mInt**2)
else
  NMCart = .false.
  External_Coor_List = .false.
  mInt = 3*nAtoms
  nLambda = 0
  call mma_allocate(BMtrx,3*nAtoms,3*nAtoms,Label='BMtrx')
  call mma_allocate(TMtrx,mInt,mInt,Label='TMtrx')
  call unitmat(BMtrx,3*nAtoms)
  call unitmat(TMtrx,mInt)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((Method(5:7) == 'SCF') .or. &
    (Method(1:6) == 'KS-DFT') .or. &
    (Method(1:6) == 'CASSCF') .or. &
    (Method(1:6) == 'RASSCF') .or. &
    (Method(1:6) == 'GASSCF') .or. &
    (Method(1:6) == 'CASPT2') .or. &
    (Method(1:5) == 'MBPT2') .or. &
    (Method(1:5) == 'CCSDT') .or. &
    (Method(1:4) == 'CHCC') .or. &
    (Method(1:6) == 'MCPDFT') .or. &
    (Method(1:6) == 'MSPDFT') .or. &
    (Method(1:4) == 'CHT3') .or. &
    (Method(1:8) == 'EXTERNAL')) then
  if (iPL_Save >= 3) then
    write(LuWr,*)
    write(LuWr,'(A,A,A)') ' Numerical_Gradient: Original ',Method,' Energies:'
    write(LuWr,'(G21.14)') Energies_Ref(:)
    write(LuWr,*)
  end if
else
  write(LuWr,'(A,A,A)') 'Numerical gradient for ',Method,' is not implemented yet.'
  call Abend()
end if

nDisp2 = 2*3*nAtoms
call mma_Allocate(EnergyArray,nRoots,nDisp2)
call FZero(EnergyArray,nRoots*nDisp2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up symmetry information

call Get_iScalar('nSym',nSym)
nIrrep = nSym
call Get_iArray('Symmetry operations',iOper,nSym)
MaxDCR = nIrrep
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up displacement vector

call Get_nAtoms_All(nAll)
call mma_allocate(AllC,3,nAll,Label='AllC')
call Get_Coord_All(AllC,nAll)
!                                                                      *
!***********************************************************************
!                                                                      *
if (External_Coor_List) then

  ! Externally define displacement list

  nDisp = mInt
  call mma_allocate(Disp,mInt,Label='Disp')
  call Get_dArray('DList',Disp,mInt)
  !call RecPrt('Dlist',' ',Disp,1,mInt)

else

  ! Cartesian displacement list

  nDisp = 3*nAtoms
  call mma_allocate(Disp,nDisp,Label='Disp')
  call FZero(Disp,nDisp)

  do i=1,nAtoms

    ! Find the stabilizer of this center

    iChxyz = iChAtm(Coor(:,i))
    call Stblz(iChxyz,nStab,jStab,MaxDCR,iCoSet)

    iDispXYZ(:) = 0
    do j=0,nStab-1
      if (btest(jStab(j),0)) then
        iDispXYZ(1) = iDispXYZ(1)-1
      else
        iDispXYZ(1) = iDispXYZ(1)+1
      end if
      if (btest(jStab(j),1)) then
        iDispXYZ(2) = iDispXYZ(2)-1
      else
        iDispXYZ(2) = iDispXYZ(2)+1
      end if
      if (btest(jStab(j),2)) then
        iDispXYZ(3) = iDispXYZ(3)-1
      else
        iDispXYZ(3) = iDispXYZ(3)+1
      end if
    end do

    ! If this is a MM atom, do not make displacements

    if (DoTinker .and. (IsMM(i) == 1)) iDispXYZ(:) = 0
    DispX = iDispXYZ(1) /= 0
    DispY = iDispXYZ(2) /= 0
    DispZ = iDispXYZ(3) /= 0
    x0 = Coor(1,i)
    y0 = Coor(2,i)
    z0 = Coor(3,i)

    ! Find the shortest distance to another atom!

    rMax = huge(rMax)
    do j=1,nAll
      x = AllC(1,j)
      y = AllC(2,j)
      z = AllC(3,j)
      rTest = sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
      if (rTest == Zero) rTest = huge(rTest)
      rMax = min(rMax,rTest)
    end do
    if (DispX) Disp((i-1)*3+1) = rDelta*rMax
    if (DispY) Disp((i-1)*3+2) = rDelta*rMax
    if (DispZ) Disp((i-1)*3+3) = rDelta*rMax
  end do
  !call RecPrt('Disp',' ',Disp,1,nDisp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Deg,3,nAtoms,Label='Deg')
Deg(:,:) = Zero
do i=1,nAtoms
  rDeg = real(iDeg(Coor(:,i)),kind=wp)
  Deg(1,i) = rDeg
  Deg(2,i) = rDeg
  Deg(3,i) = rDeg
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(XYZ,3*nAtoms,2*nDisp,Label='XYZ')
XYZ(:,:) = Zero
if (External_Coor_List) then
  mDisp = nDisp*2
  call Get_dArray('CList',XYZ,3*nAtoms*mDisp)
else
  mDisp = 0
  do iDisp=1,nDisp2
    iCoor = (iDisp+1)/2
    if (Disp(icoor) /= Zero) then
      mDisp = mDisp+1

      ! Modify the geometry

      call dcopy_(3*nAtoms,Coor,1,XYZ(:,mDisp),1)
      Sgn = One
      if (mod(iDisp,2) == 0) Sgn = -One
      XYZ(iCoor,mDisp) = XYZ(iCoor,mDisp)+Sgn*Disp(icoor)
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Save the "new geometry" field from the RunFile, if any

call qpg_dArray('GeoNew',Found,nGNew)
if (.not. Found) nGNew = 0
if (nGNew > 0) then
  call mma_allocate(GNew,nGNew)
  call Get_dArray('GeoNew',GNew,nGNew)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Save global print level

iPL_Save = iPrintLevel(-1)
!iPL_Base = 0
!If (iPL_Save >= 3) iPl_Base = iPL_Save

#ifdef _DEBUGPRINT_
call RecPrt('BMtrx',' ',BMtrx,3*nAtoms,mInt)
call RecPrt('TMtrx',' ',TMtrx,mInt,mInt)
call RecPrt('Degeneracy vector',' ',Deg,3,nAtoms)
call RecPrt('Coordinate List',' ',XYZ,3*nAtoms,mDisp)
#endif

!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Loop over displacements                                              *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (iPL >= 2) then
  write(LuWr,*) 'Root to use: ',iRoot
  write(LuWr,*) 'Number of internal degrees            ',nDisp
  write(LuWr,*) 'Number of constraints                 ',nLambda
  write(LuWr,*) 'Number of displacements               ',nDisp*2
  if ((nAtMM /= 0) .and. DoTinker .and. (.not. DoDirect)) then
    write(LuWr,*) 'Number of MM degrees (analytical grad)',nAtMM*3
  end if
  write(LuWr,*) 'Effective number of displacements     ',2*(nDisp-nLambda-3*nAtMM)
  write(LuWr,*) 'Relative displacements                ',rDelta
  write(LuWr,*)
end if

! This printout is disabled since it means different things for
! externally defined coordinates or not.

if (iPL >= 3) then
# ifdef _HIDE_
  write(LuWr,'(1x,A)') '---------------------------------------------'
  write(LuWr,'(1x,A)') '               X           Y           Z     '
  write(LuWr,'(1x,A)') '---------------------------------------------'
  do iAtom=1,nAtoms
    TempX = Disp(3*(iAtom-1)+1)
    TempY = Disp(3*(iAtom-1)+2)
    TempZ = Disp(3*(iAtom-1)+3)
    Namei = AtomLbl(iAtom)
    write(LuWr,'(2X,A,3X,3F12.6)') Namei,TempX,TempY,TempZ
  end do
# endif
  write(LuWr,'(1x,A)') '---------------------------------------------'
  write(LuWr,*)
  write(LuWr,*) 'Here we go ...'
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Change output unit

LuWr_save = LuWr
if (MyRank /= 0) then
  LuWr = 55
  LuWr = isFreeUnit(LuWr)
  call molcas_open(luwr,'Temp_OutPut')
end if

! FM 16/4/2013
! ESPF charges are set to zero so that no microiterations will be
! performed during numerical gradient
! IFG: swap files, as now NG uses a subdirectory

iSave = 15
if (Do_ESPF) then
  iSave = IsFreeUnit(iSave)
  call Molcas_Open(iSave,'ESPF.SAV')
  iData = IsFreeUnit(iSave)
  call Molcas_Open(iData,'ESPF.DATA')
  do
    Line = Get_Ln(iData)
    ESPFKey = Line(1:10)
    if (ESPFKey == 'ENDOFESPF ') then
      write(iSave,'(A132)') Line
      exit
    else
      if (ESPFKey /= 'MULTIPOLE ') then
        write(iSave,'(A132)') Line
      else
        call Get_I1(2,nMult)
        do iMlt=1,nMult
          Line = Get_Ln(iData)
        end do
      end if
    end if
  end do
  close(iSave)
  close(iData)
end if
! FM End
!                                                                      *
!***********************************************************************
!                                                                      *
! Parallel loop over the nDisp displacements.
! Reserve task on global task list and get task range in return.
! Function will be false if no more tasks to execute.
#if !defined(_GA_) && defined(_MOLCAS_MPP_)
call getenvf('MOLCAS_SSTMNGR',SSTMNGR)
if (SSTMNGR(1:1) == 'Y') then
  SSTMODE = 1
else
  SSTMODE = 0
end if
if (SSTMODE == 1) then
  call Init_Tsk(id_Tsk,mDisp-2*nLambda)
else
  call Init_Tsk_Even(id_Tsk,mDisp-2*nLambda)
end if
#else
call Init_Tsk(id_Tsk,mDisp-2*nLambda)
#endif
call mma_allocate(C,3,nAtoms,Label='C')
do
# if defined (_MOLCAS_MPP_) && !defined(_GA_)
  if (SSTMODE == 1) then
    if (.not. Rsv_Tsk(id_Tsk,iDisp)) exit
  else
    if (.not. Rsv_Tsk_Even(id_Tsk,iDisp)) exit
  end if
# else
  if (.not. Rsv_Tsk(id_Tsk,iDisp)) exit
# endif

  ! Offset for the constraints

  iDisp = iDisp+2*nLambda

  ! Get the displaced geometry

  call dcopy_(3*nAtoms,xyz(:,iDisp),1,C,1)
  !call RecPrt('C',' ',C,3*nAtoms,1)
  call Put_Coord_New(C,nAtoms)

  ! Compute integrals

  !jPL = iPrintLevel(Max(iPL_Base,0)) ! Silent
  call Set_Do_Parallel(.false.)

  ! Switch to a new directory
  ! WARNING WARNING WARNING WARNING WARNING
  ! ugly hack, do not try this at home
  call SubWorkDir()
  ! WARNING WARNING WARNING WARNING WARNING

  call StartLight('seward')
  call init_run_use()
  call init_ppu(.true.)
  call Disable_Spool()
  call Seward(ireturn)
  if (iReturn /= 0) then
    write(LuWr,*) 'Numerical_Gradient failed ...'
    write(LuWr,*) 'Seward returned with return code, rc = ',iReturn
    write(LuWr,*) 'for the perturbation iDisp = ',iDisp
    call Abend()
  end if

  ! Compute the ESPF stuff

  if (Do_ESPF) then
    call StartLight('espf')
    call init_run_use()
    call Disable_Spool()
    StandAlone = .true.
    call ESPF(ireturn,StandAlone)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'ESPF returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  ! Compute the wave function

  if ((Method(5:7) == 'SCF') .or. &
      (Method(1:6) == 'KS-DFT') .or. &
      (Method(1:5) == 'MBPT2') .or. &
      (Method(1:4) == 'CHCC') .or. &
      (Method(1:4) == 'CHT3')) then
    call StartLight('scf')
    call init_run_use()
    call Disable_Spool()
    call xml_open('module',' ',' ',0,'scf')
    call SCF(iReturn)
    call xml_close('module')
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'SCF returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  else if ((Method(1:6) == 'RASSCF') .or. &
           (Method(1:6) == 'GASSCF') .or. &
           (Method(1:6) == 'CASSCF') .or. &
           (Method(1:6) == 'MCPDFT') .or. &
           (Method(1:6) == 'MSPDFT') .or. &
           (Method(1:6) == 'CASPT2') .or. &
           (Method(1:5) == 'CCSDT')) then
    call StartLight('rasscf')
    call init_run_use()
    call Disable_Spool()
    call RASSCF(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'RASSCF returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  else if (Method(1:8) == 'EXTERNAL') then
    call StartLight('false')
    call init_run_use()
    call Disable_Spool()
    call False_program(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'FALSE returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if (Method(1:5) == 'MBPT2') then
    call StartLight('mbpt2')
    call init_run_use()
    call Disable_Spool()
    call MP2_Driver(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'MBPT2 returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if (Method(1:5) == 'CCSDT') then
    call StartLight('motra')
    call init_run_use()
    call Disable_Spool()
    call Motra(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'Motra returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if

    call StartLight('ccsdt')
    call init_run_use()
    call Disable_Spool()
    call CCSDT(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'CCSDT returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if ((Method(1:4) == 'CHCC') .or. &
      (Method(1:4) == 'CHT3')) then
    call StartLight('chcc')
    call init_run_use()
    call Disable_Spool()
    call CHCC(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'CHCC returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if (Method(1:4) == 'CHT3') then
    call StartLight('cht3')
    call init_run_use()
    call Disable_Spool()
    call CHT3(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'CHT3 returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if (Method(1:6) == 'CASPT2') then
    call StartLight('caspt2')
    call init_run_use()
    call Disable_Spool()
    call CASPT2(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'CASPT2 returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  if ((Method(1:6) == 'MCPDFT') .or. (Method(1:6) == 'MSPDFT')) then
    call StartLight('mcpdft')
    call init_run_use()
    call Disable_Spool()
    call MCPDFT(ireturn)
    if (iReturn /= 0) then
      write(LuWr,*) 'Numerical_Gradient failed ...'
      write(LuWr,*) 'MCPDFT returned with return code, rc = ',iReturn
      write(LuWr,*) 'for the perturbation iDisp = ',iDisp
      call Abend()
    end if
  end if

  !call Get_Energy(EnergyArray(iRoot,iDisp))
  if (nRoots > 1) then
    call Get_dArray('Last energies',EnergyArray(1,iDisp),nRoots)
  else
    call Get_dScalar('Last energy',EnergyArray(iRoot,iDisp))
  end if

  ! Restore directory and prgm database

  call ParentWorkDir()
  call prgmfree()
  call prgminit('numerical_gradient')

  call Set_Do_Parallel(.true.)

  if (iPL >= 2) then
    write(LuWr,200) '   * Point #',iDisp-2*nLambda,' of ',2*(nDisp-nLambda-3*nAtMM),' done.'
    write(LuWr,200) '    (Perturbation ',iDisp,')'
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# if defined (_MOLCAS_MPP_) && !defined(_GA_)
  if (King() .and. (nProcs > 1) .and. (SSTMODE == 1)) exit
# endif
end do
!_MPP End Do
#if !defined(_GA_) && defined(_MOLCAS_MPP_)
if (SSTMODE == 1) then
  call Free_Tsk(id_Tsk)
else
  call Free_Tsk_Even(id_Tsk)
end if
#else
call Free_Tsk(id_Tsk)
#endif
call GADSum(EnergyArray,nRoots*mDisp)
call mma_deallocate(C)
call mma_deallocate(XYZ)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! End of Loop                                                          *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Flush the output from the other nodes.
!
do iRank=1,nProcs-1
  call GASync()
  if (iRank == MyRank) then
    rewind(LuWr)
    do
      read(LuWr,'(A)',iostat=error) Line
      if (error /= 0) exit
      write(LuWr_Save,*) Line
    end do
  end if
  call GASync()
end do
if (MyRank /= 0) then
  close(LuWr)
  LuWr = LuWr_save
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore the "new geometry" field to the RunFile, if any

if (nGNew == 0) then
  call Put_Coord_New([Zero],0)
else
  call Put_Coord_New(GNew,nGNew/3)
  call mma_deallocate(GNew)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((nProcs >= 2) .and. (iPL >= 2)) then
  write(LuWr,*)
  write(LuWr,*) ' Points were printed only by master node'
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read old gradient(s) and convert to internal coordinates

call mma_Allocate(OldGrads,3*nAtoms,nRoots)
call FZero(OldGrads,3*nAtoms*nRoots)
call Get_lScalar('Keep old gradient',KeepOld)
if (KeepOld) then
  call Query_Grads(Found,i,j)
  if (Found .and. (i >= nRoots) .and. (j == 3*nAtoms)) then
    ! If there is a valid GRADS file, read each gradient
    do iR=1,nRoots
      rc = Read_Grad(OldGrads(1,iR),3*nAtoms,iR,0,0)
      if (rc <= 0) then
        write(LuWr,*)
        write(LuWr,'(2X,A,I4,A)') 'No gradient found for root ',iR,', using 0'
        write(LuWr,*)
      end if
    end do
  else
    ! If the GRADS file does not exist or has the wrong sizes,
    ! read a single gradient from the runfile
    call qpg_dArray('GRAD',Found,nGrad)
    if (Found) then
      call Get_dArray('GRAD',OldGrads(1,1),nGrad)
      do iR=2,nRoots
        call dCopy_(nGrad,OldGrads(1,1),1,OldGrads(1,iR),1)
      end do
      if (nRoots > 1) then
        write(LuWr,*)
        write(LuWr,'(2X,A)') 'Using stored gradient for all roots'
        write(LuWr,*)
      end if
    else
      write(LuWr,*)
      write(LuWr,'(2X,A)') 'No gradient found, using 0'
      write(LuWr,*)
    end if
  end if
end if
call mma_allocate(Tmp2,nDisp,Label='Tmp2')
do iR=1,nRoots
  call Eq_Solver('N',3*nAtoms,nDisp,1,BMtrx,.true.,rDum(1),OldGrads(1,iR),Tmp2)
  call FZero(OldGrads(1,iR),3*nAtoms)
  call Eq_Solver('N',nDisp,nDisp,1,TMtrx,.true.,rDum(1),Tmp2,OldGrads(1,iR))
end do
call mma_deallocate(Tmp2)
!                                                                      *
!***********************************************************************
!                                                                      *
!jPL = iPrintLevel(iPL_Save)
if (iPL_Save >= 3) call RecPrt('Energies','(8G16.10)',EnergyArray,nRoots,mDisp)
call mma_Allocate(GradArray,nDisp,nRoots)
call FZero(GradArray,nDisp*nRoots)
call mma_Allocate(Grad,nRoots)

iDisp = nLambda
do i=1,nDisp

  if (Disp(i) == Zero) then
    if (iPL_Save >= 3) then
      if (KeepOld) then
        write(u6,*) 'gradient set to old value'
      else
        write(u6,*) 'gradient set to zero'
      end if
    end if
    do iR=1,nRoots
      Grad(iR) = OldGrads(i,iR)
    end do
  else
    iDisp = iDisp+1
    do iR=1,nRoots
      iEp = iDisp*2-1
      EPlus = EnergyArray(iR,iEp)
      iEm = iDisp*2
      EMinus = EnergyArray(iR,iEm)
      Dsp = Disp(i)
      Grad(iR) = (EPlus-EMinus)/(Two*Dsp)

      ! If the gradient is not close to zero check that it is
      ! consistent. For CASPT2/CASSCF sometimes the active space
      ! breaks down and the computed gradient is rubbish. If we
      ! are lucky this happens only in one of the displacement
      ! directions. In that case compute the gradient with the
      ! one-point equation. The one with the lowest gradient is
      ! the one which is most likely to be correct.
      if (abs(Grad(iR)) > 0.1_wp) then
        Energy_Ref = Energies_Ref(iR)
        Grada = (EPlus-Energy_Ref)/Dsp
        Gradb = (Energy_Ref-EMinus)/Dsp
        if ((abs(Grad(iR)/Grada) > OneHalf) .or. (abs(Grad(iR)/Gradb) > OneHalf)) then
          if (abs(Grada) <= abs(Gradb)) then
            Grad(iR) = Grada
          else
            Grad(iR) = Gradb
          end if
        end if
      end if
    end do
  end if

  call dCopy_(nRoots,Grad,1,GradArray(i,1),nDisp)

end do
!call RecPrt('Grads (old)',' ',OldGrads,3*nAtoms,nRoots)
!call RecPrt('BMtrx',' ',BMtrx,3*nAtoms,nDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute the gradient in Cartesian coordinates
!
call mma_allocate(Tmp2,nDisp,Label='Tmp2')
call mma_allocate(Tmp,3,nAtoms,Label='Tmp')
do iR=1,nRoots
  Tmp2(:) = Zero
  Tmp(:,:) = Zero

  ! Transform the gradient in PCO basis to the internal coordinate
  ! format.

  call DGEMM_('N','N',nDisp,1,nDisp,One,TMtrx,nDisp,GradArray(1,iR),nDisp,Zero,Tmp2,nDisp)

  ! Transform internal coordinates to Cartesian.

  call DGEMM_('N','N',3*nAtoms,1,nDisp,One,BMtrx,3*nAtoms,Tmp2,nDisp,Zero,Tmp,3*nAtoms)
  !call RecPrt('Tmp',' ',Tmp,3,nAtoms)

  ! Modify with degeneracy factors.

  if (.not. NMCart) then
    do iAtom=1,nAtoms
      do ixyz=1,3
        Tmp(ixyz,iAtom) = Tmp(ixyz,iAtom)/Deg(ixyz,iAtom)
      end do
    end do
  end if

  ! Add the MM contribution for MM atoms

  if (nAtMM /= 0) then
    call daxpy_(3*nAtoms,One,MMGrd,1,Tmp,1)
  end if

  ! Apply Morokuma's scheme if needed

  call F_Inquire('QMMM',Exists)
  if (Exists .and. DoTinker) then
    call LA_Morok(nAtoms,Tmp,1)
  end if

  if (iR == iRoot) call Put_dArray('GRAD',Tmp,3*nAtoms)
  call Add_Info('Grad',Tmp,3*nAtoms,6)
  call Store_grad(Tmp,3*nAtoms,iR,0,0)

  if (iPL >= 2) then
    write(LuWr,*)
    write(LuWr,'(1x,A,I5)') 'Numerical gradient, root ',iR
    write(LuWr,'(2x,A)') '---------------------------------------------'
    write(LuWr,'(2x,A)') '               X           Y           Z'
    write(LuWr,'(2x,A)') '---------------------------------------------'
    do iAtom=1,nAtoms
      TempX = Tmp(1,iAtom)
      TempY = Tmp(2,iAtom)
      TempZ = Tmp(3,iAtom)
      Namei = AtomLbl(iAtom)
      write(LuWr,'(2X,A,3X,3F12.6)') Namei,TempX,TempY,TempZ
    end do
    write(LuWr,'(2x,A)') '---------------------------------------------'
  end if
end do
call mma_deallocate(Tmp2)
call mma_deallocate(Tmp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate

call mma_deallocate(Deg)
call mma_deallocate(Disp)
call mma_deallocate(AllC)
call mma_Deallocate(TMtrx)
call mma_Deallocate(BMtrx)
call mma_Deallocate(EnergyArray)
call mma_Deallocate(GradArray)
call mma_Deallocate(OldGrads)
call mma_Deallocate(Grad)
call mma_deallocate(Coor)
call mma_deallocate(Energies_Ref)
call mma_deallocate(AtomLbl)
if (allocated(Mltp)) call mma_deallocate(Mltp)
if (DoTinker) call mma_deallocate(IsMM)
if (nAtMM > 0) call mma_deallocate(MMGrd)

! Restore iRlxRoot if changed as set by the RASSCF module.

if ((Method(1:6) == 'CASSCF') .or. &
    (Method(1:6) == 'RASSCF')) then
  call Get_iScalar('Relax CASSCF root',irlxroot1)
  call Get_iScalar('Relax Original root',irlxroot2)
  if (iRlxRoot1 /= iRlxRoot2) then
    call Put_iScalar('Relax CASSCF root',irlxroot2)
    call Put_iScalar('NumGradRoot',irlxroot2)
  end if
end if

return

200 format(A,I4,A,I4,A)

end subroutine Numerical_Gradient
