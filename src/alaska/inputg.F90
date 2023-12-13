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
! Copyright (C) 1991,1992, Roland Lindh                                *
!***********************************************************************

subroutine Inputg(LuSpool)
!***********************************************************************
!                                                                      *
! Object: input module for the gradient code                           *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September '91                                            *
!                                                                      *
!             Modified to complement GetInf, January '92.              *
!***********************************************************************

use Alaska_Info, only: Am, Auto, ForceNAC
use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iChTbl, iOper, lIrrep, lBsFnc
use Gateway_global, only: Onenly, Test
use Gateway_Info, only: CutInt
use RI_glob, only: Timings_default
use Cholesky, only: timings
use OFembed, only: Do_OFemb, KEonly, OFE_first, Xsigma, dFMD, OFE_KSDFT
use pso_stuff, only: No_Nuc
use Disp, only: ChDisp, CutGrd, Dirct, Disp_Fac, HF_Force, IndDsp, IndxEQ, InxDsp, l2DI, lDisp, lEQ, Mult_Disp, nTR, TRSymm
use NAC, only: DoCSF, EDiff, isCSF, isNAC, NACStates
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: i, iCar, iCnt, iCnttp, iCo, iComp, iElem, iGroup, iIrrep, ijSym, iPL, iPrint, iRout, istatus, iSym(3), iTR, &
                     j, jIrrep, jOper, jPrint, jRout, jTR, k, kTR, ldsp, lTR, LuWr, mc, mdc, mDisp, n, nCnttp_Valence, nDisp, &
                     nElem, nGroup, nRoots, nSlct
real(kind=wp) :: alpha, Fact, ovlp
logical(kind=iwp) :: TstFnc, ltype, Slct, T_Only, No_Input_OK, Skip
character(len=80) :: KWord, Key
integer(kind=iwp), allocatable :: IndCar(:), iTemp(:)
real(kind=wp), allocatable :: Tmp(:), C(:,:), Scr(:,:), Temp(:,:)
character, parameter :: xyz(0:2) = ['x','y','z']
integer(kind=iwp), external :: iPrintLevel, iPrmt, NrOpr
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: Reduce_Prt

iRout = 99
iPrint = nPrint(iRout)
do i=1,nRout
  nPrint(i) = 5
end do
if (ForceNAC) isNAC = .true.
DoCSF = .true.
isCSF = .false.
Auto = .false.
Test = .false.
T_Only = .false.
TRSymm = .false.
lEq = .false.
Slct = .false.
l2DI = .true.
HF_Force = .false.
NO_NUC = .false.
Timings_default = Timings
Xsigma = 1.0e4_wp
dFMD = Zero
Do_OFemb = .false.
KEonly = .false.
OFE_first = .true.
Show = .true.
LuWr = u6
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = iPL-1
if (iPL == 0) then
  jPrint = 0
else if (iPL == 1) then
  jPrint = 0
else if (iPL == 2) then
  jPrint = 6
else if (iPL == 3) then
  jPrint = 6
else if (iPL == 4) then
  jPrint = 49
else
!else if (iPL == 5) then
  jPrint = 98
end if

do i=1,nRout
  nPrint(i) = jPrint
end do

call mma_allocate(iTemp,3*MxAtom,label='iTemp')

! First CutGrd cannot be more accurate than CutInt!
CutGrd = max(1.0e-7_wp,CutInt)
! Second CutInt should now locally for Alaska be reset to the value
! of CutInt/100!
CutInt = CutGrd*1.0e-2_wp
do i=1,3*MxAtom
  IndxEq(i) = i
end do
do ldsp=1,3*MxAtom
  Dirct(ldsp) = .true.
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! KeyWord directed input

rewind(LuSpool)
No_Input_OK = .true.
call RdNLst_(LuSpool,'ALASKA',No_Input_OK)
KWord = ' &ALASKA'
do
  read(LuSpool,'(A72)',iostat=istatus) Key
  if (istatus < 0) then
    exit
  else if (istatus > 0) then
    call Error()
  end if
  KWord = Key
  call UpCase(KWord)
  if (KWord(1:1) == '*') cycle
  if (KWord == '') cycle
  select case (KWord(1:4))
    case ('VERB')
      !                                                                *
      !***** VERB ******************************************************
      !                                                                *
      ! Verbose mode.

      nPrint(80) = 6
      nPrint(1) = 6
      nPrint(9) = 6
      nPrint(99) = 6
    case ('PRIN')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Print level

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*,iostat=istatus) n
      call Error()
      do i=1,n
        do
          read(LuSpool,'(A)',iostat=istatus) KWord
          call Error()
          if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
        end do
        read(KWord,*,iostat=istatus) jRout,iPrint
        call Error()
        nPrint(jRout) = iPrint
      end do
    case ('EQUI')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Equivalence option

      if (T_Only) then
        call WarningMessage(2,'Error in InputG')
        write(LuWr,*) 'EQUI option does not ork with RF calculations!'
        call Quit_OnUserError()
      end if
      lEq = .true.
      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) nGroup
      do iGroup=1,nGroup
        do
          read(LuSpool,'(A)',iostat=istatus) KWord
          call Error()
          if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
        end do
        read(KWord,*) nElem,(iTemp(iElem),iElem=1,nElem)
        do iElem=2,nElem
          IndxEq(iTemp(iElem)) = iTemp(1)
          Dirct(iTemp(iElem)) = .false.
        end do
      end do
    case ('CUTO')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Cutoff for computing primitive gradients

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*,iostat=istatus) CutGrd
      call Error()
      CutGrd = abs(CutGrd)
    case ('HF-F')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Compute Hellmann-Feynman forces

      HF_Force = .true.
    case ('NOIN')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Disable the utilization of translational and
      ! rotational invariance of the energy in the
      ! computation of the molecular gradient.

      TRSymm = .false.
    case ('SELE')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! selection option

      if (T_Only) then
        call WarningMessage(2,'Error in InputG')
        write(LuWr,*) 'SELE option does not work with RF calculations!'
        call Quit_OnUserError()
      end if
      Slct = .true.
      if (lEq) then
        call WarningMessage(2,'Error in InputG')
        write(LuWr,*) ' The Selection option must preceed the Equivalence option to work together.'
        call Quit_OnUserError()
      end if
      do i=1,3*MxAtom
        Dirct(i) = .false.
      end do
      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) nSlct

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) (iTemp(iElem),iElem=1,nSlct)
      do iElem=1,nSlct
        Dirct(iTemp(iElem)) = .true.
      end do
    case ('2DOP')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Change default for the prescreening.

      l2DI = .false.
    case ('2DIP')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Change default for the prescreening.

      l2DI = .true.
    case ('ONEO')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Do not compute two electron integrals.

      Onenly = .true.
    case ('TEST')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Process only the input.

      Test = .true.
    case ('SHOW')
      !                                                                *
      !*****************************************************************
      !                                                                *
      !-----Raise the printlevel to show gradient contributions

      if (iPL >= 2) then
        nPrint(112) = 15
        nPrint(1) = 15
        nPrint(33) = 15
      end if
    case ('PNEW')
      !                                                                *
      !***** PNEW ******************************************************
      !                                                                *
      ! Print gradient in NEW human-readable format

      nPrint(1) = 4
    case ('POLD')
      !                                                                *
      !***** POLD ******************************************************
      !                                                                *
      ! Print gradient in OLD format

      nPrint(1) = 5
    case ('NONU')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Do not compute the nuclear charge contribution

      NO_NUC = .true.
    case ('EXTR')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Put the program name and the time stamp onto the extract file

      write(LuWr,*) 'InputG: EXTRACT option is redundant and is ignored!'
    case ('CHOI')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Cholesky input section

      call Cho_alaska_rdInp(LuSpool)
    case ('OFEM')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Orbital-Free Embedding (OFE) input section

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      call UpCase(KWord)
      KWord = adjustl(KWord)
      read(KWord,'(A)') OFE_KSDFT
      Do_OFemb = .true.
    case ('KEON')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Mode "Kinetic Energy Only" for OFE input section

      KEonly = .true.
    case ('DFMD')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Mode "Kinetic Energy Only" for OFE input section

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) dFMD,Xsigma
    ! Keyword 'NUMErical' checked earlier - forces numerical gradients
    ! Keyword 'DELTa' selects the scaling factor for the displacements
    !                 in the numerical_gradient module
    ! Keyword 'KEEP' does not remove the old gradient
    ! Keyword 'INVErt' inverts the treatment of constraints
    ! Here it's only included for consistency
    case ('NUME')
    case ('DELT')
    case ('KEEP')
    case ('INVE')
    case ('ROOT')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Root keyword, now also for analytical gradient
      ! This is a dummy, the keyword is already read in chk_numerical

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) i
    case ('NAC ')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! NAC keyword: compute non-adiabatic couplings between 2 states
      ! The keyword is already read in chk_numerical

      do
        read(LuSpool,'(A)',iostat=istatus) KWord
        call Error()
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      read(KWord,*) NACstates(1),NACstates(2)
      isNAC = .true.
    case ('NOCS')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! NOCSF keyword, to neglect the CSF contribution to the NAC,
      ! which is the cause for translational variance

      DoCSF = .false.
    case ('AUTO')
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! AUTO keyword, used by SLAPAF, to signal this is an automated
      ! call to ALASKA

      Auto = .true.
    case ('END ')
      exit
    case default
      istatus = 1
      call Error()
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!                          End of input section.                       *
!                                                                      *
!***********************************************************************
!                                                                      *

! NAC could have been activated through explicit input or through
! a previous MCLR

if (isNAC) then
  No_Nuc = .true.
  ! Get the state energies
  call Get_iScalar('Number of roots',nRoots)
  call mma_Allocate(Tmp,nRoots,Label='Tmp')
  call Get_dArray('Last energies',Tmp,nRoots)
  Ediff = Tmp(NACstates(1))-Tmp(NACstates(2))
  call mma_deallocate(Tmp)
end if

nCnttp_Valence = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux) exit
  nCnttp_Valence = nCnttp_Valence+1
end do

if (lEq) TRSymm = .false.
if (Slct) TRSymm = .false.
iPrint = nPrint(iRout)

TRsymm = (TRsymm .or. T_Only) .and. (.not. Test)

! Compute number of centers and displacements. Ignore pseudo centers.
! If any pseudo centers disable use of translational and rotational
! invariance.

mDisp = 0
mdc = 0
do iCnttp=1,nCnttp_Valence
  if (dbsc(iCnttp)%pChrg) then
    TRSymm = .false.
    mdc = mdc+dbsc(iCnttp)%nCntr
  else
    if ((dbsc(iCnttp)%nFragType > 0) .or. dbsc(iCnttp)%Frag) TRSymm = .false.
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      mDisp = mDisp+3*(nIrrep/dc(mdc)%nStab)
    end do
  end if
end do

if (HF_Force .and. Show .and. (iPrint >= 6)) then
  write(LuWr,*)
  write(LuWr,'(A)') '            O B S E R V E ! '
  write(LuWr,'(A)') '            Option for computation of interstate coupling vector or'
  write(LuWr,'(A)') '            Hellmann-Feynman gradient is active.'
  write(LuWr,*)
end if
if (Show .and. (iPrint >= 6)) then
  write(LuWr,*)
  write(LuWr,'(20X,A,ES10.3)') ' Threshold for contributions to the gradient:',CutGrd
  write(LuWr,*)
end if

! Generate symmetry adapted cartesian displacements

if (Show .and. (iPrint >= 6)) then
  write(LuWr,*)
  write(LuWr,'(20X,A)') '********************************************'
  write(LuWr,'(20X,A)') '* Symmetry Adapted Cartesian Displacements *'
  write(LuWr,'(20X,A)') '********************************************'
  write(LuWr,*)
end if

IndDsp(:,:) = 0
InxDsp(:,:) = 0
Disp_Fac(:,:,:) = One
mult_Disp(:) = 1
nDisp = 0
do iIrrep=0,nIrrep-1
  lDisp(iIrrep) = 0
  ltype = .true.
  ! Loop over basis function definitions
  mdc = 0
  mc = 1
  do iCnttp=1,nCnttp_Valence
    ! Loop over unique centers associated with this basis set.
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      IndDsp(mdc,iIrrep) = nDisp
      ! Loop over the cartesian components
      do iCar=0,2
        iComp = 2**iCar
        if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab) .and. (.not. dbsc(iCnttp)%pChrg)) then
          nDisp = nDisp+1
          if (iIrrep == 0) InxDsp(mdc,iCar+1) = nDisp
          lDisp(iIrrep) = lDisp(iIrrep)+1
          mult_Disp(nDisp) = nIrrep/dc(mdc)%nStab
          if (ltype) then
            if (Show .and. (iPrint >= 6)) then
              write(LuWr,*)
              write(LuWr,'(10X,A,A)') ' Irreducible representation : ',lIrrep(iIrrep)
              write(LuWr,'(10X,2A)') ' Basis function(s) of irrep: ',lBsFnc(iIrrep)
              write(LuWr,*)
              write(LuWr,'(A)') ' Basis Label        Type   Center Phase'
            end if
            ltype = .false.
          end if
          if (iIrrep == 0) then
            do jOper=0,nIrrep-1
              Disp_Fac(iCar+1,jOper,mdc) = real(iPrmt(jOper,iComp)*iChTbl(iIrrep,jOper),kind=wp)
            end do
          end if
          if (Show .and. (iPrint >= 6)) then
            write(LuWr,'(I4,3X,A8,5X,A1,7X,8(I3,4X,I2,4X))') nDisp,dc(mdc)%LblCnt,xyz(iCar), &
                                                             (mc+iCo,iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iComp)* &
                                                              iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
                                                              iCo=0,nIrrep/dc(mdc)%nStab-1)
          end if
          write(ChDisp(nDisp),'(A,1X,A1)') dc(mdc)%LblCnt,xyz(iCar)
        end if

      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do
  end do

end do

if (nDisp /= mDisp) then
  call WarningMessage(2,'Error in InputG')
  write(LuWr,*) ' Wrong number of symmetry adapted displacements',nDisp,'=/=',mDisp
  call Abend()
end if

! Set up data for the utilization of the translational
! and rotational invariance of the energy.

if (TRSymm) then
  write(u6,*) 'Unsupported option: TRSymm'
  call Abend()
  iSym(1) = 0
  iSym(2) = 0
  iSym(3) = 0
  do i=1,min(nIrrep-1,5)
    j = i
    if (i == 3) j = 4
    do k=1,3
      if (btest(iOper(j),k-1)) iSym(k) = 2**(k-1)
    end do
  end do
  nTR = 0
  ! Translational equations
  do i=1,3
    if (iSym(i) == 0) nTR = nTR+1
  end do
  if (iPrint >= 99) write(LuWr,*) ' nTR=',nTR
  ! Rotational equations
  if (.not. T_Only) then
    do i=1,3
      j = i+1
      if (j > 3) j = j-3
      k = i+2
      if (k > 3) k = k-3
      ijSym = ieor(iSym(j),iSym(k))
      if (ijSym == 0) nTR = nTR+1
    end do
  end if
  if (nTR == 0) then
    TRSymm = .false.
    Skip = .true.
  else
    if (iPrint >= 99) write(LuWr,*) ' nTR=',nTR
    call mma_allocate(Am,nTR,lDisp(0),Label='Am')
    call mma_allocate(Temp,nTR,nTR,Label='Temp')
    call mma_allocate(C,4,lDisp(0),Label='C')
    call mma_allocate(IndCar,lDisp(0),Label='IndCar')

    Am(:,:) = Zero
    C(:,:) = Zero

    ! Generate temporary information of the symmetrical
    ! displacements.

    ldsp = 0
    mdc = 0
    iIrrep = 0
    do iCnttp=1,nCnttp_Valence
      do iCnt=1,dbsc(iCnttp)%nCntr
        mdc = mdc+1
        !call RecPrt(' Coordinates',' ',dbsc(iCnttp)%Coor(1,iCnt),1,3)
        Fact = Zero
        iComp = 0
        if (dbsc(iCnttp)%Coor(1,iCnt) /= Zero) iComp = ibset(iComp,0)
        if (dbsc(iCnttp)%Coor(2,iCnt) /= Zero) iComp = ibset(iComp,1)
        if (dbsc(iCnttp)%Coor(3,iCnt) /= Zero) iComp = ibset(iComp,2)
        do jIrrep=0,nIrrep-1
          if (TstFnc(dc(mdc)%iCoSet,jIrrep,iComp,dc(mdc)%nStab)) then
            Fact = Fact+One
          end if
        end do
        do iCar=1,3
          iComp = 2**(iCar-1)
          if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)) then
            ldsp = ldsp+1
            Dirct(lDsp) = .true.
            ! Transfer the coordinates
            C(:,ldsp) = dbsc(iCnttp)%Coor(1:3,iCnt)
            ! Transfer the multiplicity factor
            C(4,ldsp) = Fact
            IndCar(ldsp) = iCar
          end if
        end do
      end do
    end do
    if (iPrint >= 99) then
      call RecPrt(' Information',' ',C,4,lDisp(0))
      write(LuWr,*) (IndCar(i),i=1,lDisp(0))
    end if

    ! Set up coefficient for the translational equations

    iTR = 0
    do i=1,3
      if (iSym(i) == 0) then
        iTR = iTR+1
        do ldsp=1,lDisp(0)
          if (IndCar(ldsp) == i) then
            Am(iTR,ldsp) = C(4,ldsp)
          end if
        end do
      end if
    end do

    ! Set up coefficient for the rotational invariance

    if (.not. T_Only) then
      do i=1,3
        j = i+1
        if (j > 3) j = j-3
        k = i+2
        if (k > 3) k = k-3
        ijSym = ieor(iSym(j),iSym(k))
        if (ijSym == 0) then
          iTR = iTR+1
          do ldsp=1,lDisp(0)
            if (IndCar(ldsp) == j) then
              Fact = C(4,ldsp)*C(k,ldsp)
              Am(iTR,ldsp) = Fact
            else if (IndCar(ldsp) == k) then
              Fact = -C(4,ldsp)*C(j,ldsp)
              Am(iTR,ldsp) = Fact
            end if
          end do
        end if
      end do
    end if
    if (iPrint >= 99) call RecPrt(' The A matrix',' ',Am,nTR,lDisp(0))

    ! Now, transfer the coefficient of those gradients which will
    ! not be computed directly.
    ! The matrix to compute the inverse of is determined via
    ! a Gram-Schmidt procedure.

    ! Pick up the other vectors
    do iTR=1,nTR
      !write(LuWr,*) ' Looking for vector #',iTR
      ovlp = Zero
      kTR = 0
      ! Check all the remaining vectors
      do ldsp=1,lDisp(0)
        Skip = .false.
        do jTR=1,iTR-1
          if (iTemp(jTR) == ldsp) then
            Skip = .true.
            exit
          end if
        end do
        if (.not. Skip) then
          !write(LuWr,*) ' Checking vector #',ldsp
          Temp(:,iTR) = AM(:,ldsp)
          !call RecPrt(' Vector',' ',Temp(1,iTR),nTR,1)
          ! Gram-Schmidt orthonormalize against accepted vectors
          do lTR=1,iTR-1
            alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,lTR),1)
            !write(LuWr,*) ' <x|y> =',alpha
            call DaXpY_(nTR,-alpha,Temp(1,lTR),1,Temp(1,iTR),1)
          end do
          !call RecPrt(' Remainings',' ',Temp(1,iTR),nTR,1)
          alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,iTR),1)
          !write(LuWr,*) ' Remaining overlap =',alpha
          ! Check the remaining magnitude of vector after Gram-Schmidt
          if (alpha > ovlp) then
            kTR = ldsp
            ovlp = alpha
          end if
          if ((.not. Dirct(ldsp)) .and. (alpha > 1.0e-2_wp)) then
            kTR = ldsp
            ovlp = huge(ovlp)
          end if
        end if
      end do
      Skip = .false.
      if (kTR == 0) then
        call WarningMessage(2,'Error in InputG')
        write(LuWr,*) ' No Vector found!'
        call Abend()
      end if
      !write(LuWr,*) ' Selecting vector #',kTR
      ! Pick up the "best" vector
      Temp(:,iTR) = Am(:,kTR)
      do lTR=1,iTR-1
        alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,lTR),1)
        call DaXpY_(nTR,-alpha,Temp(1,lTR),1,Temp(1,iTR),1)
      end do
      alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,iTR),1)
      call DScal_(nTR,One/sqrt(alpha),Temp(1,iTR),1)
      iTemp(iTR) = kTR
    end do
    do iTR=1,nTR
      Temp(:,iTR) = Am(:,iTemp(iTR))
      Am(:,iTemp(iTR)) = Zero
    end do
    if (iPrint >= 99) then
      call RecPrt(' The A matrix',' ',Am,nTR,lDisp(0))
      call RecPrt(' The T matrix',' ',Temp,nTR,nTR)
      write(LuWr,*) (iTemp(iTR),iTR=1,nTR)
    end if

    ! Compute the inverse of the T matrix

    call MatInvert(Temp,nTR)
    if (IPrint >= 99) call RecPrt(' The T-1 matrix',' ',Temp,nTR,nTR)
    call DScal_(nTR**2,-One,Temp,1)

    ! Generate the complete matrix

    call mma_allocate(Scr,nTR,lDisp(0),Label='Scr')
    call DGEMM_('N','N',nTR,lDisp(0),nTR,One,Temp,nTR,Am,nTR,Zero,Scr,nTR)
    if (IPrint >= 99) call RecPrt(' A-1*A',' ',Scr,nTR,lDisp(0))
    call mma_deallocate(Am)
    call mma_allocate(Am,lDisp(0),lDisp(0),Label='Am')
    call unitmat(Am,lDisp(0))
    do iTR=1,nTR
      ldsp = iTemp(iTR)
      call dcopy_(lDisp(0),Scr(1,iTR),nTR,Am(1,lDisp),lDisp(0))
    end do
    if (iPrint >= 99) call RecPrt('Final A matrix',' ',Am,lDisp(0),lDisp(0))

    call mma_deallocate(Scr)
    call mma_deallocate(IndCar)
    call mma_deallocate(C)
    call mma_deallocate(Temp)
    do iTR=1,nTR
      ldsp = iTemp(iTR)
      Dirct(ldsp) = .false.
    end do

    write(LuWr,*)
    write(LuWr,'(20X,A)') ' Automatic utilization of translational and rotational invariance of the energy is employed.'
    write(LuWr,*)
    do i=1,lDisp(0)
      if (Dirct(i)) then
        write(LuWr,'(25X,A,A)') Chdisp(i),' is independent'
      else
        write(LuWr,'(25X,A,A)') Chdisp(i),' is dependent'
      end if
    end do
    write(LuWr,*)
  end if

else
  nTR = 0
  if (Show .and. (iPrint >= 6)) then
    write(LuWr,*)
    write(LuWr,'(20X,A)') ' No automatic utilization of translational and rotational invariance of the energy is employed.'
    write(LuWr,*)
  end if
end if

call mma_deallocate(iTemp)

if (Slct .and. (.not. Skip)) then
  write(LuWr,*)
  write(LuWr,'(20X,A)') ' The Selection option is used'
  write(LuWr,*)
  do i=1,lDisp(0)
    if (Dirct(i)) then
      write(LuWr,'(25X,A,A)') Chdisp(i),' is computed'
    else
      write(LuWr,'(25X,A,A)') Chdisp(i),' is set to zero'
    end if
  end do
  write(LuWr,*)
end if

Onenly = HF_Force

return

contains

subroutine Error()
  if (istatus > 0) then
    call WarningMessage(2,'Error in InputG')
    write(LuWr,*) 'Inputg: Illegal keyword'
    write(LuWr,'(A,A)') 'KWord=',KWord
    call Quit_OnUserError()
  else if (istatus < 0) then
    call Abend()
  end if
end subroutine Error

end subroutine Inputg
