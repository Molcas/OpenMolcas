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
!               1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine Inputh(Run_MCLR)
!***********************************************************************
!                                                                      *
! Object: input module for the gradient code                           *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September 1991                                           *
!                                                                      *
!             Modified to complement GetInf, January 1992              *
!***********************************************************************

use McKinley_global, only: lGrd, lHss, nFck, Nona, PreScr, sIrrep
use MckDat, only: sNew
use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nBas, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: iChTbl, iOper, lBsFnc, lIrrep, nIrrep
use Gateway_global, only: Onenly, Test
use Gateway_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u5, u6

implicit none
logical(kind=iwp), intent(out) :: Run_MCLR
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
integer(kind=iwp) :: i, iCar, iCnt, iCnttp, iCo, iComp, idum, iDummer, iElem, iIrrep, ijSym, iOpt, ipert, iprint, iRC, iRout, &
                     istatus, iSym(3), iTR, j, jIrrep, jTR, k, kIrrep, kTR, ldsp, lTR, Lu_Mck, LuRd, mc, mdc, mDisp, nd(1), nDisp, &
                     nSlct
real(kind=wp) :: alpha, Fact, ovlp
logical(kind=iwp) :: defPert, ltype, Slct !, DoCholesky
character(len=80) :: Key, KWord
character(len=32) :: Label2
character(len=8) :: Label, labelop
integer(kind=iwp), allocatable :: ATDisp(:), Car(:), DEGDisp(:), iTemp(:), TDisp(:)
logical(kind=iwp), allocatable :: lPert(:)
real(kind=wp), allocatable :: AM(:,:), C(:,:), Scr(:,:), Tmp(:,:)
character, parameter :: xyz(0:2) = ['x','y','z']
integer(kind=iwp), external :: iPrmt, NrOpr
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: TstFnc

!call DecideOnCholesky(DoCholesky)
!if (DoCholesky) then
!  write(u6,*)'** Cholesky or RI/DF not yet implemented in McKinley '
!  call abend()
!end if

iRout = 99
nPrint(:) = 5
show = .false.
Onenly = .false.
Test = .false.
TRSymm = .false.
lEq = .false.
Slct = .false.
PreScr = .true.
lGrd = .true.
lHss = .true.
Nona = .false.
Run_MCLR = .true.
CutInt = 1.0e-7_wp
ipert = 2
defPert = .true.
call mma_allocate(lPert,0,label='lPert')
lPert(:) = defPert
sIrrep = .false.
iprint = 0
do i=1,3*MxAtom
  IndxEq(i) = i
end do

! KeyWord directed input

LuRd = u5
call RdNLst(LuRd,'MCKINLEY')
do
  read(u5,'(A72)',iostat=istatus) Key
  if (istatus < 0) call Error(1)
  if (istatus > 0) call Error(2)
  KWord = Key
  call UpCase(KWord)
  if (KWord(1:1) == '*') cycle
  if (KWord == '') cycle
  select case (KWord(1:4))
    !case ('EQUI')
    !  !                                                                *
    !  !***** EQUI ******************************************************
    !  !                                                                *
    !  ! Equivalence option
    !
    !  lEq = .true.
    !  do
    !    read(u5,'(A)',iostat=istatus) KWord
    !    if (istatus > 0) call Error(2)
    !    if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
    !  end do
    !  read(KWord,*) nGroup
    !  do iGroup=1,nGroup
    !    do
    !      read(u5,'(A)',iostat=istatus) KWord
    !      if (istatus > 0) call Error(2)
    !      if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
    !    end do
    !    read(KWord,*) nElem,(iTemp(iElem),iElem=1,nElem)
    !    do iElem=2,nElem
    !      IndxEq(iTemp(iElem)) = iTemp(1)
    !      Direct(iTemp(iElem)) = .false.
    !    end do
    !  end do

    !case ('NOIN')
    !  !                                                                *
    !  !***** NOIN ******************************************************
    !  !                                                                *
    !  ! Disable the utilization of translational and
    !  ! rotational invariance of the energy in the
    !  ! computation of the molecular gradient.
    !
    !  TRSymm = .false.

    case ('SHOW')
      !                                                                *
      !***** SHOW ******************************************************
      !                                                                *
      ! Raise the printlevel to show gradient contributions

      Show = .true.

    case ('CUTO')
      !                                                                *
      !***** CUTO ******************************************************
      !                                                                *
      ! Cutoff for computing primitive gradients

      !do
      !  read(u5,'(A)',iostat=istatus) KWord
      !  if (istatus > 0) call Error(2)
      !  if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      !end do
      !read(KWord,*,iostat=istatus) CutInt
      !if (istatus > 0) call Error(2)
      read(u5,*) Cutint
      CutInt = abs(CutInt)

    case ('VERB')
      !                                                                *
      !***** VERB ******************************************************
      !                                                                *
      ! Verbose output

      nPrint(1) = 6
      nPrint(99) = 6

    case ('NOSC')
      !                                                                *
      !***** NOSC ******************************************************
      !                                                                *
      ! Change default for the prescreening.

      PreScr = .false.

    case ('ONEO')
      !                                                                *
      !***** ONEO ******************************************************
      !                                                                *
      ! Do not compute two electron integrals.

      Onenly = .true.

    case ('SELE')
      !                                                                *
      !***** SELE ******************************************************
      !                                                                *
      ! selection option

      Slct = .true.
      defPert = .false.
      lPert(:) = defPert
      !do
      !  read(u5,'(A)',iostat=istatus) KWord
      !  if (istatus > 0) call Error(2)
      !  if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      !end do
      !read(KWord,*) nSlct
      read(u5,*) nSlct

      call mma_allocate(iTemp,nSlct,label='iTemp')
      read(u5,*) (iTemp(iElem),iElem=1,nSlct)
      do iElem=1,nSlct
        call extend_lPert(iTemp(iElem))
        lPert(iTemp(iElem)) = .true.
      end do
      call mma_deallocate(iTemp)

    case ('REMO')
      !                                                                *
      !***** REMO ******************************************************
      !                                                                *
      Slct = .true.
      read(u5,*) nSlct

      call mma_allocate(iTemp,nSlct,label='iTemp')
      read(u5,*) (iTemp(iElem),iElem=1,nSlct)
      do iElem=1,nSlct
        call extend_lPert(iTemp(iElem))
        lPert(iTemp(iElem)) = .false.
      end do
      call mma_deallocate(iTemp)

    case ('PERT')
      !                                                                *
      !***** PERT ******************************************************
      !                                                                *
      ! Select which part of the Hessian will be computed.

      do
        read(u5,'(A)',iostat=istatus) KWord
        if (istatus > 0) call Error(2)
        if ((KWord(1:1) /= '*') .and. (KWord /= '')) exit
      end do
      call UpCase(KWord)
      if (KWORD(1:4) == 'HESS') then
        ipert = 2
      else if (KWORD(1:4) == 'GEOM') then
        ipert = 1
      else
        write(u6,*) 'InputH: Illegal perturbation keyword'
        write(u6,'(A,A)') 'KWord=',KWord
        call Abend()
      end if

    case ('TEST')
      !                                                                *
      !***** TEST ******************************************************
      !                                                                *
      ! Process only the input.

      Test = .true.

    case ('EXTR')
      !                                                                *
      !***** EXTR ******************************************************
      !                                                                *
      ! Put the program name and the time stamp onto the extract file

      write(u6,*) 'InputH: EXTRACT option is redundant and is ignored!'

    case ('NONA')
      !                                                                *
      !***** NONA ******************************************************
      !                                                                *
      ! Compute the anti-symmetric overlap gradient only.

      Nona = .true.
      Run_MCLR = .false.

    case ('NOMC')
      !                                                                *
      !***** NOMC ******************************************************
      !                                                                *
      ! Request no automatic run of MCLR

      Run_MCLR = .false.

    case ('END ')
      !                                                                *
      !***** END  ******************************************************
      !                                                                *
      exit

    case default
      write(u6,*) 'InputH: Illegal keyword'
      write(u6,'(A,A)') 'KWord=',KWord
      call Abend()
  end select
end do
!***********************************************************************
!                                                                      *
!                          End of input section.                       *
!                                                                      *
!***********************************************************************

iPrint = nPrint(iRout)

iOpt = ibset(0,sNew)
if (onenly) iopt = 0
iRC = -1
Lu_Mck = 35
call OpnMck(irc,iOpt,'MCKINT',Lu_Mck)
if (iRC /= 0) then
  write(u6,*) 'InputH: Error opening MCKINT'
  call Abend()
end if
irc = -1
iopt = 0
if (ipert == 1) then
  Label2 = 'Geometry'
  LabelOp = 'PERT    '
  call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
  sIrrep = .true.
else if (ipert == 2) then
  Label2 = 'Hessian'
  LabelOp = 'PERT    '
  call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
else if (ipert == 3) then
  LabelOp = 'PERT    '
  Label2 = 'Magnetic'
  call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
  write(u6,*) 'InputH: Illegal perturbation option'
  write(u6,*) 'iPert=',iPert
  call Abend()
else if (ipert == 4) then
  LabelOp = 'PERT    '
  Label2 = 'Relativistic'
  call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
  write(u6,*) 'InputH: Illegal perturbation option'
  write(u6,*) 'iPert=',iPert
  call Abend()
else
  write(u6,*) 'InputH: Illegal perturbation option'
  write(u6,*) 'iPert=',iPert
  call Abend()
end if

!if (lEq) TRSymm = .false.
!if (Slct) TRSymm = .false.

mDisp = 0
mdc = 0
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    mDisp = mDisp+3*(nIrrep/dc(mdc)%nStab)
  end do
end do

write(u6,*)
write(u6,'(20X,A,E10.3)') ' Threshold for contributions to the gradient or Hessian:',CutInt
write(u6,*)

if (Nona) then
  write(u6,*)
  write(u6,'(20X,A)') ' McKinley only is computing the antisymmetric gradient of the overlap integrals for the NonAdiabatic '// &
                      'Coupling.'
  write(u6,*)
end if

! Generate symmetry adapted cartesian displacements

if (iPrint >= 6) then
  write(u6,*)
  write(u6,'(20X,A)') '********************************************'
  write(u6,'(20X,A)') '* Symmetry Adapted Cartesian Displacements *'
  write(u6,'(20X,A)') '********************************************'
  write(u6,*)
end if
IndDsp(:,:) = 0
InxDsp(:,:) = 0
call mma_allocate(ATDisp,mDisp,Label='ATDisp')
call mma_allocate(DEGDisp,mDisp,Label='DEGDisp')
nDisp = 0
do iIrrep=0,nIrrep-1
  lDisp(iIrrep) = 0
  ltype = .true.
  ! Loop over basis function definitions
  mdc = 0
  mc = 1
  do iCnttp=1,nCnttp
    ! Loop over unique centers associated with this basis set.
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      IndDsp(mdc,iIrrep) = nDisp
      ! Loop over the cartesian components
      do iCar=0,2
        iComp = 2**iCar
        if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)) then
          nDisp = nDisp+1
          if (nDisp > mDisp) then
            write(u6,*) 'nDisp > mDisp'
            call Abend()
          end if
          if (iIrrep == 0) InxDsp(mdc,iCar+1) = nDisp
          lDisp(iIrrep) = lDisp(iIrrep)+1
          if (ltype) then
            if (iPrint >= 6) then
              write(u6,*)
              write(u6,'(10X,A,A)') ' Irreducible representation : ',lIrrep(iIrrep)
              write(u6,'(10X,2A)') ' Basis function(s) of irrep: ',lBsFnc(iIrrep)
              write(u6,*)
              write(u6,'(A)') ' Basis Label        Type   Center Phase'
            end if
            ltype = .false.
          end if
          if (iPrint >= 6) &
            write(u6,'(I4,3X,A8,5X,A1,7X,8(I3,4X,I2,4X))') nDisp,dc(mdc)%LblCnt,xyz(iCar), &
                                                           (mc+iCo,iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iComp)* &
                                                                   iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))), &
                                                            iCo=0,nIrrep/dc(mdc)%nStab-1)
          write(ChDisp(nDisp),'(A,1X,A1)') dc(mdc)%LblCnt,xyz(iCar)
          ATDisp(ndisp) = icnttp
          DEGDisp(ndisp) = nIrrep/dc(mdc)%nStab
        end if

      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do
  end do

end do

if (nDisp /= mDisp) then
  write(u6,*) 'InputH: nDisp /= mDisp'
  write(u6,*) 'nDisp,mDisp=',nDisp,mDisp
  call Abend()
end if
if (sIrrep) then
  ndisp = ldisp(0)
  lDisp(1:nIrrep-1) = 0
end if
call mma_allocate(TDisp,nDisp,Label='TDisp')
TDisp(:) = 30
iOpt = 0
iRC = -1
labelOp = 'ndisp   '
nd(1) = ndisp
call WrMck(iRC,iOpt,labelop,1,nd,iDummer)
if (iRC /= 0) then
  write(u6,*) 'InputH: Error writing to MCKINT'
  write(u6,'(A,A)') 'labelOp=',labelOp
  call Abend()
end if
LABEL = 'DEGDISP'
iRc = -1
iOpt = 0
call WrMck(iRC,iOpt,Label,idum,DEGDISP,idum)
if (iRC /= 0) then
  write(u6,*) 'InputH: Error writing to MCKINT'
  write(u6,'(A,A)') 'LABEL=',LABEL
  call Abend()
end if
call mma_deallocate(DEGDisp)
LABEL = 'NRCTDISP'
iRc = -1
iOpt = 0
call WrMck(iRC,iOpt,Label,idum,ATDisp,idum)
if (iRC /= 0) then
  write(u6,*) 'InputH: Error writing to MCKINT'
  write(u6,'(A,A)') 'LABEL=',LABEL
  call Abend()
end if
call mma_deallocate(ATDisp)
LABEL = 'TDISP'
iRc = -1
iOpt = 0
call WrMck(iRC,iOpt,Label,idum,TDisp,idum)
if (iRC /= 0) then
  write(u6,*) 'InputH: Error writing to MCKINT'
  write(u6,'(A,A)') 'LABEL=',LABEL
  call Abend()
end if
call mma_deallocate(TDisp)

! Set up data for the utilization of the translational
! and rotational invariance of the energy.

if (TRSymm) then
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
  if (iPrint >= 99) write(u6,*) ' nTR=',nTR
  ! Rotational equations
  do i=1,3
    j = i+1
    if (j > 3) j = j-3
    k = i+2
    if (k > 3) k = k-3
    ijSym = ieor(iSym(j),iSym(k))
    if (ijSym == 0) nTR = nTR+1
  end do
  if (nTR == 0) then
    TRSymm = .false.
  else
    if (iPrint >= 99) write(u6,*) ' nTR=',nTR
    call mma_allocate(AM,nTR,lDisp(0),Label='AM')
    call mma_allocate(Tmp,nTR,nTR,Label='Tmp')
    call mma_allocate(C,4,lDisp(0),Label='C')
    call mma_allocate(Car,lDisp(0),Label='Car')

    AM(:,:) = Zero
    C(:,:) = Zero

    ! Generate temporary information of the symmetrical displacements.

    ldsp = 0
    mdc = 0
    iIrrep = 0
    do iCnttp=1,nCnttp
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
        do iCar=0,2
          iComp = 2**iCar
          if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)) then
            ldsp = ldsp+1
            ! Transfer the coordinates
            C(1:3,ldsp) = dbsc(iCnttp)%Coor(:,iCnt)
            ! Transfer the multiplicity factor
            C(4,ldsp) = Fact
            Car(ldsp) = iCar+1
          end if
        end do
      end do
    end do
    if (iPrint >= 99) then
      call RecPrt(' Information',' ',C,4,lDisp(0))
      write(u6,*) (Car(i),i=1,lDisp(0))
    end if

    ! Set up coefficient for the translational equations

    iTR = 0
    do i=1,3
      if (iSym(i) /= 0) cycle
      iTR = iTR+1
      do ldsp=1,lDisp(0)
        if (Car(ldsp) == i) AM(iTR,ldsp) = C(4,ldsp)
      end do
    end do

    ! Set up coefficient for the rotational invariance

    do i=1,3
      j = i+1
      if (j > 3) j = j-3
      k = i+2
      if (k > 3) k = k-3
      ijSym = ieor(iSym(j),iSym(k))
      if (ijSym /= 0) cycle
      iTR = iTR+1
      do ldsp=1,lDisp(0)
        if (Car(ldsp) == j) then
          Fact = C(4,ldsp)*C(k,ldsp)
        else if (Car(ldsp) == k) then
          Fact = -C(4,ldsp)*C(j,ldsp)
        else
          Fact = Zero
          write(u6,*) 'Inputh: Error'
          call Abend()
        end if
        AM(iTR,ldsp) = Fact
      end do
    end do
    if (iPrint >= 99) call RecPrt(' The A matrix',' ',AM,nTR,lDisp(0))

    ! Now, transfer the coefficient of those gradients which will
    ! not be computed directly.
    ! The matrix to compute the inverse of is determined via
    ! a Gram-Schmidt procedure.

    ! Pick up the other vectors
    do iTR=1,nTR
      !write(u6,*) ' Looking for vector #',iTR
      ovlp = Zero
      kTR = 0
      ! Check all the remaining vectors
      loop1: do ldsp=1,lDisp(0)
        do jTR=1,iTR-1
          if (iTemp(jTR) == ldsp) cycle loop1
        end do
        !write(u6,*) ' Checking vector #', ldsp
        Tmp(:,iTR) = AM(:,ldsp)
        !call RecPrt(' Vector',' ',Tmp(:,iTR),nTR,1)
        ! Gram-Schmidt orthonormalize against accepted vectors
        do lTR=1,iTR-1
          alpha = DDot_(nTR,Tmp(:,iTR),1,Tmp(:,lTR),1)
          !write(u6,*) ' <x|y> =', alpha
          Tmp(:,iTR) = Tmp(:,iTR)-alpha*Tmp(:,lTR)
        end do
        !call RecPrt(' Remainings',' ',Tmp(:,iTR),nTR,1)
        alpha = DDot_(nTR,Tmp(:,iTR),1,Tmp(:,iTR),1)
        !write(u6,*) ' Remaining overlap =', alpha
        ! Check the remaining magnitude of vector after Gram-Schmidt
        if (alpha > ovlp) then
          kTR = ldsp
          ovlp = alpha
        end if
      end do loop1
      if (kTR == 0) then
        write(u6,*) ' No Vector found!'
        call Abend()
      end if
      !write(u6,*) ' Selecting vector #', kTR
      ! Pick up the "best" vector
      Tmp(:,iTR) = AM(:,kTR)
      do lTR=1,iTR-1
        alpha = DDot_(nTR,Tmp(:,iTR),1,Tmp(:,lTR),1)
        Tmp(:,iTR) = Tmp(:,iTR)-alpha*Tmp(:,lTR)
      end do
      alpha = DDot_(nTR,Tmp(:,iTR),1,Tmp(:,iTR),1)
      Tmp(:,iTR) = Tmp(:,iTR)/sqrt(alpha)
      iTemp(iTR) = kTR
    end do
    do iTR=1,nTR
      Tmp(:,iTR) = AM(:,iTemp(iTR))
      AM(:,iTemp(iTR)) = Zero
    end do
    if (iPrint >= 99) then
      call RecPrt(' The A matrix',' ',AM,nTR,lDisp(0))
      call RecPrt(' The T matrix',' ',Tmp,nTR,nTR)
      write(u6,*) (iTemp(iTR),iTR=1,nTR)
    end if

    ! Compute the inverse of the T matrix

    call MatInvert(Tmp,nTR)
    if (IPrint >= 99) call RecPrt(' The T-1 matrix',' ',Tmp,nTR,nTR)
    Tmp(:,:) = -Tmp

    ! Generate the complete matrix

    call mma_allocate(Scr,nTR,lDisp(0),Label='Scr')
    call DGEMM_('N','N',nTR,lDisp(0),nTR,One,Tmp,nTR,AM,nTR,Zero,Scr,nTR)
    if (IPrint >= 99) call RecPrt(' A-1*A',' ',Scr,nTR,lDisp(0))
    call mma_deallocate(AM)
    call mma_allocate(AM,lDisp(0),lDisp(0),Label='AM')
    call unitmat(AM,lDisp(0))
    do iTR=1,nTR
      ldsp = iTemp(iTR)
      AM(ldsp,:) = Scr(iTR,1:lDisp(0))
    end do
    if (iPrint >= 99) call RecPrt('Final A matrix',' ',AM,lDisp(0),lDisp(0))

    call mma_deallocate(Scr)
    call mma_deallocate(Car)
    call mma_deallocate(C)
    call mma_deallocate(Tmp)
    do iTR=1,nTR
      ldsp = iTemp(iTR)
      call extend_lPert(ldsp)
      lPert(ldsp) = .false.
    end do

    write(u6,*)
    write(u6,'(20X,A)') ' Automatic utilization of translational and rotational invariance of the energy is employed.'
    write(u6,*)
    call extend_lPert(lDisp(0))
    do i=1,lDisp(0)
      if (lPert(i)) then
        write(u6,'(25X,A,A)') Chdisp(i),' is independent'
      else
        write(u6,'(25X,A,A)') Chdisp(i),' is dependent'
      end if
    end do
    write(u6,*)
  end if

else
  nTR = 0
  if (iPrint >= 6) then
    write(u6,*)
    write(u6,'(20X,A)') ' No automatic utilization of translational and rotational invariance of the energy is employed.'
    write(u6,*)
  end if
end if

if (Slct) then
  write(u6,*)
  write(u6,'(20X,A)') ' The Selection option is used'
  write(u6,*)
  call extend_lPert(lDisp(0))
  do i=1,lDisp(0)
    if (lPert(i)) then
      write(u6,'(25X,A,A)') Chdisp(i),' is computed'
    else
      write(u6,'(25X,A,A)') Chdisp(i),' is set to zero'
    end if
  end do
  write(u6,*)
end if

call Datimx(KWord)
nFck(0:nIrrep-1) = 0
do iIrrep=0,nIrrep-1
  if (iIrrep /= 0) then
    do jIrrep=0,nIrrep-1
      kIrrep = NrOpr(ieor(ioper(jIrrep),ioper(iIrrep)))
      if (kIrrep < jIrrep) nFck(iIrrep) = nFck(iIrrep)+nBas(jIrrep)*nBas(kIrrep)
    end do
  else
    do jIrrep=0,nIrrep-1
      nFck(0) = nFck(0)+nTri_Elem(nBas(jIrrep))
    end do
  end if
end do

call mma_deallocate(lPert)

return

contains

subroutine Extend_lPert(n)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp) :: m
  logical(kind=iwp), allocatable :: Temp(:)

  m = size(lPert)
  if (n > m) then
    call mma_allocate(Temp,n,label='lPert')
    Temp(:m) = lPert
    Temp(m+1:) = defPert
    call mma_deallocate(lPert)
    call move_alloc(Temp,lPert)
  end if

end subroutine Extend_lPert

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      write(u6,*) 'InputH: end of input file.'
    case (2)
      write(u6,*) 'InputH: error reading input file.'
  end select
  write(u6,'(A,A)') 'Last command=',KWord
  call Abend()

end subroutine Error

end subroutine Inputh
