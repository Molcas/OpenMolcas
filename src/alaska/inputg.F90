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

use Alaska_Info, only: Am
use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep, iChTbl, iOper, lIrrep, lBsFnc
use Temporary_Parameters
use Real_Info, only: CutInt
use OFembed, only: Do_OFemb, KEonly, OFE_first, Xsigma, dFMD, OFE_KSDFT

implicit real*8(A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "disp.fh"
#include "iavec.fh"
#include "stdalloc.fh"
#include "columbus_gamma.fh"
#include "exterm.fh"
#include "nac.fh"
#include "alaska_root.fh"
logical TstFnc, type, Slct, T_Only, No_Input_OK
real*8, allocatable :: Tmp(:), C(:,:), Scr(:,:), Temp(:,:)
integer, allocatable :: IndCar(:)

character(LEN=1) :: xyz(0:2) = ['x','y','z']
character(LEN=80) KWord, Key
integer iSym(3), iTemp(3*MxAtom)
logical Reduce_Prt
external Reduce_Prt
#include "chotime.fh"

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
Xsigma = 1.0d4
dFMD = 0.0d0
Do_OFemb = .false.
KEonly = .false.
OFE_first = .true.
Show = .true.
LuWr = 6
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

! First CutGrd can not be more accurate than CutInt!
CutGrd = max(1.0D-07,CutInt)
! Second CutInt should now locally for Alaska be reset to the value
! of CutInt/100!
CutInt = CutGrd*1.0D-2
do i=1,3*MxAtom
  IndxEq(i) = i
end do
do ldsp=1,3*MxAtom
  Direct(ldsp) = .true.
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! KeyWord directed input

rewind(LuSpool)
No_Input_OK = .true.
call RdNLst_(LuSpool,'ALASKA',No_Input_OK)
KWord = ' &ALASKA'
998 read(LuSpool,'(A72)',end=997,ERR=988) Key
KWord = Key
call UpCase(KWord)
if (KWord(1:1) == '*') Go To 998
if (KWord == '') Go To 998
if (KWord(1:4) == 'VERB') Go To 912
if (KWord(1:4) == 'PRIN') Go To 930
if (KWord(1:4) == 'EQUI') Go To 935
if (KWord(1:4) == 'CUTO') Go To 942
if (KWord(1:4) == 'HF-F') Go To 993
if (KWord(1:4) == 'NOIN') Go To 953
if (KWord(1:4) == 'SELE') Go To 960
if (KWord(1:4) == '2DOP') Go To 965
if (KWord(1:4) == '2DIP') Go To 966
if (KWord(1:4) == 'ONEO') Go To 990
if (KWord(1:4) == 'TEST') Go To 991
if (KWord(1:4) == 'SHOW') Go To 992
if (KWord(1:4) == 'PNEW') Go To 994
if (KWord(1:4) == 'POLD') Go To 995
if (KWord(1:4) == 'NONU') Go To 996
if (KWord(1:4) == 'EXTR') Go To 971
if (KWord(1:4) == 'CHOI') Go To 972
if (KWord(1:4) == 'OFEM') Go To 973
if (KWord(1:4) == 'KEON') Go To 974
if (KWord(1:4) == 'DFMD') Go To 975
! Keyword 'NUMErical' checked earlier - forces numerical gradients
! Keyword 'DELTa' selects the scaling factor for the displacements
!                 in the numerical_gradient module
! Keyword 'KEEP' does not remove the old gradient
! Keyword 'INVErt' inverts the treatment of constraints
! Here it's only included for consistency
if (KWord(1:4) == 'NUME') Go To 998
if (KWord(1:4) == 'DELT') Go To 998
if (KWord(1:4) == 'KEEP') Go To 998
if (KWord(1:4) == 'INVE') Go To 998
if (KWord(1:4) == 'ROOT') Go To 976
if (KWord(1:4) == 'NAC ') Go To 977
if (KWord(1:4) == 'NOCS') Go To 978
if (KWord(1:4) == 'AUTO') Go To 979
if (KWord(1:4) == 'END ') Go To 997
call WarningMessage(2,'Error in InputG')
write(LuWr,*) 'Inputg: Illegal keyword'
write(LuWr,'(A,A)') 'KWord=',KWord
call Quit_OnUserError()

988 call WarningMessage(2,'Error in InputG')
write(LuWr,*) 'Inputg: Error reading the input'
write(LuWr,'(A,A)') 'Last read line=',KWord
call Quit_OnUserError()
!                                                                      *
!***********************************************************************
!                                                                      *
! Print level

930 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 930
if (KWord == '') Go To 930
read(KWord,*,Err=988) n
do i=1,n
9301 read(LuSpool,'(A)',Err=988) KWord
  if (KWord(1:1) == '*') Go To 9301
  if (KWord == '') Go To 9301
  read(KWord,*,Err=988) jRout,iPrint
  nPrint(jRout) = iPrint
end do
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Equivalence option

935 continue
if (T_Only) then
  call WarningMessage(2,'Error in InputG')
  write(LuWr,*) 'EQUI option does not ork with RF calculations!'
  call Quit_OnUserError()
end if
lEq = .true.
936 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 936
if (KWord == '') Go To 936
read(KWord,*) nGroup
do iGroup=1,nGroup
938 read(LuSpool,'(A)',Err=988) KWord
  if (KWord(1:1) == '*') Go To 938
  if (KWord == '') Go To 938
  read(KWord,*) nElem,(iTemp(iElem),iElem=1,nElem)
  do iElem=2,nElem
    IndxEq(iTemp(iElem)) = iTemp(1)
    Direct(iTemp(iElem)) = .false.
  end do
end do
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Cutoff for computing primitive gradients

942 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 942
if (KWord == '') Go To 942
read(KWord,*,Err=988) CutGrd
CutGrd = abs(CutGrd)
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Disable the utilization of translational and
! rotational invariance of the energy in the
! computation of the molecular gradient.

953 TRSymm = .false.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! selection option

960 continue
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
  Direct(i) = .false.
end do
962 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 962
if (KWord == '') Go To 962
read(KWord,*) nSlct

963 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 963
if (KWord == '') Go To 963
read(KWord,*)(iTemp(iElem),iElem=1,nSlct)
do iElem=1,nSlct
  Direct(iTemp(iElem)) = .true.
end do
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Change default for the prescreening.

965 l2DI = .false.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Change default for the prescreening.

966 l2DI = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Do not compute two electron integrals.

990 Onenly = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Process only the input.

991 Test = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Raise the printlevel to show gradient contributions

992 continue
if (iPL >= 2) then
  nPrint(112) = 15
  nPrint(1) = 15
  nPrint(33) = 15
end if
Go To 998
!                                                                      *
!***** PNEW ************************************************************
!                                                                      *
! Print gradient in NEW human-readable format

994 continue
nPrint(1) = 4
Go To 998
!                                                                      *
!***** POLD ************************************************************
!                                                                      *
! Print gradient in OLD format

995 continue
nPrint(1) = 5
Go To 998
!                                                                      *
!***** VERB ************************************************************
!                                                                      *
! Verbose mode.

912 continue
nPrint(80) = 6
nPrint(1) = 6
nPrint(9) = 6
nPrint(99) = 6
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute Hellmann-Feynman forces

993 HF_Force = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Do not compute the nuclear charge contribution

996 NO_NUC = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Put the program name and the time stamp onto the extract file

971 write(LuWr,*) 'InputG: EXTRACT option is redundant and is ignored!'
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Cholesky input section

972 continue
call Cho_alaska_rdInp(LuSpool)
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Orbital-Free Embedding (OFE) input section

973 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 973
if (KWord == '') Go To 973
call UpCase(KWord)
call LeftAd(KWord)
read(KWord,'(A)') OFE_KSDFT
Do_OFemb = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Mode "Kinetic Energy Only" for OFE input section

974 continue
KEonly = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Mode "Kinetic Energy Only" for OFE input section

975 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 975
if (KWord == '') Go To 975
read(KWord,*) dFMD,Xsigma
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! Root keyword, now also for analytical gradient
! This is a dummy, the keyword is already read in chk_numerical

976 read(LuSpool,'(A)',Err=988) KWord
if (KWord(1:1) == '*') Go To 976
if (KWord == '') Go To 976
read(KWord,*) i
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! NAC keyword: compute non-adiabatic couplings between 2 states
! The keyword is already read in chk_numerical

977 read(LuSpool,'(A)',Err=988) KWord
isNAC = .true.
if (KWord(1:1) == '*') Go To 977
if (KWord == '') Go To 977
read(KWord,*) NACstates(1),NACstates(2)
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! NOCSF keyword, to neglect the CSF contribution to the NAC,
! which is the cause for translational variance

978 DoCSF = .false.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
! AUTO keyword, used by SLAPAF, to signal this is an automated
! call to ALASKA

979 Auto = .true.
Go To 998
!                                                                      *
!***********************************************************************
!                                                                      *
!                          End of input section.                       *
!                                                                      *
!***********************************************************************
!                                                                      *
997 continue

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
    Go To 10
  else if ((dbsc(iCnttp)%nFragType > 0) .or. dbsc(iCnttp)%Frag) then
    TRSymm = .false.
  end if
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    mDisp = mDisp+3*(nIrrep/dc(mdc)%nStab)
  end do
10 continue
end do

if (HF_Force .and. Show .and. (iPrint >= 6)) then
  write(LuWr,*)
  write(LuWr,'(A)') '            O B S E R V E ! '
  write(LuWr,'(A)') '            Option for computation of interstate couling vector or'
  write(LuWr,'(A)') '            Hellmann-Feynman gradient is active.'

  write(LuWr,*)
end if
if (Show .and. (iPrint >= 6)) then
  write(LuWr,*)
  write(LuWr,'(20X,A,E10.3)') ' Threshold for contributions to the gradient:',CutGrd
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

call ICopy(MxAtom*8,[0],0,IndDsp,1)
call ICopy(MxAtom*3,[0],0,InxDsp,1)
call dcopy_(3*MxSym*MxAtom,[One],0,Disp_Fac,1)
call ICopy(3*MxAtom,[1],0,mult_Disp,1)
nDisp = 0
do iIrrep=0,nIrrep-1
  lDisp(iIrrep) = 0
  type = .true.
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
          if (type) then
            if (Show .and. (iPrint >= 6)) then
              write(LuWr,*)
              write(LuWr,'(10X,A,A)') ' Irreducible representation : ',lIrrep(iIrrep)
              write(LuWr,'(10X,2A)') ' Basis function(s) of irrep: ',lBsFnc(iIrrep)
              write(LuWr,*)
              write(LuWr,'(A)') ' Basis Label        Type   Center Phase'
            end if
            type = .false.
          end if
          if (iIrrep == 0) then
            do jOper=0,nIrrep-1
              Disp_Fac(iCar+1,jOper,mdc) = dble(iPrmt(jOper,iComp)*iChTbl(iIrrep,jOper))
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
  write(6,*) 'Unsupported option: TRSymm'
  call Abend()
  iSym(1) = 0
  iSym(2) = 0
  iSym(3) = 0
  do i=1,min(nIrrep-1,5)
    j = i
    if (i == 3) j = 4
    do k=1,3
      if (iand(iOper(j),2**(k-1)) /= 0) iSym(k) = 2**(k-1)
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
    Go To 9876
  end if
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
      if (dbsc(iCnttp)%Coor(1,iCnt) /= Zero) iComp = ior(iComp,1)
      if (dbsc(iCnttp)%Coor(2,iCnt) /= Zero) iComp = ior(iComp,2)
      if (dbsc(iCnttp)%Coor(3,iCnt) /= Zero) iComp = ior(iComp,4)
      do jIrrep=0,nIrrep-1
        if (TstFnc(dc(mdc)%iCoSet,jIrrep,iComp,dc(mdc)%nStab)) then
          Fact = Fact+One
        end if
      end do
      do iCar=1,3
        iComp = 2**(iCar-1)
        if (TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)) then
          ldsp = ldsp+1
          Direct(lDsp) = .true.
          ! Transfer the coordinates
          call dcopy_(3,dbsc(iCnttp)%Coor(1:3,iCnt),1,C(1,ldsp),1)
          ! Transfer the multiplicity factor
          C(4,ldsp) = Fact
          IndCar(ldsp) = iCar
        end if
      end do
    end do
  end do
  if (iPrint >= 99) then
    call RecPrt(' Information',' ',C,4,lDisp(0))
    write(LuWr,*)(IndCar(i),i=1,lDisp(0))
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
      if (ijSym /= 0) Go To 1210
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
1210  continue
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
      do jTR=1,iTR-1
        if (iTemp(jTR) == ldsp) Go To 1231
      end do
      !write(LuWr,*) ' Checking vector #',ldsp
      call dcopy_(nTR,Am(1,ldsp),1,Temp(1,iTR),1)
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
      if ((.not. Direct(ldsp)) .and. (alpha > 1.0D-2)) then
        kTR = ldsp
        ovlp = 1.0d99
      end if
1231  continue
    end do
    if (kTR == 0) then
      call WarningMessage(2,'Error in InputG')
      write(LuWr,*) ' No Vector found!'
      call Abend()
    end if
    !write(LuWr,*) ' Selecting vector #',kTR
    ! Pick up the "best" vector
    call dcopy_(nTR,Am(1,kTR),1,Temp(1,iTR),1)
    do lTR=1,iTR-1
      alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,lTR),1)
      call DaXpY_(nTR,-alpha,Temp(1,lTR),1,Temp(1,iTR),1)
    end do
    alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,iTR),1)
    call DScal_(nTR,One/sqrt(alpha),Temp(1,iTR),1)
    iTemp(iTR) = kTR
  end do
  do iTR=1,nTR
    call dcopy_(nTR,Am(1,iTemp(iTR)),1,Temp(1,iTR),1)
    Am(:,iTemp(iTR)) = Zero
  end do
  if (iPrint >= 99) then
    call RecPrt(' The A matrix',' ',Am,nTR,lDisp(0))
    call RecPrt(' The T matrix',' ',Temp,nTR,nTR)
    write(LuWr,*)(iTemp(iTR),iTR=1,nTR)
  end if

  ! Compute the inverse of the T matrix

  call MatInvert(Temp,nTR)
  if (IPrint >= 99) call RecPrt(' The T-1 matrix',' ',Temp,nTR,nTR)
  call DScal_(nTR**2,-One,Temp,1)

  ! Generate the complete matrix

  call mma_allocate(Scr,nTR,lDisp(0),Label='Scr')
  call DGEMM_('N','N',nTR,lDisp(0),nTR,1.0d0,Temp,nTR,Am,nTR,0.0d0,Scr,nTR)
  if (IPrint >= 99) call RecPrt(' A-1*A',' ',Scr,nTR,lDisp(0))
  call mma_deallocate(Am)
  call mma_allocate(Am,lDisp(0),lDisp(0),Label='Am')
  Am(:,:) = Zero
  do i=1,lDisp(0)
    Am(i,i) = One
  end do
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
    Direct(ldsp) = .false.
  end do

  write(LuWr,*)
  write(LuWr,'(20X,A)') ' Automatic utilization of translational and rotational invariance of the energy is employed.'
  write(LuWr,*)
  do i=1,lDisp(0)
    if (Direct(i)) then
      write(LuWr,'(25X,A,A)') Chdisp(i),' is independent'
    else
      write(LuWr,'(25X,A,A)') Chdisp(i),' is dependent'
    end if
  end do
  write(LuWr,*)

else
  nTR = 0
  if (Show .and. (iPrint >= 6)) then
    write(LuWr,*)
    write(LuWr,'(20X,A)') ' No automatic utilization of translational and rotational invariance of the energy is employed.'
    write(LuWr,*)
  end if
end if

if (Slct) then
  write(LuWr,*)
  write(LuWr,'(20X,A)') ' The Selection option is used'
  write(LuWr,*)
  do i=1,lDisp(0)
    if (Direct(i)) then
      write(LuWr,'(25X,A,A)') Chdisp(i),' is computed'
    else
      write(LuWr,'(25X,A,A)') Chdisp(i),' is set to zero'
    end if
  end do
  write(LuWr,*)
end if

9876 continue

! Set up the angular index vector

i = 0
do iR=0,iTabMx
  do ix=iR,0,-1
    do iy=iR-ix,0,-1
      iz = iR-ix-iy
      i = i+1
      ixyz(1,i) = ix
      ixyz(2,i) = iy
      ixyz(3,i) = iz
    end do
  end do
end do

Onenly = HF_Force

return

end subroutine Inputg
