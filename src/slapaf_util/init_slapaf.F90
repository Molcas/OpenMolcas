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

subroutine Init_SlapAf()

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: Analytic_Hessian, ANr, ApproxNADC, AtomLbl, Coor, Degen, dMass, Force_dB, Grd, Header, iCoSet, jStab, &
                       Line_Search, MaxItr, mB_Tot, mdB_Tot, mq, mTtAtm, MxItr, NADC, nDimBC, nStab, q_nuclear, RootMap, Smmtrc, &
                       ThrCons, ThrEne, ThrGrd, ThrMEP !, lRP, R12
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "print.fh"
integer(kind=iwp) :: Columbus, i, iAdd(0:7), iChxyz, iComp, iIrrep, iMAX, iMode, ind, iPL, isAtom, ISPIN1, ISPIN2, itest, jCoSet, &
                     jPrint, jTest, LSYM1, LSYM2, n, nCoSet, nHess, nRM, nRoots, nStb
real(kind=wp) :: tmp
logical(kind=iwp) :: Do_ESPF, Exist_2, Found, Same, Skip
character(len=8) :: CMAX
real(kind=wp), allocatable :: xMass(:)
integer(kind=iwp), external :: iDeg, iPrintLevel, iPrmt
logical(kind=iwp), external :: Reduce_Prt

!
!***********************************************************************
!************************* StartUp section   ***************************
!***********************************************************************
!                                                                      *
! Set the default value of iterations from MOLCAS_MAXITER if it
! has been defined.

call GetEnvf('MOLCAS_MAXITER',CMAX)
!write(u6,'(3A)') 'CMAX="',CMAX,'"'
if (CMAX /= ' ') then
  read(CMAX,'(I8)') iMAX
  MxItr = min(MaxItr,iMax)
else
  MxItr = MaxItr
end if
!                                                                      *
!***********************************************************************
!                                                                      *
jPrint = 10
!                                                                      *
!***********************************************************************
!                                                                      *
call DecideOnESPF(Do_ESPF)
if (Do_ESPF) then
  ThrGrd = 0.003_wp
  ThrEne = 1.0e-5_wp
  Line_Search = .false.
else
  ThrGrd = 0.0003_wp
  ThrEne = 1.0e-6_wp
  Line_Search = .true.
end if
ThrMEP = ThrGrd
ThrCons = 1.0e10_wp
!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (iPL == 2) then
  iPL = 5
else if (iPL == 3) then
  iPL = 6
else if (iPL == 4) then
  iPL = 99
else if (iPL == 5) then
  iPL = 99
end if
nPrint(:) = iPL

! Reduced print level of Slapaf parameters after the first iteration

if (Reduce_Prt() .and. (iPL <= 5)) nPrint(:) = iPL-1

!                                                                      *
!***********************************************************************
!                                                                      *
! Get Molecular data

! Read the title

call Get_cArray('Seward Title',Header,144)

! Read number of atoms, charges, coordinates, gradients and atom labels

call Get_Molecule()
!                                                                      *
!***********************************************************************
!                                                                      *
NADC = .false.
ApproxNADC = .false.
call Get_iScalar('Columbus',Columbus)
if (Columbus == 1) then

  ! C&M mode

  call Get_iScalar('ColGradMode',iMode)
  if (iMode == 3) NADC = .true.
else

  ! M mode

  ! ISPIN should only be found for RASSCF-based
  ! methods, so no CI mode for SCF, MP2, etc. (or that's the idea)

  !write(u6,*) 'See if CI'
  call Qpg_iScalar('ISPIN',Found)
  if (Found) then
    call Get_iScalar('ISPIN',ISPIN1)
    call Get_iScalar('STSYM',LSYM1)
  else
    ISPIN1 = 0
    LSYM1 = 0
  end if
  !write(u6,*) 'iSpin=',ISPIN1
  !write(u6,*) 'stSym=',LSYM1

  call f_Inquire('RUNFILE2',Exist_2)
  !write(u6,*) 'Exist_2=',Exist_2
  if (Exist_2) then
    call NameRun('RUNFILE2')
    call Qpg_iScalar('ISPIN',Found)
    if (Found) then
      call Get_iScalar('ISPIN',ISPIN2)
      call Get_iScalar('STSYM',LSYM2)
    else
      ISPIN2 = 0
      LSYM2 = 0
    end if
    call NameRun('#Pop')
  else
    ISPIN2 = ISPIN1
    LSYM2 = LSYM1
  end if
  !write(u6,*) 'iSpin=',ISPIN1,ISPIN2
  !write(u6,*) 'stSym=',LSYM1,LSYM2

  ! Do not add the constraint at the NumGrad stage

  if (SuperName /= 'numerical_gradient') then
    if ((ISPIN1 /= 0) .and. (LSYM1 /= 0)) NADC = (ISPIN1 == ISPIN2) .and. (LSYM1 == LSYM2)
    !NADC = .false. ! for debugging
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read or initialize the root map

call Qpg_iArray('Root Mapping',Found,nRM)
if (nRM > 0) then
  call mma_allocate(RootMap,nRM,Label='RootMap')
  call Get_iArray('Root Mapping',RootMap,nRM)
else
  call Qpg_iScalar('Number of roots',Found)
  nRoots = 1
  if (Found) call Get_iScalar('Number of roots',nRoots)
  call mma_allocate(RootMap,nRoots,Label='RootMap')
  do i=1,nRoots
    RootMap(i) = i
  end do
end if

! Check if there is an analytic Hessian
call qpg_dArray('Analytic Hessian',Analytic_Hessian,nHess)

if (.not. Analytic_Hessian) then
  call NameRun('RUNOLD')
  call qpg_dArray('Analytic Hessian',Analytic_Hessian,nHess)
  call NameRun('#Pop')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of total symmetric displacements

call mma_allocate(jStab,[0,7],[1,size(Coor,2)],Label='jStab ')
call mma_allocate(nStab,[1,size(Coor,2)],Label='nStab ')
call mma_allocate(iCoSet,[0,7],[1,size(Coor,2)],Label='iCoSet')
call mma_allocate(Smmtrc,3,size(Coor,2),Label='Smmtrc')
jStab(:,:) = 0
nStab(:) = 0
iCoSet(:,:) = 0
Smmtrc(:,:) = .false.

nDimbc = 0
! Loop over the unique atoms
do isAtom=1,size(Coor,2)
  ! Find character of center
  iChxyz = 0
  do i=1,3
    if (Coor(i,isAtom) /= Zero) then
      do iIrrep=0,nIrrep-1
        if (btest(iOper(iIrrep),i-1)) iChxyz = ibset(iChxyz,i-1)
      end do
    end if
  end do
  nStb = 0
  do iIrrep=0,nIrrep-1
    if (iand(iChxyz,iOper(iIrrep)) == 0) then
      jStab(nStb,isAtom) = iOper(iIrrep)
      nStb = nStb+1
    end if
  end do
  nStab(isAtom) = nStb
  ! Find the coset representatives
  iCoSet(0,size(Coor,2)) = 0      ! Put in the unit operator
  nCoSet = 1
  do iIrrep=1,nIrrep-1
    itest = iand(iChxyz,iOper(iIrrep))
    Same = .false.
    do jCoSet=0,nCoSet-1
      jTest = iand(iChxyz,iCoSet(jCoSet,isAtom))
      Same = jTest == iTest
      if (Same) exit
    end do
    if (.not. Same) then
      nCoSet = nCoSet+1
      iCoSet(nCoSet-1,isAtom) = iOper(iIrrep)
    end if
  end do
  if (nIrrep/nStb /= nCoSet) then
    call WarningMessage(2,' Error while doing cosets.')
    call Abend()
  end if
  do i=1,3
    iComp = ibset(0,i-1)
    iAdd(0:nCoSet-1) = 0
    do iIrrep=0,nIrrep-1
      ! find the stabilizer index
      iTest = iand(iChxyz,iOper(iIrrep))
      n = -1
      do jCoset=0,nCoset-1
        jTest = iand(iChxyz,iCoset(jCoSet,isAtom))
        if (iTest == jTest) n = jCoset
      end do
      if ((n < 0) .or. (n > nCoset-1)) then
        call WarningMessage(2,' Error finding coset element')
        call Abend()
      end if
      iAdd(n) = iAdd(n)+iPrmt(iIrrep,iComp)
    end do
    Skip = .false.
    do jCoSet=0,nCoSet-1
      if (iAdd(jCoSet) == 0) then
        Skip = .true.
        exit
      end if
    end do
    if (.not. Skip) then
      nDimbc = nDimbc+1
      Smmtrc(i,isAtom) = .true.
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Transform charges to masses (C=12)

call mma_allocate(dMass,size(Coor,2),Label='dMass')
call mma_allocate(xMass,size(Coor,2),Label='xMass')
call Get_Mass(xMass,size(Coor,2))
!call RecPrt(' Charges',' ',Q_nuclear,SIZE(Coor,2),1)
call mma_allocate(ANr,size(Coor,2),Label='ANr')
do isAtom=1,size(Coor,2)
  ind = int(Q_nuclear(isAtom))
  if (ind <= 0) then
    dMass(isAtom) = 1.0e-10_wp
  else
    dMass(isAtom) = xMass(isAtom)
  end if
  ANr(isAtom) = ind
end do
call mma_deallocate(xMass)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the multiplicities of the cartesian coordinates and the
! total number of atoms.

mTtAtm = 0
call mma_Allocate(Degen,3,size(Coor,2),Label='Degen')
do isAtom=1,size(Coor,2)
  mTtAtm = mTtAtm+iDeg(Coor(:,isAtom))
  tmp = real(iDeg(Coor(:,isAtom)),kind=wp)
  Degen(:,isAtom) = tmp
end do
#ifdef _DEBUGPRINT_
call RecPrt('Degen',' ',Degen,3,size(Coor,2))
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!call qpg_dArray('Transverse',lRP,nRP)
!if (lRP) then
!  call mma_allocate(R12,3,nRP/3,Label='R12')
!  call Get_dArray('Transverse',R12,nRP)
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute center of mass and molecular mass. The molecule is
! translated so origin and center of mass is identical.

if (jPrint >= 99) call Prlist('Symmetry Distinct Nuclear Coordinates / Bohr',AtomLbl,size(Coor,2),Coor,3,size(Coor,2))
if (jPrint >= 99) call PrList('Symmetry Distinct Nuclear Forces / au',AtomLbl,size(Coor,2),Grd,3,size(Coor,2))
!                                                                      *
!***********************************************************************
!                                                                      *
mB_Tot = 0
mdB_Tot = 0
mq = 0
Force_dB = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Init_SlapAf
