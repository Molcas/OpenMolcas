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

subroutine MpProp(iReturn)

use MPProp_globals, only: AtBoMltPl, AtBoMltPlCopy, AtMltPl, AtPol, AtBoPol, BondMat, Cor, CordMltPl, EneV, Frac, iAtomType, &
                          iAtomPar, iAtPrTab, Labe, Method, MltPl, nAtomPBas, Qnuc
use Data_Structures, only: Allocate_DT, Deallocate_DT
use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: i, iCmp, iComp, iDum(1), iErr, iMltpl, iOff1, iOff2, iopt, iPol, iPrint, irc, iSmLbl, iSym, iWarn, iWFtype, &
                     Lu_, LuYou, nAtoms, nBas(8), nCenters, nComp, n_Int, nIrrep, nMltPl, nOcc, NOCOB, nOcOb_b, nOrbi, nPrim(8), &
                     nSize, nSum, nSym, nThrs, nTM, nVec, nVec_p
real(kind=wp) :: dLimmo(2), Thrs1, Thrs2, ThrsMul
character(len=6) :: FName
character(len=8) :: Label, MemLabel
character(len=80) :: VTitle
logical(kind=iwp) :: LNearestAtom, LAllCenters, AveOrb, Diffuse(3), LFirstRun, LLumOrb, Exists
integer(kind=iwp), allocatable :: ANr(:)
real(kind=wp), allocatable :: Atype(:), CenX(:), CenY(:), CenZ(:), Coor(:,:), D(:), D_p(:), D_p_b(:), EC(:,:), Ene(:,:), MP(:,:), &
                              Occ(:,:), Ocen(:,:), Ocen_b(:,:), Ocof(:), Ocof_b(:), TM(:), TP(:), Ttot(:,:), Ttot_Inv(:,:), &
                              Vec(:,:), Vec_p(:), Vec_p_b(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
! Set some zeroes
nPrim(:) = 0
nBas(:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Set some defaults
!iPol = 1
LNearestAtom = .true.
LAllCenters = .false.
AveOrb = .false.
LLumOrb = .false.
Diffuse(1) = .false.
Diffuse(2) = .false.
Diffuse(3) = .false.
dLimmo(1) = 0.65_wp
dLimmo(2) = Two
Thrs1 = 1.0e-5_wp
Thrs2 = 1.0e-4_wp
nThrs = 3
ThrsMul = 1.0e-2_wp
iPrint = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Get some coords, nuc. charges and labels

! Memory debugging
!call Setmem('TRACE=ON')
!
! IO check
!call fastio('TRACE=ON')
!
!call Bnnr()
call Get_cArray('Relax Method',Method,8)
if (Method == 'RHF-SCF') then
  iPol = 1
else if (Method == 'UHF-SCF') then
  iPol = 1
else if (Method == 'MBPT2') then
  iPol = 0
end if
! Runfile update
!call Get_nAtoms(nAtoms)
call Get_iScalar('Unique atoms',nAtoms)

call mma_allocate(Labe,nAtoms,label='Labe')
call mma_allocate(iAtomType,nAtoms,label='iAtomType')
call mma_allocate(iAtomPar,nAtoms,label='iAtomPar')
call mma_allocate(Cor,3,nAtoms,nAtoms,label='Cor')
call mma_allocate(Frac,nAtoms,nAtoms,label='Frac')
call mma_allocate(BondMat,nAtoms,nAtoms,label='BondMat')
call mma_allocate(nAtomPBas,nAtoms,label='nAtomPBas')
call mma_allocate(Qnuc,nAtoms,label='Qnuc')
iAtomPar(:) = 1
BondMat(:,:) = .false.

call mma_allocate(Coor,3,nAtoms,label='Coor')
call mma_allocate(Atype,nAtoms,label='Atype')
nSum = nAtoms*4

! Runfile update
!call Get_Coord(Coor,nAtoms)
call Get_dArray('Unique Coordinates',Coor,nAtoms*3)
! Runfile update
!call Get_AtomLabel(Labe,nAtoms)
call Get_cArray('Unique Atom Names',Labe,nAtoms*len(Labe))
! Runfile update
!call Get_Charge(Atype,nAtoms)
call Get_dArray('Nuclear charge',Atype,nAtoms)

! Runfile update
!call Get_Charge_Eff(Qnuc,nAtoms)
call Get_dArray('Effective nuclear Charge',Qnuc,nAtoms)

do i=1,nAtoms
  iAtomType(i) = int(Atype(i))
  Cor(:,i,i) = Coor(:,i)
end do
call Wr_Cord(nAtoms)
call mma_deallocate(Atype)
nCenters = nAtoms*(nAtoms+1)/2
!                                                                      *
!***********************************************************************
!                                                                      *
! Get information from input

call Get_Mpprop_input(nAtoms,iPol,LNearestAtom,LAllCenters,AveOrb,LLumOrb,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,iPrint)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read first the size of the contracted basis using the default file
! ONEINT and COMFILE

! Runfile update
!call Get_nSym(nSym)
call Get_iScalar('nSym',nSym)
!
! Runfile update
!call Get_nBas(nBas)
nIrrep = 1
call Get_iArray('nBas',nBas,1)

nVec = 0
nOcc = 0
do iSym=1,nSym
  nVec = nVec+nBas(iSym)**2
  nOcc = nOcc+nBas(iSym)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Read first the size of the primitive basis using ONEREL and COMREL

! Runfile update
!call Get_nBas(nPrim)
call Get_iArray('nBas_Prim',nPrim,1)

nSize = 0
do iSym=1,nSym
  nSize = nSize+nPrim(iSym)*(nPrim(iSym)+1)/2
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Do a dirty trick since the one-integral file is not explicitly opened.
call Put_iArray('nBas',nPrim,1)
!call OneBas('PRIM')
!                                                                      *
!***********************************************************************
!                                                                      *
! Read overlap, dipole moment and quadrupole moment integrals

! Count multipoles

iMltpl = 0
iComp = 1
do
  write(label,'(a,i2)') 'PLTPL ',iMltpl
  irc = -1
  iopt = ibset(0,sOpSiz)
  call iRdOne(irc,iopt,label,iComp,iDum,iSmLbl)
  if (irc /= 0) exit
  iMltpl = iMltpl+1
end do
nMltpl = max(0,iMltpl-1)

call Allocate_DT(MltPl,[0,nMltPl],'MltPl')
call mma_allocate(CordMltPl,[1,3],[0,nMltPl],label='CordMltPl')
nSum = nSum+3*(nMltPl+1)

! Do it general

do iMltpl=0,nMltPl
  write(label,'(a,i2)') 'PLTPL ',iMltpl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
  do iComp=1,nComp
    iCmp = iComp
    irc = -1
    iopt = ibset(0,sOpSiz)
    !EB call RdOne(irc,iopt,label,iComp,n_Int,iSmLbl)
    call iRdOne(irc,iopt,label,iCmp,iDum,iSmLbl)
    if (irc /= 0) then
      write(u6,'(2A)') 'MPProp: Error reading label=',label
      call Abend()
    end if
    if (iComp == 1) then
      n_Int = iDum(1)
      call mma_allocate(MltPl(iMltPl)%A,n_Int+4,nComp,label=MemLabel)
      nSum = nSum+nComp*(n_Int+4)
    else if (iDum(1) /= n_Int) then
      write(u6,'(2A)') 'MPProp: Error reading iComp /= 1 label=',label
      call Abend()
    end if
    if (n_Int == 0) then
      write(u6,'(2A)') 'MPProp: Error reading n_Int=0 label=',label
      call Abend()
    end if
    irc = -1
    iopt = 0
    call RdOne(irc,iopt,label,iCmp,MltPl(iMltpl)%A(:,iComp),iSmLbl)
    if (irc /= 0) then
      write(u6,'(2A)') '2 MPProp: Error reading ',label
      call Abend()
    end if
    !???????????????????????
    call CmpInt(MltPl(iMltpl)%A(:,iComp),n_Int,nPrim,nIrrep,iSmLbl)
    do i=1,3
      CordMltPl(i,iMltpl) = MltPl(iMltpl)%A(n_Int+i,1)
    end do
  end do
end do

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate Memory For Multipoles on Atoms + Atoms and Bonds

call Allocate_DT(AtMltPl,[0,nMltPl],'AtMltPl')
call Allocate_DT(AtBoMltPl,[0,nMltPl],'AtBoMltPl')
call Allocate_DT(AtBoMltPlCopy,[0,nMltPl],'AtBoMltPlCopy')
do iMltpl=0,nMltPl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
  call mma_allocate(AtMltPl(iMltpl)%A,nComp,nCenters,label=MemLabel)
  write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
  call mma_allocate(AtBoMltPl(iMltpl)%A,nComp,nCenters,label=MemLabel)
  write(MemLabel,'(A5,i3.3)') 'MtPCp',iMltpl
  call mma_allocate(AtBoMltPlCopy(iMltpl)%A,nComp,nCenters,label=MemLabel)
  nSum = nSum+3*nComp*nCenters
  AtMltPl(iMltpl)%A(:,:) = Zero
  AtBoMltPl(iMltpl)%A(:,:) = Zero
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the P-matrix

Label = 'P_matrix'
irc = -1
iopt = ibset(0,sOpSiz)
iComp = 1
!EB call RdOne(irc,iopt,label,iComp,n_Int,iSmLbl)
call iRdOne(irc,iopt,label,iComp,iDum,iSmLbl)
n_Int = iDum(1)
if (irc /= 0) then
  write(u6,'(2A)') 'MPProp: Error getting length of ',label
  write(u6,*) 'Length of the vector',n_Int,iSmLbl
  write(u6,*) 'irc=',irc
  call Abend()
end if
if (n_Int /= nSize) then
  write(u6,*) 'MPProp: n_Int.ne.nSize'
  write(u6,*) 'n_Int=',n_Int
  write(u6,*) 'nSize=',nSize
  call Abend()
end if
irc = -1
iopt = 0

call mma_allocate(CenX,n_Int+4,label='CenX')
call mma_allocate(CenY,n_Int+4,label='CenY')
call mma_allocate(CenZ,n_Int+4,label='CenZ')
nSum = nSum+3*(n_Int+4)

iComp = 1
call RdOne(irc,iopt,label,iComp,CenX,iSmLbl)
iComp = 2
call RdOne(irc,iopt,label,iComp,CenY,iSmLbl)
iComp = 3
call RdOne(irc,iopt,label,iComp,CenZ,iSmLbl)
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore

call Put_iArray('nBas',nBas,1)
!call OneBas('CONT')
!                                                                      *
!***********************************************************************
!                                                                      *
! READ TRANSFORMATION MATRIX  nPrim(i)*nBas(i)

! The TM vectors holds the contraction coefficients of each
! basisset in the primitive base. This means that the columns in the
! TM are composed of the coefficients for each atomic basefunction
! and there are nBas(iSym) of them. The coefficients found in each
! column are exactly the same as in the basisset shifted to the
! right position in the matrix column. For a hydrogen molecule
! calculated with 2s function on each atom the first column will
! look like (c11,c21,...,0,0,...) where the zeros are to delete the
! primitive gaussians on the "second" hydrogen atom.
! The second column looks like (c12,c22,...,0,0,...) and the third
! (0,0,...,c11,c21,...) where in the third column the coefficients
! are shifted to zero out the "first" hydrogen atom.

nTM = 0
do iSym=1,nSym
  nTM = nTM+nBas(iSym)*nPrim(iSym)
end do
call mma_allocate(TM,nTM,label='TM')
nSum = nSum+nTM
! Runfile update
!call Get_TPC(TM,nTM)
call Get_dArray('NEMO TPC',TM,nTM)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the CMO's and occupation numbers

! CMO's stored as nBas(iSym)*nSO(iSym)

NOCOB = 0
nOcOb_b = 0
if (LLumOrb) then
  if (Method == 'UHF-SCF') then
    call mma_allocate(Vec,nVec,2,label='Vec')
    call mma_allocate(Occ,nOcc,2,label='Occ')
    call mma_allocate(Ene,nOcc,2,label='Ene')
    nSum = nSum+4*nOcc+2*nVec
  else
    call mma_allocate(Vec,nVec,1,label='Vec')
    call mma_allocate(Occ,nOcc,1,label='Occ')
    call mma_allocate(Ene,nOcc,2,label='Ene') ! because Ene(:,2) is later used for all cases
    ! The 3* is for technical reason
    nSum = nSum+3*nOcc+nVec
  end if
  Lu_ = 11
  FName = 'INPORB'
  if (Method == 'UHF-SCF') then
    call RdVec_(FName,Lu_,'COE',1,nSym,nBas,nBas,Vec(:,1),Vec(:,2),Occ(:,1),Occ(:,2),Ene(:,1),Ene(:,2),iDum,VTitle,iWarn,iErr, &
                iWFtype)
  else
    call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Vec,Occ,Ene,iDum,VTitle,iWarn,iErr)
  end if
  if (index(VTitle,'IVO') /= 0) then
    write(u6,*) ' MpProp not implemented for IVO orbitals!'
    call Abend()
  end if
  !if (iPol == 2) then
  !  call LauraPol()
  !end if
  do i=1,nOcc
    if (Occ(i,1) /= Zero) then
      nOcOb = nOcOb+1
    end if
  end do
  if (Method == 'UHF-SCF') then
    do i=1,nOcc
      if (Occ(i,2) /= Zero) then
        nOcOb_b = nOcOb_b+1
      end if
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Project the MO's on the primitive basis
  !
  nVec_p = 0
  do iSym=1,nSym
    nVec_p = nVec_p+nPrim(iSym)*nBas(iSym)
  end do
  call mma_allocate(Vec_p,nVec_p,label='Vec_p')
  nSum = nSum+nVec_p
  iOff1 = 1
  iOff2 = 1
  do iSym=1,nSym
    if (nPrim(iSym) > 0) then
      call DGEMM_('N','N',nPrim(iSym),nBas(iSym),nBas(iSym),One,TM(iOff2),nPrim(iSym),Vec(iOff1,1),nBas(iSym),Zero,Vec_p(iOff2), &
                  nPrim(iSym))
      iOff1 = iOff1+nBas(iSym)**2
      iOff2 = iOff2+nPrim(iSym)*nBas(iSym)
    end if
  end do

  if (Method == 'UHF-SCF') then
    call mma_allocate(Vec_p_b,nVec_p,label='Vec_p_b')
    nSum = nSum+nVec_p
    iOff1 = 1
    iOff2 = 1
    do iSym=1,nSym
      if (nPrim(iSym) > 0) then
        call DGEMM_('N','N',nPrim(iSym),nBas(iSym),nBas(iSym),One,TM(iOff2),nPrim(iSym),Vec(iOff1,2),nBas(iSym),Zero, &
                    Vec_p_b(iOff2),nPrim(iSym))
        iOff1 = iOff1+nBas(iSym)**2
        iOff2 = iOff2+nPrim(iSym)*nBas(iSym)
      end if
    end do
  end if

  call mma_allocate(Ocof,nVec_p,label='Ocof')
  nSum = nSum+nVec_p
  call Get_OCOF(nPrim(1),nBas(1),Vec_p,nVec_p,Ocof)
  if (Method == 'UHF-SCF') then
    call mma_allocate(Ocof_b,nVec_p,label='Ocofb')
    nSum = nSum+nVec_p
    call Get_OCOF(nPrim(1),nBas(1),Vec_p_b,nVec_p,Ocof_b)
  end if

  call mma_allocate(D_p,nPrim(1)*(nPrim(1)+1)/2,label='D_p')
  nSum = nSum+size(D_p)
  call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),D_p,nOcOb,Occ,Ocof)
  if (Method == 'UHF-SCF') then
    call mma_allocate(D_p_b,nPrim(1)*(nPrim(1)+1)/2,label='D_p_b')
    nSum = nSum+size(D_p_b)
    call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),D_p_b,nOcOb_b,Occ(:,2),Ocof_b)
  end if
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the densities are used for expansion then
  ! Project the densities on the primitive basis

  ! Get the orbital energy and allocate double memory due to UHF calculations
  call mma_allocate(Ene,nOcc,2,label='OrbE')
  nsum = nSum+2*nOcc
  if (Method == 'RHF-SCF') then
    call Get_OrbE(Ene,nOcc)
  else
    Ene(:,:) = Zero
  end if
  call mma_allocate(D,nBas(1)*(nBas(1)+1)/2,label='D1ao')
  call Get_Density_Matrix_mpprop(D,nBas(1),nSym)
  write(u6,*) 'No polarizability will be calculated'
  iPol = 0
  call mma_allocate(D_p,nPrim(1)*(nPrim(1)+1)/2,label='D_p')
  nSum = nSum+size(D_p)
  call Get_Prim_Density_Matrix(D,nBas(1),D_p,nPrim(1),TM)
  call mma_deallocate(D)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(iAtPrTab,nPrim(1),nAtoms,label='iAtPrTab')
call Get_Prim_Atom_Tab(nAtoms,nPrim(1),Coor,CenX,CenY,CenZ)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the SCF or RASSCF energy
! Runfile update
!call Get_Energy(EneV)
call Get_dScalar('Last energy',EneV)

write(u6,*)
write(u6,'(a,f16.8)') ' Total SCF energy ',EneV
write(u6,*)
nOrbi = nBas(1)

! If multipoles and polarizability of the bonds should go to
! the atoms they belong to or to the nearest atom
if (.not. LNearestAtom) then
  write(u6,*)
  write(u6,*) ' I will not move bonds and polarizabilities'
  write(u6,*) ' to the nearest atom, just to the bonding pair'
  write(u6,*)
else
  write(u6,*)
  write(u6,*) ' I WILL move bonds and polarizabilities'
  write(u6,*) ' to the nearest atom'
  write(u6,*)
end if

! Get multipole properties
LFirstRun = .true.

call Get_MpProp(nPrim(1),nAtoms,nMltPl,D_p,CenX,CenY,CenZ,LNearestAtom,LFirstRun,LLumOrb)
! nOcOb,nMltPl,Occ,nOcc,Ocof,

if (Method == 'UHF-SCF') then
  LFirstRun = .false.
  call Get_MpProp(nPrim(1),nAtoms,nMltPl,D_p_b,CenX,CenY,CenZ,LNearestAtom,LFirstRun,LLumOrb)
  ! nOcOb_b,nMltPl,Occ(:,2),nOcc,Ocof_b,
  LFirstRun = .true.
end if

! If the esteemed user wants to obtain diffuse local distributions,
! then proceed here. Much of the code is common with LoProp, hence
! there is first a call to a routine to order some quantities from
! MpProp in the same say as in LoProp, then Mother Goose is called.

if (Diffuse(1)) then
  call mma_allocate(ANr,nAtoms,label='ANr')
  call mma_allocate(Ttot,nBas(1),nBas(1),label='T')
  call mma_allocate(Ttot_Inv,nBas(1),nBas(1),label='Tinv')
  call mma_allocate(MP,nCenters,(nMltPl+1)*(nMltPl+2)*(nMltPl+3)/6,label='MultMom')
  call mma_allocate(EC,3,nCenters,label='ExpCent')
  call StoreMpAsLop(nAtoms,ANr,nBas(1),Ttot,Ttot_Inv,MP,nMltPl,EC)
  call mma_allocate(TP,nAtoms,label='ToPoint')
  call CoreToPoint(nAtoms,MP,TP)
  LuYou = IsFreeUnit(81)
  call OpnFl('DIFFPR',LuYou,Exists)
  call Diff_MotherGoose(Diffuse,nAtoms,nBas(1),MP,nCenters,EC,ANr,Ttot,Ttot_Inv,nMltPl,TP,dLimmo,Thrs1,Thrs2,nThrs,iPrint,ThrsMul, &
                        LuYou)
  close(LuYou)
  call mma_deallocate(ANr)
  call mma_deallocate(Ttot)
  call mma_deallocate(Ttot_Inv)
  call mma_deallocate(MP)
  call mma_deallocate(EC)
  call mma_deallocate(TP)
end if

! End of Diffuse.

if (LLumOrb) then
  ! Get center of charge for each molecular orbital
  call mma_allocate(Ocen,3,nOrbi,label='Ocen')
  nSum = nSum+3*nOrbi
  if (Method == 'UHF-SCF') then
    call mma_allocate(Ocen_b,3,nOrbi,label='Ocen_b')
    nSum = nSum+3*nOrbi
  end if
end if
! Get polarizabillities if iPol
call mma_allocate(AtPol,6,nAtoms,label='AtPol')
call mma_allocate(AtBoPol,6,nCenters,label='AtBoPol')
nSum = nSum+6*(nCenters+nAtoms)
AtPol(:,:) = Zero
AtBoPol(:,:) = Zero
if (iPol > 0) then
  !EB call Get_OrbCen(nPrim(1),nBas(1),NORBI,MltPl(0)%A(:,:))
  call Get_OrbCen(nPrim(1),NORBI,MltPl(0)%A(:,1),Ocen,CenX,CenY,CenZ,Ocof)
  if (Method == 'UHF-SCF') call Get_OrbCen(nPrim(1),NORBI,MltPl(0)%A(:,1),Ocen_b,CenX,CenY,CenZ,Ocof_b)
  if (iPol == 1) then
    if (nOcOb < nOcc) then
      call Get_Polar(nPrim(1),nBas(1),nAtoms,nOcOb,Ene,nOcc,Ocof,Ocen,LNearestAtom,LFirstRun)
      !EB  Ene,Occ,nOcc,Ocof,Ocen,
      if (Method == 'UHF-SCF') then
        LFirstRun = .false.
        call Get_Polar(nPrim(1),nBas(1),nAtoms,nOcOb_b,Ene(:,2),nOcc,Ocof_b,Ocen_b,LNearestAtom,LFirstRun)
      end if
    else
      write(u6,*)
      write(u6,*) 'I will not do an analysis of the polarizability'
      write(u6,*) 'no. of occupied orb. is to large'
      write(u6,*)
      if (Method == 'UHF-SCF') then
        write(u6,*) ' nOcOb nOcOb_b nOcc ',nOcOb,nOcOb_b,nOcc
      else
        write(u6,*) ' nOcOb nOcc ',nOcOb,nOcc
      end if
      write(u6,*)
      iPol = 0
    end if
  else if (iPol == 2) then
    call LauraPol()
  end if
end if

! Write output, the properties
call Wr_Prop(nAtoms,nCenters,nBas(1),nMltPl,NOCOB,NOCOB_b,Ene,Ene(:,2),iPol,LAllCenters)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate

write(u6,*)
write(u6,*) 'Number of allocated real numbers',nsum
write(u6,'(a,f6.2,a)') ' That is ',real(nsum*RtoB,kind=wp)/(1024.0_wp**2),' MBytes'
write(u6,*)
call mma_deallocate(D_p)
!if (iPol == 0) then
call mma_deallocate(AtPol)
call mma_deallocate(AtBoPol)
!end if
if (LLumorb) then
  if (Method == 'UHF-SCF') then
    call mma_deallocate(D_p_b)
    call mma_deallocate(Ocen_b)
    call mma_deallocate(Ocof_b)
    call mma_deallocate(Vec_p_b)
  end if
  call mma_deallocate(Ocen)
  call mma_deallocate(Ocof)
  call mma_deallocate(Vec_p)
  call mma_deallocate(Occ)
  call mma_deallocate(Vec)
end if
call mma_deallocate(Ene)
call mma_deallocate(TM)
call mma_deallocate(CenX)
call mma_deallocate(CenY)
call mma_deallocate(CenZ)
call Deallocate_DT(AtMltPl)
call Deallocate_DT(AtBoMltPl)
call Deallocate_DT(AtBoMltPlCopy)
call Deallocate_DT(MltPl)
call mma_deallocate(CordMltPl)
call mma_deallocate(Coor)
call mma_deallocate(Labe)
call mma_deallocate(iAtomType)
call mma_deallocate(iAtomPar)
call mma_deallocate(Cor)
call mma_deallocate(Frac)
call mma_deallocate(BondMat)
call mma_deallocate(nAtomPBas)
call mma_deallocate(Qnuc)
call mma_deallocate(iAtPrTab)
!                                                                      *
!***********************************************************************
!                                                                      *
iReturn = 0

return

end subroutine MpProp
