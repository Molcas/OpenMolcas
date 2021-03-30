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

use MPProp_globals, only: AtPol, AtBoPol, BondMat, Cor, CordMltPl, EneV, Frac, iAtomType, iAtomPar, iAtPrTab, Labe, Method, &
                          nAtomPBas, Qnuc
use MPProp_globals, only: AtBoMltPl, AtBoMltPlCopy, AtMltPl, MltPl
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
#include "LenIn.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iAtype, iCenX, iCenY, iCenZ, iComp, iDum(1), iEne, iErr, iMltpl, iOcen, iOcen_b, iOcof, iOcof_b, iOff1, &
                     iOff2, iOff3, iopt, ip_ANr, ip_Coor, ip_D, ip_D_p, ip_D_p_b, ip_EC, ip_Ene, ip_Occ, ip_TM, ip_Ttot, &
                     ip_Ttot_Inv, ip_Vec, ip_Vec_p, ip_vec_p_b, ipMP, iPol, iPrint, irc, iSmLbl,iSym, iTP, iWarn, iWFtype, Lu_, &
                     LuYou, nAtoms, nBas(8), nCenters, nComp, n_Int, nIrrep, nMltPl, nOcc, NOCOB, nOcOb_b, nOrbi, nPrim(8), nSize, &
                     nSize1, nSize2, nSum, nSym, nThrs, nTM, nVec, nVec_p
real(kind=wp) :: dLimmo(2), Thrs1, Thrs2, ThrsMul
character(len=6) :: FName
character(len=8) :: Label, MemLabel
character(len=80) :: VTitle
logical(kind=iwp) :: LNearestAtom, LAllCenters, AveOrb, Diffuse(3), LFirstRun, LLumOrb, Exists
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
! Set some zeroes
do i=1,8
  nPrim(i) = 0
  nBas(i) = 0
end do
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
elseif (Method == 'UHF-SCF') then
  iPol = 1
elseif (Method == 'MBPT2') then
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

call GetMem('Coord','Allo','Real',ip_Coor,3*nAtoms)
call GetMem('Atype','Allo','Real',iAtype,nAtoms)
! Runfile update
!call Get_Coord(Work(ip_Coor),nAtoms)
call Get_dArray('Unique Coordinates',Work(ip_Coor),nAtoms*3)
! Runfile update
!call Get_AtomLabel(Labe,nAtoms)
call Get_cArray('Unique Atom Names',Labe,LenIn*nAtoms)
! Runfile update
!call Get_Charge(Work(iAtype),nAtoms)
call Get_dArray('Nuclear charge',Work(iAtype),nAtoms)

! Runfile update
!call Get_Charge_Eff(Work(iQnuc),nAtoms)
call Get_dArray('Effective nuclear Charge',Qnuc,nAtoms)

do i=1,nAtoms
  iAtomType(i) = int(Work(iAtype+i-1))
  Cor(1,i,i) = Work(ip_Coor+(i-1)*3)
  Cor(2,i,i) = Work(ip_Coor+(i-1)*3+1)
  Cor(3,i,i) = Work(ip_Coor+(i-1)*3+2)
end do
call Wr_Cord(nAtoms)
call GetMem('Atype','Free','Real',iAtype,nAtoms)
nCenters = nAtoms*(nAtoms+1)/2
nSum = nAtoms*6
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
do
  write(label,'(a,i2)') 'PLTPL ',iMltpl
  irc = -1
  iopt = 1
  call iRdOne(irc,iopt,label,1,iDum,iSmLbl)
  if (irc /= 0) exit
  iMltpl = iMltpl+1
end do
nMltpl = max(0,iMltpl-1)

allocate(MltPl(0:nMltPl))!,label='MltPl')
call mma_allocate(CordMltPl,[1,3],[0,nMltPl],label='CordMltPl')

! Do it general

do iMltpl=0,nMltPl
  write(label,'(a,i2)') 'PLTPL ',iMltpl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
  nSum = nSum+nComp
  do iComp=1,nComp
    irc = -1
    iopt = 1
    !EB call RdOne(irc,iopt,label,iComp,n_Int,iSmLbl)
    call iRdOne(irc,iopt,label,iComp,iDum,iSmLbl)
    if (irc /= 0) then
      write(u6,'(2A)') 'MPProp: Error reading label=',label
      call Abend()
    end if
    if (iComp == 1) then
      n_Int = iDum(1)
      call mma_allocate(MltPl(iMltPl)%M,n_Int+4,nComp,label=MemLabel)
    else if (iDum(1) /= n_Int) then
      write(u6,'(2A)') 'MPProp: Error reading iComp /= 1 label=',label
      call Abend()
    end if
    if (n_Int == 0) then
      write(u6,'(2A)') 'MPProp: Error reading n_Int=0 label=',label
      call Abend()
    end if
    nSum = nSum+n_Int+4
    irc = -1
    iopt = 0
    call RdOne(irc,iopt,label,iComp,MltPl(iMltpl)%M(:,iComp),iSmLbl)
    if (irc /= 0) then
      write(u6,'(2A)') '2 MPProp: Error reading ',label
      call Abend()
    end if
    !???????????????????????
    call CmpInt(MltPl(iMltpl)%M(:,iComp),n_Int,nPrim,nIrrep,iSmLbl)
    do i=1,3
      CordMltPl(i,iMltpl) = MltPl(iMltpl)%M(n_Int+i,1)
    end do
  end do
end do

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate Memory For Multipoles on Atoms + Atoms and Bonds

allocate(AtMltPl(0:nMltPl))
allocate(AtBoMltPl(0:nMltPl))
allocate(AtBoMltPlCopy(0:nMltPl))
do iMltpl=0,nMltPl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
  call mma_allocate(AtMltPl(iMltpl)%M,nComp,nCenters,label=MemLabel)
  write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
  call mma_allocate(AtBoMltPl(iMltpl)%M,nComp,nCenters,label=MemLabel)
  write(MemLabel,'(A5,i3.3)') 'MtPCp',iMltpl
  call mma_allocate(AtBoMltPlCopy(iMltpl)%M,nComp,nCenters,label=MemLabel)
  nSum = nSum+nComp*(nAtoms+nCenters)
  AtMltPl(iMltpl)%M(:,:) = Zero
  AtBoMltPl(iMltpl)%M(:,:) = Zero
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the P-matrix

Label = 'P_matrix'
irc = -1
iopt = 1
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

call GetMem('CenX','Allo','Real',iCenX,n_Int+4)
call GetMem('CenY','Allo','Real',iCenY,n_Int+4)
call GetMem('CenZ','Allo','Real',iCenZ,n_Int+4)
nSum = nSum+n_Int*3

iComp = 1
call RdOne(irc,iopt,label,iComp,Work(iCenX),iSmLbl)
iComp = 2
call RdOne(irc,iopt,label,iComp,Work(iCenY),iSmLbl)
iComp = 3
call RdOne(irc,iopt,label,iComp,Work(iCenZ),iSmLbl)
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
call GetMem('TM','Allo','Real',ip_TM,nTM)
nSum = nSum+nTM
! Runfile update
!call Get_TPC(Work(ip_TM),nTM)
call Get_dArray('NEMO TPC',Work(ip_TM),nTM)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the CMO's and occupation numbers

! CMO's stored as nBas(iSym)*nSO(iSym)

NOCOB = 0
nOcOb_b = 0
if (LLumOrb) then
  if (Method == 'UHF-SCF') then
    call GetMem('Vec','Allo','Real',ip_Vec,2*nVec)
    call GetMem('Occ','Allo','Real',ip_Occ,2*nOcc)
    call GetMem('Ene','Allo','Real',ip_Ene,2*nOcc)
    nSum = nSum+4*nOcc+2*nVec
  else
    call GetMem('Vec','Allo','Real',ip_Vec,nVec)
    call GetMem('Occ','Allo','Real',ip_Occ,nOcc)
    call GetMem('Ene','Allo','Real',ip_Ene,2*nOcc)
    ! The 2* is for technical reason
    nSum = nSum+2*nOcc+nVec
  end if
  Lu_ = 11
  FName = 'INPORB'
  if (Method == 'UHF-SCF') then
    call RdVec_(FName,Lu_,'COE',1,nSym,nBas,nBas,Work(ip_Vec),Work(ip_Vec+nVec),Work(ip_Occ),Work(ip_Occ+nOcc),Work(ip_Ene), &
                Work(ip_Ene+nOcc),iDum,VTitle,iWarn,iErr,iWFtype)
  else
    call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Work(ip_Vec),Work(ip_Occ),Work(ip_Ene),iDum,VTitle,iWarn,iErr)
  end if
  if (index(VTitle,'IVO') /= 0) then
    write(u6,*) ' MpProp not implemented for IVO orbitals!'
    call Abend
  end if
  !if (iPol == 2) then
  !  call LauraPol()
  !end if
  if (Method == 'UHF-SCF') then
    do i=0,nOcc-1
      if (Work(ip_Occ+i) /= Zero) then
        nOcOb = nOcOb+1
      end if
    end do
    do i=nOcc,2*nOcc-1
      if (Work(ip_Occ+i) /= Zero) then
        nOcOb_b = nOcOb_b+1
      end if
    end do
  else
    do i=0,nOcc-1
      if (Work(ip_Occ+i) /= Zero) then
        nOcOb = nOcOb+1
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
  call GetMem('Vec_p','Allo','Real',ip_Vec_p,nVec_p)
  nSum = nSum+nVec_p
  iOff1 = ip_Vec
  iOff2 = ip_TM
  iOff3 = ip_Vec_p
  do iSym=1,nSym
    if (nPrim(iSym) > 0) then
      call DGEMM_('N','N',nPrim(iSym),nBas(iSym),nBas(iSym),One,Work(iOff2),nPrim(iSym),Work(iOff1),nBas(iSym),Zero,Work(iOff3), &
                  nPrim(iSym))
      iOff1 = iOff1+nBas(iSym)**2
      iOff2 = iOff2+nPrim(iSym)*nBas(iSym)
      iOff3 = iOff3+nPrim(iSym)*nBas(iSym)
    end if
  end do

  if (Method == 'UHF-SCF') then
    call GetMem('Vec_p_b','Allo','Real',ip_Vec_p_b,nVec_p)
    nSum = nSum+nVec_p
    iOff1 = ip_Vec+nVec
    iOff2 = ip_TM
    iOff3 = ip_Vec_p_b
    do iSym=1,nSym
      if (nPrim(iSym) > 0) then
        call DGEMM_('N','N',nPrim(iSym),nBas(iSym),nBas(iSym),One,Work(iOff2),nPrim(iSym),Work(iOff1),nBas(iSym),Zero,Work(iOff3), &
                    nPrim(iSym))
        iOff1 = iOff1+nBas(iSym)**2
        iOff2 = iOff2+nPrim(iSym)*nBas(iSym)
        iOff3 = iOff3+nPrim(iSym)*nBas(iSym)
      end if
    end do
  end if

  call GetMem('Ocof','Allo','Real',iOcof,nVec_p)
  nSum = nSum+nVec_p
  call Get_OCOF(nPrim(1),nBas(1),Work(ip_Vec_p),nVec_p,Work(iOcof))
  if (Method == 'UHF-SCF') then
    call GetMem('Ocofb','Allo','Real',iOcof_b,nVec_p)
    nSum = nSum+nVec_p
    call Get_OCOF(nPrim(1),nBas(1),Work(ip_Vec_p_b),nVec_p,Work(iOcof_b))
  end if

  call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim(1)*(nPrim(1)+1)/2)
  call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p,nOcOb,Work(ip_Occ),Work(iOcof))
  if (Method == 'UHF-SCF') then
    call GetMem('D_p','ALLO','REAL',ip_D_p_b,nPrim(1)*(nPrim(1)+1)/2)
    call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p_b,nOcOb_b,Work(ip_Occ+nOcc),Work(iOcof_b))
  end if
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the densities are used for expansion then
  ! Project the densities on the primitive basis

  ! Get the orbital energy and allocate double memory due to UHF calculations
  call GetMem('OrbE','Allo','Real',ip_Ene,2*nOcc)
  if (Method == 'RHF-SCF') then
    call Get_OrbE_mpprop(ip_Ene,nOcc)
  else
    do iEne=1,2*nOcc
      Work(ip_Ene+iEne-1) = Zero
    end do
  end if
  call GetMem('D1ao','Allo','Real',ip_D,nBas(1)*(nBas(1)+1)/2)
  call Get_Density_Matrix_mpprop(ip_D,nBas(1),nSym)
  write(u6,*) 'No polarizability will be calculated'
  iPol = 0
  call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim(1)*(nPrim(1)+1)/2)
  call Get_Prim_Density_Matrix(ip_D,nBas(1),ip_D_p,nPrim(1),Work(ip_TM))
  call Free_Work(ip_D)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(iAtPrTab,nPrim(1),nAtoms,label='iAtPrTab')
call Get_Prim_Atom_Tab(nAtoms,nPrim(1),Work(ip_Coor),Work(iCenX),Work(iCenY),Work(iCenZ))
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

call Get_MpProp(nPrim(1),nBas(1),nAtoms,nCenters,nMltPl,ip_D_p,Work(iCenX),Work(iCenY),Work(iCenZ),LNearestAtom,LFirstRun,LLumOrb)
! nOcOb,nMltPl,Work(ip_Occ),nOcc,Work(iOcof),

if (Method == 'UHF-SCF') then
  LFirstRun = .false.
  call Get_MpProp(nPrim(1),nBas(1),nAtoms,nCenters,nMltPl,ip_D_p_b,Work(iCenX),Work(iCenY),Work(iCenZ),LNearestAtom,LFirstRun, &
                  LLumOrb)
  ! nOcOb_b,nMltPl,Work(ip_Occ+nOcc),nOcc,Work(iOcof_b),
  LFirstRun = .true.
end if

! If the esteemed user wants to obtain diffuse local distributions,
! then proceed here. Much of the code is common with LoProp, hence
! there is first a call to a routine to order some quantaties from
! MpProp in the same say as in LoProp, then Mother Goose is called.

if (Diffuse(1)) then
  call Allocate_iWork(ip_ANr,nAtoms)
  call GetMem('T','Allo','Real',ip_Ttot,nBas(1)**2)
  call GetMem('Tinv','Allo','Real',ip_Ttot_Inv,nBas(1)**2)
  call GetMem('ExpCent','Allo','Real',ip_EC,3*nAtoms*(nAtoms+1)/2)
  call GetMem('MultMom','Allo','Real',ipMP,nAtoms*(nAtoms+1)/2*(nMltPl*(nMltPl**2+6*nMltPl+11)+6)/6)
  call StoreMpAsLop(nAtoms,ip_ANr,nBas(1),ip_Ttot,ip_Ttot_Inv,ipMP,nMltPl,ip_EC)
  call GetMem('ToPoint','Allo','Real',iTP,nAtoms)
  call CoreToPoint(nAtoms,Work(ipMP),Work(iTP))
  LuYou = IsFreeUnit(81)
  call OpnFl('DIFFPR',LuYou,Exists)
  call Diff_MotherGoose(Diffuse,nAtoms,nBas(1),ipMP,ip_Coor,nCenters,ip_EC,ip_ANr,ip_Ttot,ip_Ttot_Inv,nMltPl,iTP,dLimmo,Thrs1, &
                        Thrs2,nThrs,iPrint,ThrsMul,LuYou)
  close(LuYou)
  call Free_iWork(ip_ANr)
  call GetMem('T','Free','Real',ip_Ttot,nBas(1)**2)
  call GetMem('Tinv','Free','Real',ip_Ttot_Inv,nBas(1)**2)
  call GetMem('ExpCent','Free','Real',ip_EC,3*nAtoms*(nAtoms+1)/2)
  nSize1 = nAtoms*(nAtoms+1)/2
  nSize2 = (nMltPl*(nMltPl**2+6*nMltPl+11)+6)/6
  call GetMem('MultMom','Free','Real',ipMP,nSize1*nSize2)
  call GetMem('ToPoint','Free','Real',iTP,nAtoms)
end if

! End of Diffuse.


if (LLumOrb) then
  ! Get center of charge for each molecular orbital
  call GetMem('Ocen','Allo','Real',iOcen,3*nOrbi)
  if (Method == 'UHF-SCF') call GetMem('Ocen_b','Allo','Real',iOcen_b,3*nOrbi)
  nSum = nSum+3*nOrbi*2
end if
! Get polarizabillities if iPol
call mma_allocate(AtPol,6,nAtoms,label='AtPol')
call mma_allocate(AtBoPol,6,nCenters,label='AtBoPol')
nSum = nSum+6*(nCenters+nAtoms)
AtPol(:,:) = Zero
AtBoPol(:,:) = Zero
if (iPol > 0) then
  !EB call Get_OrbCen(nPrim(1),nBas(1),NORBI,MltPl(0)%M(:,:))
  call Get_OrbCen(nPrim(1),NORBI,MltPl(0)%M(:,1),Work(iOcen),Work(iCenX),Work(iCenY),Work(iCenZ),Work(iOcof))
  if (Method == 'UHF-SCF') call Get_OrbCen(nPrim(1),NORBI,MltPl(0)%M(:,1),Work(iOcen_b),Work(iCenX),Work(iCenY),Work(iCenZ), &
                                           Work(iOcof_b))
  if (iPol == 1) then
    if (nOcOb < nOcc) then
      call Get_Polar(nPrim(1),nBas(1),nAtoms,nCenters,nOcOb,Work(ip_Ene),nOcc,Work(iOcof),Work(iOcen),LNearestAtom,LFirstRun)
      !EB  Work(ip_Ene),Work(ip_Occ),nOcc,Work(iOcof),Work(iOcen),
      if (Method == 'UHF-SCF') then
        LFirstRun = .false.
        call Get_Polar(nPrim(1),nBas(1),nAtoms,nCenters,nOcOb_b,Work(ip_Ene+nOcc),nOcc,Work(iOcof_b),Work(iOcen_b),LNearestAtom, &
                       LFirstRun)
      end if
    else
      write(u6,*)
      write(u6,*) 'I will not do an analyze of the polarizability'
      write(u6,*) 'no of occupied orb. is to large'
      write(u6,*)
      if (Method == 'UHF-SCF') then
        write(u6,*) ' nOcOb nOcOb_b nOcc ',nOcOb,nOcOb_b,nOcc
      else
        write(u6,*) ' nOcOb nOcc ',nOcOb,nOcc
      end if
      write(u6,*)
      iPol = 0
    end if
  elseif (iPol == 2) then
    call LauraPol()
  end if
end if

! Write output, the properties
call Wr_Prop(nAtoms,nCenters,nBas(1),nMltPl,NOCOB,NOCOB_b,Work(ip_Ene),Work(ip_Ene+nOcc),iPol,LAllCenters)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate Work

write(u6,*)
write(u6,*) 'Number of allocated real*8 words',nsum
write(u6,'(a,f6.2,a)') ' That is ',nsum*Eight/(1024.0_wp**2),' MBytes'
write(u6,*)
call Free_Work(ip_D_p)
!if (iPol == 0) then
call mma_deallocate(AtPol)
call mma_deallocate(AtBoPol)
!end if
if (LLumorb) then
  if (Method == 'UHF-SCF') then
    call Free_Work(ip_D_p_b)
    call GetMem('Ocen_b','Free','Real',iOcen_b,3*nOrbi)
    call GetMem('Ocof','Free','Real',iOcof_b,nVec_p)
    call GetMem('Vec_p','Free','Real',ip_Vec_p_b,nVec_p)
  end if
  call GetMem('Ocen','Free','Real',iOcen,3*nOrbi)
  call GetMem('Ocof','Free','Real',iOcof,nVec_p)
  call GetMem('Vec_p','Free','Real',ip_Vec_p,nVec_p)
  call GetMem('Occ','Free','Real',ip_Occ,nOcc)
  call GetMem('Vec','Free','Real',ip_Vec,nVec)
end if
call GetMem('Ene','Free','Real',ip_Ene,nOcc)
call GetMem('TM','Free','Real',ip_TM,nTM)
call GetMem('CenZ','Free','Real',iCenZ,n_Int+4)
call GetMem('CenY','Free','Real',iCenY,n_Int+4)
call GetMem('CenX','Free','Real',iCenX,n_Int+4)
do iMltpl=0,nMltPl
  call mma_deallocate(AtMltPl(iMltpl)%M)
  call mma_deallocate(AtBoMltPl(iMltpl)%M)
  call mma_deallocate(AtBoMltPlCopy(iMltpl)%M)
  call mma_deallocate(MltPl(iMltpl)%M)
end do
deallocate(AtMltPl)
deallocate(AtBoMltPl)
deallocate(AtBoMltPlCopy)
deallocate(MltPl)
call mma_deallocate(CordMltPl)
call GetMem('Coord','Free','Real',ip_Coor,3*nAtoms)
call GetMem('Coord','Check','Real',ip_Coor,3*nAtoms)
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
