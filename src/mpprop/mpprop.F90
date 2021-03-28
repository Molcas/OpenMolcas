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

use Constants, only: Zero, One, Two, Eight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
#include "MpData.fh"
#include "WrkSpc.fh"
#include "MolProp.fh"
integer(kind=iwp) :: i, iAtype, iCenX, iCenY, iCenZ, iComp, iDum(1), iEne, iErr, iMltpl, iOcen, iOcen_b, iOcof, iOcof_b, iOff1, &
                     iOff2, iOff3, iopt, ip_ANr, ip_Coor, ip_D, ip_D_p, ip_D_p_b, ip_EC, ip_Ene, ip_Occ, ip_TM, ip_Ttot, &
                     ip_Ttot_Inv, ip_Vec, ip_Vec_p, ip_vec_p_b, ipMP, iPol, iPrint, irc, iSmLbl,iSym, iTP, iWarn, iWFtype, j, Lu_, &
                     LuYou, nAtoms, nBas(8), nCenters, nComp, nDens, n_Int, nIrrep, nMltPl, nOcc, NOCOB, nOcOb_b, nOrbi, nPrim(8), &
                     nSize, nSize1, nSize2, nSum, nSym, nThrs, nTM, nVec, nVec_p
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
do i=1,mxAtomMP
  iAtomPar(i) = 1
  nub(i) = 0
  do j=1,mxAtomMP
    nbi(i,j) = 0
  end do
end do
do i=1,mxCen
  iBondPar(i) = 1
end do
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
if (nAtoms > mxAtomMP) then
  write(u6,'(A)') 'MPProp: Too many atoms'
  call Abend()
end if

call GetMem('Coord','Allo','Real',ip_Coor,3*nAtoms)
call GetMem('Atype','Allo','Real',iAtype,nAtoms)
! Runfile update
!call Get_Coord(Work(ip_Coor),nAtoms)
call Get_dArray('Unique Coordinates',Work(ip_Coor),nAtoms*3)
! Runfile update
!call Get_AtomLabel(Labe,nAtoms)
call Get_cArray('Unique Atom Names',Labe,LENIN*nAtoms)
! Runfile update
!call Get_Charge(Work(iAtype),nAtoms)
call Get_dArray('Nuclear charge',Work(iAtype),nAtoms)

call GetMem('Qnuc','Allo','Real',iQnuc,nAtoms)
! Runfile update
!call Get_Charge_Eff(Work(iQnuc),nAtoms)
call Get_dArray('Effective nuclear Charge',Work(iQnuc),nAtoms)

do i=1,nAtoms
  iAtomType(i) = int(Work(iAtype+i-1))
  COR(1,i,i) = Work(ip_Coor+(i-1)*3)
  COR(2,i,i) = Work(ip_Coor+(i-1)*3+1)
  COR(3,i,i) = Work(ip_Coor+(i-1)*3+2)
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
! Read first the size of the primitiv basis using ONEREL and COMREL

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
! Do a dirty trick since the one-intergral file is not explicitly
! opened.
call Put_iArray('nBas',nPrim,1)
!call OneBas('PRIM')
!                                                                      *
!***********************************************************************
!                                                                      *
! Read overlap, dipole moment and quadrupole moment integrals

! Do it general

do iMltpl=0,mxMltPl
  write(label,'(a,i2)') 'PLTPL ',iMltpl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
  call GetMem(MemLabel,'Allo','Inte',iMltPlAd(iMltpl),nComp)
  nSum = nSum+nComp
  do iComp=1,nComp
    irc = -1
    iopt = 1
    !EB call RdOne(irc,iopt,label,iComp,n_Int,iSmLbl)
    call iRdOne(irc,iopt,label,iComp,iDum,iSmLbl)
    if (irc == 0) n_Int = iDum(1)
    if (irc /= 0) then
      if (iComp /= 1) then
        write(u6,'(2A)') 'MPProp: Error reading iComp.ne.0 label=',label
        call Abend()
      else
        call GetMem(MemLabel,'Free','Inte',iMltPlAd(iMltpl),nComp)
        nMltPl = iMltPl-1
        nSum = nSum-nComp
        go to 100
      end if
    end if
    if (n_Int /= 0) then
      write(MemLabel,'(i3.3,i5.5)') iMltpl,iComp
      call GetMem(MemLabel,'Allo','Real',iWork(iMltPlAd(iMltpl)+iComp-1),n_Int+4)
      nSum = nSum+n_Int+4
      irc = -1
      iopt = 0
      call RdOne(irc,iopt,label,iComp,Work(iWork(iMltPlAd(iMltpl)+iComp-1)),iSmLbl)
    else
      write(u6,'(2A)') 'MPProp: Error reading n_Int=0 label=',label
      call Abend()
    end if
    if (irc /= 0) then
      write(u6,'(2A)') '2 MPProp: Error reading ',label
      call Abend()
    end if
    !???????????????????????
    if (n_Int /= 0) call CmpInt(Work(iWork(iMltPlAd(iMltpl)+iComp-1)),n_Int,nPrim,nIrrep,iSmLbl)
    do i=1,3
      CordMltPl(i,iMltpl) = Work(iWork(iMltPlAd(iMltpl))+n_Int+i-1)
    end do
  end do
end do

100 continue

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate Memory For Multipoles on Atoms + Atoms and Bonds

do iMltpl=0,nMltPl
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
  call GetMem(MemLabel,'Allo','Real',iAtMltPlAd(iMltpl),nComp*nAtoms)
  write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
  call GetMem(MemLabel,'Allo','Real',iAtBoMltPlAd(iMltpl),nComp*nCenters)
  write(MemLabel,'(A5,i3.3)') 'MtPCp',iMltpl
  call GetMem(MemLabel,'Allo','Real',iAtBoMltPlAdCopy(iMltpl),nComp*nCenters)
  nSum = nSum+nComp*(nAtoms+nCenters)
  do i=1,nComp*nAtoms
    Work(iAtMltPlAd(iMltpl)+i-1) = Zero
  end do
  do i=1,nComp*nCenters
    Work(iAtBoMltPlAd(iMltpl)+i-1) = Zero
  end do
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

  call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p,nOcOb,Work(ip_Occ),Work(iOcof))
  if (Method == 'UHF-SCF') then
    call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p_b,nOcOb_b,Work(ip_Occ+nOcc),Work(iOcof_b))
  end if
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If the densities are used for expansion then
  ! Project the densities on the primitive basis

  ! Get the orbital energy and allocate dubble memory due to UHF calculations
  if (Method == 'RHF-SCF') then
    call Get_OrbE_mpprop(ip_Ene,nOcc)
  else
    call GetMem('OrbE','Allo','Real',ip_Ene,2*nOcc)
    do iEne=1,2*nOcc
      Work(ip_Ene+iEne-1) = Zero
    end do
  end if
  call Get_Density_Matrix_mpprop(ip_D,nDens,nBas(1),nSym)
  write(u6,*) 'No polarizability will be calculated'
  iPol = 0
  call Get_Prim_Density_Matrix(ip_D,nBas(1),ip_D_p,nPrim(1),Work(ip_TM))
  call Free_Work(ip_D)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
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
call GetMem('AtPol','Allo','Real',iAtPolAd,nAtoms*6)
call GetMem('AtBoPol','Allo','Real',iAtBoPolAd,nCenters*6)
nSum = nSum+6*(nCenters+nAtoms)
do i=0,nAtoms*6-1
  Work(iAtPolAd+i) = Zero
end do
do i=0,nCenters*6-1
  Work(iAtBoPolAd+i) = Zero
end do
if (iPol > 0) then
  !EB call Get_OrbCen(nPrim(1),nBas(1),NORBI,Work(iWork(iMltPlAd(0)))
  call Get_OrbCen(nPrim(1),NORBI,Work(iWork(iMltPlAd(0))),Work(iOcen),Work(iCenX),Work(iCenY),Work(iCenZ),Work(iOcof))
  if (Method == 'UHF-SCF') call Get_OrbCen(nPrim(1),NORBI,Work(iWork(iMltPlAd(0))),Work(iOcen_b),Work(iCenX),Work(iCenY), &
                                           Work(iCenZ),Work(iOcof_b))
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
call GetMem('AtBoPol','Free','Real',iAtBoPolAd,nCenters*6)
call GetMem('AtPol','Free','Real',iAtPolAd,nAtoms*6)
!end if
if (LLumorb) then
  if (Method == 'UHF-SCF') then
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
  nComp = (iMltpl+1)*(iMltpl+2)/2
  write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
  call GetMem(MemLabel,'Free','Real',iAtMltPlAd(iMltpl),nComp*nAtoms)
  write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
  call GetMem(MemLabel,'Free','Real',iAtBoMltPlAd(iMltpl),nComp*nCenters)
  call GetMem(MemLabel,'Free','Real',iAtBoMltPlAdCopy(iMltpl),nComp*nCenters)
  do iComp=1,nComp
    write(MemLabel,'(i3.3,i5.5)') iMltpl,iComp
    call GetMem(MemLabel,'Free','Real',iWork(iMltPlAd(iMltpl)+iComp-1),n_Int+4)
  end do
  write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
  call GetMem(MemLabel,'Free','Inte',iMltPlAd(iMltpl),nComp)
end do
call GetMem('Qnuc','Free','Real',iQnuc,nAtoms)
call GetMem('Coord','Free','Real',ip_Coor,3*nAtoms)
call GetMem('Coord','Check','Real',ip_Coor,3*nAtoms)
!                                                                      *
!***********************************************************************
!                                                                      *
iReturn = 0

return

end subroutine MpProp
