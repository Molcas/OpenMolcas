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

subroutine h1_espf(h1,RepNuc,nh1,First,Do_DFT)
! Driver for computing the ESPF one-electron modification
! of the core hamiltonian

use espf_global, only: MxExtPotComp
use Index_Functions, only: nTri_Elem1
use Basis_Info, only: nBas
use Data_Structures, only: Alloc2DArray_Type, Alloc4DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nh1
real(kind=wp), intent(inout) :: h1(nh1), RepNuc
logical(kind=iwp), intent(in) :: First, Do_DFT
#include "print.fh"
integer(kind=iwp) :: iAt, ibla, iGrdTyp, ii, iMlt, iMode, iPL, IPotFl, iQMchg, iRMax, ITkQMMM, jAt, MltOrd, nAtMM, natom, nAtQM, &
                     nChg, nGrdPt, nMult, nSym
real(kind=wp) :: DeltaR, RealDummy, rms2, rms3, rms4, sum1, sum2, sum3, sum4
logical(kind=iwp) :: DoDirect, DoGromacs, DoTinker, DynExtPot, Exists, lMorok, StandAlone, UpdateVMM
character(len=180) :: Line
character(len=10) :: ESPFKey
type(Alloc2DArray_Type) :: Grid
type(Alloc4DArray_Type) :: DGrid
integer(kind=iwp), allocatable :: IsMM(:)
real(kind=wp), allocatable :: B(:), Cord(:,:), Ext(:,:), Mltp(:), OldMltp(:), T(:,:), TT(:,:), TTT(:,:)
integer(kind=iwp), external :: iPrintLevel, IsFreeUnit
character(len=180), external :: Get_Ln

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)

RealDummy = Zero

! Recover some ESPF data

StandAlone = .false.
DoTinker = .false.
DoGromacs = .false.
lMorok = .false.
DoDirect = .false.
call F_Inquire('ESPF.DATA',Exists)
if (Exists) then
  IPotFl = IsFreeUnit(1)
  call Molcas_Open(IPotFl,'ESPF.DATA')
  do
    Line = Get_Ln(IPotFl)
    ESPFKey = Line(1:10)
    if (ESPFKey == 'MLTORD    ') then
      call Get_I1(2,MltOrd)
      ibla = 0
      do ii=0,MltOrd
        ibla = ibla+nTri_Elem1(ii)
      end do
      MltOrd = ibla
    else if (ESPFKey == 'IRMAX     ') then
      call Get_I1(2,iRMax)
    else if (ESPFKey == 'DELTAR    ') then
      call Get_F1(2,DeltaR)
    else if (ESPFKey == 'GRIDTYPE  ') then
      call Get_I1(2,iGrdTyp)
    else if (ESPFKey == 'TINKER    ') then
      DoTinker = .true.
    else if (ESPFKey == 'GROMACS   ') then
      DoGromacs = .true.
    else if (ESPFKey == 'LA_MOROK  ') then
      lMorok = .true.
    else if (ESPFKey == 'DIRECT    ') then
      DoDirect = .true.
    else if (ESPFKey == 'MULTIPOLE ') then
      call Get_I1(2,nMult)
      call mma_allocate(OldMltp,nMult,label='OldMltp')
      do iMlt=1,nMult,MltOrd
        Line = Get_Ln(IPotFl)
        call Get_I1(1,iAt)
        call Get_F(2,OldMltp(iMlt),MltOrd)
      end do
    else if (ESPFKey == 'ENDOFESPF ') then
      exit
    end if
  end do
  close(IPotFl)
else
  write(u6,*) 'No ESPF.DATA file. Abort'
  call Quit_OnUserError()
end if

! Is this a fully coupled QM/MM calculation?

DynExtPot = .false.
iMode = 0
if (DoTinker) then
  ITkQMMM = IsFreeUnit(30)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'FullCoupling') /= 0) then
      DynExtPot = .true.
      iMode = max(iMode,1)
    else if (index(Line,'MMPolar') /= 0) then
      DynExtPot = .true.
      iMode = max(iMode,2)
    end if
  end do
  close(ITkQMMM)
end if
if (.not. DynExtPot) then
  if (allocated(OldMltp)) call mma_deallocate(OldMltp)
  return
end if

call Get_iScalar('Unique atoms',natom)
call mma_allocate(Cord,3,natom,label='AtomCoord')
call Get_dArray('Unique Coordinates',Cord,3*natom)
call mma_allocate(IsMM,natom,label='IsMM for atoms')
call mma_allocate(Ext,MxExtPotComp,natom,label='ExtPot')
Ext(:,:) = Zero
call MMCount(natom,nAtMM,IsMM)
nAtQM = natom-nAtMM
nMult = MltOrd*nAtQM

! Compute the grid around the molecule

nGrdPt = 0
call StatusLine(' espf:',' Making the grid')
if (iGrdTyp == 1) call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,.false.,IsMM,-iGrdTyp,DGrid,nAtQM)
call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,.false.,IsMM,iGrdTyp,DGrid,nAtQM)

! Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
! and B=ExtPot[TtT^-1]Tt
! Tt means the transpose of T

! Warning, at this point ExtPot is not filled (only TTT is needed)

call mma_allocate(T,nMult,nGrdPt,label='CartTensor')
call mma_allocate(TT,nMult,nMult,label='TT')
call mma_allocate(TTT,nGrdPt,nMult,label='TTT')
call mma_allocate(B,nGrdPt,label='ExtPot*TTT')
call InitB(nMult,natom,nAtQM,nGrdPt,Cord,Grid%A,T,TT,TTT,Ext,B,IsMM)

! Read the old MM electrostatic potential (and derivatives) from PotFile

IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.EXTPOT')
Line = Get_Ln(IPotFl)
call Get_I1(1,nChg)
if (nChg /= 0) then
  write(u6,*) 'ESPF: nChg ne 0 in h1_espf'
  call Abend()
end if
do iAt=1,natom
  Line = Get_Ln(IPotFl)
  call Get_I1(1,jAt)
  call Get_F(2,Ext(:,jAt),MxExtPotComp)
end do
close(IPotFl)

! Compute the quantum atomic multipoles

call mma_allocate(Mltp,nMult,label='ESPFMltp')
call espf_mltp(natom,MltOrd,nMult,nGrdPt,TTT,Mltp,Grid%A,IsMM,Ext,iPL-1)

! Run Tinker

UpdateVMM = .false.
if (allocated(OldMltp)) then
  sum1 = Zero
  sum2 = Zero
  sum3 = Zero
  sum4 = Zero
  do iMlt=1,nMult,MltOrd
    sum1 = abs(Mltp(iMlt)-OldMltp(iMlt))
    UpdateVMM = UpdateVMM .or. (sum1 > 1.0e-3_wp)
    if (MltOrd == 4) then
      sum2 = sum2+(Mltp(iMlt+1)-OldMltp(iMlt+1))**2
      sum3 = sum3+(Mltp(iMlt+2)-OldMltp(iMlt+2))**2
      sum4 = sum4+(Mltp(iMlt+3)-OldMltp(iMlt+3))**2
    end if
  end do
  rms2 = sqrt(sum2/nMult)
  rms3 = sqrt(sum3/nMult)
  rms4 = sqrt(sum4/nMult)
  if (MltOrd == 4) then
    UpdateVMM = UpdateVMM .or. (rms2 > -1.0e-2_wp)
    UpdateVMM = UpdateVMM .or. (rms3 > -1.0e-2_wp)
    UpdateVMM = UpdateVMM .or. (rms4 > -1.0e-2_wp)
  end if
  call mma_deallocate(OldMltp)
else
  UpdateVMM = .true.
end if
iQMchg = 1
if (First .and. Do_DFT) UpdateVMM = .true.
if (UpdateVMM) call RunTinker(natom,Cord,Mltp,.false.,IsMM,MltOrd,DynExtPot,iQMchg,nAtMM,StandAlone,DoDirect)

! Read the MM electrostatic potential (and derivatives) from PotFile

IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.EXTPOT')
Line = Get_Ln(IPotFl)
call Get_I1(1,nChg)
if (nChg /= 0) then
  write(u6,*) 'ESPF: nChg ne 0 in h1_espf'
  call Abend()
end if
do iAt=1,natom
  Line = Get_Ln(IPotFl)
  call Get_I1(1,jAt)
  call Get_F(2,Ext(:,jAt),MxExtPotComp)
end do
close(IPotFl)

! Recompute the B matrix

call InitB(nMult,natom,nAtQM,nGrdPt,Cord,Grid%A,T,TT,TTT,Ext,B,IsMM)

! Compute the modification of the core hamiltonian

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call StatusLine(' espf:',' Computing energy components')
call espf_energy(nBas(0),natom,nGrdPt,Ext,Grid%A,B,h1,nh1,RepNuc,RealDummy,DoTinker,DoGromacs,DynExtPot)

! Save the modified h1 matrix

call Put_Temp('h1    XX',h1,nh1)
call Put_Temp('PotNucXX',[RepNuc],1)
if (.not. DynExtPot) then
  call Put_Temp('h1_raw  ',h1,nh1)
  call Put_Temp('PotNuc00',[RepNuc],1)
end if

! Save data in the ESPF.DATA file

call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,Mltp,nMult,IsMM,natom,.false.,.false.,DoDirect)

! The end

call mma_deallocate(T)
call mma_deallocate(TT)
call mma_deallocate(TTT)
call mma_deallocate(B)
call mma_deallocate(Cord)
call mma_deallocate(IsMM)
call mma_deallocate(Ext)
call mma_deallocate(Mltp)
call mma_deallocate(Grid%A)
if (allocated(DGrid%A)) call mma_deallocate(DGrid%A)

return

end subroutine h1_espf
