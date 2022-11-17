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
use Basis_Info, only: nBas
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nh1
real(kind=wp) :: h1(nh1), RepNuc
logical(kind=iwp) :: First, Do_DFT
#include "WrkSpc.fh"
#include "print.fh"
integer(kind=iwp) :: iAt, ibla, iGrdTyp, ii, iMlt, iMode, ipB, ipCord, ipDGrd, ipExt, ipGrid, ipIsMM, iPL, ipMltp, ipOldMltp, &
                     IPotFl, ipT, ipTT, ipTTT, iQMchg, iRMax, iSize1, iSize2, iSize3, ITkQMMM, jAt, MltOrd, natMM, natom, nAtQM, &
                     nChg, nGrdPt, nMult, nSym
real(kind=wp) :: DeltaR, RealDummy, rms2, rms3, rms4, sum1, sum2, sum3, sum4
logical(kind=iwp) :: DoDirect, DoGromacs, DoTinker, DynExtPot, Exists, lMorok, StandAlone, UpdateVMM
character(len=180) :: Line
character(len=10) :: ESPFKey
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
ipOldMltp = ip_Dummy
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
        ibla = ibla+(ii+2)*(ii+1)/2
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
      call Allocate_Work(ipOldMltp,nMult)
      do iMlt=1,nMult,MltOrd
        Line = Get_Ln(IPotFl)
        call Get_I1(1,iAt)
        call Get_F(2,Work(ipOldMltp+iMlt-1),MltOrd)
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
  if (ipOldMltp /= ip_Dummy) call Free_Work(ipOldMltp)
  return
end if

ipIsMM = ip_iDummy
call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
nMult = MltOrd*nAtQM

! Compute the grid around the molecule

nGrdPt = 0
ipGrid = ip_Dummy
ipDGrd = ip_Dummy
call StatusLine(' espf:',' Making the grid')
if (iGrdTyp == 1) then
  call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.false.,ipIsMM,-iGrdTyp,ipDGrd,nAtQM)
  call GetMem('ESPF_Grid','Allo','Real',ipGrid,3*nGrdPt)
  call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.false.,ipIsMM,iGrdTyp,ipDGrd,nAtQM)
else
  call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.false.,ipIsMM,iGrdTyp,ipDGrd,nAtQM)
end if

! Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
! and B=ExtPot[TtT^-1]Tt
! Tt means the transpose of T

! Warning, at this point ExtPot is not filled (only TTT is needed)

iSize1 = nMult*nGrdPt
iSize2 = nMult*nMult
iSize3 = nMult*max(nMult,nGrdPt)
call GetMem('CartTensor','Allo','Real',ipT,iSize1)
call GetMem('TT','Allo','Real',ipTT,iSize2)
call GetMem('TTT','Allo','Real',ipTTT,iSize3)
call GetMem('ExtPot*TTT','Allo','Real',ipB,nGrdPt)
call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipB,ipIsMM)

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
  call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
end do
close(IPotFl)

! Compute the quantum atomic multipoles

call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,ipIsMM,ipExt,iPL-1)

! Run Tinker

UpdateVMM = .false.
if (ipOldMltp /= ip_Dummy) then
  sum1 = Zero
  sum2 = Zero
  sum3 = Zero
  sum4 = Zero
  do iMlt=1,nMult,MltOrd
    sum1 = abs(Work(ipMltp+iMlt-1)-Work(ipOldMltp+iMlt-1))
    UpdateVMM = UpdateVMM .or. (sum1 > 1.0e-3_wp)
    if (MltOrd == 4) then
      sum2 = sum2+(Work(ipMltp+iMlt)-Work(ipOldMltp+iMlt))**2
      sum3 = sum3+(Work(ipMltp+iMlt+1)-Work(ipOldMltp+iMlt+1))**2
      sum4 = sum4+(Work(ipMltp+iMlt+2)-Work(ipOldMltp+iMlt+2))**2
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
  call Free_Work(ipOldMltp)
else
  UpdateVMM = .true.
end if
iQMchg = 1
if (First .and. Do_DFT) UpdateVMM = .true.
if (UpdateVMM) call RunTinker(natom,Work(ipCord),ipMltp,iWork(ipIsMM),MltOrd,DynExtPot,iQMchg,natMM,StandAlone,DoDirect)

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
  call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
end do
close(IPotFl)

! Recompute the B matrix

call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipB,ipIsMM)

! Compute the modification of the core hamiltonian

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call StatusLine(' espf:',' Computing energy components')
call espf_energy(nBas(0),natom,nGrdPt,ipExt,ipGrid,ipB,h1,nh1,RepNuc,RealDummy,DoTinker,DoGromacs,DynExtPot)

! Save the modified h1 matrix

call Put_Temp('h1    XX',h1,nh1)
call Put_Temp('PotNucXX',[RepNuc],1)
if (.not. DynExtPot) then
  call Put_Temp('h1_raw  ',h1,nh1)
  call Put_Temp('PotNuc00',[RepNuc],1)
end if

! Save data in the ESPF.DATA file

call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,.false.,.false.,DoDirect)

! The end

call GetMem('ExtPot*TTT','Free','Real',ipB,nGrdPt)
call GetMem('ExtPot','Free','Real',ipExt,natom*MxExtPotComp)
call GetMem('TTT','Free','Real',ipTTT,iSize3)
call GetMem('TT','Free','Real',ipTT,iSize2)
call GetMem('CartTensor','Free','Real',ipT,iSize1)
call GetMem('ESPFMltp','Free','Real',ipMltp,nMult)
call GetMem('ESPF_Grid','Free','Real',ipGrid,3*nGrdPt)
call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
call GetMem('AtomCoord','Free','Real',ipCord,3*natom)

return

end subroutine h1_espf
