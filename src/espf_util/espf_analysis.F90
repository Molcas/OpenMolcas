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

subroutine espf_analysis(lSave)

use espf_global, only: MxExtPotComp
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp) :: lSave
#include "WrkSpc.fh"
integer(kind=iwp) :: iAt, ibla, iGrdTyp, ii, ipB, ipCord, ipDGrd, ipExt, ipGrid, ipIsMM, iPL, ipMltp, IPotFl, ipT, ipTT, ipTTT, &
                     iRMax, iSize1, iSize2, iSize3, jAt, MltOrd, natom, nAtQM, nChg, nGrdPt, nMult
real(kind=wp) :: DeltaR
logical(kind=iwp) :: DoDirect, DoGromacs, DoTinker, Exists, lMorok
character(len=180) :: ESPFLine
character(len=10) :: ESPFKey
integer(kind=iwp), external :: iPL_espf, IsFreeUnit
character(len=180), external :: Get_Ln

iPL = iPL_espf()

if (iPL >= 2) then
  write(u6,*)
  call CollapseOutput(1,'   ESPF analysis')
  write(u6,'(3X,A)') '   -------------'
end if

! Recover some ESPF data

MltOrd = 0
iRMax = 0
DeltaR = Zero
iGrdTyp = 0
nGrdPt = 0
DoTinker = .false.
DoGromacs = .false.
lMorok = .false.
DoDirect = .false.
call F_Inquire('ESPF.DATA',Exists)
if (Exists) then
  IPotFl = IsFreeUnit(1)
  call Molcas_Open(IPotFl,'ESPF.DATA')
  do
    ESPFLine = Get_Ln(IPotFl)
    ESPFKey = ESPFLine(1:10)
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
    else if (ESPFKey == 'GRID      ') then
      call Get_I1(2,nGrdPt)
    else if (ESPFKey == 'TINKER    ') then
      DoTinker = .true.
    else if (ESPFKey == 'GROMACS   ') then
      DoGromacs = .true.
    else if (ESPFKey == 'LA_MOROK  ') then
      lMorok = .true.
    else if (ESPFKey == 'DIRECT    ') then
      DoDirect = .true.
    else if (ESPFKey == 'ENDOFESPF ') then
      exit
    end if
  end do
  close(IPotFl)
else
  write(u6,*) 'No ESPF.DATA file. Abort'
  call Quit_OnUserError()
end if

call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
nMult = MltOrd*nAtQM

! Read the ESPF potential (and derivatives) from PotFile

IPotFl = IsFreeUnit(33)
IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.EXTPOT')
ESPFLine = Get_Ln(IPotFl)
call Get_I1(1,nChg)
if (nChg /= 0) then
  write(u6,*) 'ESPF: nChg ne 0 in espf_analysis'
  call Abend()
end if
do iAt=1,natom
  ESPFLine = Get_Ln(IPotFl)
  call Get_I1(1,jAt)
  call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
end do
close(IPotFl)

! Compute the grid around the molecule

!nGrdPt = 0
ipGrid = ip_Dummy
ipDGrd = ip_Dummy
if (iGrdTyp == 1) call GetMem('ESPF_Grid','Allo','Real',ipGrid,3*nGrdPt)
call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.false.,ipIsMM,iGrdTyp,ipDGrd,nAtQM)

! Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
! and B=ExtPot[TtT^-1]Tt
! Tt means the transpose of T

iSize1 = nMult*nGrdPt
iSize2 = nMult*nMult
iSize3 = nMult*max(nMult,nGrdPt)
call GetMem('CartTensor','Allo','Real',ipT,iSize1)
call GetMem('TT','Allo','Real',ipTT,iSize2)
call GetMem('TTT','Allo','Real',ipTTT,iSize3)
call GetMem('ExtPot*TTT','Allo','Real',ipB,nGrdPt)
call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipB,ipIsMM)

! Now the analysis

call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,ipIsMM,ipExt,iPL+1)
call Add_Info('ESPF multipoles',Work(ipMltp),nMult,6)

! Save some data

if (lSave) call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,.false.,.false., &
                           DoDirect)

! The end

if (iPL >= 2) then
  call CollapseOutput(0,'   ESPF analysis')
  write(u6,*)
end if
call GetMem('ESPFMltp','Free','Real',ipMltp,nMult)
call GetMem('ExtPot*TTT','Free','Real',ipB,nGrdPt)
call GetMem('TTT','Free','Real',ipTTT,iSize3)
call GetMem('TT','Free','Real',ipTT,iSize2)
call GetMem('CartTensor','Free','Real',ipT,iSize1)
call GetMem('ESPF_Grid','Free','Real',ipGrid,3*nGrdPt)
call GetMem('ExtPot','Free','Real',ipExt,natom*MxExtPotComp)
call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
call GetMem('AtomCoord','Free','Real',ipCord,3*natom)

call ClsSew()

return

end subroutine espf_analysis
