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

subroutine espf_energy(nBas0,natom,nGrdPt,ipExt,ipGrid,ipB,h1,nh1,RepNuc,EnergyCl,DoTinker,DoGromacs,DynExtPot)
! Compute the integrals <mu|B/R_grid|nu>, where B weights every
! point of the grid and R_grid is the distance to one grid point.

use OneDat, only: sOpSiz

implicit real*8(A-H,O-Z)
#include "espf.fh"
character*180 Line, Get_Ln
external Get_Ln
character*8 Label
logical DoTinker, DoGromacs, DynExtPot
real*8 h1(nh1)
dimension opnuc(1), idum(1)

iPL = iPL_espf()

! Read the MM contribution to the total energy and add it
! to the Nuclear Repulsion term

if (DoTinker) then
  ITkQMMM = IsFreeUnit(30)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'MMEnergy ') /= 0) call Get_F1(2,TkE)
  end do
  close(ITkQMMM)
  TkE = TkE*ToHartree
  RepNuc_old = RepNuc
  RepNuc = RepNuc+TkE
  if (iPL >= 3) write(6,3000) RepNuc_old,TkE,RepNuc
else if (DoGromacs) then
  RepNuc_old = RepNuc
  RepNuc = RepNuc+EnergyCl
  if (iPL >= 3) write(6,3000) RepNuc_old,EnergyCl,RepNuc
end if

! Call to DrvPot to compute the integrals
! Here we don't care about opnuc (nuclear potential)

nSize = nBas0*(nBas0+1)/2+4
if (nSize /= (nh1+4)) then
  write(6,*) 'In espf_energy, nSize ne nh1',nSize,nh1+4
  call Abend()
end if
opnuc = Dum

ncmp = 1
iAddPot = 1
if (iPL >= 4) then
  do i=1,NGrdPt
    write(6,1234) i,(Work(ipGrid+(i-1)*3+j),j=0,2),Work(ipB+(i-1))
  end do
end if
call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipB),nGrdPt,iAddPot)
Label = 'Pot     '
iComp = 1
iSyLbl = 1
iRc = -1
iOpt = ibset(0,sOpSiz)
call iRdOne(iRc,iOpt,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) then
  write(6,'(A)') ' ESPF: Error reading ONEINT'
  write(6,'(A,A8)') ' Label = ',Label
  call Abend()
end if
if (nInts+4 /= nSize) then
  write(6,'(A,2I5)') ' ESPF: nInts+4 /= nSize',nInts+4,nSize
  call Abend()
end if
call GetMem('IntOnGrid','Allo','Real',ipInt,nSize)
iOpt = 0
call RdOne(iRc,iOpt,Label,iComp,Work(ipInt),iSyLbl)
if (iPL >= 4) call TriPrt(Label,' ',Work(ipInt),nBas0)

! The core Hamiltonian must be updated

call daxpy_(nInts,One,Work(ipInt),1,h1,1)
if (DynExtPot) then
  iSyLbl = 1
  iRc = -1
  iOpt = 0
  iComp = 1
  Label = 'OneHamRF'
  call WrOne(iRc,iOpt,Label,iComp,Work(ipInt),iSyLbl)
end if
call GetMem('IntOnGrid','Free','Real',ipInt,nSize)

! The electrostatic energy between the external potential
! and the nuclei is added to the nuclear energy

EQC = ExtNuc(ipExt,natom)
RepNuc = RepNuc+EQC
if (IsStructure() == 1) then
  call Add_Info('PotNuc',[RepNuc],1,6)
else
  call Add_Info('PotNuc',[RepNuc],1,12)
end if

return

1234 format('Grid point ',I4,/,4F12.6)
3000 format(/,' RepNuc + MM = ',F13.8,' + ',F13.8,' = ',F13.8)

end subroutine espf_energy
