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

subroutine espf_energy(nBas0,natom,nGrdPt,Ext,Grid,B,h1,nh1,RepNuc,EnergyCl,DoTinker,DoGromacs,DynExtPot)
! Compute the integrals <mu|B/R_grid|nu>, where B weights every
! point of the grid and R_grid is the distance to one grid point.

use espf_global, only: MxExtPotComp
use Index_Functions, only: nTri_Elem
use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auTokcalmol
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nBas0, natom, nGrdPt, nh1
real(kind=wp), intent(in) :: Ext(MxExtPotComp,natom), Grid(3,nGrdPt), EnergyCl
real(kind=wp), intent(_IN_) :: B(nGrdPt)
real(kind=wp), intent(inout) :: h1(nh1), RepNuc
logical(kind=iwp), intent(in) :: DoTinker, DoGromacs, DynExtPot
integer(kind=iwp) :: i, iAddPot, iComp, idum(1), iOpt, iPL, iRc, iSyLbl, ITkQMMM, ncmp, nInts, nSize
real(kind=wp) :: EQC, opnuc(1), RepNuc_old, TkE
character(len=180) :: Line
character(len=8) :: Label
real(kind=wp), allocatable :: IntOnGrid(:)
integer(kind=iwp), external :: iPL_espf, IsFreeUnit, IsStructure
real(kind=wp), external :: ExtNuc
character(len=180), external :: Get_Ln

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
  TkE = TkE/auTokcalmol
  RepNuc_old = RepNuc
  RepNuc = RepNuc+TkE
  if (iPL >= 3) write(u6,3000) RepNuc_old,TkE,RepNuc
else if (DoGromacs) then
  RepNuc_old = RepNuc
  RepNuc = RepNuc+EnergyCl
  if (iPL >= 3) write(u6,3000) RepNuc_old,EnergyCl,RepNuc
end if

! Call to DrvPot to compute the integrals
! Here we don't care about opnuc (nuclear potential)

nSize = nTri_Elem(nBas0)+4
if (nSize /= (nh1+4)) then
  write(u6,*) 'In espf_energy, nSize ne nh1',nSize,nh1+4
  call Abend()
end if
opnuc = Zero

ncmp = 1
iAddPot = 1
if (iPL >= 4) then
  do i=1,NGrdPt
    write(u6,1234) i,Grid(:,i),B(i)
  end do
end if
call DrvPot(Grid,opnuc,ncmp,B,nGrdPt,iAddPot)
Label = 'Pot     '
iComp = 1
iSyLbl = 1
iRc = -1
iOpt = ibset(0,sOpSiz)
call iRdOne(iRc,iOpt,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) then
  write(u6,'(A)') ' ESPF: Error reading ONEINT'
  write(u6,'(A,A8)') ' Label = ',Label
  call Abend()
end if
if (nInts+4 /= nSize) then
  write(u6,'(A,2I5)') ' ESPF: nInts+4 /= nSize',nInts+4,nSize
  call Abend()
end if
call mma_allocate(IntOnGrid,nSize,label='IntOnGrid')
iOpt = 0
call RdOne(iRc,iOpt,Label,iComp,IntOnGrid,iSyLbl)
if (iPL >= 4) call TriPrt(Label,' ',IntOnGrid,nBas0)

! The core Hamiltonian must be updated

h1(1:nInts) = h1(1:nInts)+IntOnGrid(1:nInts)
if (DynExtPot) then
  iSyLbl = 1
  iRc = -1
  iOpt = 0
  iComp = 1
  Label = 'OneHamRF'
  call WrOne(iRc,iOpt,Label,iComp,IntOnGrid,iSyLbl)
end if
call mma_deallocate(IntOnGrid)

! The electrostatic energy between the external potential
! and the nuclei is added to the nuclear energy

EQC = ExtNuc(Ext,natom)
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
