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
use Index_Functions, only: nTri_Elem1
use Data_Structures, only: Alloc2DArray_Type, Alloc4DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: lSave
integer(kind=iwp) :: iAt, ibla, iGrdTyp, ii, iPL, IPotFl, iRMax, jAt, MltOrd, natom, nAtMM, nAtQM, nChg, nGrdPt, nMult
real(kind=wp) :: DeltaR
logical(kind=iwp) :: DoDirect, DoGromacs, DoTinker, Exists, lMorok
character(len=180) :: ESPFLine
character(len=10) :: ESPFKey
type(Alloc2DArray_Type) :: Grid
type(Alloc4DArray_Type) :: DGrid
integer(kind=iwp), allocatable :: IsMM(:)
real(kind=wp), allocatable :: B(:), Cord(:,:), Ext(:,:), Mltp(:), T(:,:), TT(:,:), TTT(:,:)
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
        ibla = ibla+nTri_Elem1(ii)
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

call Get_iScalar('Unique atoms',natom)
call mma_allocate(Cord,3,natom,label='AtomCoord')
call Get_dArray('Unique Coordinates',Cord,3*natom)
call mma_allocate(IsMM,natom,label='IsMM for atoms')
call mma_allocate(Ext,MxExtPotComp,natom,label='ExtPot')
Ext(:,:) = Zero
call MMCount(natom,nAtMM,IsMM)
nAtQM = natom-nAtMM
nMult = MltOrd*nAtQM

! Read the ESPF potential (and derivatives) from PotFile

IPotFl = IsFreeUnit(33)
IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.EXTPOT')
ESPFLine = Get_Ln(IPotFl)
call Get_I1(1,nChg)
if (nChg /= 0) then
  write(u6,*) 'ESPF: nChg /= 0 in espf_analysis'
  call Abend()
end if
do iAt=1,natom
  ESPFLine = Get_Ln(IPotFl)
  call Get_I1(1,jAt)
  call Get_F(2,Ext(:,jAt),MxExtPotComp)
end do
close(IPotFl)

! Compute the grid around the molecule

!nGrdPt = 0
call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,.false.,IsMM,iGrdTyp,DGrid,nAtQM)

! Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
! and B=ExtPot[TtT^-1]Tt
! Tt means the transpose of T

call mma_allocate(T,nMult,nGrdPt,label='CartTensor')
call mma_allocate(TT,nMult,nMult,label='TT')
call mma_allocate(TTT,nGrdPt,nMult,label='TTT')
call mma_allocate(B,nGrdPt,label='ExtPot*TTT')
call InitB(nMult,natom,nAtQM,nGrdPt,Cord,Grid%A,T,TT,TTT,Ext,B,IsMM)

! Now the analysis

call mma_allocate(Mltp,nMult,label='ESPFMltp')
call espf_mltp(natom,MltOrd,nMult,nGrdPt,TTT,Mltp,Grid%A,IsMM,Ext,iPL+1)
call Add_Info('ESPF multipoles',Mltp,nMult,6)

! Save some data

if (lSave) &
  call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,Mltp,nMult,IsMM,natom,.false.,.false.,DoDirect)

! The end

if (iPL >= 2) then
  call CollapseOutput(0,'   ESPF analysis')
  write(u6,*)
end if
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

call ClsSew()

return

end subroutine espf_analysis
