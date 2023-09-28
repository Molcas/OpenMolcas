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

subroutine espf(ireturn,StandAlone)

use espf_global, only: MxExtPotComp
use Index_Functions, only: nTri_Elem
use Basis_Info, only: nBas
use OneDat, only: sOpSiz
use Symmetry_Info, only: VarR, VarT, Symmetry_Info_Dmp
use Data_Structures, only: Alloc1DArray_Type, Alloc2DArray_Type, Alloc4DArray_Type
use NAC, only: isNAC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
logical(kind=iwp), intent(in) :: StandAlone
integer(kind=iwp) :: iComp, idum(1), iGrdTyp, iOpt, iOption, iPL, iRc, iRMax, iSyLbl, MltOrd, nAtMM, natom, nAtQM, nBas0, nGrdPt, &
                     nInts, nMult, nSize, nSym
real(kind=wp) :: DeltaR, EnergyCl, RepNuc
logical(kind=iwp) :: Close_Seward, DoDirect, DoGromacs, DoTinker, DynExtPot, Forces, isNAC_tmp, lMorok, Show_espf
character(len=8) :: Label
type(Alloc1DArray_Type) :: Mltp
type(Alloc2DArray_Type) :: GradCl, Grid
type(Alloc4DArray_Type) :: DGrid
integer(kind=iwp), allocatable :: IsMM(:)
real(kind=wp), allocatable :: B(:), Cord(:,:), DB(:,:,:), Ext(:,:), H(:), T(:,:), TT(:,:), TTT(:,:)
integer(kind=iwp), external :: iPL_espf

iReturn = 99

! Print

iPL = iPL_espf()
Close_Seward = .false.

! Some warnings

call Get_iScalar('nSym',nSym)
if (nSym > 1) then
  write(u6,'(A)') ' Symmetry cannot be used together with ESPF.'
  call Quit_OnUserError()
end if

! Set on the System Bit 11

call Get_iScalar('System Bitswitch',iOption)
iOption = ior(iOption,2**11)
call Put_iScalar('System Bitswitch',iOption)

! Some initializations

call Get_iScalar('Unique atoms',natom)
call mma_allocate(Cord,3,natom,label='AtomCoord')
call Get_dArray('Unique Coordinates',Cord,3*natom)
call mma_allocate(IsMM,natom,label='IsMM for atoms')
call mma_allocate(Ext,MxExtPotComp,natom,label='ExtPot')
Ext(:,:) = Zero
call MMCount(natom,nAtMM,IsMM)
nAtQM = natom-nAtMM
nGrdPt = 0
isNAC_tmp = isNAC

! Read the input and compute the external potential

call StatusLine(' espf:',' Reading input')
call ReadIn_ESPF(natom,Cord,Ext,MltOrd,iRMax,DeltaR,Forces,Show_espf,IsMM,StandAlone,iGrdTyp,DoTinker,DoGromacs,DynExtPot,Mltp, &
                 nAtMM,lMorok,DoDirect,GradCl,EnergyCl)

! If the present calculation does not use ESPF but the Direct scheme

if (DoDirect) then

  call No_ESPF(Forces,DoTinker)

else

  nMult = MltOrd*nAtQM
  if (iPL >= 2) write(u6,'(/,A,I2,A,i4,A,i6)') ' Number of ESPF operators (nMult=',MltOrd,' * nAtQM=',nAtQM,'): ',nMult

  ! Compute the grid around the molecule

  call StatusLine(' espf:',' Making the grid')
  if (iGrdTyp == 1) then
    if (nGrdPt == 0) call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,Forces,IsMM,-iGrdTyp,DGrid,nAtQM)
    call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,Forces,IsMM,iGrdTyp,DGrid,nAtQM)
    if (iPL >= 2) then
      write(u6,'(A)') ' PNT Grid (Warning: no grid derivatives)'
      write(u6,'(A)') ' (C. Chipot and J. Angyan, Henri Poincare University, Nancy, France)'
      write(u6,'(5X,I5,A)') nGrdPt,' grid points'
    end if
  else
    call MkGrid(natom,Cord,Grid,nGrdPt,iRMax,DeltaR,Forces,IsMM,iGrdTyp,DGrid,nAtQM)
    if (iPL >= 2) then
      write(u6,'(A)') ' GEPOL Grid, using United Atoms radii'
      write(u6,'(5X,I5,A)') nGrdPt,' grid points'
    end if
  end if

  ! If this is a standalone call to &ESPF, there are 2 options:
  !    1) static external potential: compute here the ESPF contributions
  !    2) dynamic external potential: nothing more to compute here

  if (.not. (StandAlone .and. DynExtPot)) then

    ! Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
    ! and B=ExtPot[TtT^-1]Tt
    ! Tt means the transpose of T

    call mma_allocate(T,nMult,nGrdPt,label='CartTensor')
    call mma_allocate(TT,nMult,nMult,label='TT')
    call mma_allocate(TTT,nGrdPt,nMult,label='TTT')
    call mma_allocate(B,nGrdPt,label='ExtPot*TTT')
    call InitB(nMult,natom,nAtQM,nGrdPt,Cord,Grid%A,T,TT,TTT,Ext,B,IsMM)
    call mma_allocate(DB,nGrdPt,3,nAtQM,label='DerivB')
    call InitDB(nMult,natom,nAtQM,nGrdPt,Cord,Grid%A,T,TT,TTT,Ext,DB,IsMM)

    ! Here we must distinguish between an energy run and a gradient run

    if (.not. Forces) then
      call StatusLine(' espf:',' Computing energy components')
      call Get_iArray('nBas',nBas,nSym)
      nBas0 = nBas(0)
      nSize = nTri_Elem(nBas0)+4
      call mma_allocate(H,nSize,label='H')
      iComp = 1
      iSyLbl = 1
      Label = 'OneHam  '
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
      iRc = -1
      iOpt = 0
      call RdOne(iRc,iOpt,Label,iComp,H,iSyLbl)
      call Get_dScalar('PotNuc',RepNuc)
      call espf_energy(nBas0,natom,nGrdPt,Ext,Grid%A,B,H,nSize-4,RepNuc,EnergyCl,DoTinker,DoGromacs,DynExtPot)
      call Put_dScalar('PotNuc',RepNuc)
      call WrOne(iRc,iOpt,Label,iComp,H,iSyLbl)
      if (iRC /= 0) then
        write(u6,*) 'ESPF: Error writing to ONEINT'
        write(u6,'(A,A8)') 'Label=',Label
        call Abend()
      end if
      call mma_deallocate(H)
      if (iPL >= 3) write(u6,*) 'The 1-e hamiltonian is now updated.'
      if (iPL >= 2) write(u6,'(A,F16.10)') ' Nuclear energy, including Ext Pot = ',RepNuc
    else
      call StatusLine(' espf:',' Computing gradient components')
      if (.not. allocated(GradCl%A)) call mma_allocate(GradCl%A,0,0,label='GradCl')
      call espf_grad(natom,nGrdPt,nAtQM,Ext,Grid%A,B,DB,IsMM,GradCl%A,DoTinker,DoGromacs)
      if (.not. allocated(Mltp%A)) call mma_allocate(Mltp%A,nMult,label='ESPFMltp')
      call espf_mltp(natom,MltOrd,nMult,nGrdPt,TTT,Mltp%A,Grid%A,IsMM,Ext,iPL)
    end if
    Close_Seward = .true.

  end if

end if

! Save data in the ESPF.DATA file

if (.not. allocated(Mltp%A)) then
  nMult = 0
  call mma_allocate(Mltp%A,nMult,label='ESPFMltp')
end if
call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,Mltp%A,nMult,IsMM,natom,Show_espf,Forces,DoDirect)

! Exit

isNAC = isNAC_tmp
if (.not. (StandAlone .and. DynExtPot)) then
  call mma_deallocate(T)
  call mma_deallocate(TT)
  call mma_deallocate(TTT)
  call mma_deallocate(B)
  call mma_deallocate(DB)
end if
call mma_deallocate(Cord)
call mma_deallocate(IsMM)
call mma_deallocate(Ext)
call mma_deallocate(Mltp%A)
if (allocated(Grid%A)) call mma_deallocate(Grid%A)
if (allocated(DGrid%A)) call mma_deallocate(DGrid%A)
if (allocated(GradCl%A)) call mma_deallocate(GradCl%A)

! Slapaf needs to know that the gradient is NOT translational
! and rotational invariant.

if ((.not. Forces) .and. (nAtMM > 0)) then
  VarR = .true.
  VarT = .true.
  call Symmetry_Info_Dmp()
end if
if (Close_Seward) call ClsSew()

iReturn = 0
return

end subroutine espf
