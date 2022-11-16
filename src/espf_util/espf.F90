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

use Basis_Info, only: nBas
use OneDat, only: sOpSiz
use Symmetry_Info, only: VarR, VarT, Symmetry_Info_Dmp
use Definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp) :: ireturn
logical(kind=iwp) :: StandAlone
#include "espf.fh"
#include "nac.fh"
integer(kind=iwp) :: iComp, idum(1), iGrdTyp, iOpt, iOption, ipB, ipCord, ipDB, ipDGrd, ipExt, ipGradCl, ipGrid, ipH, ipIsMM, iPL, &
                     ipMltp, ipT, ipTT, ipTTT, iRc, iRMax, iSize1, iSize2, iSize3, iSyLbl, MltOrd, natMM, natom, nAtQM, nBas0, &
                     nGrdPt, nInts, nMult, nSize, nSym
real(kind=wp) :: DeltaR, EnergyCl, RepNuc
logical(kind=iwp) :: Close_Seward, DoDirect, DoGromacs, DoTinker, DynExtPot, Forces, isNAC_tmp, lMorok, Show_espf
character(len=8) :: Label
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

call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
ipMltp = ip_Dummy
nGrdPt = 0
ipGrid = ip_Dummy
ipDGrd = ip_Dummy
isNAC_tmp = isNAC

! Read the input and compute the external potential

call StatusLine(' espf:',' Reading input')
call ReadIn_ESPF(natom,ipCord,ipExt,MltOrd,iRMax,DeltaR,Forces,Show_espf,ipIsMM,StandAlone,iGrdTyp,DoTinker,DoGromacs,DynExtPot, &
                 ipMltp,natMM,lMorok,DoDirect,ipGradCl,EnergyCl)

! If the present calculation does not use ESPF but the Direct scheme

if (DoDirect) then

  call No_ESPF(Forces,DoTinker)

else

  nMult = MltOrd*nAtQM
  if (iPL >= 2) write(u6,'(/,A,I2,A,i4,A,i6)') ' Number of ESPF operators (nMult=',MltOrd,' * nAtQM=',nAtQM,'): ',nMult

  ! Compute the grid around the molecule

  call StatusLine(' espf:',' Making the grid')
  if (iGrdTyp == 1) then
    if (nGrdPt == 0) call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,Forces,ipIsMM,-iGrdTyp,ipDGrd,nAtQM)
    call GetMem('ESPF_Grid','ALLO','REAL',ipGrid,3*nGrdPt)
    call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,Forces,ipIsMM,iGrdTyp,ipDGrd,nAtQM)
    if (iPL >= 2) then
      write(u6,'(A)') ' PNT Grid (Warning: no grid derivatives)'
      write(u6,'(A)') ' (C. Chipot and J. Angyan, Henri Poincare University, Nancy, France)'
      write(u6,'(5X,I5,A)') nGrdPt,' grid points'
    end if
  else
    call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,Forces,ipIsMM,iGrdTyp,ipDGrd,nAtQM)
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

    iSize1 = nMult*nGrdPt
    iSize2 = nMult*nMult
    iSize3 = nMult*max(nMult,nGrdPt)
    call GetMem('CartTensor','Allo','Real',ipT,iSize1)
    call GetMem('TT','Allo','Real',ipTT,iSize2)
    call GetMem('TTT','Allo','Real',ipTTT,iSize3)
    call GetMem('ExtPot*TTT','Allo','Real',ipB,nGrdPt)
    call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipB,ipIsMM)
    call GetMem('DerivB','Allo','Real',ipDB,nGrdPt*3*nAtQM)
    call InitDB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,ipExt,ipDB,ipIsMM)
    if ((iGrdTyp == 2) .and. (ipDGrd /= ip_Dummy)) call GetMem('ESPF_DGrid','Free','Real',ipDGrd,3*nGrdPt*3*nAtQM)

    ! Here we must distinguish between an energy run and a gradient run

    if (.not. Forces) then
      call StatusLine(' espf:',' Computing energy components')
      call Get_iArray('nBas',nBas,nSym)
      nBas0 = nBas(0)
      nSize = nBas0*(nBas0+1)/2+4
      call Allocate_Work(ipH,nSize)
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
      call RdOne(iRc,iOpt,Label,iComp,Work(ipH),iSyLbl)
      call Get_dScalar('PotNuc',RepNuc)
      call espf_energy(nBas0,natom,nGrdPt,ipExt,ipGrid,ipB,Work(ipH),nSize-4,RepNuc,EnergyCl,DoTinker,DoGromacs,DynExtPot)
      call Put_dScalar('PotNuc',RepNuc)
      call WrOne(iRc,iOpt,Label,iComp,Work(ipH),iSyLbl)
      if (iRC /= 0) then
        write(u6,*) 'ESPF: Error writing to ONEINT'
        write(u6,'(A,A8)') 'Label=',Label
        call Abend()
      end if
      call Free_Work(ipH)
      if (iPL >= 3) write(u6,*) 'The 1-e hamiltonian is now updated.'
      if (iPL >= 2) write(u6,'(A,F16.10)') ' Nuclear energy, including Ext Pot = ',RepNuc
    else
      call StatusLine(' espf:',' Computing gradient components')
      call espf_grad(natom,nGrdPt,ipExt,ipGrid,ipB,ipDB,ipIsMM,ipGradCl,DoTinker,DoGromacs)
      if (ipMltp == ip_Dummy) call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
      call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,ipIsMM,ipExt,iPL)
    end if
    Close_Seward = .true.

  end if

end if

! Save data in the ESPF.DATA file

call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,Show_espf,Forces,DoDirect)

! Exit

isNAC = isNAC_tmp
if (.not.(StandAlone .and. DynExtPot)) then
  call GetMem('DerivB','FREE','REAL',ipDB,nGrdPt*3*nAtQM)
  call GetMem('CartTensor','FREE','REAL',ipT,iSize1)
  call GetMem('TT','FREE','REAL',ipTT,iSize2)
  call GetMem('TTT','FREE','REAL',ipTTT,iSize3)
  call GetMem('ExtPot*TTT','FREE','REAL',ipB,nGrdPt)
end if
if (ipMltp /= ip_Dummy) call GetMem('ESPFMltp','FREE','REAL',ipMltp,nMult)
if (ipGrid /= ip_Dummy) call GetMem('ESPF_Grid','FREE','REAL',ipGrid,3*nGrdPt)
call GetMem('AtomCoord','FREE','REAL',ipCord,3*natom)
call GetMem('ExtPot','FREE','REAL',ipExt,natom*MxExtPotComp)
call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
if (DoGromacs .and. Forces) then
  call GetMem('GradCl','FREE','REAL',ipGradCl,3*natom)
end if

! Slapaf needs to know that the gradient is NOT translational
! and rotational invariant.

if ((.not. Forces) .and. (natMM > 0)) then
  VarR = .true.
  VarT = .true.
  call Symmetry_Info_Dmp()
end if
if (Close_Seward) call ClsSew()

iReturn = 0
return

end subroutine espf
