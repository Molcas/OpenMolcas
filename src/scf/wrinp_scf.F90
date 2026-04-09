!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine WrInp_SCF(SIntTh)
!***********************************************************************
!                                                                      *
!     purpose: Write input                                             *
!                                                                      *
!     input: SIntTh: Threshold for Integral prescreening               *
!                                                                      *
!***********************************************************************

use Functionals, only: Print_Info
use KSDFT_Info, only: CoefR, CoefX
use InfSCF, only: Algo, Aufb, DDnoff, DelThr, DIIS, DIISTh, DltNth, dmpk, DoCholesky, DSCF, DThr, EThr, Expand, FThr, Header, &
                  iAU_ab, InVec, isHDF5, IterSO_Max, jPrint, jVOut, kIVO, kOptim_Max, KSDFT, LKOn, lpaper, MiniDn, nAufb, nBas, &
                  nCore, nD, nDel, nDIsc, nFro, nIter, nMem, nOcc, NoExchange, nOrb, nSym, nTit, One_Grid, PreSch, QNRTh, ReOrd, &
                  RFPert, RGEK, rTemp, Scrmbl, StVec, Teee, TemFac, Thize, Title, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, &
                  TStop, VTitle
use Fock_util_global, only: Deco
use RICD_Info, only: Do_DCCD
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: SIntTh
integer(kind=iwp) :: i, iCharge, iDoRI, iSym, iTit, mTmp
logical(kind=iwp) :: NonEq
character(len=72) :: Line
character(len=60) :: Frmt, FmtI, FmtR
character(len=7) :: ExpMeth
character(len=3) :: lIrrep(8)

if (jPrint >= 2) then
  call CollapseOutput(1,'   Input section:')
  write(u6,'(3X,A)') '   --------------'
  write(u6,*)
end if

! Print out header of the integral file
if (jPrint >= 2) then
  write(u6,'(6X,A)') 'Header of the integral files:'
  write(Line,'(A72)') Header(1)
  write(u6,'(6X,A)') trim(Line)
  write(Line,'(A72)') Header(2)
  write(u6,'(6X,A)') trim(Line)
  write(u6,*)
end if

! Print out title (if any)
if (nTit > 0) then
  if (jPrint >= 3) then
    call Banner(Title,nTit,lPaper-7)
    write(u6,*)
  else if (jPrint >= 2) then
    write(u6,'(6X,A)') ' Title:'
    do iTit=1,nTit
      write(u6,'(8X,A)') trim(Title(iTit))
    end do
  end if
end if

! Print the coordinates of the system
if (jPrint >= 2) call PrCoor()

! Print out Print level
!write(u6,'(1X,A,I3)') 'Print level=',iPrint
!write(u6,*)

call Get_cArray('Irreps',lIrrep,24)
do iSym=1,nSym
  lIrrep(iSym) = adjustr(lIrrep(iSym))
end do

if (jPrint >= 2) then
  call CollapseOutput(0,'   Input section:')
  write(u6,*)
end if
! Print out orbital informations
Frmt = '(6X,A,T35,8I4)'
if (jPrint >= 2) then
  call CollapseOutput(1,'   Orbital specifications:')
  write(u6,'(3X,A)') '   -----------------------'
  write(u6,*)
  write(u6,Frmt) 'Symmetry species',(i,i=1,nSym)
  write(u6,'(6X,A,T35,8(1X,A))') '                ',(lIrrep(i),i=1,nSym)
  write(u6,Frmt) 'Frozen orbitals',(nFro(i),i=1,nSym)
end if

if (Aufb) then
  if (nD == 1) then
    if (nAufb(1) == -1) then
      Tot_El_Charge = Tot_Charge-Tot_Nuc_Charge
      ! Check that Tot_El_Charge is a multiple of two!
      mtmp = int(-Tot_El_Charge+0.1_wp)
      if (mod(mtmp,2) /= 0) then
        write(u6,*) 'WrInp: Error in definition of molecular charge!'
        write(u6,*) 'Current implementation only allows double occupations.'
        write(u6,*) 'Tot_Charge    :',Tot_Charge
        write(u6,*) 'Tot_El_Charge :',Tot_El_Charge
        write(u6,*) 'Tot_Nuc_Charge:',Tot_Nuc_Charge
        call Abend()
      end if
      nAufb(1) = mtmp/2
    end if
  else
    if (nAufb(1) == -1) then
      Tot_El_Charge = Tot_Charge-Tot_Nuc_Charge
      ! Check that Tot_El_Charge is a multiple of two!
      mtmp = int(-Tot_El_Charge+0.1_wp)
      ! if ZSPIN is not set - make difference alpha-beta = 0 or 1
      nAufb(2) = (mtmp-iAu_ab)/2
      nAufb(1) = int(-Tot_El_Charge-nAufb(2))
    end if
    !write(u6,*) ' CHARGE + UHF is un'
    !call Abend()
  end if
  if ((nD == 1) .and. (jPrint >= 2)) then
    write(u6,Frmt) 'Aufbau',nAufb(1)
  else if (jPrint >= 3) then
    write(u6,Frmt) 'Aufbau alpha',nAufb(1)
    write(u6,Frmt) 'Aufbau beta ',nAufb(2)
  end if
  if (Teee .and. (jPrint >= 2)) then
    write(u6,'(a,f6.3)') '      Start temperature =',RTemp
    write(u6,'(a,f6.3)') '      End temperature   =',TStop
    write(u6,'(a,f6.3)') '      Temperature Factor=',TemFac
  end if
else
  if ((nD == 1) .and. (jPrint >= 2)) then
    write(u6,Frmt) 'Occupied orbitals',(nOcc(i,1),i=1,nSym)
    write(u6,Frmt) 'Secondary orbitals',(nOrb(i)-nOcc(i,1),i=1,nSym)
  else if (jPrint >= 2) then
    write(u6,Frmt) 'Occupied orbitals alpha',(nOcc(i,1),i=1,nSym)
    write(u6,Frmt) 'Occupied orbitals beta ',(nOcc(i,2),i=1,nSym)
    write(u6,Frmt) 'Secondary orbitals alpha',(nOrb(i)-nOcc(i,1),i=1,nSym)
    write(u6,Frmt) 'Secondary orbitals beta',(nOrb(i)-nOcc(i,2),i=1,nSym)

  end if
end if

Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
call put_dscalar('Total Charge    ',Tot_Charge)
if (jPrint >= 2) then
  write(u6,Frmt) 'Deleted orbitals         ',(nDel(i),i=1,nSym)
  write(u6,Frmt) 'Total number of orbitals ',(nOrb(i),i=1,nSym)
  write(u6,Frmt) 'Number of basis functions',(nBas(i),i=1,nSym)
  call CollapseOutput(0,'   Orbital specifications:')
  write(u6,*)
  write(u6,'(6X,A,T45,F10.3)') 'Molecular charge',Tot_Charge
  write(u6,*)
end if

! Print out reaction field specifications

iCharge = int(Tot_Charge)
NonEq = .false.
call PrRF(.false.,NonEq,iCharge,jPrint)
if (RFpert .and. (jPrint >= 2)) then
  write(u6,'(6X,A)') 'Reaction field specifications:'
  write(u6,'(6X,A)') '------------------------------'
  write(u6,*)
  write(u6,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previous calculation'
  write(u6,*)
end if

! Print out grid information in case of DFT

if (KSDFT /= 'SCF') then
  call Put_dScalar('DFT exch coeff',CoefX)
  call Put_dScalar('DFT corr coeff',CoefR)
  call Put_dScalar('EThr',EThr)
  call Funi_Print()
  if (jPrint >= 2) then
    if (One_Grid) then
      write(u6,'(6X,A)') 'The same grid will be used for all iterations.'
    else
      write(u6,'(6X,A)') 'A smaller intermediate grid will be used the first few iterations.'
    end if
    write(u6,*)
    write(u6,'(6X,A)') 'DFT functional specifications'
    write(u6,'(6X,A)') '-----------------------------'
    call libxc_version()
    call Print_Info()
    write(u6,*)
  end if
end if

! Print out informations concerning Direct/Conventional scheme
FmtI = '(6X,A,T50,I9)'
FmtR = '(6X,A,T50,ES9.2)'
if (jPrint >= 2) then
  call CollapseOutput(1,'   Optimization specifications:')
  write(u6,'(3X,A)') '   ----------------------------'
  write(u6,*)
end if
if (DSCF .and. (.not. Do_DCCD) .and. (jPrint >= 2)) then
  if ((nDisc == 0) .and. (nCore == 0)) then
    write(u6,'(6X,A)') 'SCF Algorithm: Direct'
  else if (nDisc*1024 >= nCore) then
    write(u6,'(6X,A)') 'SCF Algorithm: Semi-direct'

    ! The threshold to be used in the direct SCF procedure is
    ! defined by the energy threshold.
    write(u6,FmtI) 'Max MByte of integrals on disk/process:',nDisc
    write(u6,FmtR) 'Threshold for saving integrals on disc',Thize
  else
    write(u6,'(6X,A)') 'SCF Algorithm: Semi-direct in-core'

    ! The threshold to be used in the direct SCF procedure is
    ! defined by the energy threshold.
    write(u6,FmtI) 'Max kByte of integrals in memory/process:',nCore
    write(u6,FmtR) 'Threshold for saving integrals in memory',Thize
  end if
  if (PreSch) then
    write(u6,'(6X,A)') 'Prescreening Scheme: Only Integral value'
  else
    write(u6,'(6X,A)') 'Prescreening Scheme: Integral*Density value'
  end if
else if (jPrint >= 2) then
  if (nD == 1) then
    if (.not. DoCholesky) then
      write(u6,'(6X,A)') 'SCF Algorithm: Conventional'
    else
      call Get_iScalar('System BitSwitch',iDoRI)
      if (iand(iDoRI,1024) == 1024) then
        if (LKon) then
          write(u6,'(6X,A)') 'SCF Algorithm: LK-RI/DF'
          write(u6,FmtR) 'LK screening threshold:',dmpk
        else
          write(u6,'(6X,A)') 'SCF Algorithm: RI/DF'
        end if
      else
        if (LKon) then
          write(u6,'(6X,A)') 'SCF Algorithm: LK-Cholesky'
          write(u6,FmtR) 'LK screening threshold:',dmpk
        else
          write(u6,'(6X,A,I1)') 'SCF Algorithm: Cholesky'
        end if
      end if

      if (ALGO == 0) then
        if (iand(iDoRI,1024) == 1024) then
          write(u6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
        else
          write(u6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
        end if
      else if (ALGO == 1) then
        write(u6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
      else if (ALGO == 2) then
        write(u6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
      else if (ALGO == 3) then
        write(u6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
      else if (ALGO == 4) then
        write(u6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
      end if

      if (Do_DCCD) write(u6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
      if (ReOrd) write(u6,'(6X,A)') ' - the Cholesky vectors are reordered'
      if (DeCo) write(u6,'(6X,A)') ' - the density matrix is decomposed'
      if (DSCF .and. (nDisc == 0) .and. (nCore == 0)) write(u6,'(6X,A)') ' - SCF Algorithm: Direct'
    end if
  else
    if (.not. DoCholesky) then
      write(u6,'(6X,A)') 'SCF Algorithm: Conventional USCF'
    else
      call Get_iScalar('System BitSwitch',iDoRI)
      if (iand(iDoRI,1024) == 1024) then
        if (LKon) then
          write(u6,'(6X,A)') 'SCF Algorithm: LK-RI/DF USCF'
          write(u6,FmtR) 'LK screening threshold:',dmpk
        else
          write(u6,'(6X,A)') 'SCF Algorithm: RI/DF USCF'
        end if
      else
        if (LKon) then
          write(u6,'(6X,A)') 'SCF Algorithm: LK-Cholesky USCF'
          write(u6,FmtR) 'LK screening threshold:',dmpk
        else
          write(u6,'(6X,A)') 'SCF Algorithm: Cholesky USCF'
        end if
      end if

      if (ALGO == 0) then
        if (iand(iDoRI,1024) == 1024) then
          write(u6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
        else
          write(u6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
        end if
      else if (ALGO == 1) then
        write(u6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
      else if (ALGO == 2) then
        write(u6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
      else if (ALGO == 3) then
        write(u6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
      else if (ALGO == 4) then
        write(u6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
      end if

      if (Do_DCCD) write(u6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
      if (ReOrd) write(u6,'(6X,A)') ' - the Cholesky vectors are reordered'
      if (DeCo) write(u6,'(6X,A)') ' - the density matrix is decomposed'
    end if
  end if
end if
if (NoExchange) write(u6,'(6X,A)') 'NOTE: exchange contributions will not be computed'

if (jPrint >= 2) then

  ! Print out informations concerning difference scheme used
  if (MiniDn) then
    write(u6,'(6X,A)') 'Minimized density differences are used'
  else
    if (.not. DDnOFF) then
      write(u6,'(6X,A)') 'D(i)-D(i-1) density differences are used'
    else
      write(u6,'(6X,A)') 'The actual AO density is used'
    end if
  end if
  write(u6,FmtI) 'Number of density matrices in core',nMem

  ! Print out number of iterations
  write(u6,FmtI) 'Maximum number of NDDO SCF iterations',nIter(0)
  write(u6,FmtI) 'Maximum number of HF SCF iterations',nIter(1)

  ! Print out thresholds for SCF
  write(u6,FmtR) 'Threshold for SCF energy change',EThr
  write(u6,FmtR) 'Threshold for density matrix',DThr
  write(u6,FmtR) 'Threshold for Fock matrix',FThr
  write(u6,FmtR) 'Threshold for linear dependence',DelThr
  if (Diis) then
    write(u6,FmtR) 'Threshold at which DIIS is turned on',DiisTh
    write(u6,FmtI) 'Maximum depth in the DIIS procedure',kOptim_Max
    write(u6,FmtI) 'Maximum depth in the BFGS Hessian update',IterSO_Max
    write(u6,FmtR) 'Threshold at which QNR/C2DIIS is turned on',QNRTh
    write(u6,FmtR) 'Threshold for Norm(delta) (QNR/C2DIIS)',DltNTh
  end if
  if (RGEK) then
    write(u6,'(6X,A)') 'RVO optimization with a subspace GEK surrogate model'
    select case(Expand)
      case (1)
        ExpMeth = 'DIIS'
      case (2)
        ExpMeth = 'BFGS'
      case (3)
        ExpMeth = 'RS-RFO'
      case default
        ExpMeth = 'Unknown'
    end select
    write(u6,'(6X,A,T50,A)') 'Subspace expansion method:',trim(ExpMeth)
  end if
  if (DSCF) write(u6,FmtR) 'Threshold for contribution to Fock matrix',SIntTh

  Frmt = '(6x,A,A)'

  ! Print out IVO information (if any)
  if (kIvo /= 0) write(u6,Frmt) 'Improved virtual orbitals.'

  ! Print out information about orbitals punched on the file
  if (jVOut <= 0) then
    write(u6,Frmt) 'No vectors punched'
  else if (jVOut == 1) then
    if (nD == 1) then
      write(u6,Frmt) 'All non deleted orbitals punched on: SCFORB'
    else
      write(u6,Frmt) 'All non deleted orbitals punched on: UHFORB'
    end if
  else
    if (nD == 1) then
      write(u6,Frmt) 'All orbitals punched on: SCFORB'
    else
      write(u6,Frmt) 'All orbitals punched on: UHFORB'
    end if
  end if
  call CollapseOutput(0,'   Optimization specifications:')
  write(u6,*)

  ! Print out
  if (InVec == 0) then
    write(u6,Frmt) 'Starting vectors from core diagonalization'
  else if (InVec == 1) then
    write(u6,Frmt) 'NDDO MOs are generated before actual HF SCF computation'
  else if (InVec == 2) then
    if (isHDF5) then
      write(u6,Frmt) 'Input vectors read from HDF5 file'
    else
      write(u6,Frmt) 'Input vectors read from INPORB'
      write(u6,Frmt) 'Orbital file label: ',trim(VTitle)
    end if
  else if (InVec == 3) then
    write(u6,Frmt) 'Input density matrix read from RUNFILE'
  else if (InVec == 4) then
    write(u6,Frmt) 'Restart...'
  else if (InVec == 5) then
    write(u6,Frmt) 'Input vectors from NDDO calculation'
  else
    write(u6,Frmt) StVec
  end if
  if (Scrmbl) write(u6,Frmt) 'Start orbitals are scrambled in order to introduce symmetry breaking'
  write(u6,*)

end if

call XFlush(u6)

end subroutine WrInp_SCF
