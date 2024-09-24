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
use InfSO, only: DltNth, QNRTh, IterSO_Max
use InfSCF, only: Aufb, DDnoff, DelThr, DIIS, DIISTh, DoCholesky, DSCF, DThr, EThr, FThr, iAU_ab, InVec, isHDF5, nD, jPrint, &
                  jVOut, kIVO, kOptim_Max, KSDFT, LKOn, lpaper, MiniDn, nCore, nDIsc, nMem, NoExchange, nSym, nTit, One_Grid, &
                  PreSch, RFPert, rTemp, Scrmbl, StVec, Teee, TemFac, Thize, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge, TStop, &
                  VTitle, Header, Title, nFro, nAufb, nOcc, nOrb, nBas, nIter, nDel
use ChoSCF, only: dmpk, Algo, ReOrd
use Fock_util_global, only: Deco
use RICD_Info, only: Do_DCCD

implicit none
real*8 SIntTh
! Define local variables
integer i, iCharge, iDoRI, iSym, iTit, mTmp
character(len=60) Fmt, FmtR, FmtI
character(len=72) Line
character(len=3) lIrrep(8)
logical NonEq

if (jPrint >= 2) then
  call CollapseOutput(1,'   Input section:')
  write(6,'(3X,A)') '   --------------'
  write(6,*)
end if

! Print out header of the integral file
if (jPrint >= 2) then
  write(6,'(6X,A)') 'Header of the integral files:'
  write(Line,'(A72)') Header(1)
  write(6,'(6X,A)') trim(Line)
  write(Line,'(A72)') Header(2)
  write(6,'(6X,A)') trim(Line)
  write(6,*)
end if

! Print out title (if any)
if (nTit > 0) then
  if (jPrint >= 3) then
    call Banner(Title,nTit,lPaper-7)
    write(6,*)
  else if (jPrint >= 2) then
    write(6,'(6X,A)') ' Title:'
    do iTit=1,nTit
      write(6,'(8X,A)') trim(Title(iTit))
    end do
  end if
end if

! Print the coordinates of the system
if (jPrint >= 2) call PrCoor()

! Print out Print level
!write(6,'(1X,A,I3)') 'Print level=',iPrint
!write(6,*)

call Get_cArray('Irreps',lIrrep,24)
do iSym=1,nSym
  lIrrep(iSym) = adjustr(lIrrep(iSym))
end do

if (jPrint >= 2) then
  call CollapseOutput(0,'   Input section:')
  write(6,*)
end if
! Print out orbital informations
Fmt = '(6X,A,T35,8I4)'
if (jPrint >= 2) then
  call CollapseOutput(1,'   Orbital specifications:')
  write(6,'(3X,A)') '   -----------------------'
  write(6,*)
  write(6,Fmt) 'Symmetry species',(i,i=1,nSym)
  write(6,'(6X,A,T35,8(1X,A))') '                ',(lIrrep(i),i=1,nSym)
  write(6,Fmt) 'Frozen orbitals',(nFro(i),i=1,nSym)
end if

if (Aufb) then
  if (nD == 1) then
    if (nAufb(1) == -1) then
      Tot_El_Charge = Tot_Charge-Tot_Nuc_Charge
      ! Check that Tot_El_Charge is a multiple of two!
      mtmp = int(-Tot_El_Charge+0.1d0)
      if (mod(mtmp,2) /= 0) then
        write(6,*) 'WrInp: Error in definition of molecular charge!'
        write(6,*) 'Current implementation only allows double occupations.'
        write(6,*) 'Tot_Charge    :',Tot_Charge
        write(6,*) 'Tot_El_Charge :',Tot_El_Charge
        write(6,*) 'Tot_Nuc_Charge:',Tot_Nuc_Charge
        call Abend()
      end if
      nAufb(1) = mtmp/2
    end if
  else
    if (nAufb(1) == -1) then
      Tot_El_Charge = Tot_Charge-Tot_Nuc_Charge
      ! Check that Tot_El_Charge is a multiple of two!
      mtmp = int(-Tot_El_Charge+0.1d0)
      ! if ZSPIN is not set - make difference alpha-beta = 0 or 1
      nAufb(2) = (mtmp-iAu_ab)/2
      nAufb(1) = int(-Tot_El_Charge-nAufb(2))
    end if
    !write(6,*) ' CHARGE + UHF is un'
    !call Abend()
  end if
  if ((nD == 1) .and. (jPrint >= 2)) then
    write(6,Fmt) 'Aufbau',nAufb(1)
  else if (jPrint >= 3) then
    write(6,Fmt) 'Aufbau alpha',nAufb(1)
    write(6,Fmt) 'Aufbau beta ',nAufb(2)
  end if
  if (Teee .and. (jPrint >= 2)) then
    write(6,'(a,f6.3)') '      Start temperature =',RTemp
    write(6,'(a,f6.3)') '      End temperature   =',TStop
    write(6,'(a,f6.3)') '      Temperature Factor=',TemFac
  end if
else
  if ((nD == 1) .and. (jPrint >= 2)) then
    write(6,Fmt) 'Occupied orbitals',(nOcc(i,1),i=1,nSym)
    write(6,Fmt) 'Secondary orbitals',(nOrb(i)-nOcc(i,1),i=1,nSym)
  else if (jPrint >= 2) then
    write(6,Fmt) 'Occupied orbitals alpha',(nOcc(i,1),i=1,nSym)
    write(6,Fmt) 'Occupied orbitals beta ',(nOcc(i,2),i=1,nSym)
    write(6,Fmt) 'Secondary orbitals alpha',(nOrb(i)-nOcc(i,1),i=1,nSym)
    write(6,Fmt) 'Secondary orbitals beta',(nOrb(i)-nOcc(i,2),i=1,nSym)

  end if
end if

Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
call put_dscalar('Total Charge    ',Tot_Charge)
if (jPrint >= 2) then
  write(6,Fmt) 'Deleted orbitals         ',(nDel(i),i=1,nSym)
  write(6,Fmt) 'Total number of orbitals ',(nOrb(i),i=1,nSym)
  write(6,Fmt) 'Number of basis functions',(nBas(i),i=1,nSym)
  call CollapseOutput(0,'   Orbital specifications:')
  write(6,*)
  write(6,'(6X,A,T45,F10.3)') 'Molecular charge',Tot_Charge
  write(6,*)
end if

! Print out reaction field specifications

iCharge = int(Tot_Charge)
NonEq = .false.
call PrRF(.false.,NonEq,iCharge,jPrint)
if (RFpert .and. (jPrint >= 2)) then
  write(6,'(6X,A)') 'Reaction field specifications:'
  write(6,'(6X,A)') '------------------------------'
  write(6,*)
  write(6,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previous calculation'
  write(6,*)
end if

! Print out grid information in case of DFT

if (KSDFT /= 'SCF') then
  call Put_dScalar('DFT exch coeff',CoefX)
  call Put_dScalar('DFT corr coeff',CoefR)
  call Put_dScalar('EThr',EThr)
  call Funi_Print()
  if (jPrint >= 2) then
    if (One_Grid) then
      write(6,'(6X,A)') 'The same grid will be used for all iterations.'
    else
      write(6,'(6X,A)') 'A smaller intermediate grid will be used the first few iterations.'
    end if
    write(6,*)
    write(6,'(6X,A)') 'DFT functional specifications'
    write(6,'(6X,A)') '-----------------------------'
    call libxc_version()
    call Print_Info()
    write(6,*)
  end if
end if

! Print out informations concerning Direct/Conventional scheme
FmtI = '(6X,A,T50,I9)'
FmtR = '(6X,A,T50,ES9.2)'
if (jPrint >= 2) then
  call CollapseOutput(1,'   Optimization specifications:')
  write(6,'(3X,A)') '   ----------------------------'
  write(6,*)
end if
if (DSCF .and. (.not. Do_DCCD) .and. (jPrint >= 2)) then
  if ((nDisc == 0) .and. (nCore == 0)) then
    write(6,'(6X,A)') 'SCF Algorithm: Direct'
  else if (nDisc*1024 >= nCore) then
    write(6,'(6X,A)') 'SCF Algorithm: Semi-direct'

    ! The threshold to be used in the direct SCF procedure is
    ! defined by the energy threshold.
    write(6,FmtI) 'Max MByte of integrals on disk/process:',nDisc
    write(6,FmtR) 'Threshold for saving integrals on disc',Thize
  else
    write(6,'(6X,A)') 'SCF Algorithm: Semi-direct in-core'

    ! The threshold to be used in the direct SCF procedure is
    ! defined by the energy threshold.
    write(6,FmtI) 'Max kByte of integrals in memory/process:',nCore
    write(6,FmtR) 'Threshold for saving integrals in memory',Thize
  end if
  if (PreSch) then
    write(6,'(6X,A)') 'Prescreening Scheme: Only Integral value'
  else
    write(6,'(6X,A)') 'Prescreening Scheme: Integral*Density value'
  end if
else if (jPrint >= 2) then
  if (nD == 1) then
    if (.not. DoCholesky) then
      write(6,'(6X,A)') 'SCF Algorithm: Conventional'
    else
      call Get_iScalar('System BitSwitch',iDoRI)
      if (iand(iDoRI,1024) == 1024) then
        if (LKon) then
          write(6,'(6X,A)') 'SCF Algorithm: LK-RI/DF'
          write(6,FmtR) 'LK screening threshold:',dmpk
        else
          write(6,'(6X,A)') 'SCF Algorithm: RI/DF'
        end if
      else
        if (LKon) then
          write(6,'(6X,A)') 'SCF Algorithm: LK-Cholesky'
          write(6,FmtR) 'LK screening threshold:',dmpk
        else
          write(6,'(6X,A,I1)') 'SCF Algorithm: Cholesky'
        end if
      end if

      if (ALGO == 0) then
        if (iand(iDoRI,1024) == 1024) then
          write(6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
        else
          write(6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
        end if
      else if (ALGO == 1) then
        write(6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
      else if (ALGO == 2) then
        write(6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
      else if (ALGO == 3) then
        write(6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
      else if (ALGO == 4) then
        write(6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
      end if

      if (Do_DCCD) write(6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
      if (ReOrd) write(6,'(6X,A)') ' - the Cholesky vectors are reordered'
      if (DeCo) write(6,'(6X,A)') ' - the density matrix is decomposed'
      if (DSCF .and. (nDisc == 0) .and. (nCore == 0)) write(6,'(6X,A)') ' - SCF Algorithm: Direct'
    end if
  else
    if (.not. DoCholesky) then
      write(6,'(6X,A)') 'SCF Algorithm: Conventional USCF'
    else
      call Get_iScalar('System BitSwitch',iDoRI)
      if (iand(iDoRI,1024) == 1024) then
        if (LKon) then
          write(6,'(6X,A)') 'SCF Algorithm: LK-RI/DF USCF'
          write(6,FmtR) 'LK screening threshold:',dmpk
        else
          write(6,'(6X,A)') 'SCF Algorithm: RI/DF USCF'
        end if
      else
        if (LKon) then
          write(6,'(6X,A)') 'SCF Algorithm: LK-Cholesky USCF'
          write(6,FmtR) 'LK screening threshold:',dmpk
        else
          write(6,'(6X,A)') 'SCF Algorithm: Cholesky USCF'
        end if
      end if

      if (ALGO == 0) then
        if (iand(iDoRI,1024) == 1024) then
          write(6,'(6X,A)') 'Integral regeneration from RI vectors reordered on disk'
        else
          write(6,'(6X,A)') 'Integral regeneration from Cholesky vectors reordered on disk'
        end if
      else if (ALGO == 1) then
        write(6,'(6X,A)') 'Density-based Cholesky. Default reorder: on the fly'
      else if (ALGO == 2) then
        write(6,'(6X,A)') 'MO-based-Exchange Cholesky. Default reorder: on the fly'
      else if (ALGO == 3) then
        write(6,'(6X,A)') 'MO-based-Exchange Cholesky. MO-transformation in reduced sets'
      else if (ALGO == 4) then
        write(6,'(6X,A)') 'Local-Exchange (LK) algorithm.'
      end if

      if (Do_DCCD) write(6,'(6X,A)') ' - Corrected with exact 1-center two-electron integrals'
      if (ReOrd) write(6,'(6X,A)') ' - the Cholesky vectors are reordered'
      if (DeCo) write(6,'(6X,A)') ' - the density matrix is decomposed'
    end if
  end if
end if
if (NoExchange) write(6,'(6X,A)') 'NOTE: exchange contributions will not be computed'

if (jPrint >= 2) then

  ! Print out informations concerning difference scheme used
  if (MiniDn) then
    write(6,'(6X,A)') 'Minimized density differences are used'
  else
    if (.not. DDnOFF) then
      write(6,'(6X,A)') 'D(i)-D(i-1) density differences are used'
    else
      write(6,'(6X,A)') 'The actual AO density is used'
    end if
  end if
  write(6,FmtI) 'Number of density matrices in core',nMem

  ! Print out number of iterations
  write(6,FmtI) 'Maximum number of NDDO SCF iterations',nIter(0)
  write(6,FmtI) 'Maximum number of HF SCF iterations',nIter(1)

  ! Print out thresholds for SCF
  write(6,FmtR) 'Threshold for SCF energy change',EThr
  write(6,FmtR) 'Threshold for density matrix',DThr
  write(6,FmtR) 'Threshold for Fock matrix',FThr
  write(6,FmtR) 'Threshold for linear dependence',DelThr
  if (Diis) then
    write(6,FmtR) 'Threshold at which DIIS is turned on',DiisTh
    write(6,FmtI) 'Maximum depth in the DIIS procedure',kOptim_Max
    write(6,FmtI) 'Maximum depth in the BFGS Hessian update',IterSO_Max
    write(6,FmtR) 'Threshold at which QNR/C2DIIS is turned on',QNRTh
    write(6,FmtR) 'Threshold for Norm(delta) (QNR/C2DIIS)',DltNTh
  end if
  if (DSCF) write(6,FmtR) 'Threshold for contribution to Fock matrix',SIntTh

  Fmt = '(6x,A,A)'

  ! Print out IVO information (if any)
  if (kIvo /= 0) write(6,Fmt) 'Improved virtual orbitals.'

  ! Print out information about orbitals punched on the file
  if (jVOut <= 0) then
    write(6,Fmt) 'No vectors punched'
  else if (jVOut == 1) then
    if (nD == 1) then
      write(6,Fmt) 'All non deleted orbitals punched on: SCFORB'
    else
      write(6,Fmt) 'All non deleted orbitals punched on: UHFORB'
    end if
  else
    if (nD == 1) then
      write(6,Fmt) 'All orbitals punched on: SCFORB'
    else
      write(6,Fmt) 'All orbitals punched on: UHFORB'
    end if
  end if
  call CollapseOutput(0,'   Optimization specifications:')
  write(6,*)

  ! Print out
  if (InVec == 0) then
    write(6,Fmt) 'Starting vectors from core diagonalization'
  else if (InVec == 1) then
    write(6,Fmt) 'NDDO MOs are generated before actual HF SCF computation'
  else if (InVec == 2) then
    if (isHDF5) then
      write(6,Fmt) 'Input vectors read from HDF5 file'
    else
      write(6,Fmt) 'Input vectors read from INPORB'
      write(6,Fmt) 'Orbital file label: ',trim(VTitle)
    end if
  else if (InVec == 3) then
    write(6,Fmt) 'Input density matrix read from RUNFILE'
  else if (InVec == 4) then
    write(6,Fmt) 'Restart...'
  else if (InVec == 5) then
    write(6,Fmt) 'Input vectors from NDDO calculation'
  else
    write(6,Fmt) StVec
  end if
  if (Scrmbl) write(6,Fmt) 'Start orbitals are scrambled in order to introduce symmetry breaking'
  write(6,*)

end if

call XFlush(6)

end subroutine WrInp_SCF
