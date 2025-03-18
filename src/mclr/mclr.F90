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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine MCLR(ireturn)
!***********************************************************************
!                                                                      *
!               #     #  #####     #       ######                      *
!               ##   ## #     #    #       #     #                     *
!               # # # # #          #       #     #                     *
!               #  #  # #          #       ######                      *
!               #     # #          #       #   #                       *
!               #     # #     #    #       #    #                      *
!               #     #  #####     ####### #     #                     *
!                                                                      *
!     A linear response function program for general RASSCF states     *
!                                                                      *
!      OK right now we can just handle CASSCF, but the only thing      *
!      we need is yet another stupid PhD student who wants to code     *
!      the preconditioning elements between active and active orbitals *
!                                                                      *
!      The spin correction part of the code is a experimental          *
!      test platform, it is not finished and is not working,           *
!      If someone wants to help me finishing it please contact me      *
!                                                                      *
!                                                                      *
!   "Real programmers don't comment their codes.                       *
!    It was hard to write, and it must be hard to understand..."       *
!                                                                      *
!***********************************************************************

use Basis_Info, only: Basis_Info_Free
use Center_Info, only: Center_Info_Free
use External_Centers, only: External_Centers_Free
use Symmetry_Info, only: Symmetry_Info_Free
use Arrays, only: Hss, FAMO, FAMO_SpinP, FAMO_SpinM, SFock, G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m, CMO_Inv, CMO, Int1, pINT1, INT2, &
                  pINT2, G2t, G2sq, G1t, FIMO, F0SQMO
use Str_Info, only: DFTP, CFTP, DTOC, CNSM
use negpre, only: SS
use PDFT_Util, only: Do_Hybrid, WF_Ratio, PDFT_Ratio
! Added for CMS NACs
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6, wp
use MCLR_Data, only: nA, nNA, nAcPar, nAcPr2
use MCLR_Data, only: nrec
use MCLR_Data, only: iAllo
use MCLR_Data, only: SA, NACSTATES
use MCLR_Data, only: LuPT2
use DetDim, only: MXCNSM
use dmrginfo, only: DoDMRG, RGRAS2, DoMCLR
use input_mclr, only: ntAsh, ntAtri, ntASqr, nSym, iMethod, SpinPol, iMCPD, iMSPD, PT2, TimeDep, TwoStep, StepType, McKinley, &
                      RASSI, NewCho, Fail, double, LuAChoVec, LuChoInt, LuIChoVec, nAsh, nDisp, nRS2

implicit none
#include "warnings.h"
#include "SysDef.fh"
integer, allocatable :: ifpK(:), ifpS(:), ifpRHS(:), ifpCI(:), ifpSC(:), ifpRHSCI(:)
logical Reduce_Prt
external Reduce_Prt
integer get_MBl_wa
external get_MBl_wa
! first a few things for making Markus Happy
logical converged(8)
logical DoCholesky
! Additional things for CMS-NACs Optimization
character(len=8) Method
integer(kind=iwp) :: LuInput, istatus, LuSpool2
character(len=16) :: StdIn
character(len=180) :: Line
character(len=128) :: FileName
logical(kind=iwp) :: Exists
integer(kind=iwp), external :: isFreeUnit
logical :: CalcNAC_Opt = .false., MECI_via_SLAPAF = .false.
integer(kind=iwp) :: iPL, nSymX, iSym, nISP, I, iRC, iReturn
integer(kind=iwp) :: CMSNACStates(2)
integer(kind=iwp), external :: ipClose, iPrintLevel
real(kind=wp) :: TCPU1, TCPU2, TCPU3, TWall1, TWall2, TWall3

! This used to be after the CWTIME() functional call
!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = iPL-1
!                                                                      *
!***********************************************************************
!                                                                      *
! This is where you should put the information for CMS
call Get_cArray('Relax Method',Method,8)
if (Method == 'MSPDFT') then
  call Get_iArray('cmsNACstates',cmsNACstates,2)
  call Get_iArray('NACstatesOpt',NACstates,2)
  call Get_lscalar('CalcNAC_Opt',CalcNAC_Opt)
  call Get_lscalar('MECI_via_SLAPAF',MECI_via_SLAPAF)

  if (.not. MECI_via_SLAPAF) then
    if ((cmsNACstates(1) /= NACstates(1)) .or. (cmsNACstates(2) /= NACstates(2))) NACstates(:) = cmsNACstates(:)
  end if

  if ((cmsNACstates(1) /= NACstates(1)) .or. (cmsNACstates(2) /= NACstates(2))) then
    if (iPL >= 3) then
      write(u6,*)
      write(u6,*) 'MS-PDFT Potentials for root(s)'
      write(u6,*) cmsNACstates(1),cmsNACstates(2)
      write(u6,*) 'MCLR Lag. Mult. for root'
      write(u6,*) NACstates(1),NACstates(2)
      write(u6,*) 'MS-PDFT roots and MCLR roots do not match'
      write(u6,*) 'MCLR requests MCPDFT to be run before'
      write(u6,*) 'it states again'
      write(u6,*)
    end if

    LuInput = 11
    LuInput = IsFreeUnit(LuInput)
    call StdIn_Name(StdIn)
    call Molcas_open(LuInput,StdIn)

    write(LuInput,'(A)') '>ECHO OFF'
    write(LuInput,'(A)') '>export MCLR_OLD_TRAP=$MOLCAS_TRAP'
    write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

    write(LuInput,'(A)') ' &MCPDFT &END'
    write(LuInput,'(A)') ' KSDFT=T:PBE'
    write(LuInput,'(A)') ' MSPDft'
    write(LuInput,'(A)') ' GRAD'
    if (NACstates(2) /= 0) then
      write(LuInput,'(A)') 'NAC'
      write(LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
      if (CalcNAC_Opt) write(LuInput,'(A)') 'MECI'
    end if
    write(LuInput,'(A)') 'End of Input'
    write(LuInput,'(A)') ''

    FileName = 'MCLRINP'
    call f_inquire(Filename,Exists)

    if (Exists) then
      LuSpool2 = 77
      LuSpool2 = IsFreeUnit(LuSpool2)
      call Molcas_Open(LuSpool2,Filename)
      do
        read(LuSpool2,'(A)',iostat=istatus) Line
        if (istatus > 0) call Abend()
        if (istatus < 0) exit
        write(LuInput,'(A)') Line
      end do
      close(LuSpool2)
    else
      write(LuInput,'(A)') ' &MCLR &End'
    end if

    write(LuInput,'(A)') '>export MOLCAS_TRAP=$MCLR_OLD_TRAP'
    write(LuInput,'(A)') '>ECHO ON'
    close(LuInput)
    call Finish(_RC_INVOKED_OTHER_MODULE_)
  end if
end if

! check the status; if quantities needed for MCLR have not been
! computed, call CASPT2
if (Method == 'CASPT2') call check_caspt2(0)
!                                                                      *
!***********************************************************************
!                                                                      *
!call MCLR_banner()
call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
iAllo = 0
!idp = rtoi
nrec = get_MBl_wa()/rtob

call DecideOnCholesky(DoCholesky)
call get_iScalar('nSym',nSymX)

if (DoCholesky .and. (nSymX > 1)) then
  write(6,*) '** Cholesky or RI/DF not implemented with symmetry **'
  call Quit(_RC_INPUT_ERROR_)
end if

call Init_Data()

doDMRG = .false.
doMCLR = .false.
!#ifdef _DMRG_
! read dmrg parameters for mclr part
!call read_dmrg_parameter_for_mclr()
!if (doDMRG) then
!  doMCLR = .true.
!  write(6,*) 'ndets_RGLR : ',ndets_RGLR
!  write(6,*) 'nstates_RGLR ',nstates_RGLR
!  write(6,*) 'RGras2 : ',RGras2
!  write(6,*) 'LRras2 : ',LRras2
!  open(unit=117,file='mclr_dets.initial')
!end if
!#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! open files

call OpnFls_MCLR(iPL)
call IpInit()
!                                                                      *
!***********************************************************************
!                                                                      *
! Read input

call InpCtl_MCLR(iPL)
!                                                                      *
!***********************************************************************
!                                                                      *
! Transform integrals and calculate fock matrixes etc.

if (doDMRG) then  ! yma
  call dmrg_spc_change_mclr(RGras2(1:8),nAsh)
  call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
end if
ntAsh = 0
ntAtri = 0
ntAsqr = 0
nnA = 0
do iSym=1,nSym
  ntAsh = ntAsh+nAsh(iSym)
  ntAtri = ntAtri+nAsh(iSym)*(nAsh(iSym)+1)/2
  ntAsqr = ntAsqr+nAsh(iSym)*nAsh(iSym)
  nA(iSym) = nna
  nnA = nnA+nAsh(isym)
end do
nacpar = (nnA+1)*nnA/2
nacpr2 = (nacpar+1)*nacpar/2

call Start_MCLR()

nisp = max(8,nDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
! File pointers

call mma_allocate(ifpK,nisp,Label='ifpK')
ifpK(:) = -1
call mma_allocate(ifpS,nisp,Label='ifpS')
ifpS(:) = -1
call mma_allocate(ifpRHS,nisp,'ifpRHS')
ifpRHS(:) = -1
if (iMethod == 2) then
  call mma_allocate(ifpCI,nisp,Label='ifpCI')
  call mma_allocate(ifpSC,nisp,Label='ifpSC')
  call mma_allocate(ifpRHSCI,nisp,Label='ifpRHSCI')
else
  call mma_allocate(ifpCI,1,Label='ifpCI')
  call mma_allocate(ifpSC,1,Label='ifpSC')
  call mma_allocate(ifpRHSCI,1,Label='ifpRHSCI')
end if
ifpCI(:) = -1
ifpSC(:) = -1
ifpRHSCI(:) = -1
!                                                                      *
!***********************************************************************
!                                                                      *
! Calculate response

! Output is stored on disk

if (SPINPOL) then
  call WfCtl_SP(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI)
!else if (ELECHESS) then
!  call WfCtl_PCG(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI)
!  call Abend()
else if (iMCPD) then !pdft

  Do_Hybrid = .false.
  call qpg_DScalar('R_WF_HMC',Do_Hybrid)
  if (Do_Hybrid) then
    call Get_DScalar('R_WF_HMC',WF_Ratio)
    PDFT_Ratio = 1.0d0-WF_Ratio
  end if

  if (iMSPD) then
    if (Do_Hybrid) then
      call WarningMessage(2,'Hybrid MS-PDFT gradient not supported yet')
      call Quit(_RC_EXIT_EXPECTED_)
    end if
    call WfCtl_MSPD(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
  else
    call WfCtl_PDFT(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
  end if
else if (SA .or. PT2) then
  call WfCtl_SA(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
else if (TimeDep) then
  call WfCtl_td(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,converged)
else
  call WfCtl_Hess(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,converged)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Contract response to hessian etc

if (.not. (TwoStep .and. (StepType == 'RUN1'))) then
  if (PT2 .or. SA .or. iMCPD) then
    call Out_PT2(ifpK,ifpCI)
    if (PT2) close(LUPT2) !! this file is opend in wfctl_sa
  else if (TimeDep) then
    call Output_td(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,converged)
  else
    call Output_mclr(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,converged)
    if (mckinley) call isoloop(double)
  end if

  if (RASSI) call OutRAS(ifpK,ifpCI)

  if (TimeDep) call OutRAS_td(ifpK,ifpCI)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Basis_Info_Free()
call Center_Info_Free()
call Symmetry_Info_Free()
call External_Centers_Free()
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

! Arrays in csf.f
if (iMethod == 2) then
  call mma_deallocate(DTOC)
  call mma_deallocate(CFTP)
  call mma_deallocate(DFTP)
end if
do i=1,MXCNSM
  call mma_deallocate(CNSM(i)%ICONF,safe='*')
  call mma_deallocate(CNSM(i)%ICTS,safe='*')
end do

! Free arrays allocated by memstr.f
if (iMethod == 2) call FreeStr()
! Arrays in inpone.f
call mma_deallocate(INT1)
! Arrays in detctl.f
if (iMethod == 2) then
  call mma_deallocate(pINT2)
  call mma_deallocate(pINT1)
end if

! Array in rdab.f
if (iMethod == 1) call mma_deallocate(CMO)
! Arrays in rdjobip_td.f and rdjobiph.f
if (iMethod == 2) then
  call mma_deallocate(G2t)
  if (Timedep) call mma_deallocate(G2sq)
  call mma_deallocate(G1t)
  call mma_deallocate(CMO)
end if

! Arrays allocated in stpert.f
call mma_deallocate(Hss)
if (SPINPOL) then
  call mma_deallocate(FAMO_SpinP)
  call mma_deallocate(FAMO_SpinM)
  call mma_deallocate(G2mm)
  call mma_deallocate(G2mp)
  call mma_deallocate(G2pp)
  call mma_deallocate(Fp)
  call mma_deallocate(Fm)
  call mma_deallocate(G1p)
  call mma_deallocate(G1m)
  call mma_deallocate(SFock)
end if
! Arrays allocated in fckmat.f
call mma_deallocate(FAMO)
call mma_deallocate(Int2)
call mma_deallocate(FIMO)
call mma_deallocate(F0SQMO)

call mma_deallocate(ifpRHSCI)
call mma_deallocate(ifpSC)
call mma_deallocate(ifpCI)
call mma_deallocate(ifpRHS)
call mma_deallocate(ifpS)
call mma_deallocate(ifpK)
call mma_deallocate(SS,safe='*')

! Close files

call ClsFls_MCLR()

if (NewCho) then
  call Cho_X_Final(irc)
  do i=1,nSym
    call DACLOS(LuAChoVec(i))
    call DACLOS(LuIChoVec(i))
  end do
  call DACLOS(LuChoInt(1))
  call DACLOS(LuChoInt(2))
  call mma_deallocate(CMO_Inv)
end if

if (TwoStep .and. (StepType == 'RUN1')) irc = ipclose(-1)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. fail) then
  if (iPL >= 2) then
    write(6,*)
    write(6,'(6X,A)') 'The response parameters are written to the file RESP.'
  end if
  ireturn = _RC_ALL_IS_WELL_
else
  ireturn = _RC_NOT_CONVERGED_
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu3,TWall3)
if (iPL >= 3) then
  write(6,*)
  write(6,'(2X,A)') 'Timings'
  write(6,'(2X,A)') '-------'
  write(6,*)
  write(6,'(2X,A)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,'(2X,A,T44,A,A,A)') ' ',' ','    CPU time','     elapsed'
  write(6,'(2X,A,T44,A,2F12.2)') '1) Initialization',':',TCpu2-TCpu1,TWall2-TWall1
  write(6,'(2X,A,T44,A,2F12.2)') '2) Response calculation',':',TCpu3-TCpu2,TWall3-TWall2
  write(6,'(2X,A)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
  write(6,'(2X,A,T44,A,2F12.2)') 'Total',':',TCpu3-TCpu1,TWall3-TWall1
  write(6,'(2X,A)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPL >= 3) call FastIO('STATUS')
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine MCLR
