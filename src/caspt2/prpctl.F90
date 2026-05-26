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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine PRPCTL(MODE,UEFF,U0,nState)

use Index_Functions, only: nTri_Elem
use PT2WFN, only: PT2WFN_DENSSTORE
use OneDat, only: sNoNuc, sNoOri
use PrintLevel, only: USUAL, VERBOSE
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use caspt2_global, only: do_grad
#endif
use EQSOLV, only: IVECX, NLSTOT
use caspt2_global, only: CMO, CMO_Internal, CMOPT2, do_nac, DPT2_tot, DPT2C_tot, iPrGlb, iRoot1, iRoot2, LISTS, LUONEM, NCMO, &
                         SLag, TORB
use caspt2_module, only: BNAME, Energy, IAD1M, IFMSCOUP, IFPROP, irlxroot, ISCF, JSTATE, MSTATE, MSTATE, NASH, NASHT, NBAS, NBAST, &
                         NCONF, NDEL, NFRO, NISH, NORB, NRAS1, NRAS2, NRAS3, NSYM, OUTFMT, PRORB, THRENE, THROCC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Five, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Mode, nState
real(kind=wp), intent(in) :: UEFF(NSTATE,NSTATE), U0(nState,nState)
integer(kind=iwp) :: iComp, IDISK, IDMAT, IDMOFF, IERR, II, II2, IJ, IJ2, IndType(7,8), iOpt, iRc, ISTATE, iSyLbl, ISYM, iUHF, &
                     KSTATE, LUTMP, N, nDens, NDMAT, NFROSAV(NSYM), NO, NOCC, NORBSAV(NSYM)
real(kind=wp) :: Dummy(2), DUM(1), SCAL
logical(kind=iwp) :: Do_ESPF, FullMlk, lSave
character(len=128) :: FILENAME, MDNAME
character(len=80) :: Note
character(len=8) :: Label
real(kind=wp), allocatable :: CI1(:), CI2(:), CNAT(:), DMAT(:), OCC(:), Scr(:), SGM(:), TG1(:,:)
integer(kind=iwp), external :: IsFreeUnit

if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,'(A)') repeat('*',80)
  write(u6,*) ' CASPT2 PROPERTY SECTION'
end if

#ifdef _MOLCAS_MPP_
if (Is_Real_Par() .and. (IPRGLB >= USUAL) .and. (.not. do_grad)) then
  write(u6,'(1X,A)') ' ====================================='
  write(u6,'(1X,A)') ' CASPT2 properties were requested, but'
  write(u6,'(1X,A)') ' these are not efficiently implemented'
  write(u6,'(1X,A)') ' in parallel. If you do not need them,'
  write(u6,'(1X,A)') ' not using the PROP keyword could lead'
  write(u6,'(1X,A)') ' to a significant speed up.'
  write(u6,'(1X,A)') ' ====================================='
end if
#endif

! PAM2008 When this subroutine is called, theMSTATE calculation has been done
! for the (individual) state nr JSTATE in 1,2,..,NSTATE.
! The corresponding CI-root from rasscf, the root state for this PT2,
! is number MSTATE(JSTATE) on the input JOBIPH file.
! JSTATE,NSTATE and MSTATE() are in the caspt2_module file.
IERR = 0
if (NSTATE > 1) then
  N = MSTATE(JSTATE)
  if ((N <= 0) .or. (N > 999)) then
    write(u6,*) ' Subroutine PRPCTL fails -- It seems to get data'
    write(u6,*) ' computed for a root nr ',N
    write(u6,*) ' which is surely wrong.'
    write(u6,*) ' PRPCTL gives up, there will be no calculations'
    write(u6,*) ' done of orbitals, properties, etc for this state.'
    write(u6,*) ' This was state nr JSTATE=',JSTATE
    write(u6,*) ' in the MS-CASPT2 calculation.'
    IERR = 1
  end if
end if
if (IERR > 0) return

FullMlk = .true.
if (.not. PRORB) FullMlk = .false.

if (IPRGLB >= USUAL) write(u6,'(A)') repeat('-',80)

! Compute density matrix, output orbitals, and properties.

! Compute density matrix of CASPT2 wave function, in MO basis,
! to produce output orbitals.
! This density matrix may be approximated in several ways, see DENS.
NDMAT = 0
NOCC = sum(NBAS(1:NSYM))
if (MODE == 0) then
  !! Density matrix for each state (single-state)
  NDMAT = sum(nTri_Elem(NORB(1:NSYM)))
  call mma_allocate(DMAT,NDMAT,Label='DMAT')
  DMAT(:) = Zero
  call mma_allocate(LISTS,NLSTOT,LABEL='LISTS')
  call MKLIST(LISTS,NLSTOT)
  call DENS(IVECX,NDMAT,NSTATE,DMAT,UEFF,U0)
  call mma_deallocate(LISTS)
else
  !! Density matrix for the target adiabatic state
  !! MODE = 1 is called after gradient stuff only for MS, so
  !! the density computed here is what we call a correct unrelaxed
  !! correlated CASPT2 density
  !! Called after gradient calculations, only, from GrdCls
  NDMAT = sum(nTri_Elem(NBAS(1:NSYM)))
  call mma_allocate(DMAT,NDMAT,Label='DMAT')
  DMAT(:) = Zero
  !! Copy the unrelaxed density matrix to triangular
  !! The basis of DPT2_tot is natural (CASSCF)
  IDMAT = 0
  IDMOFF = 0
  do ISYM=1,NSYM
    NO = NBAS(ISYM)
    do II=1,NO
      do IJ=1,II
        !! second-order (DPT2) and first-order (DPT2C)
        DMAT(1+IDMAT) = DPT2_TOT(IDMOFF+II+NO*(IJ-1))+DPT2C_TOT(IDMOFF+II+NO*(IJ-1))*Quart
        if (.not. DO_NAC) then
          !! Add the reference density matrix (inactive)
          if ((II == IJ) .and. (II <= NFRO(ISYM)+NISH(ISYM))) DMAT(1+IDMAT) = DMAT(1+IDMAT)+Two
        end if
        IDMAT = IDMAT+1
      end do
    end do
    IDMOFF = IDMOFF+NO*NO
  end do
  !! Add the reference density matrix (active)
  call mma_allocate(CI1,NCONF,Label='CI1')
  call mma_allocate(CI2,NCONF,Label='CI2')
  call mma_allocate(SGM,NCONF,Label='SGM')
  call mma_allocate(TG1,NASHT,NASHT,Label='TG1')
  do ISTATE=1,NSTATE
    if (ISCF /= 0) then
      CI1(1) = One
    else
      call LOADCI_XMS('N',1,NCONF,NSTATE,CI1,ISTATE,U0)
    end if
    do KSTATE=1,NSTATE
      SCAL = SLag(ISTATE,KSTATE)
      if (.not. DO_NAC) then
        if ((ISTATE == IROOT1) .and. (KSTATE == IROOT2)) SCAL = SCAL+One
      end if
      if (abs(SCAL) <= 1.0e-9_wp) cycle
      if (ISCF /= 0) then
        CI2(1) = One
      else
        call LOADCI_XMS('N',1,NCONF,NSTATE,CI2,KSTATE,U0)
      end if
      call Dens1T_RPT2(CI1,CI2,SGM,TG1,nAshT)
      TG1(:,:) = SCAL*TG1(:,:)
      do II=1,NASH(1)
        II2 = II+NFRO(1)+NISH(1)
        IJ2 = nTri_Elem(II2-1)+NFRO(1)+NISH(1)
        DMAT(IJ2+1:IJ2+II) = DMAT(IJ2+1:IJ2+II)+(TG1(II,1:II)+TG1(1:II,II))*Half
      end do
    end do
  end do
  call mma_deallocate(CI1)
  call mma_deallocate(CI2)
  call mma_deallocate(SGM)
  call mma_deallocate(TG1)

  !! The density matrix above has frozen orbitals.
  !! Accordingly, the natural orbitals should be generated
  !! by diagonalizing all orbitals (including frozen), so
  !! make the number of frozen orbitals zero for the moment.
  NFROSAV(1:NSYM) = NFRO(1:NSYM)
  NORBSAV(1:NSYM) = NORB(1:NSYM)
  NFRO(1:NSYM) = 0
  NORB(1:NSYM) = NBAS(1:NSYM)
end if

! Compute natural orbitals of CASPT2 wave function.
call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
CMO => CMO_Internal
CMO(:) = CMOPT2(:)
if (MODE == 1) then
  !! Read the CASSCF orbital (really?)
  IDISK = IAD1M(1)
  call ddafile(LUONEM,2,CMO,NCMO,IDISK)
end if
call mma_allocate(CNAT,NCMO,Label='CNAT')
call mma_allocate(OCC,NOCC,Label='OCC')
call NATORB_CASPT2(DMAT,nDMAT,CMO,nCMO,OCC,nOCC,CNAT,nCMO)
call mma_deallocate(CMO_Internal)
nullify(CMO)
! Backtransform density matrix to original MO basis before storing
call TRANSFOCK(TORB,size(TORB),DMAT,NDMAT,-1)
call PT2WFN_DENSSTORE(DMAT,NDMAT)
call mma_deallocate(DMAT)
if (MODE == 1) then
  !! Restore with frozen orbitals
  NFRO(1:NSYM) = NFROSAV(1:NSYM)
  NORB(1:NSYM) = NORBSAV(1:NSYM)
end if

if (IFMSCOUP .and. (MODE == 0) .and. (.not. IFPROP)) then
  !! Do not show properties, if PROP calculation is not requested
  call mma_deallocate(CNAT)
  call mma_deallocate(OCC)
  return
end if

! Write natural orbitals as standard orbital file on PT2ORB
! PAM2008: Separate PT2ORB files for each state:
FILENAME = 'PT2ORB'
MDNAME = 'MD_PT2'
if ((NSTATE > 1) .and. (MODE == 0)) then
  FILENAME = 'PT2ORB.x'
  MDNAME = 'MD_PT2.x'
  N = MSTATE(JSTATE)
  if (N <= 9) then
    write(FILENAME(8:8),'(I1)') N
    write(MDNAME(8:8),'(I1)') N
  else if (N <= 99) then
    write(FILENAME(8:9),'(I2)') N
    write(MDNAME(8:9),'(I2)') N
  else if (N <= 999) then
    write(FILENAME(8:10),'(I3)') N
    write(MDNAME(8:10),'(I3)') N
  end if
end if
! PAM2008: For MS-CASPT2 with more than one state, orbital files will
! now be numbered PT2ORB.1 ... PT2ORB.999
! depending on which CI-root of the rasscf calculation that was
! the root state of the PT2.
LUTMP = 19
LUTMP = ISFREEUNIT(LUTMP)
! PAM 2008: Add typeindex information:
!----------------------------------------------------------------------*
!     Make typeindex information                                       *
!----------------------------------------------------------------------*
IndType(1,1:NSYM) = NFRO(1:NSYM)
IndType(2,1:NSYM) = NISH(1:NSYM)
IndType(3,1:NSYM) = NRAS1(1:NSYM)
IndType(4,1:NSYM) = NRAS2(1:NSYM)
IndType(5,1:NSYM) = NRAS3(1:NSYM)
IndType(6,1:NSYM) = NDEL(1:NSYM)
IndType(7,1:NSYM) = NBAS(1:NSYM)-sum(IndType(1:6,1:NSYM),dim=1)
if (NSTATE > 1) then
  write(Note,'(A41,I3,A3,f22.12)') '* CASPT2 natural orbitals for root number',N,' E=',Energy(JSTATE)
else
  Note = '* CASPT2 natural orbitals'
end if

call WRVEC(FILENAME,LUTMP,'COI',NSYM,NBAS,NBAS,CNAT,OCC,Dummy,IndType,Note)
iUHF = 0
call Molden_Interface(iUHF,FILENAME,MDNAME)

! Write natural orbitals to standard output.
if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,'(A)') '  The CASPT2 orbitals are computed as natural orbitals of a density matrix'
  write(u6,'(A)') '  defined as:'
  write(u6,'(A)') '   D = (D0 + D1 + D2)/<PSI|PSI>'
  write(u6,'(A)') ' where D0..D2 are zeroth..2nd order contributions'
  write(u6,'(A)') ' and |PSI> is the total wave function.'
  write(u6,'(A)') ' A new RasOrb file named PT2ORB is prepared.'
  if (PRORB) then
    if (OUTFMT == 'LONG') then
      THRENE = Two**31
      THROCC = -Two**31
    else if (OUTFMT == 'DEFAULT ') then
      THRENE = Five
      THROCC = 5.0e-4_wp
    end if
    call PRIMO('Output orbitals from CASPT2',.true.,.false.,THROCC,THRENE,NSYM,NBAS,NBAS,BNAME,DUM,OCC,CNAT,-1)
  end if
end if

! compute Mulliken's orbital populations

if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') 'Mulliken population Analysis:'
  write(u6,'(6X,A)') '-----------------------------'

  call mma_allocate(Scr,NBAST**2,Label='Scr')
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSyLbl = 1
  Label = 'Mltpl  0'
  call RdOne(iRc,iOpt,Label,iComp,Scr,iSyLbl)
  if (iRc == 0) then
    lSave = MSTATE(JSTATE) == irlxroot
    call Charge(nSym,nBas,bName,CNAT,OCC,Scr,2,FullMlk,lSave)
  end if
  call mma_deallocate(Scr)
end if

! compute one-electron properties

if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,'(6X,A)') 'Expectation values of various properties:'
  write(u6,'(6X,A)') '-----------------------------------------'
end if

nDens = sum(nTri_Elem(nBas(1:nSym)))
call mma_allocate(Scr,NDENS,Label='Scr')

call DOne_Caspt2(CNAT,nCMO,OCC,nOcc,Scr,nDENS)
call Put_dArray('D1ao',Scr,nDens)

Note = 'Temporary orbital file used by prpt.'
if (mode == 1) Note = 'var'
LuTmp = 50
LuTmp = IsFreeUnit(LuTmp)
call WrVec('TMPORB',LuTmp,'CO',nSym,nBas,nBas,CNAT,OCC,Dummy,IndType,Note)
call Prpt()

!nf
!------- ESPF analysis
call DecideOnESPF(Do_ESPF)
lSave = MSTATE(JSTATE) == irlxroot
if (Do_ESPF) call espf_analysis(lSave)
!nf

! On return from PrPt the 1-particle matrix is stored
! in the beginning of the scratch array.

call mma_deallocate(Scr)
call mma_deallocate(CNAT)
call mma_deallocate(OCC)

end subroutine PRPCTL
