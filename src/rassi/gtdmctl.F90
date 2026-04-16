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

subroutine GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,NZ,IDISK)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: MUL, nIrrep
use frenkel_global_vars, only: DoCoul
use gugx, only: CIStruct, EXStruct, SGStruct
use mspt2_eigenvectors, only: Heff_evc_pc, Heff_evc_sc, prpdata_mspt2_eigenvectors
use rasdef, only: NRAS, NRASEL, NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T, NRSPRT
use rasscf_global, only: DoDMRG
use rassi_aux, only: AO_Mode, iDisk_TDM, ipglob, jDisk_TDM
use rassi_data, only: ENUC, NASH, NASHT, NCMO, NDEL, NFRO, NISH, NISHT, NOSH, NSSH, NTDMAB, NTDMZZ, NTRA
use rassi_global_arrays, only: CnfTab1, CnfTab2, FSBTAB1, FSBTAB2, HAM, OrbTab, PART, REST1, REST2, SFDYS, SPNTAB1, SPNTAB2, &
                               SSTAB, TRANS1, TRANS2
use Cntrl, only: CITHR, DCHO, DCHS, DOGSOR, DYSO, ERFNUC, IFEJOB, IFHAM, IFHEFF, IFHEXT, IFNTO, IFTRD1, IFTRD2, IRREP, ISTAT, &
                 iToc15, JBNAME, LSYM1, LSYM2, LuIph, LuTDM, MLTPLT, NACTE, NATO, NCONF1, NDET, NELE3, NHOLE1, NPROP, NSTAT, &
                 NSTATE, PRCI, QDPT2EV, QDPT2SC, RASTYP, SAVEDENS, SECOND_TIME, TDYS, sonatnstate
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use Cntrl, only: CIH5
use RASSIWfn, only: wfn_CMO, wfn_CMO_OR, wfn_DetCoeff, wfn_DetCoeff_OR, wfn_DetOcc, wfn_DetOcc_OR
#endif
#ifdef _DMRG_
use qcmaquis_info, only: qcm_prefixes
use qcmaquis_interface_cfg, only: dmrg_external
use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_overlap
use qcmaquis_interface_utility_routines, only: pretty_print_util
use rassi_global_arrays, only: LROOT
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auToEV
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: PROP(NSTATE,NSTATE,NPROP), OVLP(NSTATE,NSTATE), DYSAMPS(NSTATE,NSTATE)
integer(kind=iwp), intent(in) :: JOB1, JOB2, NZ
integer(kind=iwp), intent(inout) :: IDISK
integer(kind=iwp) :: AUGSPIN, DCHIJ, I, IAD, IE1, IE1MN, IE1MX, IE2, IE2MN, IE2MX, IE3, IE3MN, IE3MX, IEMPTY, IERR, IFORM, IGO, &
                     IJ, IO, IOP1, IOP2, IOP3, IOPT, IRC, ISOIND, ISORB, ISPART, IST, ISTATE, ISUM, ISY, ISY12, ISYM, IT, ITABS, &
                     J, JST, JSTATE, JSY, KOINFO, KSPART, LSY, LUCITH, LUIPHn, MAXOP, MINOP, MPLET1, MPLET2, MSPROJ1, MSPROJ2, &
                     nActE1, NACTE2, NASHES(8), NASORB, NASPRT, NCONF2, NDCHSM, NDET1, NDET2, NDYSAB, NDYSZZ, NELE31, NELE32, &
                     NGAS, NGASLIM(2,10), NGASORB(100), NGL11, NGL12, NGL13, NGL21, NGL22, NGL23, NHOL11, NHOL12, NI, NJ, NL, NO, &
                     NPART, NRT2M, NRT2MAB, NTDM1, NTDM2, NTRAD
real(kind=wp) :: BEi, BEij, BEj, Dot_Prod, DYNORM, DYSAMP, ECORE, Energies(1:20), fac1, fac2, HII, HIJ, HJJ, HONE, HTWO, HZERO, &
                 Norm_fac, OVERLAP_RASSI, SIJ
real(kind=wp), pointer :: DET1(:), DET2(:)
logical(kind=iwp) :: DoNTO, IF00, IF01, IF02, IF10, IF11, IF12, IF20, IF21, IF22, IFTWO, mstate_dens, TRORB
character(len=48) :: STLNE2
character(len=8) :: WFTP1, WFTP2
type(CIStruct) :: CIS(2)
type(EXStruct) :: EXS(2)
type(SGStruct) :: SGS(2)
integer, allocatable :: OMAP(:)
real(kind=wp), allocatable :: CI1(:), CI2(:), CI2_o(:), CMO1(:), CMO2(:), DCHSM(:), detcoeff1(:), detcoeff2(:), DYSAB(:), &
                              DYSCOF(:), DYSZZ(:), FMO(:), mixed_1p_overlap(:,:), mixed_1p_rtdm(:,:,:), mixed_1p_stdm(:,:,:), &
                              mixed_1p_wtdm(:,:,:), RT2M(:), RT2MAB(:), TDM2(:), TDMAB(:), TDMZZ(:), Theta1(:), ThetaM(:), &
                              ThetaN(:), TRA1(:), TRA2(:), TRAD(:), TRASD(:), TSDMAB(:), TSDMZZ(:), TUVX(:), WDMAB(:), WDMZZ(:), &
                              WERD(:)
real(kind=wp), allocatable, target :: DETTOT1(:,:), DETTOT2(:,:)
character(len=NASHT+1), allocatable :: detocc(:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: DDot_

#define _TIME_GTDM
#ifdef _TIME_GTDM_
call CWTime(TCpu1,TWall1)
#endif

! WF parameters for ISTATE and JSTATE

NACTE1 = NACTE(JOB1)
MPLET1 = MLTPLT(JOB1)
LSYM1 = IRREP(JOB1)
NHOL11 = NHOLE1(JOB1)
NELE31 = NELE3(JOB1)
WFTP1 = RASTYP(JOB1)
NACTE2 = NACTE(JOB2)
MPLET2 = MLTPLT(JOB2)
LSYM2 = IRREP(JOB2)
NHOL12 = NHOLE1(JOB2)
NELE32 = NELE3(JOB2)
WFTP2 = RASTYP(JOB2)
SGS(1)%IFRAS = 1
SGS(2)%IFRAS = 1
if (IPGLOB >= 4) then
  write(u6,*) ' Entered GTDMCTL.'
  write(u6,'(1X,A,I3,A,I3)') '  JOB1:',JOB1,'        JOB2:',JOB2
  write(u6,'(1X,A,I3,A,I3)') 'NACTE1:',NACTE1,'      NACTE2:',NACTE2
  write(u6,'(1X,A,I3,A,I3)') 'MPLET1:',MPLET1,'      MPLET2:',MPLET2
  write(u6,'(1X,A,I3,A,I3)') ' LSYM1:',LSYM1,'       LSYM2:',LSYM2
  write(u6,'(1X,A,A8,A,A8)') ' WFTP1:',WFTP1,'       WFTP2:',WFTP2
  write(u6,'(1X,A,I3,A,I3)') ' NPROP:',NPROP
end if

if (doDMRG .and. (nacte1 /= nacte2)) then
  call WarningMessage(2,'Problem in gtdmctl for MPS-SI: no match for #e- in bra/ket')
  call abend()
end if

!> Logical variables, controlling which GTDM's to compute

!>  Overlap
IF00 = (NACTE1 == NACTE2) .and. (MPLET1 == MPLET2)
IF00 = IF00 .and. (LSYM1 == LSYM2)

!> Dyson amplitudes
IF10 = (NACTE1-NACTE2) == 1
IF10 = IF10 .and. (abs(MPLET1-MPLET2) == 1)
IF01 = (NACTE1-NACTE2) == -1
IF01 = IF01 .and. (abs(MPLET1-MPLET2) == 1)

!> Pair amplitudes:
IF20 = (NACTE1-NACTE2) == 2
IF02 = (NACTE1-NACTE2) == -2

!> 1-TDMs and transition spin densities
IF11 = (NACTE1 == NACTE2) .and. (NACTE1 >= 1)
IF11 = IF11 .and. (abs(MPLET1-MPLET2) <= 2)

!> 2h1p and 1h2p amplitudes:
IF21 = IF10 .and. (NACTE2 >= 1)
IF21 = IF21 .and. (abs(MPLET1-MPLET2) <= 3)
IF12 = IF01 .and. (NACTE1 >= 1)
IF12 = IF12 .and. (abs(MPLET1-MPLET2) <= 3)

!> 2-TDMs and transition spin densities
IF22 = (NACTE1 == NACTE2) .and. (NACTE1 >= 2)
IF22 = IF22 .and. (abs(MPLET1-MPLET2) <= 4)

!> check if they are needed at all:
!> It may be that the Hamiltonian matrix should be used in
!> diagonalization (IFHAM is .TRUE.), but it does not have
!> to be computed (because IFHEXT or IFHEFF or IFEJOB are true).
IFTWO = IFHAM .and. (.not. (IFHEXT .or. IFHEFF .or. IFEJOB))

!> For the moment, we have no use for the two-electron density
!> except when used for the scalar two-body Hamiltonian matrix:
IF22 = IF22 .and. IFTWO .and. (MPLET1 == MPLET2) .and. (LSYM1 == LSYM2)

if (IPGLOB >= 4) then
  if (IF00) write(u6,*) ' Overlap will be computed.'
  if (IF10 .or. IF01) write(u6,*) ' Dyson orbital will be computed.'
  if (IF20 .or. IF02) write(u6,*) ' Pair amplitudes will be computed.'
  if (IF11) write(u6,*) ' Density 1-matrix will be computed.'
  if (IF21 .or. IF12) write(u6,*) ' 2h1p amplitudes will be computed.'
  if (IF22) write(u6,*) ' Density 2-matrix will be computed.'
end if

! Pick up orbitals of ket and bra states.
call mma_allocate(CMO1,nCMO,Label='CMO1')
call mma_allocate(CMO2,nCMO,Label='CMO2')
call RDCMO_RASSI(JOB1,CMO1)
call RDCMO_RASSI(JOB2,CMO2)

! Nr of active spin-orbitals
NASORB = 2*NASHT
NTDM1 = NASHT**2
NTDM2 = nTri_Elem(NTDM1)

! Size of some data sets of reduced-2TDM in terms of active
! orbitals NASHT (For Auger matrix elements):
NRT2M = NASHT**3
! Size of Symmetry blocks
ISY12 = MUL(LSYM1,LSYM2)
NRT2MAB = 0
do ISY=1,nIrrep
  NI = NOSH(ISY)
  if (NI == 0) cycle
  do JSY=1,nIrrep
    NJ = NOSH(JSY)
    if (NJ == 0) cycle
    do LSY=1,nIrrep
      NL = NOSH(LSY)
      if (NL == 0) cycle
      if (MUL(ISY,MUL(JSY,LSY)) == ISY12) NRT2MAB = NRT2MAB+NI*NJ*NL
    end do
  end do
end do
if (TDYS .and. (.not. DYSO)) then
  write(u6,*)
  write(u6,*) 'Auger (TDYS) requires Dyson calculation.'
  write(u6,*) 'Make sure to activate Dyson in your RASSI input.'
  write(u6,*) 'For now, Auger computation will be skipped.'
end if
! evaluation of DCH
if (DCHS) NDCHSM = NASHT**2

! +++ J. Norell 13/7 - 2018
! 1D arrays for Dyson orbital coefficients
! COF = active biorthonormal orbital base
! AB  = inactive+active biorthonormal orbital base
! ZZ  = atomic (basis function) base
if ((IF10 .or. IF01) .and. DYSO) then
  call mma_allocate(DYSCOF,NASORB,Label='DYSCOF')
  ! Number of inactive+active orbitals
  NDYSAB = NASHT+NISHT
  call mma_allocate(DYSAB,nDYSAB,Label='DYSAB')
  ! Number of atomic / basis functions
  NDYSZZ = NZ
  call mma_allocate(DYSZZ,nDYSZZ,Label='DYSZZ')
  DYSZZ(:) = Zero
end if
! +++

! Transition density matrices, TDMAB is for active biorthonormal
! orbitals only, while TDMZZ is in the fixed AO basis.
! WDMAB, WDMZZ similar, but WE-reduced 'triplet' densities.
if (IF11 .and. (NATO .or. (NPROP > 0))) then
  call mma_allocate(TDMAB,nTDMAB,Label='TDMAB')
  call mma_allocate(TSDMAB,nTDMAB,Label='TSDMAB')
  call mma_allocate(WDMAB,nTDMAB,Label='WDMAB')
  call mma_allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
  call mma_allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
  call mma_allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
end if

if (IF11) then
  NTRAD = NASHT**2 ! NTRAD == NWERD == NTRASD
  call mma_allocate(TRAD,nTRAD+1,Label='TRAD')
  call mma_allocate(TRASD,nTRAD+1,Label='TRASD')
  call mma_allocate(WERD,nTRAD+1,Label='WERD')
end if
if (IF22) then
  call mma_allocate(TDM2,nTDM2,Label='TDM2')
else
  ! To avoid passing an unallocated argument
  call mma_allocate(TDM2,0,Label='TDM2')
end if

if (JOB1 /= JOB2) then
  ! Transform to biorthonormal orbital system
  if (DoGSOR) then
    call FCopy(trim(JBNAME(JOB2)),'JOBGS',ierr)
    call DANAME(LUIPH,'JOBGS')
    IAD = 0
    call IDAFile(LUIPH,2,ITOC15,30,IAD)
    IAD = ITOC15(2)
    call DDAFile(LUIPH,1,CMO2,nCMO,IAD)
    call DACLOS(LUIPH)
  end if !DoGSOR

  call mma_allocate(TRA1,nTRA,Label='TRA1')
  call mma_allocate(TRA2,nTRA,Label='TRA2')
  call FINDT(CMO1,CMO2,TRA1,TRA2)
# ifdef _HDF5_
  ! put the pair of transformed orbitals to h5
  if (CIH5) then
    call mh5_put_dset(wfn_cmo,CMO1,[nCMO,1],[0,JOB1-1])
    call mh5_put_dset(wfn_cmo,CMO2,[nCMO,1],[0,JOB2-1])
  end if
# endif
  TRORB = .true.
else
# ifdef _HDF5_
  ! put original orbitals to hdf5 file
  if (CIH5) call mh5_put_dset(wfn_cmo_or,CMO1,[nCMO,1],[0,JOB1-1])
# endif
  TRORB = .false.
end if

!> check whether we do RASSI with an effective multi-state PT2 Hamiltonian
!> whose eigenvectors are stored in Heff_evc_pc/Heff_evc_sc
!> i.e., we do not use mixed CI coefficients / MPS wave functions but rather mix the TDMs

mstate_dens = (job1 == job2) .and. (allocated(Heff_evc_pc) .or. allocated(Heff_evc_sc))
mstate_dens = mstate_dens .and. if11
mstate_dens = mstate_dens .and. qdpt2ev

if (mstate_dens) then
  call mma_allocate(mixed_1p_overlap,nstat(job1),nstat(job1),Label='mixed_1p_overlap')
  call mma_allocate(mixed_1p_rtdm,NTDMZZ,nstat(job1),nstat(job1),Label='mixed_1p_rtdm')
  call mma_allocate(mixed_1p_stdm,NTDMZZ,nstat(job1),nstat(job1),Label='mixed_1p_stdm')
  call mma_allocate(mixed_1p_wtdm,NTDMZZ,nstat(job1),nstat(job1),Label='mixed_1p_wtdm')
  mixed_1p_rtdm(:,:,:) = 0
  mixed_1p_stdm(:,:,:) = 0
  mixed_1p_wtdm(:,:,:) = 0
  mixed_1p_overlap(:,:) = 0
end if

#ifdef _DMRG_
dmrg_external%MPSrotated = trorb
#endif

! OBTAIN CORE ENERGY, FOCK MATRIX, AND TWO-ELECTRON INTEGRALS
! IN THE MIXED ACTIVE MO BASIS:
ECORE = Zero
if (IFTWO .and. (MPLET1 == MPLET2)) then
  call mma_allocate(FMO,nTDM1,Label='FMO')
  call mma_allocate(TUVX,nTDM2,Label='TUVX')
  TUVX(:) = Zero
  !TEST  write(u6,*) 'GTDMCTL calling TRINT.'
  call TRINT(CMO1,CMO2,ECORE,nTDM1,FMO,nTDM2,TUVX)
  ECORE = ENUC+ERFNUC+ECORE
  !TEST  write(u6,*) 'GTDMCTL back from TRINT.'
  !TEST  write(u6,*) 'ENUC  =',ENUC
  !TEST  write(u6,*) 'ERFNUC=',ERFNUC
  !TEST  write(u6,*) 'ECORE =',ECORE
end if

! In the calculation of matrix elements ( S1, S2 ), we will use
! the same Ms quantum numbers, if S1 and S2 differ by 0 or an int,
! else we will use Ms quantum numbers that differ by 1/2:
if (mod(abs(MPLET1-MPLET2),2) == 0) then
  MSPROJ1 = min(MPLET1,MPLET2)-1
  MSPROJ2 = MSPROJ1
else
  if (MPLET1 > MPLET2) then
    MSPROJ2 = MPLET2-1
    MSPROJ1 = MSPROJ2+1
  else
    MSPROJ1 = MPLET1-1
    MSPROJ2 = MSPROJ1+1
  end if
end if

#ifdef _DMRG_
!> set spin-up/spin-down # of electrons for target state(s)
if (doDMRG) then
  dmrg_external%nalpha = (nacte1+msproj1)/2
  dmrg_external%nbeta = (nacte1-msproj1)/2
end if
#endif

!---------------  For all wave functions: ---------------------
! Define structures ('tables') pertinent all jobs.
! (Later, move this up before the GTDMCTL calls).
! These are at:
! PART
! ORBTAB
! SSTAB
call NEWPRTTAB(nIrrep,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
if (IPGLOB >= 4) call PRPRTTAB(PART)

call NEWORBTAB(PART)
if (IPGLOB >= 4) call PRORBTAB(ORBTAB)

call NEWSSTAB(ORBTAB)
if (IPGLOB >= 4) call PRSSTAB(SSTAB)

! Mapping from active spin-orbital to active orbital in external order.
! Note that these differ, not just because of the existence of two
! spin-orbitals for each orbital, but also because the active orbitals
! (external order) are grouped by symmetry and then RAS space, but the
! spin orbitals are grouped by subpartition.
call mma_allocate(OMAP,NASORB,Label='OMAP')
NASPRT = ORBTAB(9)
KSPART = ORBTAB(10)
KOINFO = 19
ISUM = 0
do ISYM=1,nIrrep
  NASHES(ISYM) = ISUM
  ISUM = ISUM+NASH(ISYM)
end do
ISORB = 0
do ISPART=1,NASPRT
  NO = OrbTab(KSPART-1+ISPART)
  do IO=1,NO
    ISORB = ISORB+1
    ! Orbital symmetry:
    ISYM = OrbTab(KOINFO+1+(ISORB-1)*8)
    ! In-Symmetry orbital index:
    ISOIND = OrbTab(KOINFO+2+(ISORB-1)*8)
    ! Subtract nr of inactive orbitals in that symmetry:
    IT = ISOIND-NISH(ISYM)
    ! Add nr of actives in earlier symmetries:
    ITABS = NASHES(ISYM)+IT
    OMAP(ISORB) = ITABS
  end do
end do

!---------------    JOB1 wave functions: ---------------------
! Initialize SGUGA tables for JOB1 functions.
! These are structures stored in user defined types:
! SGS(1),CIS(1) and EXS(1).

! Set variables in /RASDEF/, used by SGUGA codes, which define
! the SGUGA space of JOB1. General RAS:
if (WFTP1 == 'GENERAL') then
  NRSPRT = 3
  NRAS(:,1) = NRS1(:)
  NRAS(:,2) = NRS2(:)
  NRAS(:,3) = NRS3(:)
  NRASEL(1) = 2*NRS1T-NHOL11
  NRASEL(2) = NACTE1-NELE31
  NRASEL(3) = NACTE1

  if (.not. doDMRG) then
    call SGINIT(nIrrep,NACTE1,MPLET1,SGS(1),CIS(1))
    if (IPGLOB > 4) then
      write(u6,*) 'Split-graph structure for JOB1=',JOB1
      call SGPRINT(SGS(1))
    end if
    call CXINIT(SGS(1),CIS(1),EXS(1))
    ! CI sizes, as function of symmetry, are now known.
    NCONF1 = CIS(1)%NCSF(LSYM1)
  else
    NCONF1 = 1
  end if
else
  ! Presently, the only other cases are HISPIN, CLOSED or EMPTY.
  ! Note: the HISPIN case may be buggy and is not used presently.
  NCONF1 = 1
end if
call mma_allocate(CI1,NCONF1,Label='CI1')

! Still JOB1, define structures ('tables') pertinent to JOB1
! These are at:
! REST1
! CNFTAB1
! FSBTAB1
! SPNTAB1

NPART = 3
NGAS = NPART
do ISYM=1,nIrrep
  NGASORB(1:nIrrep) = NRS1(1:nIrrep)
  NGASORB(nIrrep+1:2*nIrrep) = NRS2(1:nIrrep)
  NGASORB(2*nIrrep+1:3*nIrrep) = NRS3(1:nIrrep)
end do

!PAM2008: The old MAXOP was far too generous:
!MAXOP = NASHT
!PAM2008: MAXOP is determined by RAS restrictions:
! Preliminary ranges of nr of electrons:
IE1MN = max(0,2*NRS1T-NHOL11)
IE1MX = 2*NRS1T
IE3MN = max(0,NACTE1-IE1MX-2*NRS2T)
IE3MX = min(2*NRS3T,NELE31)
IE2MN = max(0,NACTE1-IE1MX-IE3MX)
IE2MX = min(2*NRS2T,NACTE1-IE1MN-IE3MN)
! Preliminary NGASLIM:
NGL11 = IE1MX
NGL21 = IE1MN
NGL12 = IE2MX
NGL22 = IE2MN
NGL13 = IE3MX
NGL23 = IE3MN
! Start with MAXOP=0, then increase:
MAXOP = 0
! Loop over possible ranges:
do IE1=IE1MN,IE1MX
  IOP1 = min(IE1,(2*NRS1T-IE1))
  if (IOP1 >= 0) then
    do IE3=IE3MN,IE3MX
      IOP3 = min(IE3,(2*NRS3T-IE3))
      if (IOP3 >= 0) then
        IE2 = NACTE1-IE1-IE3
        IOP2 = min(IE2,(2*NRS2T-IE2))
        if (IOP2 >= 0) then
          ! Actually possible combination:
          MAXOP = max(MAXOP,IOP1+IOP2+IOP3)
          NGL11 = min(NGL11,IE1)
          NGL21 = max(NGL21,IE1)
          NGL12 = min(NGL12,IE2)
          NGL22 = max(NGL22,IE2)
          NGL13 = min(NGL13,IE3)
          NGL23 = max(NGL23,IE3)
        end if
      end if
    end do
  end if
end do
NGASLIM(1,1) = NGL11
NGASLIM(2,1) = NGL21
NGASLIM(1,2) = NGL12
NGASLIM(2,2) = NGL22
NGASLIM(1,3) = NGL13
NGASLIM(2,3) = NGL23

if (.not. doDMRG) then
  IFORM = 1
  MINOP = 0
  call NEWGASTAB(nIrrep,NGAS,NGASORB,NGASLIM,1)
  if (IPGLOB >= 4) call PRGASTAB(REST1)

  ! At present, we will only annihilate, at most 2 electrons will
  ! be removed. This limits the possible MAXOP:
  MAXOP = min(MAXOP+1,NACTE1,NASHT)
  call NEWCNFTAB(NACTE1,NASHT,MINOP,MAXOP,LSYM1,NGAS,NGASORB,NGASLIM,IFORM,1)
  if (IPGLOB >= 4) call PRCNFTAB(CNFTAB1,100)

  call NEWFSBTAB(NACTE1,MSPROJ1,LSYM1,REST1,SSTAB,1)
  if (IPGLOB >= 4) call PRFSBTAB(FSBTAB1)
  NDET1 = FSBTAB1(5)
  if (ndet1 /= ndet(job1)) ndet(job1) = ndet1
  call NEWSCTAB(MINOP,MAXOP,MPLET1,MSPROJ1,1)
  if (IPGLOB > 4) then
    !PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
    ! do not get their knickers into a twist.
    call PRSCTAB(SPNTAB1,TRANS1)
  end if
else
  NDET1 = 1 ! minimum to avoid runtime error
end if
!---------------    JOB2 wave functions: ---------------------
! Initialize SGUGA tables for JOB2 functions.
! These are structures stored in arrays:
! SGS(2),CIS(2) and EXS(2).

! Set variables in /RASDEF/, used by SGUGA codes, which define
! the SGUGA space of JOB1. General RAS:
if (WFTP2 == 'GENERAL') then
  NRSPRT = 3
  NRAS(:,1) = NRS1(:)
  NRAS(:,2) = NRS2(:)
  NRAS(:,3) = NRS3(:)
  NRASEL(1) = 2*NRS1T-NHOL12
  NRASEL(2) = NACTE2-NELE32
  NRASEL(3) = NACTE2

  if (.not. doDMRG) then
    call SGINIT(nIrrep,NACTE2,MPLET2,SGS(2),CIS(2))
    if (IPGLOB > 4) then
      write(u6,*) 'Split-graph structure for JOB2=',JOB2
      call SGPRINT(SGS(2))
    end if
    call CXINIT(SGS(2),CIS(2),EXS(2))
    ! CI sizes, as function of symmetry, are now known.
    NCONF2 = CIS(2)%NCSF(LSYM2)
  else
    NCONF2 = 1
  end if
else
  ! Presently, the only other cases are HISPIN, CLOSED or EMPTY.
  ! Note: the HISPIN case may be buggy and is not used presently.
  NCONF2 = 1
end if
call mma_allocate(CI2,NCONF2,Label='CI2')
if (DoGSOR) call mma_allocate(CI2_o,NCONF2,Label='CI2_o')

NPART = 3
NGAS = NPART
do ISYM=1,nIrrep
  NGASORB(1:nIrrep) = NRS1(1:nIrrep)
  NGASORB(nIrrep+1:2*nIrrep) = NRS2(1:nIrrep)
  NGASORB(2*nIrrep+1:3*nIrrep) = NRS3(1:nIrrep)
end do
!PAM2008: The old MAXOP was far too generous:
!MAXOP = NASHT
!PAM2008: MAXOP is determined by RAS restrictions:
! Preliminary ranges of nr of electrons:
IE1MN = max(0,2*NRS1T-NHOL12)
IE1MX = 2*NRS1T
IE3MN = max(0,NACTE2-IE1MX-2*NRS2T)
IE3MX = min(2*NRS3T,NELE32)
IE2MN = max(0,NACTE2-IE1MX-IE3MX)
IE2MX = min(2*NRS2T,NACTE2-IE1MN-IE3MN)
! Preliminary NGASLIM:
NGL11 = IE1MX
NGL21 = IE1MN
NGL12 = IE2MX
NGL22 = IE2MN
NGL13 = IE3MX
NGL23 = IE3MN
! Start with MAXOP=0, then increase:
MAXOP = 0
! Loop over possible ranges:
do IE1=IE1MN,IE1MX
  IOP1 = min(IE1,(2*NRS1T-IE1))
  if (IOP1 >= 0) then
    do IE3=IE3MN,IE3MX
      IOP3 = min(IE3,(2*NRS3T-IE3))
      if (IOP3 >= 0) then
        IE2 = NACTE2-IE1-IE3
        IOP2 = min(IE2,(2*NRS2T-IE2))
        if (IOP2 >= 0) then
          ! Actually possible combination:
          MAXOP = max(MAXOP,IOP1+IOP2+IOP3)
          NGL11 = min(NGL11,IE1)
          NGL21 = max(NGL21,IE1)
          NGL12 = min(NGL12,IE2)
          NGL22 = max(NGL22,IE2)
          NGL13 = min(NGL13,IE3)
          NGL23 = max(NGL23,IE3)
        end if
      end if
    end do
  end if
end do
NGASLIM(1,1) = NGL11
NGASLIM(2,1) = NGL21
NGASLIM(1,2) = NGL12
NGASLIM(2,2) = NGL22
NGASLIM(1,3) = NGL13
NGASLIM(2,3) = NGL23

if (.not. dodmrg) then
  call NEWGASTAB(nIrrep,NGAS,NGASORB,NGASLIM,2)
  if (IPGLOB >= 4) call PRGASTAB(REST2)

  IFORM = 1
  MINOP = 0
  ! At present, we will only annihilate. This limits the possible MAXOP:
  MAXOP = min(MAXOP+1,NACTE2,NASHT)
  call NEWCNFTAB(NACTE2,NASHT,MINOP,MAXOP,LSYM2,NGAS,NGASORB,NGASLIM,IFORM,2)
  if (IPGLOB >= 4) call PRCNFTAB(CNFTAB2,100)

  call NEWFSBTAB(NACTE2,MSPROJ2,LSYM2,REST2,SSTAB,2)
  if (IPGLOB >= 4) call PRFSBTAB(FSBTAB2)
  NDET2 = FSBTAB2(5)
  if (ndet2 /= ndet(job2)) ndet(job2) = ndet2
  call NEWSCTAB(MINOP,MAXOP,MPLET2,MSPROJ2,2)
  !PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
  ! do not get their knickers into a twist.
  if (IPGLOB > 4) call PRSCTAB(SPNTAB2,TRANS2)
else
  NDET2 = 1 ! minimum to avoid runtime error
end if
!-------------------------------------------------------------
call mma_allocate(DETTOT1,NDET1,NSTAT(JOB1),Label='DETTOT1')
call mma_allocate(DETTOT2,NDET2,NSTAT(JOB2),Label='DETTOT2')
call mma_allocate(detocc,max(nDet1,nDet2),label='detocc')

! Loop over the states of JOBIPH nr JOB1
do IST=1,NSTAT(JOB1)
  DET1 => DETTOT1(1:NDET1,IST)
  ISTATE = ISTAT(JOB1)-1+IST

  if (.not. doDMRG) then
    ! Read ISTATE wave function
    if (WFTP1 == 'GENERAL') then
      call READCI(ISTATE,SGS(1),CIS(1),NCONF1,CI1)
    else
      CI1(1) = One
    end if
    DET1(:) = Zero
    ! Transform to bion basis, Split-Guga format
    if (TrOrb) call CITRA(WFTP1,SGS(1),CIS(1),EXS(1),LSYM1,TRA1,NCONF1,CI1)
    call mma_allocate(detcoeff1,nDet1,label='detcoeff1')
    call PREPSD(WFTP1,SGS(1),CIS(1),LSYM1,CNFTAB1,SPNTAB1,SSTAB,FSBTAB1,NCONF1,CI1,DET1,detocc,detcoeff1,TRANS1)

    ! print transformed ci expansion
    if (JOB1 /= JOB2) then
      if (PRCI) call prwf_biorth(istate,job1,nconf1,ndet1,nasht,detocc,detcoeff1,cithr)
#     ifdef _HDF5_
      ! put transformed ci coefficients for JOB1 to h5
      if (CIH5) then
        call mh5_put_dset(wfn_detcoeff,detcoeff1,[nDet1,1],[0,istate-1])
        call mh5_put_dset(wfn_detocc,detocc(:nDet1),[nDet1,1],[0,(JOB1-1)])
      end if
    else
      ! JOB1=JOB2, put original ci coefficients for JOB1 to h5
      if (CIH5) then
        call mh5_put_dset(wfn_detcoeff_or,detcoeff1,[nDet1,1],[0,istate-1])
        call mh5_put_dset(wfn_detocc_or,detocc(:nDet1),[nDet1,1],[0,(JOB1-1)])
      end if
#     endif
    end if

    call mma_deallocate(detcoeff1)

# ifdef _DMRG_
  else ! doDMRG
    call prepMPS(TRORB,LROOT(ISTATE),LSYM1,MPLET1,MSPROJ1,NACTE1,TRA1,NTRA,NISH,NASH,NOSH,nIrrep,6,job1,ist)
# endif
  end if
end do

if (DoGSOR) then
  call mma_allocate(Theta1,NCONF2,Label='Theta1')
  Theta1(:) = Zero
end if

!-------------------------------------------------------------

do JST=1,NSTAT(JOB2)
  DET2 => DETTOT2(1:NDET2,JST)
  JSTATE = ISTAT(JOB2)-1+JST
  if (.not. doDMRG) then
    ! Read JSTATE wave function
    if (WFTP2 == 'GENERAL') then
      call READCI(JSTATE,SGS(2),CIS(2),NCONF2,CI2)
    else
      CI2(1) = One
    end if
    if (DoGSOR) CI2_o(:) = CI2(:)
    DET2(:) = Zero
    ! Transform to bion basis, Split-Guga format
    if (TrOrb) call CITRA(WFTP2,SGS(2),CIS(2),EXS(2),LSYM2,TRA2,NCONF2,CI2)
    call mma_allocate(detcoeff2,nDet2,label='detcoeff2')
    call PREPSD(WFTP2,SGS(2),CIS(2),LSYM2,CNFTAB2,SPNTAB2,SSTAB,FSBTAB2,NCONF2,CI2,DET2,detocc,detcoeff2,TRANS2)

    ! print transformed ci expansion
    if (JOB1 /= JOB2) then
      if (PRCI) call prwf_biorth(jstate,job2,nconf2,ndet2,nasht,detocc,detcoeff2,cithr)
#     ifdef _HDF5_
      ! put ci coefficients for JOB2 to h5
      if (CIH5) then
        call mh5_put_dset(wfn_detcoeff,detcoeff2,[nDet2,1],[0,jstate-1])
        call mh5_put_dset(wfn_detocc,detocc(:nDet2),[nDet2,1],[0,(JOB2-1)])
      end if
#     endif
    end if
    call mma_deallocate(detcoeff2)

  else
#   ifdef _DMRG_
    call prepMPS(TRORB,lroot(JSTATE),LSYM2,MPLET2,MSPROJ2,NACTE2,TRA2,NTRA,NISH,NASH,NOSH,nIrrep,6,job2,jst)
#   endif
  end if
end do

! Loop over the states of JOBIPH nr JOB2
job2_loop: do JST=1,NSTAT(JOB2)
  JSTATE = ISTAT(JOB2)-1+JST
  ! Loop over the states of JOBIPH nr JOB1
  job1_loop: do IST=1,NSTAT(JOB1)
    ISTATE = ISTAT(JOB1)-1+IST
    if (ISTATE < JSTATE) cycle
    !-------------------------------------------------------------------

    ! Entry into monitor: Status line
    write(STLNE2,'(A33,I5,A5,I5)') 'Trans. dens. matrices for states ',ISTATE,' and ',JSTATE
    call StatusLine('RASSI: ',STLNE2)

    ! Read ISTATE WF from TOTDET1 and JSTATE WF from TOTDET2
#   ifdef _DMRG_
    if (.not. doDMRG) then
#   endif
      DET1 => DETTOT1(1:NDET1,IST)
      DET2 => DETTOT2(1:NDET2,JST)
#   ifdef _DMRG_
    end if
#   endif

    if (doGSOR) then
      if (JOB1 /= JOB2) then
        Dot_prod = DDOT_(NCONF2,CI1,1,CI2,1)
        THETA1(:) = THETA1(:)+Dot_prod*CI2(:)
      end if
    end if

    ! Calculate whatever type of GTDM that was requested, unless
    ! it is known to be zero.
    HZERO = Zero
    HONE = Zero
    HTWO = Zero

    SIJ = Zero
    DYSAMP = Zero
    ! +++ J. Norell 12/7 - 2018
    ! +++ Modified by Bruno Tenorio, 2020
    ! Dyson amplitudes:
    ! DYSAMP = D_ij for states i and j
    ! DYSCOF = Active orbital coefficents of the DO
    if ((IF10 .or. IF01) .and. DYSO) then
      call MKDYSORB(OrbTab,SSTAB,FSBTAB1,FSBTAB2,DET1,DET2,IF10,IF01,DYSAMP,DYSCOF)

      ! Write Dyson orbital coefficients in AO basis to disk.
      ! In full biorthonormal basis:
      call MKDYSAB(DYSCOF,DYSAB)
      ! Correct Dyson norms, for a biorth. basis. Add by Bruno
      call DYSNORM(CMO2,DYSAB,DYNORM) !do not change CMO2
      if (DYNORM > 1.0e-5_wp) then
        ! In AO basis:
        call MKDYSZZ(CMO2,DYSAB,DYSZZ)  !do not change CMO2
        if (DYSO) then
          SFDYS(:,JSTATE,ISTATE) = DYSZZ(:)
          SFDYS(:,ISTATE,JSTATE) = DYSZZ(:)
        end if
        DYSZZ(:) = Zero
        ! DYSAMPS corresponds to the Dyson norms corrected
        ! for a MO biorth. basis
        DYSAMPS(ISTATE,JSTATE) = sqrt(DYNORM)
        DYSAMPS(JSTATE,ISTATE) = sqrt(DYNORM)
      end if ! AMP THRS
    end if ! IF01 IF10

    ! ------------------------------------------------------------
    ! This part computes the needed densities for Auger.
    ! (DOI:10.1021/acs.jctc.2c00252)
    if ((IF21 .or. IF12) .and. TDYS .and. DYSO) then
      call mma_allocate(RT2M,nRT2M,Label='RT2M')
      RT2M(:) = Zero
      call mma_allocate(RT2MAB,nRT2MAB,Label='RT2MAB')
      RT2MAB(:) = Zero
      ! Defining the Binding energy Ei-Ej
      BEi = HAM(ISTATE,ISTATE)
      BEj = HAM(JSTATE,JSTATE)
      BEij = abs(BEi-BEj)*auToEV

      if ((MPLET1-MPLET2) == int(1)) then
        ! evaluate K-2V spin+1 density
        AUGSPIN = 1
        call MKRTDM2(FSBTAB1,FSBTAB2,SSTAB,OMAP,DET1,DET2,IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        call RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,RT2M,CMO1,CMO2,AUGSPIN)

      else if ((MPLET1-MPLET2) == int(-1)) then
        ! evaluate K-2V spin-1 density
        AUGSPIN = -1
        call MKRTDM2(FSBTAB1,FSBTAB2,SSTAB,OMAP,DET1,DET2,IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        call RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,RT2M,CMO1,CMO2,AUGSPIN)
      else ! write then both
        AUGSPIN = 1
        call MKRTDM2(FSBTAB1,FSBTAB2,SSTAB,OMAP,DET1,DET2,IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        call RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,RT2M,CMO1,CMO2,AUGSPIN)

        AUGSPIN = -1
        call MKRTDM2(FSBTAB1,FSBTAB2,SSTAB,OMAP,DET1,DET2,IF21,IF12,NRT2M,RT2M,AUGSPIN,OrbTab)
        call RTDM2_PRINT(ISTATE,JSTATE,BEij,NDYSAB,DYSAB,NRT2MAB,RT2M,CMO1,CMO2,AUGSPIN)
      end if
      call mma_deallocate(RT2M)
      call mma_deallocate(RT2MAB)
    end if
    ! ------------------------------------------------------------

    ! evaluation of DCH shake-up intensities (DOI:10.1063/5.0062130)
    if ((IF20 .or. IF02) .and. DCHS) then
      DCHIJ = DCHO+NASHT*(DCHO-1)
      ! Defining the Binding energy Ei-Ej
      BEi = HAM(ISTATE,ISTATE)
      BEj = HAM(JSTATE,JSTATE)
      BEij = abs(BEi-BEj)*auToEV
      call mma_allocate(DCHSM,nDCHSM,Label='DCHSM')
      DCHSM(:) = Zero
      call MKDCHS(FSBTAB1,FSBTAB2,SSTAB,OMAP,DET1,DET2,IF20,IF02,NDCHSM,DCHSM,OrbTab)
      write(u6,'(A,I5,I5,A,F14.5,ES23.14)') '  RASSI Pair States:',JSTATE,ISTATE,'  ssDCH BE(eV) and Norm:  ',BEij,DCHSM(DCHIJ)
      call mma_deallocate(DCHSM)
    end if
    ! ------------------------------------------------------------

    ! General 1-particle transition density matrix:
    if (IF11) then
      call MKTDM1(LSYM1,MPLET1,MSPROJ1,FSBTAB1,LSYM2,MPLET2,MSPROJ2,FSBTAB2,SSTAB,OMAP,DET1,DET2,SIJ,NASHT,TRAD,TRASD,WERD,ISTATE, &
                  JSTATE,job1,job2,OrbTab)
      ! Calculate Natural Transition Orbital (NTO):
      if (IFNTO) then
        DoNTO = job1 /= job2
        if (DoNTO) then
          call NTOCalc(job1,job2,ISTATE,JSTATE,TRAD,TRASD,MPLET1)
          write(u6,*) 'ntocalculation finished'
        end if
      end if
      ! End of Calculating NTO

      ! Compute 1-electron contribution to Hamiltonian matrix element:
      if (IFTWO .and. (MPLET1 == MPLET2)) HONE = DDOT_(NTRAD,TRAD,1,FMO,1)

      ! BEGIN MODIFIED by Aquilante, Segatta and Kaiser (2022)
      if (DoCoul) call EXCTDM(SIJ,TRAD,TDMAB,iRC,CMO1,CMO2,TDMZZ,TRASD,TSDMAB,TSDMZZ,ISTATE,JSTATE)
      ! END MODIFIED by Aquilante, Segatta and Kaiser(2022)

      ! Write density 1-matrices in AO basis to disk.
      if (NATO .or. (NPROP > 0)) then

        iEmpty = 0
        !> regular-TDM
        call MKTDAB(SIJ,TRAD,TDMAB,iRC)
        !> transform to AO basis
        call MKTDZZ(CMO1,CMO2,TDMAB,TDMZZ,iRC)
        if (iRC == 1) iEmpty = 1

        !> spin-TDM
        call MKTDAB(Zero,TRASD,TSDMAB,iRC)
        !> transform to AO basis
        call MKTDZZ(CMO1,CMO2,TSDMAB,TSDMZZ,iRC)
        if (iRC == 1) iEmpty = iEmpty+2

        !> WE-reduced TDM's of triplet type:
        call MKTDAB(Zero,WERD,WDMAB,iRC)
        !> transform to AO basis
        call MKTDZZ(CMO1,CMO2,WDMAB,WDMZZ,iRC)
        if (iRC == 1) iEmpty = iEmpty+4

        if (.not. mstate_dens) then

          if (SaveDens) then
            ! Transition density matrices, TDMZZ, in AO or MO basis.
            ! WDMZZ similar, but WE-reduced 'triplet' densities.
            ij = nTri_Elem(iSTATE-1)+JSTATE
            jDisk_TDM(1,ij) = IDISK
            jDisk_TDM(2,ij) = iEmpty
            iOpt = 1
            iGo = ibset(ibset(ibset(0,0),1),2)
            if (AO_Mode) then
              call dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,IState,jState)
            else
              iEmpty = 0
              TRAD(nTrad+1) = SIJ
              iRC = 0
              if (DDot_(nTRAD+1,TRAD,1,TRAD,1) > Zero) iRC = 1
              if (iRC == 1) iEmpty = 1

              TRASD(nTrad+1) = Zero
              iRC = 0
              if (DDot_(nTRAD+1,TRASD,1,TRASD,1) > Zero) iRC = 1
              if (iRC == 1) iEmpty = iEmpty+2

              WERD(nTrad+1) = Zero
              iRC = 0
              if (DDot_(nTRAD+1,WERD,1,WERD,1) > Zero) iRC = 1
              if (iRC == 1) iEmpty = iEmpty+4

              call dens2file(TRAD,TRASD,WERD,nTRAD+1,LUTDM,IDISK,iEmpty,iOpt,iGo,ISTATE,JSTATE)
            end if
          end if
          !> calculate property matrix elements
          call PROPER(PROP,ISTATE,JSTATE,TDMZZ,WDMZZ)
        else

          !> scale rdm elements with eigenvector coefficients of Heff of a multi-state (PT2) Hamiltonian
          !> accumulate data first and run PROPER and other utility routines later
          do i=1,nstat(job1)
            do j=1,nstat(job2)

              if (i < j) cycle

              if (qdpt2sc) then
                fac1 = Heff_evc_sc(ist,i,job1)
                fac2 = Heff_evc_sc(jst,j,job2)
              else
                fac1 = Heff_evc_pc(ist,i,job1)
                fac2 = Heff_evc_pc(jst,j,job2)
              end if

              !> regular-TDM
              mixed_1p_rtdm(:,i,j) = mixed_1p_rtdm(:,i,j)+fac1*fac2*tdmzz(:)
              !> spin-TDM
              mixed_1p_stdm(:,i,j) = mixed_1p_stdm(:,i,j)+fac1*fac2*tsdmzz(:)
              !> WE-reduced TDM's of triplet type:
              mixed_1p_wtdm(:,i,j) = mixed_1p_wtdm(:,i,j)+fac1*fac2*wdmzz(:)
              !> overlap
              mixed_1p_overlap(i,j) = mixed_1p_overlap(i,j)+fac1*fac2*sij
            end do
          end do
        end if
      end if

    else ! IF11

      !> overlap
      if (IF00) then
#       ifdef _DMRG_
        if (.not. doDMRG) then
#       endif
          SIJ = OVERLAP_RASSI(FSBTAB1,FSBTAB2,DET1,DET2)
#       ifdef _DMRG_
        else
          sij = qcmaquis_mpssi_overlap(qcm_prefixes(job1),ist,qcm_prefixes(job2),jst,.true.)
        end if !doDMRG
#       endif
      end if ! IF00
    end if ! IF11

    if (.not. mstate_dens) then
      OVLP(ISTATE,JSTATE) = SIJ
      OVLP(JSTATE,ISTATE) = SIJ
    end if

    !> General 2-particle transition density matrix:
    if (IF22) then
      call MKTDM2(LSYM1,MSPROJ1,FSBTAB1,LSYM2,MSPROJ2,FSBTAB2,SSTAB,OMAP,DET1,DET2,NTDM2,TDM2,OrbTab)

      !> Compute 2-electron contribution to Hamiltonian matrix element:
      if (IFTWO .and. (MPLET1 == MPLET2)) HTWO = DDOT_(NTDM2,TDM2,1,TUVX,1)

    end if ! IF22

    !> PAM 2011 Nov 3, writing transition matrices if requested
    if ((IFTRD1 .or. IFTRD2) .and. (.not. mstate_dens)) call trd_print(ISTATE,JSTATE,IFTRD2 .and. IF22,TDMAB,TDM2,CMO1,CMO2,SIJ)

    ! Store SIJ temporarily
    if (IFEJOB .and. (ISTATE /= JSTATE)) then
      HAM(ISTATE,JSTATE) = SIJ
      HAM(JSTATE,ISTATE) = SIJ
    end if
    if (IFHAM .and. (.not. (IFHEXT .or. IFHEFF .or. IFEJOB))) then
      HZERO = ECORE*SIJ
      HIJ = HZERO+HONE+HTWO
      HAM(ISTATE,JSTATE) = HIJ
      HAM(JSTATE,ISTATE) = HIJ

      ! SI-PDFT related code for "second_time" case
      if (second_time) then
        Energies(:) = Zero
        call DANAME(LUIPH,'JOBGS')
        IAD = 0
        call IDAFILE(LUIPH,2,ITOC15,30,IAD)
        IAD = ITOC15(6)
        call DDAFILE(LUIPH,2,Energies,NSTAT(JOB1),IAD)
        do i=1,NSTAT(JOB1)
          HAM(i,i) = Energies(i)
        end do
        call DACLOS(LUIPH)
      end if

      if (IPGLOB >= 4) then
        write(u6,'(1x,a,2I5)') ' ISTATE, JSTATE:',ISTATE,JSTATE
        write(u6,'(1x,a,f16.8)') ' HZERO=',HZERO
        write(u6,'(1x,a,f16.8)') ' HONE =',HONE
        write(u6,'(1x,a,f16.8)') ' HTWO =',HTWO
        write(u6,'(1x,a,f16.8)') ' HIJ  =',HIJ
      end if
    end if
  end do job1_loop

end do job2_loop

! For ejob, create an approximate off-diagonal based on the overlap (temporarily stored in HIJ)

if (IFEJOB) then
  do JST=1,NSTAT(JOB2)
    JSTATE = ISTAT(JOB2)-1+JST
    do IST=1,NSTAT(JOB1)
      ISTATE = ISTAT(JOB1)-1+IST
      if (ISTATE <= JSTATE) cycle
      SIJ = HAM(ISTATE,JSTATE)
      HII = HAM(ISTATE,ISTATE)
      HJJ = HAM(JSTATE,JSTATE)
      HAM(ISTATE,JSTATE) = SIJ*(HII+HJJ)*Half
      HAM(JSTATE,ISTATE) = SIJ*(HII+HJJ)*Half
    end do
  end do
end if

if (DoGSOR) then
  if (job1 /= job2) then
    dot_prod = DDOT_(NCONF2,THETA1,1,THETA1,1)
    Norm_Fac = One/sqrt(dot_prod)
    THETA1(:) = Norm_Fac*THETA1(:)

    !Write theta1 to file.
    LUCITH = IsFreeUnit(87)
    !Open(unit=87,file='CI_THETA', action='write',iostat=ios)
    call Molcas_Open(LUCITH,'CI_THETA')
    do i=1,NCONF2
      write(LUCITH,*) Theta1(i)
    end do
    close(LUCITH)

    ! Now we need to build the other states.
    call mma_allocate(detcoeff2,nDet2,label='detcoeff2')
    do JST=2,NSTAT(JOB2)
      JSTATE = ISTAT(JOB2)-1+JST
      call READCI(JSTATE,SGS(2),CIS(2),NCONF2,CI2)
      CI2_o(:) = CI2(:)
      DET2(:) = Zero
      if (TrOrb) call CITRA(WFTP2,SGS(2),CIS(2),EXS(2),LSYM2,TRA2,NCONF2,CI2)
      call PREPSD(WFTP2,SGS(2),CIS(2),LSYM2,CNFTAB2,SPNTAB2,SSTAB,FSBTAB2,NCONF2,CI2,DET2,detocc,detcoeff2,TRANS2)

      call mma_allocate(ThetaN,NCONF2,Label='ThetaN')
      ThetaN(:) = Zero
      Norm_Fac = DDOT_(NCONF2,THETA1,1,CI2_o,1)
      ThetaN(:) = ThetaN(:)-Norm_Fac*THETA1(:)

      LUCITH = IsFreeUnit(LUCITH)
      call Molcas_Open(LUCITH,'CI_THETA')
      !open(unit=87,file='CI_THETA',action='read',iostat=ios)
      if (JST-1 >= 2) then
        do i=1,NCONF2
          read(LUCITH,*) dot_prod ! dummy
        end do
      end if
      call mma_allocate(ThetaM,NCONF2,Label='ThetaM')
      do IST=2,JST-1
        ! Read in previous theta vectors
        do i=1,NCONF2
          read(LUCITH,*) ThetaM(i)
        end do
        Dot_prod = DDOT_(NCONF2,ThetaM,1,CI2_o,1)
        ThetaN(:) = ThetaN(:)-Dot_prod*ThetaM(:)

      end do
      call mma_deallocate(detcoeff2)
      close(LUCITH)
      ! Normalize
      dot_prod = DDOT_(NCONF2,ThetaN,1,ThetaN,1)
      Norm_Fac = One/sqrt(dot_prod)
      ThetaN(:) = Norm_Fac*ThetaN(:)

      !dot_prod = DDOT_(NCONF2,THETA1,1,THETA1,1)
      !dot_prod = DDOT_(NCONF2,ThetaN,1,THETA1,1)
      !dot_prod = DDOT_(NCONF2,ThetaN,1,ThetaN,1)

      ! Write to file
      LUCITH = IsFreeUnit(LUCITH)
      call Molcas_Open(LUCITH,'CI_THETA')
      call Append_file(LUCITH)
      !open(unit=87,file='CI_THETA',position='append',iostat=ios,action='write')
      do i=1,nConf2
        write(LUCITH,*) ThetaN(i)
      end do
      close(LUCITH)
      ! Deallocate
      call mma_deallocate(ThetaN)
    end do
    ! Copy to new IPH file
    LUCITH = IsFreeUnit(LUCITH)
    call Molcas_Open(LUCITH,'CI_THETA')
    !open(unit=87,file='CI_THETA',iostat=ios,action='read')
    call DANAME(LUIPHn,'JOBGS')
    IAD = 0
    call IDAFILE(LUIPHn,2,ITOC15,30,IAD)
    IAD = ITOC15(4)
    do i=1,ISTAT(JOB1)-1
      do j=1,nCONF2
        read(LUCITH,*) ThetaM(i)
      end do
      call DDafile(LUIPHn,1,ThetaM,nCONF2,IAD)
    end do

    IAD = ITOC15(4)
    call DDAFILE(LUIPHn,2,ThetaM,nCONF2,IAD)
    call DDAFILE(LUIPHn,2,ThetaM,nCONF2,IAD)

    close(LUCITH)
    call DACLOS(LUIPHn)
    call mma_deallocate(ThetaM)
  end if
  call mma_deallocate(Theta1)
end if !DoGSOR

#ifdef _DMRG_
if (IPGLOB >= 4) then
  write(u6,*) 'full SF-HAMILTONIAN'
  write(u6,*) 'dimension: ',nstate**2
  call pretty_print_util(HAM,1,nstate,1,nstate,nstate,nstate,1,u6)
end if
#endif

!> create actual property data and put everything to file (if requested)
!> in case of using eigenvectors of a multi-state (PT2) Hamiltonian
if (mstate_dens) then
  do JST=1,NSTAT(JOB2)
    JSTATE = ISTAT(JOB2)-1+JST
    do IST=1,NSTAT(JOB1)
      ISTATE = ISTAT(JOB1)-1+IST
      if (istate < jstate) cycle

      ovlp(istate,jstate) = mixed_1p_overlap(ist,jst)
      ovlp(jstate,istate) = mixed_1p_overlap(ist,jst)

      call prpdata_mspt2_eigenvectors(mixed_1p_rtdm(:,ist,jst),mixed_1p_stdm(:,ist,jst),mixed_1p_wtdm(:,ist,jst),prop,nprop, &
                                      nstate,istate,jstate,ntdmzz,iDisk_TDM(JSTATE,ISTATE,1),iDisk_TDM(JSTATE,ISTATE,2),lutdm, &
                                      sonatnstate > 0,if11 .and. (lsym1 == lsym2))
    end do
  end do
end if

if (WFTP1 == 'GENERAL') then
  if (.not. doDMRG) call MkGUGA_Free(SGS(1),CIS(1),EXS(1))
end if
if (WFTP2 == 'GENERAL') then
  if (.not. doDMRG) call MkGUGA_Free(SGS(2),CIS(2),EXS(2))
end if

if (JOB1 /= JOB2) then
  call mma_deallocate(TRA1)
  call mma_deallocate(TRA2)
end if
nullify(DET1,DET2)
call mma_deallocate(DETTOT1)
call mma_deallocate(DETTOT2)
call mma_deallocate(detocc)
call mma_deallocate(CI2)
if (DoGSOR) call mma_deallocate(CI2_o)
call mma_deallocate(CI1)
if (.not. doDMRG) then
  call mma_deallocate(TRANS2)
  call mma_deallocate(TRANS1)
  call mma_deallocate(SPNTAB1)
  call mma_deallocate(SPNTAB2)
end if
if ((IF10 .or. IF01) .and. DYSO) then
  call mma_deallocate(DYSCOF)
  call mma_deallocate(DYSAB)
  call mma_deallocate(DYSZZ)
end if
if (IF11) then
  call mma_deallocate(TRAD)
  call mma_deallocate(TRASD)
  call mma_deallocate(WERD)
  if (NATO .or. (NPROP > 0)) then
    call mma_deallocate(TDMAB)
    call mma_deallocate(TSDMAB)
    call mma_deallocate(WDMAB)
    call mma_deallocate(TDMZZ)
    call mma_deallocate(TSDMZZ)
    call mma_deallocate(WDMZZ)
  end if
end if
call mma_deallocate(TDM2)

if (IFTWO .and. (MPLET1 == MPLET2)) then
  call mma_deallocate(FMO)
  call mma_deallocate(TUVX)
end if

call mma_deallocate(CMO2)
call mma_deallocate(CMO1)
call mma_deallocate(PART)
call mma_deallocate(OrbTab)
call mma_deallocate(SSTAB)
if (.not. doDMRG) then
  call mma_deallocate(REST2)
  call mma_deallocate(REST1)
  call mma_deallocate(CNFTAB2)
  call mma_deallocate(CNFTAB1)
  call mma_deallocate(FSBTAB2)
  call mma_deallocate(FSBTAB1)
end if
call mma_deallocate(OMAP)

!> release memory
if (mstate_dens) then
  call mma_deallocate(mixed_1p_overlap)
  call mma_deallocate(mixed_1p_rtdm)
  call mma_deallocate(mixed_1p_stdm)
  call mma_deallocate(mixed_1p_wtdm)
end if

#ifdef _TIME_GTDM_
call CWTime(TCpu2,TWall2)
write(u6,*) 'Time for GTDM : ',TCpu2-TCpu1,TWall2-TWall1
#endif

end subroutine GTDMCTL
