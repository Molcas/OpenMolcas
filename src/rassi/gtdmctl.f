************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,NZ,IDDET1,IDISK)

#ifdef _DMRG_
      use rassi_global_arrays, only: HAM, SFDYS, LROOT
#else
      use rassi_global_arrays, only: HAM, SFDYS
#endif
      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
      use qcmaquis_info
#endif
      use mspt2_eigenvectors
      use rassi_aux, only : jDisk_TDM, iDisk_TDM, AO_Mode
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='GTDMCTL')
#include "rasdim.fh"
#include "rasdef.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "rassiwfn.fh"
#include "stdalloc.fh"
      DIMENSION ISGSTR1(NSGSIZE), ISGSTR2(NSGSIZE)
      DIMENSION ICISTR1(NCISIZE), ICISTR2(NCISIZE)
      DIMENSION IXSTR1(NXSIZE), IXSTR2(NXSIZE)
      DIMENSION PROP(NSTATE,NSTATE,NPROP)
      DIMENSION NGASORB(100),NGASLIM(2,10)
      DIMENSION NASHES(8)
      DIMENSION OVLP(NSTATE,NSTATE)
      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION IDDET1(NSTATE)
      LOGICAL IF00, IF10,IF01,IF20,IF11,IF02,IF21,IF12,IF22
      LOGICAL IFTWO,TRORB
      CHARACTER*8 WFTP1,WFTP2
      CHARACTER*6 STLNE1
      CHARACTER*48 STLNE2
      Real*8 Energies(1:20)
      Integer IAD,LUIPHn,lThetaM,LUCITH
      Real*8 Norm_fac
CC    NTO section
      Logical DoNTO
CC    NTO section
      External IsFreeUnit

      type mixed_1pdensities
        real*8              :: overlap
        real*8, allocatable :: rtdm(:)
        real*8, allocatable :: stdm(:)
        real*8, allocatable :: wtdm(:)
      end type

      type(mixed_1pdensities), allocatable :: mstate_1pdens(:,:)

      logical               :: mstate_dens
      real*8                :: fac1, fac2
      real*8, Allocatable:: CMO1(:), CMO2(:)
      real*8, Allocatable:: TRAD(:), TRASD(:), WERD(:)
      real*8, Allocatable:: TDMAB(:), TSDMAB(:), WDMAB(:)
      real*8, Allocatable:: TDMZZ(:), TSDMZZ(:), WDMZZ(:)
      real*8, Allocatable:: TDM2(:), TRA1(:), TRA2(:), FMO(:), TUVX(:)
      real*8, Allocatable:: DYSCOF(:), DYSAB(:), DYSZZ(:)

#ifdef _DMRG_
!     strings for conversion of the qcmaquis h5 checkpoint names from 2u1 to su2u1
      character(len=3) :: mplet1s, msproj1s
      ! new checkpoint names
      character(len=2300) :: checkpoint1_2u1,checkpoint2_2u1
#else
      logical             :: doDMRG = .false.
#endif
#include "SysDef.fh"

      CALL QENTER(ROUTINE)
#define _TIME_GTDM
#ifdef _TIME_GTDM_
      Call CWTime(TCpu1,TWall1)
#endif
* Avoid compiler warnings about possibly unitialised mstate_1pdens
* The below can be removed if the file is compiled with
* -Wno-error=maybe-uninitialized
      allocate(mstate_1pdens(0,0))
      deallocate(mstate_1pdens)
C WF parameters for ISTATE and JSTATE
      NACTE1=NACTE(JOB1)
      MPLET1=MLTPLT(JOB1)
      LSYM1=IRREP(JOB1)
      NHOL11=NHOLE1(JOB1)
      NELE31=NELE3(JOB1)
      WFTP1=RASTYP(JOB1)
      NACTE2=NACTE(JOB2)
      MPLET2=MLTPLT(JOB2)
      LSYM2=IRREP(JOB2)
      NHOL12=NHOLE1(JOB2)
      NELE32=NELE3(JOB2)
      WFTP2=RASTYP(JOB2)
      IF(IPGLOB.GE.DEBUG) THEN
        WRITE(6,*)' Entered GTDMCTL.'
        WRITE(6,'(1X,A,I3,A,I3)')'  JOB1:',  JOB1,'        JOB2:',  JOB2
        WRITE(6,'(1X,A,I3,A,I3)')'NACTE1:',NACTE1,'      NACTE2:',NACTE2
        WRITE(6,'(1X,A,I3,A,I3)')'MPLET1:',MPLET1,'      MPLET2:',MPLET2
        WRITE(6,'(1X,A,I3,A,I3)')' LSYM1:', LSYM1,'       LSYM2:', LSYM2
        WRITE(6,'(1X,A,A8,A,A8)')' WFTP1:',WFTP1, '       WFTP2:',WFTP2
        WRITE(6,'(1X,A,I3,A,I3)')' NPROP:', NPROP
      END IF

      if(doDMRG .and. (nacte1 /= nacte2 ))then
        call WarningMessage(2,'Problem in gtdmctl for MPS-SI: no '//
     &   ' match for #e- in bra/ket')
        call abend()
      end if

      !> Logical variables, controlling which GTDM''s to compute

      !>  Overlap
      IF00 = NACTE1.EQ.NACTE2.AND. MPLET1.EQ.MPLET2
      IF00 = IF00.AND.LSYM1.EQ.LSYM2

      !> Dyson amplitudes
      IF10 = (NACTE1-NACTE2).EQ. 1
      IF10 = IF10.AND.( ABS(MPLET1-MPLET2).EQ.1)
      IF01 = (NACTE1-NACTE2).EQ.-1
      IF01 = IF01.AND.( ABS(MPLET1-MPLET2).EQ.1)

      !> Pair amplitudes:
      IF20 = (NACTE1-NACTE2).EQ. 2
      IF02 = (NACTE1-NACTE2).EQ.-2

      !> 1-TDMs and transition spin densities
      IF11 = ( NACTE1.EQ.NACTE2.AND.NACTE1.GE.1)
      IF11 = IF11.AND.( ABS(MPLET1-MPLET2).LE.2)

      !> 2h1p and 1h2p amplitudes:
      IF21 = IF10.AND.NACTE2.GE.1
      IF21 = IF21.AND.( ABS(MPLET1-MPLET2).LE.3)
      IF12 = IF01.AND.NACTE1.GE.1
      IF12 = IF12.AND.( ABS(MPLET1-MPLET2).LE.3)

      !> 2-TDMs and transition spin densities
      IF22 = (NACTE1.EQ.NACTE2.AND.NACTE1.GE.2)
      IF22 = IF22.AND.(ABS(MPLET1-MPLET2).LE.4)

      !> check if they are needed at all:
        !> It may be that the Hamiltonian matrix should be used in
        !> diagonalization (IFHAM is .TRUE.), but it does not have
        !> to be computed (because IFHEXT or IFHEFF or IFEJOB are true).
      IFTWO = IFHAM .AND..NOT.(IFHEXT.OR.IFHEFF.OR.IFEJOB)

      !> For the moment, we have no use for the two-electron density
      !> except when used for the scalar two-body Hamiltonian matrix:
      IF22 = IF22.AND.IFTWO.AND.(MPLET1.EQ.MPLET2).AND.(LSYM1.EQ.LSYM2)

      IF(IPGLOB.GE.DEBUG) THEN
        IF(IF00) WRITE(6,*)' Overlap will be computed.'
        IF(IF10.or.IF01) WRITE(6,*)' Dyson orbital will be computed.'
        IF(IF20.or.IF02) WRITE(6,*)' Pair amplitudes will be computed.'
        IF(IF11) WRITE(6,*)' Density 1-matrix will be computed.'
        IF(IF21.or.IF12) WRITE(6,*)' 2h1p amplitudes will be computed.'
        IF(IF22) WRITE(6,*)' Density 2-matrix will be computed.'
      END IF

C Pick up orbitals of ket and bra states.
      Call mma_allocate(CMO1,nCMO,Label='CMO1')
      Call mma_allocate(CMO2,nCMO,Label='CMO2')
      CALL RDCMO_RASSI(JOB1,CMO1)
      CALL RDCMO_RASSI(JOB2,CMO2)

C Nr of active spin-orbitals
      NASORB=2*NASHT
      NTDM1=NASHT**2
      NTSDM1=NASHT**2
      NTDM2=(NTDM1*(NTDM1+1))/2


! +++ J. Norell 13/7 - 2018
C 1D arrays for Dyson orbital coefficients
C COF = active biorthonormal orbital base
C AB  = inactive+active biorthonormal orbital base
C ZZ  = atomic (basis function) base
      IF ((IF10.or.IF01).and.DYSO) THEN
        Call mma_allocate(DYSCOF,NASORB,Label='DYSCOF')
        ! Number of inactive+active orbitals
        NDYSAB = NASHT+NISHT
        Call mma_allocate(DYSAB,nDYSAB,Label='DYSAB')
        ! Number of atomic / basis functions
        NDYSZZ = NZ
        Call mma_allocate(DYSZZ,nDYSZZ,Label='DYSZZ')
        DYSZZ(:)=0.0D0
      END IF
! +++

C Transition density matrices, TDMAB is for active biorthonormal
C orbitals only, while TDMZZ is in the fixed AO basis.
C WDMAB, WDMZZ similar, but WE-reduced 'triplet' densities.
      IF(IF11.AND.(NATO.OR.NPROP.GT.0)) THEN
        Call mma_allocate(TDMAB,nTDMAB,Label='TDMAB')
        Call mma_allocate(TSDMAB,nTDMAB,Label='TSDMAB')
        Call mma_allocate(WDMAB,nTDMAB,Label='WDMAB')
        Call mma_allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
        Call mma_allocate(TSDMZZ,nTDMZZ,Label='TSDMZZ')
        Call mma_allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
      END IF

      IF (IF11) THEN
        NTRAD=NASHT**2
        NTRASD=NASHT**2
        NWERD=NASHT**2
        Call mma_allocate(TRAD,nTRAD+1,Label='TRAD')
        Call mma_allocate(TRASD,nTRAD+1,Label='TRASD')
        Call mma_allocate(WERD,nTRAD+1,Label='WERD')
      END IF
      IF (IF22) Call mma_allocate(TDM2,nTDM2,Label='TDM2')

      IF(JOB1.NE.JOB2) THEN
C Transform to biorthonormal orbital system
        IF (DoGSOR) Then
          Call FCopy(Trim(JBNAME(JOB2)),'JOBGS',ierr)
          Call DANAME(LUIPH,'JOBGS')
          IAD = 0
          Call IDAFile(LUIPH,2,ITOC15,30,IAD)
          IAD=ITOC15(2)
          Call DDAFile(LUIPH,1,CMO2,nCMO,IAD)
          Call DACLOS(LUIPH)
        End if !DoGSOR


        Call mma_allocate(TRA1,nTRA,Label='TRA1')
        Call mma_allocate(TRA2,nTRA,Label='TRA2')
        CALL FINDT(CMO1,CMO2,TRA1,TRA2)
        TRORB = .true.
      else
        TRORB = .false.
      end if

!     > check whether we do RASSI with an effective multi-state PT2 Hamiltonian
!     > whose eigenvectors are stored in Heff_evc
!     > i.e., we do not use mixed CI coefficients / MPS wave functions but rather mix the TDMs

      mstate_dens = job1.eq.job2.and.
     &              (allocated(Heff_evc(job1)%pc).or.
     &               allocated(Heff_evc(job1)%sc))
      mstate_dens = mstate_dens.and.if11
      mstate_dens = mstate_dens.and.qdpt2ev

      if(mstate_dens)then
        allocate(mstate_1pdens(nstat(job1),nstat(job1)))
        do i = 1, nstat(job1)
          do j = 1, nstat(job1)
            call mma_allocate(mstate_1pdens(i,j)%rtdm, NTDMZZ)
            call mma_allocate(mstate_1pdens(i,j)%stdm, NTDMZZ)
            call mma_allocate(mstate_1pdens(i,j)%wtdm,NTDMZZ)
            mstate_1pdens(i,j)%rtdm = 0; mstate_1pdens(i,j)%stdm    = 0
            mstate_1pdens(i,j)%wtdm = 0; mstate_1pdens(i,j)%overlap = 0
          end do
        end do
      end if

#ifdef _DMRG_
      dmrg_external%MPSrotated = trorb
#endif

C OBTAIN CORE ENERGY, FOCK MATRIX, AND TWO-ELECTRON INTEGRALS
C IN THE MIXED ACTIVE MO BASIS:
      ECORE=0.0D0
      IF (IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
       Call mma_allocate(FMO,nTDM1,Label='FMO')
       Call mma_allocate(TUVX,nTDM2,Label='TUVX')
       TUVX(:)=0.0D0
CTEST       write(*,*)'GTDMCTL calling TRINT.'
       CALL TRINT(CMO1,CMO2,ECORE,nTDM1,FMO,nTDM2,TUVX)
       ECORE=ENUC+ERFNUC+ECORE
CTEST       write(*,*)'GTDMCTL back from TRINT.'
CTEST       write(*,*)'ENUC  =',ENUC
CTEST       write(*,*)'ERFNUC=',ERFNUC
CTEST       write(*,*)'ECORE =',ECORE
      END IF

C In the calculation of matrix elements ( S1, S2 ), we will use
C the same Ms quantum numbers, if S1 and S2 differ by 0 or an int,
C else we will use Ms quantum numbers that differ by 1/2:
      IF( MOD(ABS(MPLET1-MPLET2),2).EQ.0 ) THEN
        MSPROJ1=MIN(MPLET1,MPLET2)-1
        MSPROJ2=MSPROJ1
      ELSE
        IF(MPLET1.GT.MPLET2) THEN
          MSPROJ2=MPLET2-1
          MSPROJ1=MSPROJ2+1
        ELSE
          MSPROJ1=MPLET1-1
          MSPROJ2=MSPROJ1+1
        END IF
      END IF

#ifdef _DMRG_
      !> set spin-up/spin-down # of electrons for target state(s)
      if(doDMRG)then
        dmrg_external%nalpha = (nacte1 + msproj1) / 2
        dmrg_external%nbeta  = (nacte1 - msproj1) / 2
      end if
#endif

C---------------  For all wave functions: ---------------------
C Define structures ('tables') pertinent all jobs.
C (Later, move this up before the GTDMCTL calls).
C These are at:
C IWORK(LPART)
C IWORK(LORBTAB)
C IWORK(LSSTAB)
      LPART=NEWPRTTAB(NSYM,NFRO,NISH,NRS1,NRS2,NRS3,NSSH,NDEL)
      IF(IPGLOB.GE.DEBUG) CALL PRPRTTAB(IWORK(LPART))

      LORBTAB=NEWORBTAB(IWORK(LPART))
      IF(IPGLOB.GE.DEBUG) CALL PRORBTAB(LORBTAB)

      LSSTAB=NEWSSTAB(LORBTAB)
      IF(IPGLOB.GE.DEBUG) CALL PRSSTAB(LSSTAB)

C Mapping from active spin-orbital to active orbital in external order.
C Note that these differ, not just because of the existence of two
C spin-orbitals for each orbital, but also because the active orbitals
C (external order) are grouped by symmetry and then RAS space, but the
C spin orbitals are grouped by subpartition.
      CALL GETMEM('OrbMap','Allo','Inte',LOMAP,NASORB)
      NASPRT=IWORK(LORBTAB+8)
      KSPART=IWORK(LORBTAB+9)
      LSPART=LORBTAB-1+KSPART
      KOINFO=19
      LOINFO=LORBTAB-1+KOINFO
      ISUM=0
      DO ISYM=1,NSYM
        NASHES(ISYM)=ISUM
        ISUM=ISUM+NASH(ISYM)
      END DO
      ISORB=0
      DO ISPART=1,NASPRT
        NO=IWORK(LSPART-1+ISPART)
        DO IO=1,NO
          ISORB=ISORB+1
C Orbital symmetry:
          ISYM=IWORK(LOINFO+1+(ISORB-1)*8)
C In-Symmetry orbital index:
          ISOIND=IWORK(LOINFO+2+(ISORB-1)*8)
C Subtract nr of inactive orbitals in that symmetry:
          IT=ISOIND-NISH(ISYM)
C Add nr of actives in earlier symmetries:
          ITABS=NASHES(ISYM)+IT
          IWORK(LOMAP-1+ISORB)=ITABS
        END DO
      END DO

C---------------    JOB1 wave functions: ---------------------
C Initialize SGUGA tables for JOB1 functions.
C These are structures stored in arrays:
C ISGSTR1,ICISTR1 and IXSTR1.

C Set variables in /RASDEF/, used by SGUGA codes, which define
C the SGUGA space of JOB1. General RAS:
      IF(WFTP1.EQ.'GENERAL ') THEN
        NRSPRT=3
        DO I=1,8
          NRAS(I,1)=NRS1(I)
          NRAS(I,2)=NRS2(I)
          NRAS(I,3)=NRS3(I)
        END DO
        NRASEL(1)=2*NRS1T-NHOL11
        NRASEL(2)=NACTE1-NELE31
        NRASEL(3)=NACTE1

        if(.not.doDMRG)then
          CALL SGINIT(NSYM,NACTE1,MPLET1,NRSPRT,NRAS,NRASEL,ISGSTR1)
          IF(IPGLOB.GT.DEBUG) THEN
            WRITE(6,*)'Split-graph structure for JOB1=',JOB1
            CALL SGPRINT(ISGSTR1)
          END IF
          CALL SGSVAL(ISGSTR1,NSYM,NASHT,LISM,NVERT,LDRT,
     &                LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
          CALL CXINIT(ISGSTR1,ICISTR1,IXSTR1)
          CALL CXSVAL(ICISTR1,IXSTR1,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &                LNOCSF,LIOCSF,NWALK,LICASE,
     &                MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &                LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C CI sizes, as function of symmetry, are now known.
          NCONF1=IWORK(LNCSF-1+LSYM1)
        else
          NCONF1=1
        end if
      ELSE
C Presently, the only other cases are HISPIN, CLOSED or EMPTY.
* Note: the HISPIN case may be buggy and is not used presently.
        NCONF1=1
      END IF
      CALL GETMEM('GTDMCI1','ALLO','REAL',LCI1,NCONF1)

C Still JOB1, define structures ('tables') pertinent to JOB1
C These are at:
C IWORK(LREST1)
C IWORK(LCNFTAB1)
C IWORK(LFSBTAB1)
C IWORK(LSPNTAB1)

      NPART=3
      NGAS=NPART
      DO ISYM=1,NSYM
        NGASORB(ISYM)=NRS1(ISYM)
        NGASORB(ISYM+NSYM)=NRS2(ISYM)
        NGASORB(ISYM+2*NSYM)=NRS3(ISYM)
      END DO

*PAM2008: The old MAXOP was far too generous:
*      MAXOP=NASHT
*PAM2008: MAXOP is determined by RAS restrictions:
* Preliminary ranges of nr of electrons:
      IE1MN=MAX(0,2*NRS1T-NHOL11)
      IE1MX=2*NRS1T
      IE3MN=MAX(0,NACTE1-IE1MX-2*NRS2T)
      IE3MX=MIN(2*NRS3T,NELE31)
      IE2MN=MAX(0,NACTE1-IE1MX-IE3MX)
      IE2MX=MIN(2*NRS2T,NACTE1-IE1MN-IE3MN)
* Preliminary NGASLIM:
      NGL11=IE1MX
      NGL21=IE1MN
      NGL12=IE2MX
      NGL22=IE2MN
      NGL13=IE3MX
      NGL23=IE3MN
* Start with MAXOP=0, then increase:
      MAXOP=0
* Loop over possible ranges:
      DO IE1=IE1MN,IE1MX
       IOP1=MIN(IE1,(2*NRS1T-IE1))
       IF(IOP1.GE.0) THEN
        DO IE3=IE3MN,IE3MX
         IOP3=MIN(IE3,(2*NRS3T-IE3))
         IF(IOP3.GE.0) THEN
          IE2=NACTE1-IE1-IE3
          IOP2=MIN(IE2,(2*NRS2T-IE2))
          IF(IOP2.GE.0) THEN
* Actually possible combination:
           MAXOP=MAX(MAXOP,IOP1+IOP2+IOP3)
           NGL11=MIN(NGL11,IE1)
           NGL21=MAX(NGL21,IE1)
           NGL12=MIN(NGL12,IE2)
           NGL22=MAX(NGL22,IE2)
           NGL13=MIN(NGL13,IE3)
           NGL23=MAX(NGL23,IE3)
          END IF
         END IF
        END DO
       END IF
      END DO
      NGASLIM(1,1)=NGL11
      NGASLIM(2,1)=NGL21
      NGASLIM(1,2)=NGL12
      NGASLIM(2,2)=NGL22
      NGASLIM(1,3)=NGL13
      NGASLIM(2,3)=NGL23

      if(.not.doDMRG)then
        IFORM=1
        MINOP=0
        LREST1=NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM)
        IF(IPGLOB.GE.DEBUG) CALL PRGASTAB(LREST1)

C At present, we will only annihilate, at most 2 electrons will
C be removed. This limits the possible MAXOP:
        MAXOP=MIN(MAXOP+1,NACTE1,NASHT)
        LCNFTAB1=NEWCNFTAB(NACTE1,NASHT,MINOP,MAXOP,LSYM1,NGAS,
     &                     NGASORB,NGASLIM,IFORM)
        IF(IPGLOB.GE.DEBUG) CALL PRCNFTAB(LCNFTAB1,100)

        LFSBTAB1=NEWFSBTAB(NACTE1,MSPROJ1,LSYM1,LREST1,LSSTAB)
        IF(IPGLOB.GE.DEBUG) CALL PRFSBTAB(IWORK(LFSBTAB1))
        NDET1=IWORK(LFSBTAB1+4)
        LSPNTAB1=NEWSCTAB(MINOP,MAXOP,MPLET1,MSPROJ1)
        IF (IPGLOB.GT.DEBUG) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
          CALL PRSCTAB(LSPNTAB1)
        END IF
      else
        NDET1 = 1 ! minimum to avoid runtime error
      end if
C---------------    JOB2 wave functions: ---------------------
C Initialize SGUGA tables for JOB2 functions.
C These are structures stored in arrays:
C ISGSTR2,ICISTR2 and IXSTR2.

C Set variables in /RASDEF/, used by SGUGA codes, which define
C the SGUGA space of JOB1. General RAS:
      IF(WFTP2.EQ.'GENERAL ') THEN
        NRSPRT=3
        DO I=1,8
          NRAS(I,1)=NRS1(I)
          NRAS(I,2)=NRS2(I)
          NRAS(I,3)=NRS3(I)
        END DO
        NRASEL(1)=2*NRS1T-NHOL12
        NRASEL(2)=NACTE2-NELE32
        NRASEL(3)=NACTE2

        IF(.not.doDMRG)then
          CALL SGINIT(NSYM,NACTE2,MPLET2,NRSPRT,NRAS,NRASEL,ISGSTR2)
          IF(IPGLOB.GT.DEBUG) THEN
            WRITE(6,*)'Split-graph structure for JOB2=',JOB2
            CALL SGPRINT(ISGSTR2)
          END IF
          CALL SGSVAL(ISGSTR2,NSYM,NASHT,LISM,NVERT,LDRT,
     &                LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
          CALL CXINIT(ISGSTR2,ICISTR2,IXSTR2)
          CALL CXSVAL(ICISTR2,IXSTR2,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &                LNOCSF,LIOCSF,NWALK,LICASE,
     &                MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &                LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C CI sizes, as function of symmetry, are now known.
          NCONF2=IWORK(LNCSF-1+LSYM2)
        else
          NCONF2=1
        end if
      ELSE
C Presently, the only other cases are HISPIN, CLOSED or EMPTY.
* Note: the HISPIN case may be buggy and is not used presently.
        NCONF2=1
      END IF
      CALL GETMEM('GTDMCI2','ALLO','REAL',LCI2,NCONF2)
      If (DoGSOR) Then
        CALL GETMEM('GTDMCI2_o','ALLO','REAL',LCI2_o,NCONF2)
      end if!DoGSOR

      NPART=3
      NGAS=NPART
      DO ISYM=1,NSYM
        NGASORB(ISYM)=NRS1(ISYM)
        NGASORB(ISYM+NSYM)=NRS2(ISYM)
        NGASORB(ISYM+2*NSYM)=NRS3(ISYM)
      END DO
*PAM2008: The old MAXOP was far too generous:
*      MAXOP=NASHT
*PAM2008: MAXOP is determined by RAS restrictions:
* Preliminary ranges of nr of electrons:
      IE1MN=MAX(0,2*NRS1T-NHOL12)
      IE1MX=2*NRS1T
      IE3MN=MAX(0,NACTE2-IE1MX-2*NRS2T)
      IE3MX=MIN(2*NRS3T,NELE32)
      IE2MN=MAX(0,NACTE2-IE1MX-IE3MX)
      IE2MX=MIN(2*NRS2T,NACTE2-IE1MN-IE3MN)
* Preliminary NGASLIM:
      NGL11=IE1MX
      NGL21=IE1MN
      NGL12=IE2MX
      NGL22=IE2MN
      NGL13=IE3MX
      NGL23=IE3MN
* Start with MAXOP=0, then increase:
      MAXOP=0
* Loop over possible ranges:
      DO IE1=IE1MN,IE1MX
       IOP1=MIN(IE1,(2*NRS1T-IE1))
       IF(IOP1.GE.0) THEN
        DO IE3=IE3MN,IE3MX
         IOP3=MIN(IE3,(2*NRS3T-IE3))
         IF(IOP3.GE.0) THEN
          IE2=NACTE2-IE1-IE3
          IOP2=MIN(IE2,(2*NRS2T-IE2))
          IF(IOP2.GE.0) THEN
* Actually possible combination:
           MAXOP=MAX(MAXOP,IOP1+IOP2+IOP3)
           NGL11=MIN(NGL11,IE1)
           NGL21=MAX(NGL21,IE1)
           NGL12=MIN(NGL12,IE2)
           NGL22=MAX(NGL22,IE2)
           NGL13=MIN(NGL13,IE3)
           NGL23=MAX(NGL23,IE3)
          END IF
         END IF
        END DO
       END IF
      END DO
      NGASLIM(1,1)=NGL11
      NGASLIM(2,1)=NGL21
      NGASLIM(1,2)=NGL12
      NGASLIM(2,2)=NGL22
      NGASLIM(1,3)=NGL13
      NGASLIM(2,3)=NGL23

      if(.not.dodmrg)then
        LREST2=NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM)
        IF(IPGLOB.GE.DEBUG) CALL PRGASTAB(LREST2)

        IFORM=1
        MINOP=0
C At present, we will only annihilate. This limits the possible MAXOP:
        MAXOP=MIN(MAXOP+1,NACTE2,NASHT)
        LCNFTAB2=NEWCNFTAB(NACTE2,NASHT,MINOP,MAXOP,LSYM2,NGAS,
     &                     NGASORB,NGASLIM,IFORM)
        IF(IPGLOB.GE.DEBUG) CALL PRCNFTAB(LCNFTAB2,100)

        LFSBTAB2=NEWFSBTAB(NACTE2,MSPROJ2,LSYM2,LREST2,LSSTAB)
        IF(IPGLOB.GE.DEBUG) CALL PRFSBTAB(IWORK(LFSBTAB2))
        NDET2=IWORK(LFSBTAB2+4)
        LSPNTAB2=NEWSCTAB(MINOP,MAXOP,MPLET2,MSPROJ2)
        IF (IPGLOB.GT.DEBUG) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
          CALL PRSCTAB(LSPNTAB2)
        END IF
      else
        NDET2 = 1 ! minimum to avoid runtime error
      end if
C-------------------------------------------------------------
      CALL GETMEM('GTDMDET1','ALLO','REAL',LDET1,NDET1)
      CALL GETMEM('GTDMDET2','ALLO','REAL',LDET2,NDET2)

C Loop over the states of JOBIPH nr JOB1
C Disk address for writing to scratch file is IDWSCR.
      IDWSCR=0
      DO IST=1,NSTAT(JOB1)
        ISTATE=ISTAT(JOB1)-1+IST

        if(.not.doDMRG)then
C Read ISTATE wave function
          IF(WFTP1.EQ.'GENERAL ') THEN
            CALL READCI(ISTATE,ISGSTR1,ICISTR1,NCONF1,WORK(LCI1))
          ELSE
            WORK(LCI1)=1.0D0
          END IF
          CALL DCOPY_(NDET1,[0.0D0],0,WORK(LDET1),1)
C         Transform to bion basis, Split-Guga format
          If (TrOrb) CALL CITRA (WFTP1,ISGSTR1,ICISTR1,IXSTR1,LSYM1,
     &                           TRA1,NCONF1,Work(LCI1))
          CALL PREPSD(WFTP1,ISGSTR1,ICISTR1,LSYM1,
     &                IWORK(LCNFTAB1),IWORK(LSPNTAB1),
     &                IWORK(LSSTAB),IWORK(LFSBTAB1),NCONF1,WORK(LCI1),
     &                WORK(LDET1))

C Write out the determinant expansion to disk.
          IDDET1(ISTATE)=IDWSCR
          CALL  DDAFILE(LUSCR,1,WORK(LDET1),NDET1,IDWSCR)
        else
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 LROOT(ISTATE),
     &                 LSYM1,
     &                 MPLET1,
     &                 MSPROJ1,
     &                 NACTE1,
     &                 TRA1,
     &                 NTRA,
     &                 NISH,
     &                 NASH,
     &                 NOSH,
     &                 NSYM,
     &                 6,
     &                 ISTATE,
     &                 job1,
     &                 ist
     &                )
#endif
        end if
      END DO

      If (DoGSOR) Then
        CALL GETMEM('Theta1','ALLO','REAL',LTheta1,NCONF2)
        CALL DCOPY_(NCONF2,[0.0D0],0,WORK(LTheta1),1)
      End If

C-------------------------------------------------------------

C Loop over the states of JOBIPH nr JOB2
      job2_loop: DO JST=1,NSTAT(JOB2)

        JSTATE=ISTAT(JOB2)-1+JST

        if(.not.doDMRG)then
C Read JSTATE wave function
          IF(WFTP2.EQ.'GENERAL ') THEN
            CALL READCI(JSTATE,ISGSTR2,ICISTR2,NCONF2,WORK(LCI2))
          ELSE
            WORK(LCI2)=1.0D0
          END IF
          If(DoGSOR) Then
            CALL DCOPY_(NCONF2,Work(LCI2),1,WORK(LCI2_o),1)
          End If
          CALL DCOPY_(NDET2,[0.0D0],0,WORK(LDET2),1)
C         Transform to bion basis, Split-Guga format
          If (TrOrb) CALL CITRA (WFTP2,ISGSTR2,ICISTR2,IXSTR2,LSYM2,
     &                           TRA2,NCONF2,Work(LCI2))
          CALL PREPSD(WFTP2,ISGSTR2,ICISTR2,LSYM2,
     &                IWORK(LCNFTAB2),IWORK(LSPNTAB2),
     &                IWORK(LSSTAB),IWORK(LFSBTAB2),NCONF2,WORK(LCI2),
     &                WORK(LDET2))

        else
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 lroot(JSTATE),
     &                 LSYM2,
     &                 MPLET2,
     &                 MSPROJ2,
     &                 NACTE2,
     &                 TRA2,
     &                 NTRA,
     &                 NISH,
     &                 NASH,
     &                 NOSH,
     &                 NSYM,
     &                 6,
     &                 JSTATE,
     &                 job2,
     &                 jst
     &                )
#endif
        end if

C Loop over the states of JOBIPH nr JOB1
        job1_loop: DO IST=1,NSTAT(JOB1)
          ISTATE=ISTAT(JOB1)-1+IST
        IF(ISTATE.LT.JSTATE) cycle
C Entry into monitor: Status line
        WRITE(STLNE1,'(A6)') 'RASSI:'
        WRITE(STLNE2,'(A33,I5,A5,I5)')
     &      'Trans. dens. matrices for states ',ISTATE,' and ',JSTATE
        Call StatusLine(STLNE1,STLNE2)

C Read ISTATE wave function from disk
#ifdef _DMRG_
      if(.not.doDMRG)then
#endif
      IDRSCR=IDDET1(ISTATE)
      CALL  DDAFILE(LUSCR,2,WORK(LDET1),NDET1,IDRSCR)
#ifdef _DMRG_
      end if
#endif

       if(doGSOR) then
         if(JOB1.ne.JOB2) then
           ST_TOT = NSTAT(JOB2)
           Dot_prod = 0
           Dot_prod = DDOT_(NCONF2,Work(LCI1),1,Work(LCI2),1)
           Call DAXPY_(NCONF2,Dot_prod,Work(LCI2_o),1,Work(LTHETA1),1)
         end if
       end if

C Calculate whatever type of GTDM that was requested, unless
C it is known to be zero.
      HZERO=0.0D0
      HONE =0.0D0
      HTWO =0.0D0

      SIJ=0.0D0
      DYSAMP=0.0D0

! +++ J. Norell 12/7 - 2018
C Dyson amplitudes:
C DYSAMP = D_ij for states i and j
C DYSCOF = Active orbital coefficents of the DO
      IF ((IF10.or.IF01).and.DYSO) THEN
        CALL DYSON(IWORK(LFSBTAB1),
     &            IWORK(LFSBTAB2),IWORK(LSSTAB),
     &            WORK(LDET1),WORK(LDET2),
     &            IF10,IF01,
     &            DYSAMP,DYSCOF)

C Write Dyson orbital coefficients in AO basis to disk.
        IF (DYSAMP.GT.1.0D-6) THEN
C In full biorthonormal basis:
         CALL MKDYSAB(DYSCOF,DYSAB)
C In AO basis:
         CALL MKDYSZZ(CMO1,DYSAB,DYSZZ)
        IF (DYSO) THEN
          SFDYS(:,JSTATE,ISTATE)=DYSZZ(:)
          SFDYS(:,ISTATE,JSTATE)=DYSZZ(:)
        END IF
        DYSZZ(:)=0.0D0
       END IF ! AMP THRS
      END IF ! IF01 IF10
      DYSAMPS(ISTATE,JSTATE)=DYSAMP
      DYSAMPS(JSTATE,ISTATE)=DYSAMP
! +++

C General 1-particle transition density matrix:
      IF (IF11) THEN
        CALL MKTDM1(LSYM1,MPLET1,MSPROJ1,IWORK(LFSBTAB1),
     &            LSYM2,MPLET2,MSPROJ2,IWORK(LFSBTAB2),IWORK(LSSTAB),
     &            IWORK(LOMAP),WORK(LDET1),WORK(LDET2),SIJ,NASHT,
     &            TRAD,TRASD,WERD,ISTATE,
     &            JSTATE,job1,job2,ist,jst)
C Calculate Natural Transition Orbital (NTO):
        IF (IFNTO) THEN
         IF (job1.ne.job2) THEN
           DoNTO=.true.
         Else
           DoNTO=.false.
         End If
         IF (DoNTO) Then
          Call NTOCalc(job1,job2,ISTATE,JSTATE,TRAD,TRASD,MPLET1)
          write(6,*) 'ntocalculation finished'
         End If
        End If
C End of Calculating NTO

        IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
C Compute 1-electron contribution to Hamiltonian matrix element:
        HONE=DDOT_(NTRAD,TRAD,1,FMO,1)
        END IF


C             Write density 1-matrices in AO basis to disk.
            IF(NATO.OR.(NPROP.GT.0))THEN

              iEmpty=0
              !> regular-TDM
              CALL MKTDAB(SIJ,TRAD,TDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,TDMAB,TDMZZ,iRC)
              If (iRC.eq.1) iEmpty=1

              !> spin-TDM
              CALL MKTDAB(0.0D0,TRASD,TSDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,TSDMAB,TSDMZZ,iRC)
              If (iRC.eq.1) iEmpty=iEmpty+2

              !> WE-reduced TDM''s of triplet type:
              CALL MKTDAB(0.0D0,WERD,WDMAB,iRC)
              !> transform to AO basis
              CALL MKTDZZ(CMO1,CMO2,WDMAB,WDMZZ,iRC)
              If (iRC.eq.1) iEmpty=iEmpty+4

              if(.not.mstate_dens)then

                IF(SaveDens) THEN
*C Transition density matrices, TDMZZ, in AO or MO basis.
*C WDMZZ similar, but WE-reduced 'triplet' densities.
                  ij=ISTATE*(iSTATE-1)/2 + JSTATE
                  jDisk_TDM(1,ij)=IDISK
                  jDisk_TDM(2,ij)=iEmpty
                  iOpt=1
                  iGo=7
                  If (AO_Mode) Then
                     CALL dens2file(TDMZZ,TSDMZZ,WDMZZ,nTDMZZ,LUTDM,
     &                              IDISK,iEmpty,iOpt,iGo,IState,jState)
                  Else
                     iEmpty=0
                     TRAD(nTrad+1)=SIJ
                     iRC=0
                     If(DDot_(nTRAD+1,TRAD,1,TRAD,1).gt.0.0D0) iRC=1
                     If (iRC.eq.1) iEmpty=1
*
                     TRASD(nTrad+1)=0.0D0
                     iRC=0
                     If(DDot_(nTRAD+1,TRASD,1,TRASD,1).gt.0.0D0) iRC=1
                     If (iRC.eq.1) iEmpty=iEmpty+2
*
                     WERD(nTrad+1)=0.0D0
                     iRC=0
                     If(DDot_(nTRAD+1,WERD,1,WERD,1).gt.0.0D0) iRC=1
                     If (iRC.eq.1) iEmpty=iEmpty+4
*
                     CALL dens2file(TRAD,TRASD,WERD,nTRAD+1,LUTDM,
     &                              IDISK,iEmpty,iOpt,iGo,ISTATE,JSTATE)
                  End If
                END IF
                !> calculate property matrix elements
                CALL PROPER(PROP,ISTATE,JSTATE,TDMZZ,WDMZZ)
              else

!               > scale rdm elements with eigenvector coefficients of Heff of a multi-state (PT2) Hamiltonian
!               > accumulate data first and run PROPER and other utility routines later
                do i = 1, nstat(job1)
                  do j = 1, nstat(job2)

                    if(i < j) cycle

                    if(qdpt2sc)then
                      fac1 = Heff_evc(job1)%sc(ist,i)
                      fac2 = Heff_evc(job2)%sc(jst,j)
                    else
                      fac1 = Heff_evc(job1)%pc(ist,i)
                      fac2 = Heff_evc(job2)%pc(jst,j)
                    end if

                    !> regular-TDM
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          tdmzz,1,
     &                          mstate_1pdens(i,j)%rtdm,1
     &                         )
                    !> spin-TDM
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          tsdmzz,1,
     &                          mstate_1pdens(i,j)%stdm,1
     &                         )
                    !> WE-reduced TDM''s of triplet type:
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          wdmzz,1,
     &                          mstate_1pdens(i,j)%wtdm,1
     &                         )
                    !> overlap
                    mstate_1pdens(i,j)%overlap=
     &                mstate_1pdens(i,j)%overlap+fac1*fac2*sij
                  end do
                end do
              end if
            END IF

          ELSE ! IF11

            !> overlap
            IF (IF00) THEN
#ifdef _DMRG_
              if(.not.doDMRG)then
#endif
                SIJ=OVERLAP_RASSI(IWORK(LFSBTAB1),IWORK(LFSBTAB2),
     &                            WORK(LDET1),WORK(LDET2))
#ifdef _DMRG_
              else
                if (doMPSSICheckpoints) then

                  if (dmrg_external%MPSrotated) then
                    write(mplet1s,'(I3)')  MPLET1-1
                    write(msproj1s,'(I3)')  MSPROJ1

                    checkpoint1_2u1 = qcm_group_names(job1)%states(ist)
     &(1:len_trim(qcm_group_names(job1)%states(ist))-3)
     &//"."//trim(adjustl(mplet1s))//"."//trim(adjustl(msproj1s))//".h5"
                    checkpoint2_2u1 = qcm_group_names(job2)%states(jst)
     &(1:len_trim(qcm_group_names(job2)%states(jst))-3)
     &//"."//trim(adjustl(mplet1s))//"."//trim(adjustl(msproj1s))//".h5"

                    call dmrg_interface_ctl(
     &                                      task       = 'overlapU',
     &                                      energy     = sij,
     &                                      checkpoint1=checkpoint1_2u1,
     &                                      checkpoint2=checkpoint2_2u1
     &                                     )
                  else
                    call dmrg_interface_ctl(
     &                                      task        = 'overlap ',
     &                                      energy      = sij,
     &                                      checkpoint1 =
     &                               qcm_group_names(job1)%states(ist),
     &                                      checkpoint2 =
     &                               qcm_group_names(job2)%states(jst)
     &                                     )
                  end if
                else
!                 > Leon: TODO: Add possibility to calculate overlap of rotated MPS without using checkpoint names
                  call dmrg_interface_ctl(
     &                               task   = 'overlap ',
     &                               energy = sij,
     &                               state  = LROOT(istate),
     &                               stateL = LROOT(jstate)
     &                              )
                end if
              end if !doDMRG
#endif
            END IF ! IF00
          END IF ! IF11

          if(.not.mstate_dens)then
            OVLP(ISTATE,JSTATE)=SIJ
            OVLP(JSTATE,ISTATE)=SIJ
          end if

          !> General 2-particle transition density matrix:
          IF (IF22) THEN
            CALL MKTDM2(LSYM1,MPLET1,MSPROJ1,IWORK(LFSBTAB1),
     &                  LSYM2,MPLET2,MSPROJ2,IWORK(LFSBTAB2),
     &                  IWORK(LSSTAB),IWORK(LOMAP),
     &                  WORK(LDET1),WORK(LDET2),NTDM2,TDM2,
     &                  ISTATE,JSTATE,job1,job2,ist,jst)

!           > Compute 2-electron contribution to Hamiltonian matrix element:
            IF(IFTWO.AND.(MPLET1.EQ.MPLET2))
     &      HTWO=DDOT_(NTDM2,TDM2,1,TUVX,1)

          END IF ! IF22

          !> PAM 2011 Nov 3, writing transition matrices if requested
          IF ((IFTRD1.or.IFTRD2).and..not.mstate_dens) THEN
            call trd_print(ISTATE, JSTATE, IFTRD2.AND.IF22,
     &                    TDMAB,TDM2,CMO1,CMO2)
          END IF

          IF(IFHAM.AND..NOT.(IFHEXT.or.IFHEFF.or.IFEJOB))THEN
            HZERO              = ECORE*SIJ
            HIJ                = HZERO+HONE+HTWO
            HAM(ISTATE,JSTATE) = HIJ
            HAM(JSTATE,ISTATE) = HIJ

         !SI-PDFT related code for "second_time" case
          if(second_time) then
            Energies(:) =0.0d0
            CALL DANAME(LUIPH,'JOBGS')
            IAD = 0
            Call IDAFILE(LUIPH,2,ITOC15,30,IAD)
            IAD=ITOC15(6)
            Call DDAFILE(LUIPH,2,Energies,NSTAT(JOB1),IAD)
            do i=1,NSTAT(JOB1)
              HAM(i,i) = Energies(i)
            end do
            Call DACLOS(LUIPH)
          end if


            IF(IPGLOB.GE.DEBUG) THEN
              WRITE(6,'(1x,a,2I5)')' ISTATE, JSTATE:',ISTATE,JSTATE
              WRITE(6,'(1x,a,f16.8)')' HZERO=',HZERO
              WRITE(6,'(1x,a,f16.8)')' HONE =',HONE
              WRITE(6,'(1x,a,f16.8)')' HTWO =',HTWO
              WRITE(6,'(1x,a,f16.8)')' HIJ  =',HIJ
            END IF
          END IF

        END DO job1_loop

      END DO job2_loop

      IF(DoGSOR) then
        if(job1.ne.job2) then
        Norm_Fac = 0.0d0
        dot_prod = DDOT_(NCONF2,Work(LTheta1),1,Work(LTheta1),1)
        Norm_Fac = 1d0/sqrt(dot_prod)
        Call DSCAL_(NCONF2,Norm_Fac,Work(LTheta1),1)

      !Write theta1 to file.
        LUCITH=87
        LUCITH=IsFreeUnit(LUCITH)
        !Open(unit=87,file='CI_THETA', action='write',iostat=ios)
        Call Molcas_Open(LUCITH,'CI_THETA')
        do i=1,NCONF2
          write(LUCITH,*) Work(LTheta1-1+i)
        end do
        Close(LUCITH)

       !Now we need to build the other states.
        DO JST=2,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          CALL READCI(JSTATE,ISGSTR2,ICISTR2,NCONF2,WORK(LCI2))
          Call DCOPY_(NCONF2,Work(LCI2),1,Work(LCI2_o),1)
          CALL DCOPY_(NDET2,[0.0D0],0,WORK(LDET2),1)
          If (TrOrb) CALL CITRA (WFTP2,ISGSTR2,ICISTR2,IXSTR2,LSYM2,
     &                           TRA2,NCONF2,Work(LCI2))
          CALL PREPSD(WFTP2,ISGSTR2,ICISTR2,LSYM2,
     &                IWORK(LCNFTAB2),IWORK(LSPNTAB2),
     &                IWORK(LSSTAB),IWORK(LFSBTAB2),NCONF2,WORK(LCI2),
     &                WORK(LDET2))

          CALL GETMEM('ThetaN','ALLO','REAL',LThetaN,NCONF2)
          CALL DCOPY_(NCONF2,Work(LCI2_o),1,WORK(LThetaN),1)
          Norm_Fac = DDOT_(NCONF2,Work(LTheta1),1,Work(LCI2_o),1)
          Call DAXPY_(NCONF2,-Norm_Fac,Work(LTHETA1),1,Work(LThetaN),1)

          LUCITH=IsFreeUnit(LUCITH)
          Call Molcas_Open(LUCITH,'CI_THETA')
          !Open(unit=87,file='CI_THETA', action='read',iostat=ios)
          if(JST-1.ge.2) then
            do i=1,NCONF2
              Read(LUCITH,*) dummy
            end do
          end if
          CALL GETMEM('ThetaM','ALLO','REAL',LThetaM,NCONF2)
          DO IST=2,JST-1
            CALL DCOPY_(NCONF2,[0.0D0],0,WORK(LThetaM),1)
            !Read in previous theta vectors
            do i=1,NCONF2
              Read(LUCITH,*) Work(LThetaM-1+i)
            end do
            Dot_prod = DDOT_(NCONF2,Work(LThetaM),1,Work(LCI2_o),1)
           Call DAXPY_(NCONF2,-Dot_prod,Work(LThetaM),1,Work(LThetaN),1)


          END DO
          Close(LUCITH)
          !Normalize
          dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LThetaN),1)
          Norm_Fac = 1d0/sqrt(dot_prod)
          Call DSCAL_(NCONF2,Norm_Fac,Work(LThetaN),1)

        !dot_prod = DDOT_(NCONF2,Work(LTheta1),1,Work(LTheta1),1)
        !dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LTheta1),1)
        !dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LThetaN),1)

          !Write to file
          LUCITH=IsFreeUnit(LUCITH)
          Call Molcas_Open(LUCITH,'CI_THETA')
          Call Append_file(LUCITH)
          !Open(unit=87,file='CI_THETA', position='append',iostat=ios,
!    &    action='write')
          do i=1,nConf2
            write(LUCITH,*) Work(LThetaN-1+i)
          end do
          close(LUCITH)
          !Deallocate
          CALL GETMEM('ThetaN','FREE','REAL',LThetaN,NCONF2)
        END DO
!Copy to new IPH file
        LUCITH=IsFreeUnit(LUCITH)
        Call Molcas_Open(LUCITH,'CI_THETA')
!       Open(unit=87,file='CI_THETA',iostat=ios,
!    &    action='read')
        CALL DANAME(LUIPHn,'JOBGS')
        IAD = 0
        Call IDAFILE(LUIPHn,2,ITOC15,30,IAD)
        IAD=ITOC15(4)
        do i=1,ISTAT(JOB1)-1
         CALL DCOPY_(NCONF2,[0.0D0],0,WORK(LThetaM),1)
         do j=1,nCONF2
           read(LUCITH,*) Work(LThetaM-1+i)
         end do
         Call DDafile(LUIPHn,1,Work(LThetaM),nCONF2,IAD)
        end do

       IAD = ITOC15(4)
       CALL DCOPY_(NCONF2,[0.0D0],0,WORK(LThetaM),1)
       Call DDAFILE(LUIPHn,2,Work(LThetaM),nCONF2,IAD)
       Call DDAFILE(LUIPHn,2,Work(LThetaM),nCONF2,IAD)

       Close(LUCITH)
       Call DACLOS(LUIPHn)
       CALL GETMEM('ThetaM','FREE','REAL',LThetaM,NCONF2)
       end if
       CALL GETMEM('Theta1','FREE','REAL',LTheta1,NCONF1)
      end if!DoGSOR



#ifdef _DMRG_
      IF(IPGLOB.GE.DEBUG) THEN
         write(6,*) 'full SF-HAMILTONIAN '
         write(6,*) 'dimension: ',nstate**2
         call pretty_print_util(HAM,1,nstate,1,nstate,
     &                          nstate,nstate,1,6)
      END IF
#endif

!     > create actual property data and put everything to file (if requested) in case of using eigenvectors of a multi-state (PT2) Hamiltonian
      if(mstate_dens)then
        DO JST=1,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          DO IST=1,NSTAT(JOB1)
            ISTATE=ISTAT(JOB1)-1+IST
            if(istate < jstate) cycle

            ovlp(istate,jstate) = mstate_1pdens(ist,jst)%overlap
            ovlp(jstate,istate) = mstate_1pdens(ist,jst)%overlap

            call prpdata_mspt2_eigenvectors(
     &                                      mstate_1pdens(ist,jst)%rtdm,
     &                                      mstate_1pdens(ist,jst)%stdm,
     &                                      mstate_1pdens(ist,jst)%wtdm,
     &                                      prop,
     &                                      nprop,
     &                                      nstate,
     &                                      istate,
     &                                      jstate,
     &                                      ntdmzz,
     &                                      iDisk_TDM(JSTATE,ISTATE,1),
     &                                      iDisk_TDM(JSTATE,ISTATE,2),
     &                                      lutdm,
     &                                      (sonatnstate.gt.0),
     &                                      if11.and.(lsym1.eq.lsym2)
     &                                     )
          end do
        end do
      end if

      IF(WFTP1.EQ.'GENERAL ') THEN
        if(.not.doDMRG)then
          CALL CXCLOSE(ISGSTR1,ICISTR1,IXSTR1)
          CALL SGCLOSE(ISGSTR1)
        end if
      END IF
      IF(WFTP2.EQ.'GENERAL ') THEN
        if(.not.doDMRG)then
          CALL CXCLOSE(ISGSTR2,ICISTR2,IXSTR2)
          CALL SGCLOSE(ISGSTR2)
        end if
      END IF

      IF(JOB1.NE.JOB2) THEN
        Call mma_deallocate(TRA1)
        Call mma_deallocate(TRA2)
      END IF
      CALL GETMEM('GTDMDET1','FREE','REAL',LDET1,NDET1)
      CALL GETMEM('GTDMDET2','FREE','REAL',LDET2,NDET2)
      CALL GETMEM('GTDMCI2','FREE','REAL',LCI2,NCONF2)
      If (DoGSOR) CALL GETMEM('GTDMCI2_o','FREE','REAL',LCI2_o,NCONF2)
      IF (.NOT.NONA) CALL GETMEM('GTDMCI1','FREE','REAL',LCI1,NCONF1)
      if(.not.doDMRG)then
        CALL KILLSCTAB(LSPNTAB1)
        CALL KILLSCTAB(LSPNTAB2)
      end if
      IF ((IF10.or.IF01).and.DYSO) THEN
        Call mma_deallocate(DYSCOF)
        Call mma_deallocate(DYSAB)
        Call mma_deallocate(DYSZZ)
      END IF
      IF (IF11) THEN
        Call mma_deallocate(TRAD)
        Call mma_deallocate(TRASD)
        Call mma_deallocate(WERD)
        IF(NATO.OR.NPROP.GT.0) THEN
          Call mma_deallocate(TDMAB)
          Call mma_deallocate(TSDMAB)
          Call mma_deallocate(WDMAB)
          Call mma_deallocate(TDMZZ)
          Call mma_deallocate(TSDMZZ)
          Call mma_deallocate(WDMZZ)
        END IF
      END IF
      IF (IF22) Call mma_deallocate(TDM2)

      IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
        Call mma_deallocate(FMO)
        Call mma_deallocate(TUVX)
      END IF

      Call mma_deallocate(CMO2)
      Call mma_deallocate(CMO1)
      CALL KILLOBJ(LPART)
      CALL KILLOBJ(LORBTAB)
      CALL KILLOBJ(LSSTAB)
      if(.not.doDMRG)then
        CALL KILLOBJ(LREST1)
        CALL KILLOBJ(LREST2)
        CALL KILLOBJ(LCNFTAB1)
        CALL KILLOBJ(LCNFTAB2)
        CALL KILLOBJ(LFSBTAB1)
        CALL KILLOBJ(LFSBTAB2)
      end if
      CALL GETMEM('OrbMap','Free','Inte',LOMAP,NASORB)

      !> release memory
      if(mstate_dens)then
        do i = 1, nstat(job1)
          do j = 1, nstat(job1)
            if(allocated(mstate_1pdens(i,j)%rtdm))
     &      call mma_deallocate(mstate_1pdens(i,j)%rtdm)
            if(allocated(mstate_1pdens(i,j)%stdm))
     &      call mma_deallocate(mstate_1pdens(i,j)%stdm)
            if(allocated(mstate_1pdens(i,j)%wtdm))
     &      call mma_deallocate(mstate_1pdens(i,j)%wtdm)
          end do
        end do
        if(allocated(mstate_1pdens)) deallocate(mstate_1pdens)
      end if
#ifdef _TIME_GTDM_
      Call CWTime(TCpu2,TWall2)
      write(6,*) 'Time for GTDM : ',TCpu2-TCpu1,TWall2-TWall1
#endif

      CALL QEXIT(ROUTINE)
      RETURN
      END
