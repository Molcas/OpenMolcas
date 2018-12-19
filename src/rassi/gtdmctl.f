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
      SUBROUTINE GTDMCTL(PROP,JOB1,JOB2,OVLP,DYSAMPS,SFDYS,NZ,
     &     HAM,IDDET1)

      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
      use qcmaquis_info
#endif
      use mspt2_eigenvectors

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
      DIMENSION OVLP(NSTATE,NSTATE),HAM(NSTATE,NSTATE)
      DIMENSION DYSAMPS(NSTATE,NSTATE)
      DIMENSION SFDYS(NZ,NSTATE,NSTATE)
      DIMENSION IDDET1(NSTATE)
      LOGICAL IF00, IF10,IF01,IF20,IF11,IF02,IF21,IF12,IF22
      LOGICAL IFTWO,TRORB
      CHARACTER*8 WFTP1,WFTP2
      CHARACTER*6 STLNE1
      CHARACTER*3 NUM1,NUM2
      CHARACTER*12 FNM
      CHARACTER*48 STLNE2
* PAM 2011 Nov 3, added write buffer WBUF:
      DIMENSION WBUF(5)
      Real*8 Energies(1:20)
      Integer IAD,LUIPHn,lThetaM
      Real*8 Norm_fac

      type mixed_1pdensities
        real*8              :: overlap
        real*8, allocatable :: rtdm(:)
        real*8, allocatable :: stdm(:)
        real*8, allocatable :: wtdm(:)
      end type

      type(mixed_1pdensities), allocatable :: mstate_1pdens(:,:)

      logical               :: mstate_dens
      real*8                :: fac1, fac2

#ifdef _DMRG_
      ! strings for conversion of the qcmaquis h5 checkpoint names from 2u1 to su2u1
      character(len=3) :: mplet1s, msproj1s
      ! new checkpoint names
      character(len=2300) :: checkpoint1_2u1,checkpoint2_2u1
#else
      logical             :: doDMRG = .false.
#endif
#include "SysDef.fh"

      CALL QENTER(ROUTINE)
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
      CALL GETMEM('GTDMCMO1','ALLO','REAL',LCMO1,NCMO)
      CALL RDCMO(JOB1,WORK(LCMO1))
      CALL GETMEM('GTDMCMO2','ALLO','REAL',LCMO2,NCMO)


      CALL RDCMO(JOB2,WORK(LCMO2))
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
        CALL GETMEM('DYSCOF','Allo','Real',LDYSCOF,NASORB)
        ! Number of inactive+active orbitals
        NDYSAB = NASHT+NISHT
        CALL GETMEM('DYSAB','Allo','Real',LDYSAB,NDYSAB)
        ! Number of atomic / basis functions
        NDYSZZ = NZ
        CALL GETMEM('DYSZZ','Allo','Real',LDYSZZ,NDYSZZ)
        DO NDUM=1,NDYSZZ
         WORK(LDYSZZ+NDUM-1)=0.0D0
        END DO
      END IF
! +++

C Transition density matrices, TDMAB is for active biorthonormal
C orbitals only, while TDMZZ is in the fixed AO basis.
C WDMAB, WDMZZ similar, but WE-reduced 'triplet' densities.
      LTDMAB=LNILPT
      LTDMZZ=LNILPT
      LTSDMAB=LNILPT
      LTSDMZZ=LNILPT
      LWDMAB=LNILPT
      LWDMZZ=LNILPT
      IF(IF11.AND.(NATO.OR.NPROP.GT.0)) THEN
        CALL GETMEM('TDMAB','Allo','Real',LTDMAB,NTDMAB)
        CALL GETMEM('TDMZZ','Allo','Real',LTDMZZ,NTDMZZ)
        CALL GETMEM('TSDMAB','Allo','Real',LTSDMAB,NTDMAB)
        CALL GETMEM('TSDMZZ','Allo','Real',LTSDMZZ,NTDMZZ)
        CALL GETMEM('WDMAB','Allo','Real',LWDMAB,NTDMAB)
        CALL GETMEM('WDMZZ','Allo','Real',LWDMZZ,NTDMZZ)
      END IF

      LSPD1=LNILPT
      LTRAD=LNILPT
      LTRASD=LNILPT
      LWERD=LNILPT
      IF (IF11) THEN
        NSPD1=NASORB**2
        CALL GETMEM('SPD1','Allo','Real',LSPD1,NSPD1)
        NTRAD=NASHT**2
        CALL GETMEM('TRAD','Allo','Real',LTRAD,NTRAD)
        NTRASD=NASHT**2
        CALL GETMEM('TRASD','Allo','Real',LTRASD,NTRASD)
        NWERD=NASHT**2
        CALL GETMEM('WERD','Allo','Real',LWERD,NWERD)
      END IF
      LSPD2=LNILPT
      LTDM2=LNILPT
      IF (IF22) THEN
        NASGEM=(NASORB*(NASORB-1))/2
        NSPD2=NASGEM**2
        CALL GETMEM('SPD2','Allo','Real',LSPD2,NSPD2)
        CALL GETMEM('TDM2','Allo','Real',LTDM2,NTDM2)
      END IF

      LTRA1=LNILPT
      LTRA2=LNILPT
      IF(JOB1.NE.JOB2) THEN
C Transform to biorthonormal orbital system
        IF (DoGSOR) Then
          Call FCopy(Trim(JBNAME(JOB2)),'JOBGS',ierr)
          Call DANAME(LUIPH,'JOBGS')
          IAD = 0
          Call IDAFile(LUIPH,2,ITOC15,30,IAD)
          IAD=ITOC15(2)
          Call DDAFile(LUIPH,1,Work(LCMO2),nCMO,IAD)
          Call DACLOS(LUIPH)
        End if !DoGSOR


        CALL GETMEM('GTDMTRA1','ALLO','REAL',LTRA1,NTRA)
        CALL GETMEM('GTDMTRA2','ALLO','REAL',LTRA2,NTRA)
        CALL FINDT(WORK(LCMO1),WORK(LCMO2),WORK(LTRA1),WORK(LTRA2))
        TRORB = .true.
      else
        TRORB = .false.
      end if

      !> check whether we do RASSI with an effective multi-state PT2 Hamiltonian
      !> whose eigenvectors are stored in Heff_evc
      !> i.e., we do not use mixed CI coefficients / MPS wave functions but rather mix the TDMs

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
      LFMO =LNILPT
      LTUVX=LNILPT
      IF (IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
       CALL GETMEM('FMO','ALLO','REAL',LFMO,NTDM1)
       CALL GETMEM('TUVX  ','ALLO','REAL',LTUVX,NTDM2)
       CALL FZERO(WORK(LTUVX),NTDM2)
CTEST       write(*,*)'GTDMCTL calling TRINT.'
       CALL TRINT(WORK(LCMO1),WORK(LCMO2),ECORE,
     &              NTDM1,WORK(LFMO),NTDM2,WORK(LTUVX))
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
          CALL DCOPY_(NDET1,0.0D0,0,WORK(LDET1),1)
          CALL PREPSD(WFTP1,TRORB,ISGSTR1,ICISTR1,IXSTR1,LSYM1,
     &                WORK(LTRA1),IWORK(LCNFTAB1),IWORK(LSPNTAB1),
     &                IWORK(LSSTAB),IWORK(LFSBTAB1),NCONF1,WORK(LCI1),
     &                WORK(LDET1))

C Write out the determinant expansion to disk.
          IDDET1(ISTATE)=IDWSCR
          CALL  DDAFILE(LUSCR,1,WORK(LDET1),NDET1,IDWSCR)
        else
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 iWork(lLROOT+ISTATE-1),
     &                 LSYM1,
     &                 MPLET1,
     &                 MSPROJ1,
     &                 NACTE1,
     &                 WORK(LTRA1),
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
        CALL DCOPY_(NCONF2,0.0D0,0,WORK(LTheta1),1)
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
          CALL DCOPY_(NDET2,0.0D0,0,WORK(LDET2),1)
          CALL PREPSD(WFTP2,TRORB,ISGSTR2,ICISTR2,IXSTR2,LSYM2,
     &                WORK(LTRA2),IWORK(LCNFTAB2),IWORK(LSPNTAB2),
     &                IWORK(LSSTAB),IWORK(LFSBTAB2),NCONF2,WORK(LCI2),
     &                WORK(LDET2))

        else
#ifdef _DMRG_
          call prepMPS(
     &                 TRORB,
     &                 iWork(llroot+JSTATE-1),
     &                 LSYM2,
     &                 MPLET2,
     &                 MSPROJ2,
     &                 NACTE2,
     &                 WORK(LTRA2),
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
     &            DYSAMP,WORK(LDYSCOF))

C Write Dyson orbital coefficients in AO basis to disk.
        IF (DYSAMP.GT.1.0D-6) THEN
C In full biorthonormal basis:
         CALL MKDYSAB(WORK(LDYSCOF),WORK(LDYSAB))
!         WRITE(*,*)"NDYSAB=",NDYSAB
!         WRITE(*,*)"DYSAB="
!         DO NDUM=1,NDYSAB
!          WRITE(*,*)WORK(LDYSAB+NDUM-1)
!         END DO
C In AO basis:
         CALL MKDYSZZ(WORK(LCMO1),WORK(LDYSAB),
     &               WORK(LDYSZZ))
!        WRITE(*,*)"NDYSZZ=",NDYSZZ
!        WRITE(*,*)"DYSZZ="
!        DO NDUM=1,NDYSZZ
!         WRITE(*,*)WORK(LDYSZZ+NDUM-1)
!        END DO
        IF (DYSO) THEN
         DO NDUM=1,NDYSZZ
          SFDYS(NDUM,JSTATE,ISTATE)=WORK(LDYSZZ+NDUM-1)
          SFDYS(NDUM,ISTATE,JSTATE)=WORK(LDYSZZ+NDUM-1)
         END DO
        END IF
        DO NDUM=1,NDYSZZ
         WORK(LDYSZZ+NDUM-1)=0.0D0
        END DO
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
     &            WORK(LTRAD),WORK(LTRASD),WORK(LWERD),ISTATE,
     &            JSTATE,lLROOT,job1,job2,ist,jst)

        IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
C Compute 1-electron contribution to Hamiltonian matrix element:
        HONE=DDOT_(NTRAD,WORK(LTRAD),1,WORK(LFMO),1)
        END IF


C             Write density 1-matrices in AO basis to disk.
            IF(NATO.OR.(NPROP.GT.0))THEN

              !> regular-TDM
              CALL MKTDAB(SIJ,WORK(LTRAD),WORK(LTDMAB))
              !> transform to AO basis
              CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LTDMAB),
     &                    WORK(LTDMZZ))

              !> spin-TDM
              CALL MKTDAB(SIJ,WORK(LTRASD),WORK(LTSDMAB))
              !> transform to AO basis
              CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LTSDMAB),
     &                    WORK(LTSDMZZ) )

              !> WE-reduced TDM''s of triplet type:
              CALL MKTDAB(0.0D0,WORK(LWERD),WORK(LWDMAB))
              !> transform to AO basis
              CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LWDMAB),
     &                    WORK(LWDMZZ))

              if(.not.mstate_dens)then

                IF((SONATNSTATE.GT.0).OR.NATO) THEN
*C Transition density matrices, TDMZZ, in AO basis.
*C WDMZZ similar, but WE-reduced 'triplet' densities.
                  IDISK=iWork(lIDTDM+(ISTATE-1)*NSTATE+JSTATE-1)
                  CALL dens2file(WORK(LTDMZZ),WORK(LTSDMZZ),
     &                           WORK(LWDMZZ),NTDMZZ,LUTDM,IDISK)
                END IF
                !> calculate property matrix elements
                CALL PROPER(PROP,ISTATE,JSTATE,
     &                      WORK(LTDMZZ),WORK(LWDMZZ))
              else

                !> scale rdm elements with eigenvector coefficients of Heff of a multi-state (PT2) Hamiltonian
                !> accumulate data first and run PROPER and other utility routines later
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
     &                          work(ltdmzz),1,
     &                          mstate_1pdens(i,j)%rtdm,1
     &                         )
                    !> spin-TDM
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          work(ltsdmzz),1,
     &                          mstate_1pdens(i,j)%stdm,1
     &                         )
                    !> WE-reduced TDM''s of triplet type:
                    call daxpy_(ntdmzz,
     &                          fac1*fac2,
     &                          work(lwdmzz),1,
     &                          mstate_1pdens(i,j)%wtdm,1
     &                         )
                    !> overlap
                    call daxpy_(1,
     &                          fac1*fac2,
     &                          sij,0,
     &                          mstate_1pdens(i,j)%overlap,0
     &                         )
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
                  !> Leon: TODO: Add possibility to calculate overlap of rotated MPS without using checkpoint names
                  call dmrg_interface_ctl(
     &                               task   = 'overlap ',
     &                               energy = sij,
     &                               state  = iWork(lLROOT+istate-1),
     &                               stateL = iWork(lLROOT+jstate-1)
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
     &                  WORK(LDET1),WORK(LDET2),NTDM2,WORK(LTDM2),
     &                  ISTATE,JSTATE,lLROOT,job1,job2,ist,jst)

            !> Compute 2-electron contribution to Hamiltonian matrix element:
            IF(IFTWO.AND.(MPLET1.EQ.MPLET2))
     &      HTWO=DDOT_(NTDM2,WORK(LTDM2),1,WORK(LTUVX),1)

          END IF ! IF22

          !> PAM 2011 Nov 3, writing transition matrices if requested
          !> The following *long* section should later be moved to its own
          !> subroutine...
          IF ((IFTRD1.or.IFTRD2).and..not.mstate_dens) THEN
            LU=50
            LU=IsFreeUnit(LU)
            WRITE(NUM1,'(I3.3)') ISTATE
            WRITE(NUM2,'(I3.3)') JSTATE
            FNM='TRD2_'//NUM1//'_'//NUM2
            CALL Molcas_Open(LU,FNM)
            WRITE(LU,*)'#Transition density file from RASSI.'
            WRITE(LU,*)'#  States:'
            WRITE(LU,*) ISTATE, JSTATE
            WRITE(LU,*)'#  Nr of irreps:'
            WRITE(LU,*) NSYM
            WRITE(LU,*)'#  Basis functions:'
            WRITE(LU,'(8I5)') (NBASF(ISYM),ISYM=1,NSYM)
            WRITE(LU,*)'#  Frozen orbitals:'
            WRITE(LU,'(8I5)') (NFRO(ISYM),ISYM=1,NSYM)
            WRITE(LU,*)'#  Inactive orbitals:'
            WRITE(LU,'(8I5)') (NISH(ISYM),ISYM=1,NSYM)
            WRITE(LU,*)'#  Active orbitals:'
            WRITE(LU,'(8I5)') (NASH(ISYM),ISYM=1,NSYM)
            WRITE(LU,*)'#  State ',ISTATE,'    CMO coefficients:'
            LPOS=LCMO1
            DO ISYM=1,NSYM
              NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
              NB=NBASF(ISYM)
              DO IO=1,NO
                WRITE(LU,*)'#  Symm ',ISYM,'   Orbital ',IO
                WRITE(LU,'(5D19.12)')(WORK(LPOS+NB*(IO-1)+i),i=0,NB-1)
              END DO
              LPOS=LPOS+NB**2
            END DO
            WRITE(LU,*)'#  State ',JSTATE,'    CMO coefficients:'
            LPOS=LCMO2
            DO ISYM=1,NSYM
              NO=NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
              NB=NBASF(ISYM)
              DO IO=1,NO
                WRITE(LU,*)'#  Symm ',ISYM,'   Orbital ',IO
                WRITE(LU,'(5D19.12)')(WORK(LPOS+NB*(IO-1)+i),i=0,NB-1)
              END DO
              LPOS=LPOS+NB**2
            END DO
            WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Overlap:'
            WRITE(LU,'(5D19.12)') SIJ
            WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Active TRD1:'
            LSYM12=MUL(LSYM1,LSYM2)
            LPOS=LTDMAB
            DO ISYM1=1,NSYM
              NO1=NOSH(ISYM1)
              ISYM2=MUL(ISYM1,LSYM12)
              NO2=NOSH(ISYM2)
              IF (NO1*NO2 .gt. 0) THEN
                NA1=NASH(ISYM1)
                NA2=NASH(ISYM2)
                IF (NA1*NA2 .gt. 0) THEN
                  NI1=NISH(ISYM1)
                  NI2=NISH(ISYM2)
                  WRITE(LU,*)'#  Symmetries ',ISYM1,ISYM2
                  WRITE(LU,'(5D19.12)')((WORK(LPOS-1+II+NO1*(JJ-1)),
     &                                  JJ=NI2+1,NO2),II=NI1+1,NO1)
                END IF
                LPOS=LPOS+NO1*NO2
              END IF
            END DO

            IF (IFTRD2.AND.IF22) THEN
              WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Active TRD2:'
              DO ISYT=1,NSYM
                DO ISYU=1,NSYM
                  ISYTU=ISYT+NSYM*(ISYU-1)
                  DO ISYV=1,ISYT
                    LIMX=ISYV
                    IF(ISYV.EQ.ISYT) LIMX=ISYU
                    DO ISYX=1,LIMX
                      ISYVX=ISYV+NSYM*(ISYX-1)
                      !> Write out one symmetry block (4 indices!) of two-electron
                      !> transition density matrix elements.
                      !> Write a full 'rectangular' array, even if it could be made
                      !> smaller by permutation symmetry.
                      WRITE(LU,*)'#  Orbital symm:',ISYT,ISYU,ISYV,ISYX
                      IWBUF=0
                      DO IT=1,NASH(ISYT)
                        ITABS=NAES(ISYT)+IT
                        DO IU=1,NASH(ISYU)
                          IUABS=NAES(ISYU)+IU
                          ITU=ITABS+NASHT*(IUABS-1)
                          DO IV=1,NASH(ISYV)
                            IVABS=NAES(ISYV)+IV
                            DO IX=1,NASH(ISYX)
                              IXABS=NAES(ISYX)+IX
                              IVX=IVABS+NASHT*(IXABS-1)
                              IF(ITU.GE.IVX) THEN
                                ITUVX=(ITU*(ITU-1))/2+IVX
                              ELSE
                                ITUVX=(IVX*(IVX-1))/2+ITU
                              END IF
                              IWBUF=IWBUF+1
                              WBUF(IWBUF)=WORK(LTDM2-1+ITUVX)
                              IF(IWBUF.EQ.5) THEN
                                WRITE(LU,'(5D19.12)')(WBUF(I),I=1,IWBUF)
                                IWBUF=0
                              END IF
                            END DO
                          END DO
                        END DO
                      END DO
                      IF(IWBUF.GT.0) THEN
                        WRITE(LU,'(5D19.12)')(WBUF(I),I=1,IWBUF)
                        IWBUF=0
                      END IF
* End of writing a symmetry block.
                    END DO
                  END DO
                END DO
              END DO
            END IF
            CLOSE (LU)
          END IF ! TRD1/2

#ifdef _HDF5_
          if(.not.mstate_dens)then
            IF(IF11.AND.(LSYM1.EQ.LSYM2).and.
     &        ((SONATNSTATE.GT.0).OR.NATO))THEN
              call mh5_put_dset_array_real(wfn_sfs_tdm,
     $        WORK(LTDMZZ),[NTDMZZ,1,1], [0,ISTATE-1,JSTATE-1])
              call mh5_put_dset_array_real(wfn_sfs_tsdm,
     $        WORK(LTSDMZZ),[NTDMZZ,1,1], [0,ISTATE-1,JSTATE-1])
            END IF
            IF(SONATNSTATE.GT.0.OR.NATO)THEN
              call mh5_put_dset_array_real(wfn_sfs_wetdm,
     $        WORK(LWDMZZ),[NTDMZZ,1,1], [0,ISTATE-1,JSTATE-1])
            END IF
          end if
#endif
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
        Open(unit=87,file='CI_THETA', action='write',iostat=ios)
        do i=1,NCONF2
          write(87,*) Work(LTheta1-1+i)
        end do
        Close(87)

       !Now we need to build the other states.
        DO JST=2,NSTAT(JOB2)
          JSTATE=ISTAT(JOB2)-1+JST
          CALL READCI(JSTATE,ISGSTR2,ICISTR2,NCONF2,WORK(LCI2))
          Call DCOPY_(NCONF2,Work(LCI2),1,Work(LCI2_o),1)
          CALL DCOPY_(NDET2,0.0D0,0,WORK(LDET2),1)
          CALL PREPSD(WFTP2,TRORB,ISGSTR2,ICISTR2,IXSTR2,LSYM2,
     &              WORK(LTRA2),IWORK(LCNFTAB2),IWORK(LSPNTAB2),
     &          IWORK(LSSTAB),IWORK(LFSBTAB2),NCONF2,WORK(LCI2),
     &          WORK(LDET2))

          CALL GETMEM('ThetaN','ALLO','REAL',LThetaN,NCONF2)
          CALL DCOPY_(NCONF2,Work(LCI2_o),1,WORK(LThetaN),1)
          Norm_Fac = DDOT_(NCONF2,Work(LTheta1),1,Work(LCI2_o),1)
          Call DAXPY_(NCONF2,-Norm_Fac,Work(LTHETA1),1,Work(LThetaN),1)

          Open(unit=87,file='CI_THETA', action='read',iostat=ios)
          if(JST-1.ge.2) then
            do i=1,NCONF2
              Read(87,*) dummy
            end do
          end if
          CALL GETMEM('ThetaM','ALLO','REAL',LThetaM,NCONF2)
          DO IST=2,JST-1
            CALL DCOPY_(NCONF2,0.0D0,0,WORK(LThetaM),1)
            !Read in previous theta vectors
            do i=1,NCONF2
              Read(87,*) Work(LThetaM-1+i)
            end do
            Dot_prod = DDOT_(NCONF2,Work(LThetaM),1,Work(LCI2_o),1)
           Call DAXPY_(NCONF2,-Dot_prod,Work(LThetaM),1,Work(LThetaN),1)


          END DO
          Close(87)
          !Normalize
          dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LThetaN),1)
          Norm_Fac = 1d0/sqrt(dot_prod)
          Call DSCAL_(NCONF2,Norm_Fac,Work(LThetaN),1)

        !dot_prod = DDOT_(NCONF2,Work(LTheta1),1,Work(LTheta1),1)
        !dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LTheta1),1)
        !dot_prod = DDOT_(NCONF2,Work(LThetaN),1,Work(LThetaN),1)

          !Write to file
          Open(unit=87,file='CI_THETA', position='append',iostat=ios,
     &    action='write')
          do i=1,nConf2
            write(87,*) Work(LThetaN-1+i)
          end do
          close(87)
          !Deallocate
        END DO
!Copy to new IPH file
        Open(unit=87,file='CI_THETA',iostat=ios,
     &    action='read')
        CALL DANAME(LUIPHn,'JOBGS')
        IAD = 0
        Call IDAFILE(LUIPHn,2,ITOC15,30,IAD)
        IAD=ITOC15(4)
        do i=1,ISTAT(JOB1)-1
         CALL DCOPY_(NCONF2,0.0D0,0,WORK(LThetaM),1)
         do j=1,nCONF2
           read(87,*) Work(LThetaM-1+i)
         end do
         Call DDafile(LUIPHn,1,Work(LThetaM),nCONF2,IAD)
       end do

       IAD = ITOC15(4)
       CALL DCOPY_(NCONF2,0.0D0,0,WORK(LThetaM),1)
       Call DDAFILE(LUIPHn,2,Work(LThetaM),nCONF2,IAD)
       Call DDAFILE(LUIPHn,2,Work(LThetaM),nCONF2,IAD)

       Close(87)
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

      !> create actual property data and put everything to file (if requested) in case of using eigenvectors of a multi-state (PT2) Hamiltonian
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
     &                           iWork(lIDTDM+(ISTATE-1)*NSTATE+JSTATE),
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
        CALL GETMEM('GTDMTRA1','FREE','REAL',LTRA1,NTRA)
        CALL GETMEM('GTDMTRA2','FREE','REAL',LTRA2,NTRA)
      END IF
      CALL GETMEM('GTDMDET1','FREE','REAL',LDET1,NDET1)
      CALL GETMEM('GTDMDET2','FREE','REAL',LDET2,NDET2)
      CALL GETMEM('GTDMCI2','FREE','REAL',LCI2,NCONF2)
      IF (.NOT.NONA) CALL GETMEM('GTDMCI1','FREE','REAL',LCI1,NCONF1)
      if(.not.doDMRG)then
        CALL KILLSCTAB(LSPNTAB1)
        CALL KILLSCTAB(LSPNTAB2)
      end if
! +++ J. Norell 13/7 - 2018
      IF ((IF10.or.IF01).and.DYSO) THEN
        CALL GETMEM('DYSCOF','Free','Real',LDYSCOF,NASORB)
        CALL GETMEM('DYSAB','Free','Real',LDYSAB,NDYSAB)
        CALL GETMEM('DYSZZ','Free','Real',LDYSZZ,NDYSZZ)
      END IF
! +++
      IF (IF11) THEN
        CALL GETMEM('SPD1','Free','Real',LSPD1,NSPD1)
        CALL GETMEM('TRAD','Free','Real',LTRAD,NTRAD)
        CALL GETMEM('TRASD','Free','Real',LTRASD,NTRASD)
        CALL GETMEM('WERD','Free','Real',LWERD,NWERD)
        IF(NATO.OR.NPROP.GT.0) THEN
          CALL GETMEM('TDMAB','Free','Real',LTDMAB,NTDMAB)
          CALL GETMEM('TDMZZ','Free','Real',LTDMZZ,NTDMZZ)
          CALL GETMEM('TSDMAB','Free','Real',LTSDMAB,NTDMAB)
          CALL GETMEM('TSDMZZ','Free','Real',LTSDMZZ,NTDMZZ)
          CALL GETMEM('WDMAB','Free','Real',LWDMAB,NTDMAB)
          CALL GETMEM('WDMZZ','Free','Real',LWDMZZ,NTDMZZ)
        END IF
      END IF
      IF (IF22) THEN
        CALL GETMEM('SPD2','Free','Real',LSPD2,NSPD2)
        CALL GETMEM('TDM2','Free','Real',LTDM2,NTDM2)
      END IF

      IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
        CALL GETMEM('FMO','Free','REAL',LFMO,NTDM1)
        CALL GETMEM('TUVX  ','Free','REAL',LTUVX,NTDM2)
      END IF

      CALL GETMEM('GTDMCMO1','FREE','REAL',LCMO1,NCMO)
      CALL GETMEM('GTDMCMO2','FREE','REAL',LCMO2,NCMO)
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

      CALL QEXIT(ROUTINE)
      RETURN
      END
