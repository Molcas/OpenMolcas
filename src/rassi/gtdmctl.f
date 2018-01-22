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
      SUBROUTINE GTDMCTL(PROP,JOB1,JOB2,OVLP,HAM,IDDET1)

      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
#endif

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
      DIMENSION ISGSTR1(NSGSIZE), ISGSTR2(NSGSIZE)
      DIMENSION ICISTR1(NCISIZE), ICISTR2(NCISIZE)
      DIMENSION IXSTR1(NXSIZE), IXSTR2(NXSIZE)
      DIMENSION PROP(NSTATE,NSTATE,NPROP)
      DIMENSION NGASORB(100),NGASLIM(2,10)
      DIMENSION NASHES(8)
      DIMENSION OVLP(NSTATE,NSTATE),HAM(NSTATE,NSTATE)
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
#include "SysDef.fh"

      CALL QENTER(ROUTINE)

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

C Logical variables, controlling which GTDM''s to compute:

C First, see which type of GTDM''s that could be non-zero:
C Overlap:
      IF00 = NACTE1.EQ.NACTE2 .AND.  MPLET1.EQ.MPLET2
      IF00 = IF00 .AND.  LSYM1.EQ. LSYM2
C Dyson amplitudes:
      IF10 = (NACTE1-NACTE2).EQ. 1
      IF10 = IF10 .AND. ( ABS(MPLET1-MPLET2).EQ.1 )
      IF01 = (NACTE1-NACTE2).EQ.-1
      IF01 = IF01 .AND. ( ABS(MPLET1-MPLET2).EQ.1 )
C Pair amplitudes:
      IF20 = (NACTE1-NACTE2).EQ. 2
      IF02 = (NACTE1-NACTE2).EQ.-2
C Density 1-matrices:
      IF11 = ( NACTE1.EQ.NACTE2 .AND. NACTE1.GE.1 )
      IF11 = IF11 .AND. ( ABS(MPLET1-MPLET2).LE.2 )
C 2h1p and 1h2p amplitudes:
      IF21 = IF10 .AND. NACTE2.GE.1
      IF21 = IF21 .AND. ( ABS(MPLET1-MPLET2).LE.3 )
      IF12 = IF01 .AND. NACTE1.GE.1
      IF12 = IF12 .AND. ( ABS(MPLET1-MPLET2).LE.3 )
C Two-electron transition density and transition spin density
      IF22 = ( NACTE1.EQ.NACTE2 .AND. NACTE1.GE.2 )
      IF22 = IF22 .AND. ( ABS(MPLET1-MPLET2).LE.4 )

C Then, check if they are needed at all:
C It may be that the Hamiltonian matrix should be used in
C diagonalization (IFHAM is .TRUE.), but it does not have
C to be computed (because IFHEXT or IFHEFF or IFEJOB are true).
      IFTWO = IFHAM .AND.
     &            .NOT.(IFHEXT.OR.IFHEFF.OR.IFEJOB)

C For the moment, we have no use for the two-electron density
C except when used for the scalar two-body Hamiltonian matrix:
      IF22 = IF22 .AND. IFTWO .AND. (MPLET1.EQ.MPLET2) .AND.
     &        (LSYM1.EQ.LSYM2)


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
        CALL GETMEM('GTDMTRA1','ALLO','REAL',LTRA1,NTRA)
        CALL GETMEM('GTDMTRA2','ALLO','REAL',LTRA2,NTRA)
        CALL FINDT(WORK(LCMO1),WORK(LCMO2),WORK(LTRA1),WORK(LTRA2))
      END IF
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

#ifdef _DMRG_
        IF(.not.doDMRG)then
#endif
        CALL SGINIT(NSYM,NACTE1,MPLET1,NRSPRT,NRAS,NRASEL,ISGSTR1)
        IF(IPGLOB.GT.DEBUG) THEN
          WRITE(6,*)'Split-graph structure for JOB1=',JOB1
          CALL SGPRINT(ISGSTR1)
        END IF
        CALL SGSVAL(ISGSTR1,NSYM,NASHT,LISM,NVERT,LDRT,
     &              LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
        CALL CXINIT(ISGSTR1,ICISTR1,IXSTR1)
        CALL CXSVAL(ICISTR1,IXSTR1,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &                 LNOCSF,LIOCSF,NWALK,LICASE,
     &                 MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &           LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C CI sizes, as function of symmetry, are now known.
        NCONF1=IWORK(LNCSF-1+LSYM1)
#ifdef _DMRG_
        else
          NCONF1=1
        end if
#endif
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

      IFORM=1
      MINOP=0
      LREST1=NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM)
      IF(IPGLOB.GE.DEBUG) CALL PRGASTAB(LREST1)

C At present, we will only annihilate, at most 2 electrons will
C be removed. This limits the possible MAXOP:
      MAXOP=MIN(MAXOP+1,NACTE1,NASHT)
      LCNFTAB1=NEWCNFTAB(NACTE1,NASHT,MINOP,MAXOP,LSYM1,NGAS,
     &                           NGASORB,NGASLIM,IFORM)
      IF(IPGLOB.GE.DEBUG) CALL PRCNFTAB(LCNFTAB1,100)

      LFSBTAB1=NEWFSBTAB(NACTE1,MSPROJ1,LSYM1,LREST1,LSSTAB)
      IF(IPGLOB.GE.DEBUG) CALL PRFSBTAB(IWORK(LFSBTAB1))
      NDET1=IWORK(LFSBTAB1+4)
#ifdef _DMRG_
      if(doDMRG) NDET1 = 1 ! minimum to avoid runtime error
#endif
      LSPNTAB1=NEWSCTAB(MINOP,MAXOP,MPLET1,MSPROJ1)
      IF (IPGLOB.GT.DEBUG) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
        CALL PRSCTAB(LSPNTAB1)
      END IF
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

#ifdef _DMRG_
        IF(.not.doDMRG)then
#endif
        CALL SGINIT(NSYM,NACTE2,MPLET2,NRSPRT,NRAS,NRASEL,ISGSTR2)
        IF(IPGLOB.GT.DEBUG) THEN
          WRITE(6,*)'Split-graph structure for JOB2=',JOB2
          CALL SGPRINT(ISGSTR2)
        END IF
        CALL SGSVAL(ISGSTR2,NSYM,NASHT,LISM,NVERT,LDRT,
     &              LDOWN,LUP,MIDLEV,MVSTA,MVEND,LMAW,LLTV)
        CALL CXINIT(ISGSTR2,ICISTR2,IXSTR2)
        CALL CXSVAL(ICISTR2,IXSTR2,NMIDV,NIPWLK,LNOW,LIOW,LNCSF,
     &                 LNOCSF,LIOCSF,NWALK,LICASE,
     &                 MXEO,LNOCP,LIOCP,NICOUP,LICOUP,NVTAB,
     &           LVTAB,LMVL,LMVR,NT1MX,NT2MX,NT3MX,NT4MX,NT5MX)
C CI sizes, as function of symmetry, are now known.
        NCONF2=IWORK(LNCSF-1+LSYM2)
#ifdef _DMRG_
        else
          NCONF2=1
        end if
#endif
      ELSE
C Presently, the only other cases are HISPIN, CLOSED or EMPTY.
* Note: the HISPIN case may be buggy and is not used presently.
        NCONF2=1
      END IF
      CALL GETMEM('GTDMCI2','ALLO','REAL',LCI2,NCONF2)

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

      LREST2=NEWGASTAB(NSYM,NGAS,NGASORB,NGASLIM)
      IF(IPGLOB.GE.DEBUG) CALL PRGASTAB(LREST2)

      IFORM=1
      MINOP=0
C At present, we will only annihilate. This limits the possible MAXOP:
      MAXOP=MIN(MAXOP+1,NACTE2,NASHT)
      LCNFTAB2=NEWCNFTAB(NACTE2,NASHT,MINOP,MAXOP,LSYM2,NGAS,
     &                           NGASORB,NGASLIM,IFORM)
      IF(IPGLOB.GE.DEBUG) CALL PRCNFTAB(LCNFTAB2,100)

      LFSBTAB2=NEWFSBTAB(NACTE2,MSPROJ2,LSYM2,LREST2,LSSTAB)
      IF(IPGLOB.GE.DEBUG) CALL PRFSBTAB(IWORK(LFSBTAB2))
      NDET2=IWORK(LFSBTAB2+4)
#ifdef _DMRG_
      if(doDMRG) NDET2 = 1 ! minimum to avoid runtime error
#endif
      LSPNTAB2=NEWSCTAB(MINOP,MAXOP,MPLET2,MSPROJ2)
      IF (IPGLOB.GT.DEBUG) THEN
*PAM2009: Put in impossible call to PRSCTAB, just so code analyzers
* do not get their knickers into a twist.
        CALL PRSCTAB(LSPNTAB2)
      END IF
C-------------------------------------------------------------
      CALL GETMEM('GTDMDET1','ALLO','REAL',LDET1,NDET1)
      CALL GETMEM('GTDMDET2','ALLO','REAL',LDET2,NDET2)

#ifdef _DMRG_
      if(.not.doDMRG)then
#endif
C Loop over the states of JOBIPH nr JOB1
C Disk address for writing to scratch file is IDWSCR.
      IDWSCR=0
      DO IST=1,NSTAT(JOB1)
        ISTATE=ISTAT(JOB1)-1+IST

C Read ISTATE wave function
        IF(WFTP1.EQ.'GENERAL ') THEN
          CALL READCI(ISTATE,ISGSTR1,ICISTR1,NCONF1,WORK(LCI1))
        ELSE
          WORK(LCI1)=1.0D0
        END IF
        TRORB=(JOB1.NE.JOB2)
        CALL DCOPY_(NDET1,0.0D0,0,WORK(LDET1),1)
        CALL PREPSD(WFTP1,TRORB,ISGSTR1,ICISTR1,IXSTR1,LSYM1,
     &              WORK(LTRA1),IWORK(LCNFTAB1),IWORK(LSPNTAB1),
     &              IWORK(LSSTAB),IWORK(LFSBTAB1),NCONF1,WORK(LCI1),
     &              WORK(LDET1))

C Write out the determinant expansion to disk.
        IDDET1(ISTATE)=IDWSCR
        CALL  DDAFILE(LUSCR,1,WORK(LDET1),NDET1,IDWSCR)

      END DO

#ifdef _DMRG_
      end if
#endif


C-------------------------------------------------------------

C Loop over the states of JOBIPH nr JOB2
      DO JST=1,NSTAT(JOB2)
        JSTATE=ISTAT(JOB2)-1+JST

C Read JSTATE wave function
#ifdef _DMRG_
        if(doDMRG)then
!         write(6,*) 'skip determinant-based wave functions for JOB2'
        else
#endif
        IF(WFTP2.EQ.'GENERAL ') THEN
          CALL READCI(JSTATE,ISGSTR2,ICISTR2,NCONF2,WORK(LCI2))
        ELSE
          WORK(LCI2)=1.0D0
        END IF
        TRORB=(JOB1.NE.JOB2)
        CALL DCOPY_(NDET2,0.0D0,0,WORK(LDET2),1)
        CALL PREPSD(WFTP2,TRORB,ISGSTR2,ICISTR2,IXSTR2,LSYM2,
     &              WORK(LTRA2),IWORK(LCNFTAB2),IWORK(LSPNTAB2),
     &          IWORK(LSSTAB),IWORK(LFSBTAB2),NCONF2,WORK(LCI2),
     &          WORK(LDET2))

#ifdef _DMRG_
        end if
#endif

C Loop over the states of JOBIPH nr JOB1
      DO IST=1,NSTAT(JOB1)
        ISTATE=ISTAT(JOB1)-1+IST
        IF(ISTATE.LT.JSTATE) GOTO 100
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

C Calculate whatever type of GTDM that was requested, unless
C it is known to be zero.
      HZERO=0.0D0
      HONE =0.0D0
      HTWO =0.0D0

      SIJ=0.0D0
C General 1-particle transition density matrix:
      IF (IF11) THEN
        CALL MKTDM1(LSYM1,MPLET1,MSPROJ1,IWORK(LFSBTAB1),
     &            LSYM2,MPLET2,MSPROJ2,IWORK(LFSBTAB2),IWORK(LSSTAB),
     &            IWORK(LOMAP),WORK(LDET1),WORK(LDET2),SIJ,NASHT,
     &            WORK(LTRAD),WORK(LTRASD),WORK(LWERD),ISTATE,
     &            JSTATE)

        IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
C Compute 1-electron contribution to Hamiltonian matrix element:
        HONE=DDOT_(NTRAD,WORK(LTRAD),1,WORK(LFMO),1)
        END IF

C Write density 1-matrices in AO basis to disk.
        IF(NATO.OR.(NPROP.GT.0)) THEN
C First, conventional 'spin-scalar' TDM''s
          CALL MKTDAB(SIJ,WORK(LTRAD),WORK(LTDMAB))
          CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LTDMAB),
     *                        WORK(LTDMZZ) )

#ifdef _DMRG_
          !> debug print
          if(doDMRG)then
          write(6,*) ' blubb: 1-TDM in AO basis'
          call pretty_print_util(WORK(LTDMZZ),1,1,1,ntdmzz,
     &                       1,ntdmzz,1,6)
          end if
#endif

          CALL MKTDAB(SIJ,WORK(LTRASD),WORK(LTSDMAB))
          CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LTSDMAB),
     *                        WORK(LTSDMZZ) )

#ifdef _DMRG_
          !> debug print
          if(doDMRG)then
            write(6,*) ' blubb: 1-SDM in AO basis'
            call pretty_print_util(WORK(LTSDMZZ),1,1,1,ntdmzz,
     &                             1,ntdmzz,1,6)
          end if
#endif

          ISY12=MUL(LSYM1,LSYM2)
C Then, WE-reduced TDM''s of triplet type:
          CALL MKTDAB(0.0D0,WORK(LWERD),WORK(LWDMAB))
          CALL MKTDZZ(WORK(LCMO1),WORK(LCMO2),WORK(LWDMAB),WORK(LWDMZZ))

#ifdef _DMRG_
          !> debug print
          if(doDMRG)then
          write(6,*) ' blubb: WE-reduced-TDM in AO basis'
          call pretty_print_util(WORK(LWDMZZ),1,1,1,ntdmzz,
     &                       1,ntdmzz,1,6)
          end if
#endif


! jochen 02/15: sonatorb needs these files
!        we'll make the I/O conditional upon the keyword
          IF(SONATNSTATE.GT.0) THEN
            IDISK=iWork(lIDTDM+(ISTATE-1)*NSTATE+JSTATE)
            CALL DDAFILE(LUTDM,1,WORK(LTDMZZ),NTDMZZ,IDISK)
            CALL DDAFILE(LUTDM,1,WORK(LTSDMZZ),NTDMZZ,IDISK)
            CALL DDAFILE(LUTDM,1,WORK(LWDMZZ),NTDMZZ,IDISK)
          END IF
*C Transition density matrices, TDMZZ, in AO basis.
*C WDMZZ similar, but WE-reduced 'triplet' densities.
          CALL PROPER (PROP,ISTATE,JSTATE,WORK(LTDMZZ),WORK(LWDMZZ))
        END IF
      ELSE
C Overlap:
        IF (IF00) THEN
#ifdef _DMRG_
          if(.not.doDMRG)then
#endif
          SIJ=OVERLAP_RASSI(IWORK(LFSBTAB1),IWORK(LFSBTAB2),WORK(LDET1),
     &                WORK(LDET2))
#ifdef _DMRG_
          else
            call dmrg_interface_ctl(
     &                              task   = 'overlap ',
     &                              energy = sij,
     &                              state  = istate,
     &                              stateL = jstate
     &                             )

          write(6,*)"SIJ in RASSI, dmrg  gtdmctl",SIJ

          end if
#endif
        END IF
      END IF

      OVLP(ISTATE,JSTATE)=SIJ
      OVLP(JSTATE,ISTATE)=SIJ

C General 2-particle transition density matrix:
      IF (IF22) THEN
#ifdef _DMRG_
        if(.not.doDMRG)then
#endif
          CALL MKTDM2(LSYM1,MSPROJ1,IWORK(LFSBTAB1),LSYM2,MSPROJ2,
     &                IWORK(LFSBTAB2),IWORK(LSSTAB),IWORK(LOMAP),
     &                WORK(LDET1),WORK(LDET2),NTDM2,WORK(LTDM2))
#ifdef _DMRG_
        else
          write(6,*) ' DMRG-RASSI: skipped 2e-TDMs for now...      '
          write(6,*) ' contact stknecht(at)ethz.ch for help!       '
        end if
#endif

C F. Plasser: I moved this block up to allow for proper printing
C    of the HDF5 files
        IF(IFTWO.AND.(MPLET1.EQ.MPLET2)) THEN
C Compute 2-electron contribution to Hamiltonian matrix element:
          HTWO=DDOT_(NTDM2,WORK(LTDM2),1,WORK(LTUVX),1)
        END IF

C End of IF22 case
      END IF
* PAM 2011 Nov 3, writing transition matrices if requested
* The following *long* section should later be moved to its own
* subroutine...
        IF (IFTRD1.or.IFTRD2) THEN
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
            WRITE(LU,'(5D18.12)')(WORK(LPOS+NB*(IO-1)+i),i=0,NB-1)
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
            WRITE(LU,'(5D18.12)')(WORK(LPOS+NB*(IO-1)+i),i=0,NB-1)
           END DO
           LPOS=LPOS+NB**2
          END DO
          WRITE(LU,*)'#  States ',ISTATE,JSTATE,' Overlap:'
          WRITE(LU,'(5D18.12)') SIJ
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
              WRITE(LU,'(5D18.12)')((WORK(LPOS-1+II+NO1*(JJ-1)),
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
* Write out one symmetry block (4 indices!) of two-electron
* transition density matrix elements.
* Write a full 'rectangular' array, even if it could be made
* smaller by permutation symmetry.
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
                 WRITE(LU,'(5D18.12)')(WBUF(I),I=1,IWBUF)
                 IWBUF=0
                END IF
               END DO
              END DO
             END DO
            END DO
            IF (IWBUF.GT.0) THEN
             WRITE(LU,'(5D18.12)')(WBUF(I),I=1,IWBUF)
             IWBUF=0
            END IF
* End of writing a symmetry block.
               END DO
              END DO
             END DO
            END DO
          END IF
          CLOSE (LU)
* PAM 2011 Nov 3: End of added section.

#ifdef _HDF5_
c         write(6,*) 'NTDMZZ=',NTDMZZ
c         write(6,*) 'ISTATE,JSTATE=',ISTATE,JSTATE
c         write(6,*)(WORK(LTDMZZ+i-1),i=1,10)
          IF(IF11.AND.(LSYM1.EQ.LSYM2))THEN
          call mh5_put_dset_array_real(wfn_sfs_tdm,
     $      WORK(LTDMZZ),[NTDMZZ,1,1], [0,ISTATE-1,JSTATE-1])
            call mh5_put_dset_array_real(wfn_sfs_tsdm,
     $      WORK(LTSDMZZ),[NTDMZZ,1,1], [0,ISTATE-1,JSTATE-1])
        END IF
#endif

C End of TRD1/2 case
      END IF

      IF (IFHAM .AND.
     &      .NOT.(IFHEXT.or.IFHEFF.or.IFEJOB)) THEN
        HZERO=ECORE*SIJ
        HIJ = HZERO+HONE+HTWO
#ifdef _DMRG_
!       !> temporary fix for DMRG-RASSI until we also compute the 2-TDMs by default
!          anyways, the states ought to be orthogonal...
        if(doDMRG)then
          write(6,*) ' ecore is ... ',ecore
          write(6,*) ' state energy is ...  ',istate,
     &    dmrg_external%dmrg_state_specific(istate)
          write(6,*) ' state energy is ...  ',jstate,
     &    dmrg_external%dmrg_state_specific(jstate)

          hij = dmrg_external%dmrg_state_specific(istate)

          if(istate /= jstate) hij = 0.0d0

          HAM(ISTATE,JSTATE) = HIJ
          HAM(JSTATE,ISTATE) = HIJ
        else
#endif
!         if(istate /= jstate) hij = 0.0d0 ! debug test

          HAM(ISTATE,JSTATE) = HIJ
          HAM(JSTATE,ISTATE) = HIJ
#ifdef _DMRG_
        end if
#endif
        IF(IPGLOB.GE.DEBUG) THEN
          WRITE(6,'(1x,a,2I5)')' ISTATE, JSTATE:',ISTATE,JSTATE
          WRITE(6,'(1x,a,f16.8)')' HZERO=',HZERO
          WRITE(6,'(1x,a,f16.8)')' HONE =',HONE
          WRITE(6,'(1x,a,f16.8)')' HTWO =',HTWO
          WRITE(6,'(1x,a,f16.8)')' HIJ  =',HIJ
        END IF
      END IF
C End of loops over states.

 100  CONTINUE
      END DO
      END DO

#ifdef _DMRG_
      IF(IPGLOB.GE.DEBUG) THEN
         write(6,*) 'full SF-HAMILTONIAN '
         write(6,*) 'dimension: ',nstate**2
         call pretty_print_util(HAM,1,nstate,1,nstate,
     &                          nstate,nstate,1,6)
      END IF
#endif

      IF(WFTP1.EQ.'GENERAL ') THEN
#ifdef _DMRG_
        if(.not.dodmrg)then
#endif
           CALL CXCLOSE(ISGSTR1,ICISTR1,IXSTR1)
           CALL SGCLOSE(ISGSTR1)
#ifdef _DMRG_
        end if
#endif
      END IF
      IF(WFTP2.EQ.'GENERAL ') THEN
#ifdef _DMRG_
        if(.not.dodmrg)then
#endif
          CALL CXCLOSE(ISGSTR2,ICISTR2,IXSTR2)
          CALL SGCLOSE(ISGSTR2)
#ifdef _DMRG_
        end if
#endif
      END IF

      IF(JOB1.NE.JOB2) THEN
        CALL GETMEM('GTDMTRA1','FREE','REAL',LTRA1,NTRA)
        CALL GETMEM('GTDMTRA2','FREE','REAL',LTRA2,NTRA)
      END IF
      CALL GETMEM('GTDMDET1','FREE','REAL',LDET1,NDET1)
      CALL GETMEM('GTDMDET2','FREE','REAL',LDET2,NDET2)
      CALL GETMEM('GTDMCI2','FREE','REAL',LCI2,NCONF2)
      IF (.NOT.NONA) CALL GETMEM('GTDMCI1','FREE','REAL',LCI1,NCONF1)
      CALL KILLSCTAB(LSPNTAB1)
      CALL KILLSCTAB(LSPNTAB2)

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
      CALL KILLOBJ(LREST1)
      CALL KILLOBJ(LREST2)
      CALL KILLOBJ(LCNFTAB1)
      CALL KILLOBJ(LCNFTAB2)
      CALL KILLOBJ(LFSBTAB1)
      CALL KILLOBJ(LFSBTAB2)
      CALL GETMEM('OrbMap','Free','Inte',LOMAP,NASORB)

      CALL QEXIT(ROUTINE)
      RETURN
      END
