************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1990, Bjorn O. Roos                                    *
************************************************************************
      SUBROUTINE FCKPT2(CMOO,CMON,FI,FP,FTR,VEC,WO,SQ,CMOX)
*
* Purpose: To diagonalize the inactive, active,
* and external parts of the Fock matrix (FI+FA) for CASPT2
* and order the eigenvalues and eigenvectors after energy.
* This ordering is of value in subsequent CI and CASPT2
* calculations. All diagonal elements of
* the transformed Fock matrix are collected in FDIAG for
* later printing. The new MO's are written onto JOBIPH
* in address IADR15(9). Also the inactive Fock matrix FI
* is saved on Jobiph for use in CASPT2.
* Note:these orbitals leave the CI expansion
* invariant only for CAS wave functions.
* Called from SXCTL if IFINAL=1 (after last MC iteration)
*
* Modifications: - Produce full matrix FP and write it to JOBIPH.
* B. Roos, Lund, June 1990
*
* ********** IBM-3090 Release 88 09 07 **********
*

      IMPLICIT REAL*8 (A-H,O-Z)

#ifdef _ENABLE_CHEMPS2_DMRG_
      Integer iChMolpro(8)
#endif


#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
#include "raswfn.fh"
      Parameter (ROUTINE='FCKPT2  ')

      DIMENSION CMOO(*),CMON(*),FI(*),FP(*),FTR(*),VEC(*),
     &          WO(*),SQ(*),CMOX(*)

      CALL QENTER(ROUTINE)
* Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*
      IB=0
      ISTMO1=1
      ISTFCK=0
      ID=0

#ifdef _ENABLE_CHEMPS2_DMRG_
      ifock=1
      norbtot = 0
      do iiash=1,nsym
        norbtot = norbtot + nAsh(iiash)
      enddo

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab_BIS(nSym,Label,iChMolpro)

* Convert orbital symmetry into MOLPRO format
      Call Getmem('OrbSym','Allo','Inte',lOrbSym,NAC)
      iOrb=1
      Do iSym=1,nSym
        Do jOrb=1,NASH(iSym)
          iWork(lOrbSym+iOrb-1)=iChMolpro(iSym)
          iOrb=iOrb+1
        End Do
      End Do
      lSymMolpro=iChMolpro(lSym)

      LuFCK=isFreeUnit(27)
*      open ( unit = LuFCK, file = "FOCK_CHEMPS2",
*     &        action="write", status="replace" )
      call molcas_open(LuFCK,'FOCK_CHEMPS2')
      write(LuFCK,'(1X,A12,I2,A1)') '&FOCK NACT= ', norbtot,','
      write(LuFCK,'(2X,A7)',ADVANCE = "NO") 'ORBSYM='
      do iOrb=1,norbtot
        write(LuFCK,'(I1,A1)',ADVANCE = "NO") iWork(lOrbSym+iOrb-1),','
      enddo
      write(LuFCK,*)
      write(LuFCK,*) '/'
      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)
#endif

      DO ISYM=1,NSYM
       NBF=NBAS(ISYM)
       NFO=NFRO(ISYM)
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NR1=NRS1(ISYM)
       NR2=NRS2(ISYM)
       NR3=NRS3(ISYM)
       NEO=NSSH(ISYM)
       NOO=NFO+NIO+NAO
       NOT=NIO+NAO+NEO
       NOC=NIO+NAO
       ISTMO=ISTMO1+NFO*NBF
************************************************************************
* Frozen orbitals (move MO's to CMON and set FDIAG to zero)
************************************************************************
       IF(NFO.NE.0) THEN
        NFNB=NBF*NFO
        CALL DCOPY_(NFNB,CMOO(ISTMO1),1,CMON(ISTMO1),1)
        DO  NF=1,NFO
         FDIAG(IB+NF)=0.0D0
        END DO
       ENDIF
*
* Clear the MO transformation matrix CMOX
*
       CALL FZERO(CMOX,NOT*NOT)
*
************************************************************************
* Inactive part of the Fock matrix
************************************************************************
*
       IF(NIO.NE.0) THEN
* MOVE FP TO TRIANGULAR FORM
        NIJ=0
        DO NI=1,NIO
         DO NJ=1,NI
          NIJ=NIJ+1
          FTR(NIJ)=FP(NIJ+ISTFCK)
          IF(IXSYM(IB+NFO+NI).NE.IXSYM(IB+NFO+NJ)) FTR(NIJ)=0.0D0
         END DO
        END DO
* DIAGONALIZE
        NIO2=NIO**2
        CALL FZERO(VEC,NIO2)
        II=1
        DO NI=1,NIO
         VEC(II)=1.0D0
         II=II+NIO+1
        END DO
        CALL Jacob(FTR,VEC,NIO,NIO)
* MOVE EIGENVALUES TO FDIAG.
*
        II=0
        NO1=IB+NFO
        DO NI=1,NIO
         II=II+NI
         FDIAG(NO1+NI)=FTR(II)
        END DO
*
* Sort eigenvalues and orbitals after energy
*
        IF(NIO.GT.1) THEN
         NIO1=NIO-1
         DO NI=1,NIO1
          NI1=NI+1
          MIN=NI
          DO NJ=NI1,NIO
           IF(FDIAG(NO1+NJ).LT.FDIAG(NO1+MIN)) MIN=NJ
          END DO
          IF(MIN.EQ.NI) GO TO 20
          FMIN=FDIAG(NO1+MIN)
          FDIAG(NO1+MIN)=FDIAG(NO1+NI)
          FDIAG(NO1+NI)=FMIN
          CALL DSWAP_(NIO,VEC(1+NIO*(NI-1)),1,VEC(1+NIO*(MIN-1)),1)
20       CONTINUE
         END DO
        ENDIF
        CALL DGEADD(CMOX,NOT,'N',
     *              VEC,NIO,'N',CMOX,NOT,NIO,NIO)
       ENDIF
*
************************************************************************
* Active part of the Fock matrix
************************************************************************
*
       IF(NAO.NE.0) THEN
************************************************************************
* RAS1 part of the Fock matrix
************************************************************************
        IF(NR1.NE.0) THEN
* MOVE FP TO TRIANGULAR FORM
         NTU=0
         DO NT=1,NR1
          DO NU=1,NT
           NTU=NTU+1
           NTT=NT+NIO
           NUT=NU+NIO
           NTUT=ISTFCK+(NTT**2-NTT)/2+NUT
           FTR(NTU)= FP(NTUT)
           IF(IXSYM(IB+NFO+NTT).NE.IXSYM(IB+NFO+NUT)) FTR(NTU)=0.0D0
          END DO
         END DO
* DIAGONALIZE
         NR12=NR1**2
         CALL FZERO(VEC,NR12)
         II=1
         DO NT=1,NR1
          VEC(II)=1.0D0
          II=II+NR1+1
         END DO
         CALL Jacob(FTR,VEC,NR1,NR1)
*
* Move eigenvalues to FDIAG.
*
         II=0
         NO1=IB+NFO+NIO
         DO NT=1,NR1
          II=II+NT
          FDIAG(NO1+NT)=FTR(II)
         END DO
*
* Sort eigenvalues and orbitals after energy
*
         IF(NR1.GT.1) THEN
          NR11=NR1-1
          DO NT=1,NR11
           NT1=NT+1
           MIN=NT
           DO NU=NT1,NR1
            IF(FDIAG(NO1+NU).LT.FDIAG(NO1+MIN)) MIN=NU
           END DO
           IF(MIN.EQ.NT) GO TO 41
           FMIN=FDIAG(NO1+MIN)
           FDIAG(NO1+MIN)=FDIAG(NO1+NT)
           FDIAG(NO1+NT)=FMIN
           CALL DSWAP_(NR1,VEC(1+NR1*(NT-1)),1,VEC(1+NR1*(MIN-1)),1)
41         CONTINUE
          END DO
         ENDIF
         CALL DGEADD(CMOX(1+NOT*NIO+NIO),NOT,'N',
     *               VEC,NR1,'N',CMOX(1+NOT*NIO+NIO),NOT,NR1,NR1)

#ifdef _ENABLE_CHEMPS2_DMRG_
         II=0
         NO1=IB+NFO+NIO
         DO NT=1,NR1
           write(LuFCK,'(1X,E23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
           ifock = ifock + 1
         END DO
#endif
        ENDIF ! NR1




************************************************************************
* RAS2 part of the Fock matrix
************************************************************************
        IF(NR2.NE.0) THEN
* MOVE FP TO TRIANGULAR FORM
         NTU=0
         DO NT=1,NR2
          DO NU=1,NT
           NTU=NTU+1
           NTT=NT+NIO+NR1
           NUT=NU+NIO+NR1
           NTUT=ISTFCK+(NTT**2-NTT)/2+NUT
           FTR(NTU)= FP(NTUT)
           IF(IXSYM(IB+NFO+NTT).NE.IXSYM(IB+NFO+NUT)) FTR(NTU)=0.0D0
          END DO
         END DO
* DIAGONALIZE
         NR22=NR2**2
         CALL FZERO(VEC,NR22)
         II=1
         DO NT=1,NR2
          VEC(II)=1.0D0
          II=II+NR2+1
         END DO
         CALL Jacob(FTR,VEC,NR2,NR2)
*
* Move eigenvalues to FDIAG.
*
         II=0
         NO1=IB+NFO+NIO+NR1
         DO NT=1,NR2
          II=II+NT
          FDIAG(NO1+NT)=FTR(II)
         END DO
*
* Sort eigenvalues and orbitals after energy
*
         IF(NR2.GT.1) THEN
          NR21=NR2-1
          DO NT=1,NR21
           NT1=NT+1
           MIN=NT
           DO NU=NT1,NR2
            IF(FDIAG(NO1+NU).LT.FDIAG(NO1+MIN)) MIN=NU
           END DO
           IF(MIN.EQ.NT) GO TO 42
           FMIN=FDIAG(NO1+MIN)
           FDIAG(NO1+MIN)=FDIAG(NO1+NT)
           FDIAG(NO1+NT)=FMIN
           CALL DSWAP_(NR2,VEC(1+NR2*(NT-1)),1,VEC(1+NR2*(MIN-1)),1)
42         CONTINUE
          END DO
         ENDIF
         CALL DGEADD(CMOX(1+NOT*(NIO+NR1)+NIO+NR1),NOT,'N',
     *               VEC,NR2,'N',
     *               CMOX(1+NOT*(NIO+NR1)+NIO+NR1),NOT,NR2,NR2)

#ifdef _ENABLE_CHEMPS2_DMRG_
         II=0
         NO1=IB+NFO+NIO+NR1
         DO NT=1,NR2
           write(LuFCK,'(1X,E23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
           ifock = ifock + 1
         END DO
#endif
        ENDIF ! NR2

************************************************************************
* RAS3 part of the Fock matrix
************************************************************************
        IF(NR3.NE.0) THEN
* MOVE FP TO TRIANGULAR FORM
         NTU=0
         DO NT=1,NR3
          DO NU=1,NT
           NTU=NTU+1
           NTT=NT+NIO+NR1+NR2
           NUT=NU+NIO+NR1+NR2
           NTUT=ISTFCK+(NTT**2-NTT)/2+NUT
           FTR(NTU)= FP(NTUT)
           IF(IXSYM(IB+NFO+NTT).NE.IXSYM(IB+NFO+NUT)) FTR(NTU)=0.0D0
          END DO
         END DO
* DIAGONALIZE
         NR32=NR3**2
         CALL FZERO(VEC,NR32)
         II=1
         DO NT=1,NR3
          VEC(II)=1.0D0
          II=II+NR3+1
         END DO
         CALL Jacob(FTR,VEC,NR3,NR3)
*
* Move eigenvalues to FDIAG.
*
         II=0
         NO1=IB+NFO+NIO+NR1+NR2
         DO NT=1,NR3
          II=II+NT
          FDIAG(NO1+NT)=FTR(II)
         END DO
*
* Sort eigenvalues and orbitals after energy
*
         IF(NR3.GT.1) THEN
          NR31=NR3-1
          DO NT=1,NR31
           NT1=NT+1
           MIN=NT
           DO NU=NT1,NR3
            IF(FDIAG(NO1+NU).LT.FDIAG(NO1+MIN)) MIN=NU
           END DO
           IF(MIN.EQ.NT) GO TO 43
           FMIN=FDIAG(NO1+MIN)
           FDIAG(NO1+MIN)=FDIAG(NO1+NT)
           FDIAG(NO1+NT)=FMIN
           CALL DSWAP_(NR3,VEC(1+NR3*(NT-1)),1,VEC(1+NR3*(MIN-1)),1)
43         CONTINUE
          END DO
         ENDIF
         CALL DGEADD(CMOX(1+NOT*(NIO+NR1+NR2)+NIO+NR1+NR2),NOT,'N',
     *               VEC,NR3,'N',
     *               CMOX(1+NOT*(NIO+NR1+NR2)+NIO+NR1+NR2),NOT,NR3,NR3)

#ifdef _ENABLE_CHEMPS2_DMRG_
         II=0
         NO1=IB+NFO+NIO+NR1+NR2
         DO NT=1,NR3
           write(LuFCK,'(1X,E23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
           ifock = ifock + 1
         END DO
#endif
        ENDIF ! NR3

       ENDIF ! NAO
*
************************************************************************
* external part of the Fock matrix
************************************************************************
       IF(NEO.NE.0) THEN
* MOVE FP TO TRIANGULAR FORM
        NAB=0
        DO NA=1,NEO
         DO NB=1,NA
          NAB=NAB+1
          NAT=NA+NIO+NAO
          NBT=NB+NIO+NAO
          NABT=ISTFCK+(NAT**2-NAT)/2+NBT
          FTR(NAB)=FP(NABT)
          IF(IXSYM(IB+NFO+NAT).NE.IXSYM(IB+NFO+NBT)) FTR(NAB)=0.0D0
         END DO
        END DO
* DIAGONALIZE
        NEO2=NEO**2
        CALL FZERO(VEC,NEO2)
        II=1
        DO NA=1,NEO
         VEC(II)=1.0D0
         II=II+NEO+1
        END DO
        CALL Jacob(FTR,VEC,NEO,NEO)
*
* Move eigenvalues to FDIAG.
*
        II=0
        NO1=IB+NFO+NIO+NAO
        DO NA=1,NEO
         II=II+NA
         FDIAG(NO1+NA)=FTR(II)
        END DO
*
* Sort eigenvalues and orbitals after energy
*
        IF(NEO.GT.1) THEN
         NEO1=NEO-1
         DO NA=1,NEO1
          NA1=NA+1
          MIN=NA
          DO NB=NA1,NEO
           IF(FDIAG(NO1+NB).LT.FDIAG(NO1+MIN)) MIN=NB
          END DO
          IF(MIN.EQ.NA) GO TO 60
          FMIN=FDIAG(NO1+MIN)
          FDIAG(NO1+MIN)=FDIAG(NO1+NA)
          FDIAG(NO1+NA)=FMIN
          CALL DSWAP_(NEO,VEC(1+NEO*(NA-1)),1,VEC(1+NEO*(MIN-1)),1)
60       CONTINUE
         END DO
        ENDIF
        CALL DGEADD(CMOX(1+NOT*NOC+NOC),NOT,'N',
     *              VEC,NEO,'N',CMOX(1+NOT*NOC+NOC),NOT,NEO,NEO)
       ENDIF
*
* Transform molecular orbitals
*
        IF ( NBF*NTOT.GT.0 )
     &  CALL DGEMM_('N','N',
     &              NBF,NOT,NOT,
     &              1.0d0,CMOO(ISTMO),NBF,
     &              CMOX,NOT,
     &              0.0d0,CMON(ISTMO),NBF)
*
************************************************************************
* Deleted orbitals (move MO's and set zero to FDIAG)
************************************************************************
       NDO=NDEL(ISYM)
       IF(NDO.NE.0) THEN
        NDNB=NDO*NBF
        IST=ISTMO1+NBF*(NOO+NEO)
        CALL DCOPY_(NDNB,CMOO(IST),1,CMON(IST),1)
        DO ND=1,NDO
         FDIAG(IB+NBF-NDO+ND)=0.0D0
        END DO
       ENDIF
*
* Transform inactive Fock matrix FI and the CASPT2 matrix FP
*
       IF(NOT.GT.0) THEN
        CALL SQUARE(FI(ISTFCK+1),SQ,1,NOT,NOT)
        CALL DGEMM_('N','N',
     &              NOT,NOT,NOT,
     &              1.0d0,SQ,NOT,
     &              CMOX,NOT,
     &              0.0d0,VEC,NOT)
        CALL DGEMM_('T','N',
     &              NOT,NOT,NOT,
     &              1.0d0,CMOX,NOT,
     &              VEC,NOT,
     &              0.0d0,SQ,NOT)
*
* Move transformed Fock matrix back to FI
*
        NPQ=ISTFCK
        DO NP=1,NOT
         DO NQ=1,NP
          NPQ=NPQ+1
          FI(NPQ)=SQ(NOT*(NP-1)+NQ)
         END DO
        END DO
*
* The FP matrix
*
        CALL SQUARE(FP(ISTFCK+1),SQ,1,NOT,NOT)
        CALL DGEMM_('N','N',
     &              NOT,NOT,NOT,
     &              1.0d0,SQ,NOT,
     &              CMOX,NOT,
     &              0.0d0,VEC,NOT)
        CALL DGEMM_('T','N',
     &              NOT,NOT,NOT,
     &              1.0d0,CMOX,NOT,
     &              VEC,NOT,
     &              0.0d0,SQ,NOT)
*
* Move transformed Fock matrix back to FP
*
        NPQ=ISTFCK
        DO NP=1,NOT
         DO NQ=1,NP
          NPQ=NPQ+1
          FP(NPQ)=SQ(NOT*(NP-1)+NQ)
         END DO
        END DO
       ENDIF
*
       IB=IB+NBF
       ISTFCK=ISTFCK+(NOT**2+NOT)/2
       ISTMO1=ISTMO1+NBF**2
       ID=ID+(NAO**2+NAO)/2
      END DO
*
#ifdef _ENABLE_CHEMPS2_DMRG_
*      close(27)
      close(LuFCK)
#endif
      IF(IPRLEV.GE.VERBOSE) THEN
       Write(LF,*)' Diagonal elements of the Fock matrix in FCKPT2:'
       Write(LF,'(1X,10F11.6)') (FDIAG(I),I=1,NTOT)
      END IF
*
************************************************************************
* Orthogonalise new orbitals
************************************************************************
*
      CALL SUPSCH(WO,CMOO,CMON)
      CALL ORTHO_RASSCF(WO,CMOX,CMON,SQ)
*
************************************************************************
* Write new orbitals to JOBIPH/rasscf.h5
************************************************************************
*
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in FCKPT2 after diag and orthog'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=0
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If

      IAD15=IADR15(9)
      CALL DDAFILE(JOBIPH,1,CMON,NTOT2,IAD15)
#ifdef _HDF5_
        call mh5_put_dset(wfn_mocoef,CMON)
#endif
*
* Write FI, FP and FDIAG to JOBIPH
* First remove frozen and deleted part of FDIAG
*
      IF=0
      IFD=0
      DO ISYM=1,NSYM
       NBF=NBAS(ISYM)
       DO NB=1,NBF
        IFD=IFD+1
        IF(NB.GT.NFRO(ISYM).AND.NB.LE.NBF-NDEL(ISYM)) THEN
         IF=IF+1
         SQ(IF)=FDIAG(IFD)
        ENDIF
       END DO
      END DO
*
      IAD15=IADR15(10)
      CALL DDAFILE(JOBIPH,1,FI,NTOT3,IAD15)
      CALL DDAFILE(JOBIPH,1,FP,NTOT3,IAD15)
      CALL DDAFILE(JOBIPH,1,SQ,NORBT,IAD15)
*
      CALL QEXIT(ROUTINE)
      RETURN
      END
