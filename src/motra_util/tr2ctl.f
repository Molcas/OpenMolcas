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
      SUBROUTINE TR2CTL(CMO)
*
*     Two-electron integral transformation program: control section
*
*     Purpose: Set up of memory locations (decide if out of core is
*              needed. Loop over symmetry blocks of AO integrals.
*              The transformation routine TRAMO is called for each
*              symmetry block of integrals.
*
      IMPLICIT REAL*8 (A-H,O-Z)
*

#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "WrkSpc.fh"
*
      Real*8 CMO(*)
      DIMENSION NBSX(8),KEEP(8),ISTSQ(8)
      Logical FoundTwoEls,DoDirect,DoCholesky
*
      Call qEnter('Tr2Ctl')
*
*     Set time at start of transformation
*
      CALL SETTIM
*
*     Initiate unit LUTWOMO (two-electron integrals in MO basis)
*
      CALL DANAME_MF(LUTWOMO,FNTWOMO)
      IAD13=0
      Call iCopy(nTraToc,0,0,iTraToc,1)
      CALL iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)
*
*     Initiate unit LUTWOAO (two-electron integrals in AO basis)
*
      Call f_Inquire(FnTwoAo,FoundTwoEls)
      Call DecideOnDirect(.False.,FoundTwoEls,DoDirect,DoCholesky)
      If (.not.DoCholesky) CALL OPNORD(IRC,0,FNTWOAO,LUTWOAO)

      CALL GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)

*
*     Compare content of 1el and 2el integral file
*
      IF ( NSYM2.NE.NSYM ) THEN
         Write (6,*) 'Tr2Ctl: NSYM2.NE.NSYM'
         Write (6,*) 'NSYM2=',NSYM2
         Write (6,*) 'NSYM=',NSYM
         Call QTrace()
         Call Abend()
      END IF
      DO ISYM=1,NSYM
        NB1=NBAS(ISYM)
        NB2=NBSX(ISYM)
        IF ( NB1.NE.NB2 ) THEN
           Write (6,*) 'Tr2Ctl: NB1.NE.NB2'
           Write (6,*) 'NB1=',NB1
           Write (6,*) 'NB2=',NB2
           Call QTrace()
           Call Abend()
        END IF
      END DO
*
*     Precompute start points
*
      ISTBS=1
      DO 20 ISYM=1,NSYM
        ISTSQ(ISYM)=ISTBS+NBAS(ISYM)*NFRO(ISYM)
        ISTBS=ISTBS+NBAS(ISYM)*NBAS(ISYM)
20    CONTINUE
*
*     Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
*     Note that the integrals on LUTWOAO have to be sorted in the
*     same order as the loop structure below.
*
      If (iPrint.GE.0) WRITE(6,2000)
2000  FORMAT(/7X,'SYMMETRY',2X,'BASIS FUNCTIONS',6X,' ORBITALS',6X,
     *       'INTEGRALS   CPU(SEC)  I/O(SEC)')
      IBATCH=0
      DO 104 NSP=1,NSYM
        NBP=NBAS(NSP)
        NOP=NORB(NSP)
        KEEPP=KEEP(NSP)
        LMOP=ISTSQ(NSP)
        ISP=NSP
        DO 103 NSQ=1,NSP
          NBQ=NBAS(NSQ)
          NOQ=NORB(NSQ)
          KEEPQ=KEEP(NSQ)
          LMOQ=ISTSQ(NSQ)
          NSPQ=IEOR(NSP-1,NSQ-1)+1
          ISQ=NSQ
          DO 102 NSR=1,NSP
            NBR=NBAS(NSR)
            NOR=NORB(NSR)
            KEEPR=KEEP(NSR)
            LMOR=ISTSQ(NSR)
            NSPQR=IEOR(NSPQ-1,NSR-1)+1
            ISR=NSR
            NSSM=NSR
            IF(NSR.EQ.NSP) NSSM=NSQ
            DO 101 NSS=1,NSSM
              NBS=NBAS(NSS)
              NOS=NORB(NSS)
              KEEPS=KEEP(NSS)
              LMOS=ISTSQ(NSS)
              NSPQRS=IEOR(NSPQR-1,NSS-1)+1
              IF( NSPQRS.NE.1 ) GOTO 101
              IBATCH=IBATCH+1
              ISS=NSS
*
*             Check the loop conditions and skip transformation step
*             if possible
*
              KEEPT=KEEPP+KEEPQ+KEEPR+KEEPS
              NORBP=NOP*NOQ*NOR*NOS
              IF( NORBP.EQ.0 ) GOTO 101
              IF( KEEPT.NE.0 ) THEN
                Write (6,*) 'Tr2Ctl: NORBP.NE.0.AND.KEEPT.NE.0'
                Write (6,*) 'NORBP=',NORBP
                Write (6,*) 'KEEPT=',KEEPT
                Call QTrace()
                Call Abend()
              End If
*
*             Define matrix sizes
*
              IF(ISR.EQ.ISS) THEN
                NBPQ=(NBP+NBP**2)/2
                NBRS=(NBR+NBR**2)/2
                NOVX=(NOR+NOR**2)/2
              ELSE
                NBPQ=NBP*NBQ
                NBRS=NBR*NBS
                NOVX=NOR*NOS
              ENDIF
*
*             Set up dynamic memory
*
* INTBUF=Size of output buffer.
*             INTBUF=256*256
              INTBUF=16 * 256*256
              NW1=2*KBUF
              NW2=MAX(INTBUF,NBP*NOQ,NBQ*NOP)
* NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
* LBUF-1 must be at least 'klB':
              NW2=MAX(NW2,NBRS+1,NBPQ+1)
              NW3=MAX(NBR**2,NBP**2,NOQ*NBP,NOVX)
              NW4=MAX(NBR*NOS,NBQ*NOP)
              CALL GETMEM('OUTBUF','ALLO','REAL',LW1,NW1)
              CALL GETMEM('X1','ALLO','REAL',LW2,NW2)
              CALL GETMEM('X2','ALLO','REAL',LW3,NW3)
              CALL GETMEM('X3','ALLO','REAL',LW4,NW4)
              Call GetMem('iDsk','Allo','Inte',ipiDsk,3*nOVX)
              CALL GETMEM('VXPQ','MAX','REAL',LW5,MEMX)

              If (DoCholesky) Then ! save some space for GenInt
                 MEMX = MAX(MEMX-MEMX/10,0)
              End If

              CALL GETMEM('VXPQ','ALLO','REAL',LW5,MEMX)
*
              IF(MEMX.LT.NOVX) then
                Write (6,*) 'Tr2Ctl: MEMX.LT.NOVX'
                Write (6,*) 'MEMX=',MEMX
                Write (6,*) 'NOVX=',NOVX
                Call QTrace()
                Call Abend()
              End If
*
*             transform the symmetry block (ISP,ISQ|ISR,ISS)
*
              iTraToc(IBATCH)=IAD13
* NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
              CALL TRAMO(NW2   ,WORK(LW1),nW1,
     &                          WORK(LW2),nW2,
     &                          WORK(LW3),nW3,
     &                          WORK(LW4),nW4,
     &                          WORK(LW5),MEMX,CMO,
     &                          iWork(ipiDsk),nOVX)
              CALL TIMING(CPT,CPE,TIOT,TIOE)
              If (iPrint.GE.0) WRITE(6,2100)
     &    ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NOP,NOQ,NOR,NOS,LTUVX,CPE,TIOE
2100          FORMAT(7X,4I2,1X,4I4,2X,4I4,3X,I9,F11.2,F10.2)
              Call Xflush(6)

*
*             Deallocate work space
*
              CALL GETMEM('VXPQ','FREE','REAL',LW5,MEMX)
              Call GetMem('iDsk','Free','Inte',ipiDsk,3*nOVX)
              CALL GETMEM('X3','FREE','REAL',LW4,NW4)
              CALL GETMEM('X2','FREE','REAL',LW3,NW3)
              CALL GETMEM('X1','FREE','REAL',LW2,NW2)
              CALL GETMEM('OUTBUF','FREE','REAL',LW1,NW1)
*
*             End of loop over quadruples of symmetries
*
101         CONTINUE
102       CONTINUE
103     CONTINUE
104   CONTINUE
      CALL TIMING(CPT,CPE,TIOT,TIOE)
      If (iPrint.GE.0) WRITE(6,2200) CPT,TIOT
2200  FORMAT(/6X,' TOTAL CPU TIME(SEC)',F8.2,'TOTAL I/O TIME(SEC)',F8.2)
*
*     Closed LUTWOAO
*
      If (.not.DoCholesky)then
         CALL CLSORD(IRC,0)
      End If
*
*     Write address block to LUTWOMO
*
      IAD13=0
      CALL iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)
      IF ( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,2A)')'DISK ADRESSES FOR SYMMETRY BLOCKS OF ',
     *                     'TRANSFORMED TWO-ELECTRON INTEGRALS'
        WRITE(6,'(6X,10I8)')(iTraToc(I),I=1,nTraToc)
      ENDIF

      CALL DACLOS(LUTWOMO)
      Call qExit('Tr2Ctl')
*
*     exits
*
      RETURN
      END
