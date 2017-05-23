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
      SUBROUTINE TRCTL_MCLR
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

#include "Input.fh"
      PARAMETER (LIOTAB=512*512)
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "toc.fh"
#include "Files_mclr.fh"
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*     Call qEnter('TrCtl')
*
      CALL DANAME_wa(LUTRI1,FNTRI1)
      CALL DANAME_wa(LUTRI2,FNTRI2)
      CALL DANAME_wa(LUTRI3,FNTRI3)
      CALL DANAME_wa(LUTRI4,FNTRI4)
      CALL DANAME_wa(LUTRI5,FNTRI5)
      Call iCopy(8**3,-1,0,toca,1)
      Call iCopy(8*36,-1,0,tocb,1)
      Call iCopy(8*36,-1,0,tocc,1)
      IAD14=0
      IAD13=0
      IAD23=0
      IAD24=0
      IAD34=0
      Call GetMem('LIOTAB','Allo','Inte',ip_Hlf1,4*LIOTAB)
      ip_Hlf2 = ip_Hlf1+LIOTAB
      ip_Hlf3 = ip_Hlf2+LIOTAB
      ip_Hlf4 = ip_Hlf3+LIOTAB

*
*     Precompute start points
*
*
*     Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
*     Note that the integrals on LUTWOAO have to be sorted in the
*     same order as the loop structure below.
*
      CALL GETMEM('Buffer','MAX','REAL',LW5,MEMX)
      MEMX = MAX(MEMX-MEMX/10,0)
      CALL GETMEM('Buffer','ALLO','REAL',ipb,MEMX)
      IBATCH=0
      DO 104 iSP=1,NSYM
        NBP=NBAS(iSP)
        NAP=NASH(iSP)+nish(isp)
        nDP=nFro(iSP)!+nDel(iSP)
        DO 103 iSQ=1,iSP
          NBQ=NBAS(iSQ)
          NAQ=NASH(iSQ)+nish(isq)
          nDQ=nFro(iSQ)!+nDel(iSQ)
          NSPQ=IEOR(iSP-1,iSQ-1)+1
          DO 102 iSR=1,NSYM
            NBR=NBAS(iSR)
            NAR=NASH(iSR)+nish(isr)
            nDR=nFro(iSR)!+nDel(iSR)
            NSPQR=IEOR(NSPQ-1,iSR-1)+1
            DO 101 iSS=1,iSR
              NBS=NBAS(iSS)
              NAS=NASH(iSS)+nIsh(iSS)
              nDS=nFro(iSS)!+nDel(iSS)
              NSPQRS=IEOR(NSPQR-1,iSS-1)+1
*
*             Check the loop conditions and skip transformation step
*             if possible
*
              NORBP=NAP*NAQ+NAR*NAS+NAP*NAR+NAP*NAS+NAQ*NAR+NAQ*NAS
              NBPQRS=NBP*NBQ*NBR*NBS
              IF( NSPQRS.NE.1 ) GOTO 101
              IBATCH=IBATCH+1
              IF( NORBP.EQ.0 ) GOTO 101
              IF( NBPQRS.EQ.0 ) GOTO 101
              IF ( NAR+NAS.eq.0) Goto 101
*
*             Set up dynamic memory
*
              INTBUF=256*256
              ipi=ipB
              LW1=ipi
              NW1=MAX(NBP*NBQ,NBR*NBS)
              ipi=ipi+nw1
              LW2=ipi
              NW2=MAX(nBR*nAS,nAQ*nBP)
              ipi=ipi+nw2
              lw3=ipi
              NW3=MAX(nAR*nBS,nBQ*nAP)
              ipi=ipi+nw3
              lw4=ipi
              NW4=NAR*NAS
              ipi=ipi+nw4
              lw5=ipi
              NW5=MEMX-nw1-nw2-nw3-nw4
              iSPQ=iTri(iSR,iSS)
              iSRS=Max(iSP,iSQ)
              If (nAR*nAS.ne.0)
     &        TocB(iSPQ,iSRS)=iAD34
              If (nAQ*nAS.ne.0)
     &        TOCA(iSP,iSQ,iSR)=iAD24
              If (iSR.ne.iSS.and.nAQ*nAR.ne.0)
     &         TOCA(iSP,iSQ,iSS)=iAD23
              If (iSP.ne.iSQ.and.nAP*nAS.ne.0)
     &         TOCA(iSQ,iSP,iSR)=iAD14
              If (iSP.ne.iSQ.and.iSR.ne.iSS.and.nAP*nAR.ne.0)
     &         TOCA(iSQ,iSP,iSS)=iAD13
*
*             transform the symmetry block (ISP,ISQ|ISR,ISS)
*
              CALL TRAMO_MCLR
     &                   (INTBUF,Work(LW1),NW1,Work(LW2),NW2,
     &                   Work(LW3),NW3,
     &                   Work(LW4),NW4,Work(LW5),NW5,
     &                   nBP,nBQ,nBR,nBS,iSP,iSQ,iSR,iSS,
     &                   nAP,nAQ,nAR,nAS,
     &                   Work(ipCMO+ipCM(iSP)-1+nBP*nDP),
     &                   Work(ipCMO+ipCM(iSQ)-1+nBQ*nDQ),
     &                   Work(ipCMO+ipCM(iSR)-1+nBR*nDR),
     &                   Work(ipCMO+ipCM(iSS)-1+nBS*nDS),
     &                   iAD13,iAD14,iAD23,iAD24,iAD34,
     &                   TocC(1,iSPQ,iSRS),
     &                   iWork(ip_Hlf1),iWork(ip_Hlf2),
     &                   iWork(ip_Hlf3),iWork(ip_Hlf4), LIOTAB)

*
*             End of loop over quadruples of symmetries
*
101         CONTINUE
102       CONTINUE
103     CONTINUE
104   CONTINUE
*
      CALL GETMEM('Buffer','FREE','REAL',ipb,MEMX)
      Call GetMem('LIOTAB','Free','Inte',ip_Hlf1,4*LIOTAB)
      Call DAName_wa(LUHLF2,FNHLF2)
      Call DAName_wa(LUHLF3,FNHLF3)
      Call DaEras(LUHLF2)
      Call DaEras(LUHLF3)
*
      RETURN
      END
