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
      SUBROUTINE TRCTL_MCLR()
      use Arrays, only: CMO
      use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
      use MCLR_Data, only: ipCM
      use MCLR_Data, only: FnHlf2,FnHlf3,FnTri1,FnTri2,FnTri3,FnTri4,
     &                     FnTri5
      use MCLR_Data, only: LuHlf2,LuHlf3,LuTri1,LuTri2,LuTri3,LuTri4,
     &                     LuTri5
      use input_mclr, only: nSym,nAsh,nBas,nFro,nIsh
*
*     Two-electron integral transformation program: control section
*
*     Purpose: Set up of memory locations (decide if out of core is
*              needed. Loop over symmetry blocks of AO integrals.
*              The transformation routine TRAMO is called for each
*              symmetry block of integrals.
*
      IMPLICIT None
*

      Integer, PARAMETER :: LIOTAB=512*512
      integer,save:: toca(8,8,8),tocb(36,8),tocc(5,36,8)
      Integer, Allocatable:: Hlf1(:,:)
      Real*8, Allocatable:: Buffer(:)
      Integer iAD14,iAd13,iAd23,iAd24,iAd34
      Integer MemX,ipB,iBatch,iSP,nBP,nDP,iSQ,nBQ,nAQ,nDQ,nSPQ,iSR,nBR,
     &        nAR,nDR,nSPQR,iSS,nBS,nAS,nDS,nSPQRS,nORBP,nBPQRS,IntBuf,
     &        ipi,lW1,nW1,lW2,nW2,lW3,nW3,lW4,nW4,lW5,nW5,iSPQ,iSRS,nAP
*                                                                      *
************************************************************************
*                                                                      *
      integer i,j,itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
      CALL DANAME_wa(LUTRI1,FNTRI1)
      CALL DANAME_wa(LUTRI2,FNTRI2)
      CALL DANAME_wa(LUTRI3,FNTRI3)
      CALL DANAME_wa(LUTRI4,FNTRI4)
      CALL DANAME_wa(LUTRI5,FNTRI5)
      Call iCopy(8**3,[-1],0,toca,1)
      Call iCopy(8*36,[-1],0,tocb,1)
      Call iCopy(8*36,[-1],0,tocc,1)
      IAD14=0
      IAD13=0
      IAD23=0
      IAD24=0
      IAD34=0
      Call mma_allocate(Hlf1,LIOTAB,4,Label='Hlf1')
*
*     Precompute start points
*
*
*     Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
*     Note that the integrals on LUTWOAO have to be sorted in the
*     same order as the loop structure below.
*
      Call mma_maxDBLE(MEMX)
      MEMX = MAX(MEMX-MEMX/10,0)
      Call mma_allocate(Buffer,MEMX,Label='Buffer')
      ipB=1
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
     &                   (INTBUF,
     &                   Buffer(LW1:LW1+NW1-1),NW1,
     &                   Buffer(LW2:LW2+NW2-1),NW2,
     &                   Buffer(LW3:LW3+NW3-1),NW3,
     &                   Buffer(LW4:LW4+NW4-1),NW4,
     &                   Buffer(LW5:LW5+NW5-1),NW5,
     &                   nBP,nBQ,nBR,nBS,iSP,iSQ,iSR,iSS,
     &                   nAP,nAQ,nAR,nAS,
     &                   CMO(ipCM(iSP)+nBP*nDP),
     &                   CMO(ipCM(iSQ)+nBQ*nDQ),
     &                   CMO(ipCM(iSR)+nBR*nDR),
     &                   CMO(ipCM(iSS)+nBS*nDS),
     &                   iAD13,iAD14,iAD23,iAD24,iAD34,
     &                   TocC(1,iSPQ,iSRS),
     &                   Hlf1(:,1),Hlf1(:,2),
     &                   Hlf1(:,3),Hlf1(:,4), LIOTAB)

*
*             End of loop over quadruples of symmetries
*
101         CONTINUE
102       CONTINUE
103     CONTINUE
104   CONTINUE
*
      CALL mma_deallocate(Buffer)
      Call mma_deallocate(Hlf1)
      Call DAName_wa(LUHLF2,FNHLF2)
      Call DAName_wa(LUHLF3,FNHLF3)
      Call DaEras(LUHLF2)
      Call DaEras(LUHLF3)
*
      END SUBROUTINE TRCTL_MCLR
