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
      Subroutine Gamma_Blocks(iTable,nBlocks,nIrrep)
      Integer idid(8), iTable(6,nBlocks)
*
      iBlock=0
*
*----- Observe the loop ordering straight from GAMSRT in Aces 2.
*      THIS LOOP ORDER DOES NOT COMPLY WITH THE MANUAL!
*
*----- AAAA
*
c     Write (*,*) 'AAAA'
      IND=0
      DO 11 IRREP=1,NIRREP
       IND=IND+1
       iBlock=iBlock+1
       iTable(1,iBlock)=1
       iTable(2,iBlock)=IRREP-1
       iTable(3,iBlock)=IRREP-1
       iTable(4,iBlock)=IRREP-1
       iTable(5,iBlock)=IRREP-1
       iTable(6,iBlock)=IND
c      Write (*,*) IND,IRREP-1,IRREP-1,IRREP-1,IRREP-1
 11   CONTINUE
*
*----- AABB
*
c     Write (*,*) 'AABB'
      IND=0
      DO 12 IRREP1=2,NIRREP
      DO 120 IRREP2=1,IRREP1-1
       IND=IND+1
       iBlock=iBlock+1
       iTable(1,iBlock)=2
       iTable(2,iBlock)=IRREP2-1
       iTable(3,iBlock)=IRREP2-1
       iTable(4,iBlock)=IRREP1-1
       iTable(5,iBlock)=IRREP1-1
       iTable(6,iBlock)=IND
c      Write (*,*) IND,IRREP2-1,IRREP2-1,IRREP1-1,IRREP1-1
120    CONTINUE
12    CONTINUE
*
*----- ABAB
*
*     Write (*,*) 'ABAB'
      IND=0
      DO 13 IRREP=2,NIRREP
      DO 130 IRREP1=1,NIRREP
       IRREP2=IEor(IRREP-1,IRREP1-1)+1
       IF(IRREP2.GT.IRREP1) THEN
        IND=IND+1
        iBlock=iBlock+1
        iTable(1,iBlock)=3
        iTable(2,iBlock)=IRREP1-1
        iTable(3,iBlock)=IRREP2-1
        iTable(4,iBlock)=IRREP1-1
        iTable(5,iBlock)=IRREP2-1
        iTable(6,iBlock)=IND
*       Write (*,*) IND,IRREP1-1,IRREP2-1,IRREP1-1,IRREP2-1
       ENDIF
130    CONTINUE
13    CONTINUE
*
*----- ABCD
*
c     Write (*,*) 'ABCD'
      IND=0
      DO 313 IRREP=2,NIRREP
       DO 314 IRREP1=1,NIRREP
        IRREP2=iEor(IRREP-1,IRREP1-1)+1
        IF(IRREP2.LT.IRREP1)GOTO 314
        IBOT=MAX(IRREP1,IRREP2)+1
        CALL IZERO(IDID,8)
        DO 315 ITMP=IBOT,NIRREP
         IRREP4=iEor(ITMP-1,IRREP-1)+1
         IRREP3=MIN(ITMP,IRREP4)
         IRREP4=MAX(ITMP,IRREP4)
         IF(MAX(IDID(IRREP4),IDID(IRREP3)).NE.0)GOTO 315
         IDID(IRREP3)=1
         IDID(IRREP4)=1
         IND=IND+1
         iBlock=iBlock+1
         iTable(1,iBlock)=4
         iTable(2,iBlock)=IRREP3-1
         iTable(3,iBlock)=IRREP4-1
         iTable(4,iBlock)=IRREP1-1
         iTable(5,iBlock)=IRREP2-1
         iTable(6,iBlock)=IND
c        Write (*,*) IND,IRREP3-1,IRREP4-1,IRREP1-1,IRREP2-1
315     CONTINUE
314    CONTINUE
313   CONTINUE
*
      End
