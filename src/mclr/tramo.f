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
* Copyright (C) 1996, Anders Bernhardsson                              *
*               1996, Bjorn O. Roos                                    *
************************************************************************
      SUBROUTINE TRAMO_MCLR(LBUF,X1,n1,X2,n2,X3,n3,X4,n4,
     &                      Buffer,MEMX,
     &                      NBP,NBQ,NBR,NBS,iSP,iSQ,iSR,iSS,
     &                      nAP,nAQ,nAR,nAS,CMP,CMQ,CMR,CMS,
     &                      iAD13,iAD14,iAD23,iAD24,iAD34,len,
     &                      IDAHLF2,IRLHLF2,
     &                      IDAHLF3,IRLHLF3,LIOTAB)
*
*****************************************************************
*                                                               *
*     Purpose: two-electron transformation routine.             *
*              Transformation is made in core if all            *
*              half-transformed integrals fit.                  *
*              Otherwise sorted integrals are written onto      *
*              unit LUHALF.                                     *
*                                                               *
*****************************************************************
*                                                               *
*  Stolen from MOTRA and destroyed by Anders to make it         *
*  possible to generate all sort of combination of transformed  *
*  and untransformed indexes                                    *
*                                                               *
*****************************************************************
*

      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "Files_mclr.fh"
#include "SysDef.fh"
*
      Integer len(5)
      DIMENSION X1(n1),X2(n2),X3(n3),X4(n4),
     &          Buffer(MemX),
     &          CMP(nBP,nAP),CMQ(nBQ,nAQ),
     &          CMR(nBR,nAR),CMS(nBS,nAS)
*
      Integer IDAHLF2(LIOTAB),IRLHLF2(LIOTAB),
     &        IDAHLF3(LIOTAB),IRLHLF3(LIOTAB)
*
      Call TRAMO_MCLR_INTERNAL(Buffer)
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(len)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine TRAMO_MCLR_INTERNAL(Buffer)
      Use Iso_C_Binding
      Real*8, Target :: Buffer(*)
      Integer, Pointer :: iBuffer(:)
      ione=1
      if (nofile) ione=0
*
*     Set some constants
*
      MemXX=MemX-lBuf-1010
      If (ISS.eq.isr) Then
          NBPQ=(NBP+NBP**2)/2
          NBRS=(NBR+NBR**2)/2
          NARS=(NAR+NAR**2)/2
      ELSE
          NBPQ=NBP*NBQ
          NBRS=NBR*NBS
          NARS=NAR*NAS
      ENDIF
      iBUF2=0
      iBUF3=0
      iAD2=0                        ! Disk address for (kL|ij)
      iAD3=0                        ! Disk address for (Kl|ij)
      iPQUT2=0                      ! Index(ij)  in buffer.
      iPQUT3=0                      ! Just one should be enough but...
      iPQUT4=0                      ! ..easy is better than intelligent
      If (iSP.eq.iSQ) Then
       iMax=nBP*(nBP+1)/2
      Else
       iMax=nBP*nBQ               ! Size of buffer
      End If
      Incore=0                    ! First we try to  transform in core
      If (iss.ne.isr) Then
      nBuf=nBR*nAS+nAR*nBS
      Else
      nbuf=nbr*nas
      End if
*
* Buffer needed in the final step in the calculation of the exchange type
* integrals. The factor ten is arbitrary, however a 10 writes to disk per
* batch is acceptable before paging out integrals to disk to increase
* memory, however we dont want to write in to small batches 512 dWords (4k) is
* the lower limit
*
* /*
*    __________________________________
*    |          |           |         |
*    |          |           |         |
*    ----------------------------------
*    |          |           |         |
*    |         -|-----------|---------|
*    ---------------------------------|\
*    |        | |           |         | \
*    | ORDINT | |  (IjkL)   | (Ij|Kl) |  \
*    |----------|           |         |   \ N
*    ----------------------------------   /
*    |          |           |         |  /
*    | (ij|KL)  |  (iJ|kL)  | (iJ|Kl) | /
*    |          |           |         |/
*    ----------------------------------\
*    |          |           |         | \
*    | (ij|Kl)  |  (ij|Kl)  | (ij|Kl) |  \
*    |          |           |         |   \
*    ----------------------------------    nBuf*iMax
*    |          |           |         |   /
*    | (ij|kL)  |  (ij|kL)  | (ij|kL) |  /
*    |          |           |         | /
*    ----------------------------------/
* */
      i1=nAS*(nAP*nBQ*nBR+nAQ*nBP*nBR)
      If (ISP.eq.ISQ) i1=nAS*nAP*nBQ*nBR
      i2=nAR*(nAP*nBQ*nBS+nAQ*nBP*nBS)
      If (ISP.eq.ISQ) i2=nAR*nAP*nBQ*nBS
*     Size of output buffer if all data should fit in memory
      NOUT=max(i1,i2,NARS*nBPQ)
*     NOUT=max(nAS*(nAP*nBQ*nBR+nAQ*nBP*nBR),
*    &      nAR*(nAP*nBQ*nBS+nAQ*nBP*nBS),
*    &      NAR*nAS*nBQ*nBP)
*
      If (nBuf*iMax+NOUT.le.MemXX) Then
         Incore=0 ! We can do the calculation totally in core
      Else
*       But we dont have enough memory
        Incore=1
        MemXXX=MemXX-nOut
*       So we give each buffer same amount of memory
        iMax=MemXXX/nBuf
*       But we need space to unpack the integrals
        iMax=(MemXXX-iMax)/nBuf
        iMax=(MemXXX-iMax)/nBuf
        iMax=(MemXXX-iMax)/nBuf
        If (iMax.lt.1) Then
           Write (6,*) 'TraMO_MCLR: iMax.lt.1'
           Write (6,*) 'iMax=',iMax
           Call Abend()
        End If
        Call DAName(LUHLF2,FNHLF2)
        Call DAName(LUHLF3,FNHLF3)
      End If
*
*****************************************************************
*
* OK to start with with we need buffers but in the final step
* the only memory we need is a place to load (pq) integrals plus
* unpacking area plus a area of the size n to put the final result
* if we do not have enough memory, we will have a very long nose
*
*****************************************************************
*
      ip2=1                  ! (ij|kL)
      ip3=ip2+iMax*nBR*nAS+1   ! (ij|Kl)
      Buffer(ip3-1)=-99999.0d0
      If (iss.ne.isr) Then
      ip4=ip3+iMax*nBS*nAR+1   ! (ij|KL)
      else
      ip4=ip3+1
      End if
      Buffer(ip4-1)=-99999.0d0
      ipB=ip4+nBPQ*nARS+1      ! Read from ORDINT
      Buffer(ipB-1)=-99999.0d0
*
*****************************************************************
*
*     Start loop over sorted AO-integrals: npq pq-pairs in each buffer
*
*****************************************************************
*
      IPQ=0
      LPQ=0
      NPQ=0
      IRSST=1-NBRS    ! startpoint for present batch of integrals
      iopt=1
*
      If (.not.NOFILE) Then
      DO 10 NP=1,NBP
        NQM=NBQ
        IF(ISP.EQ.ISQ) NQM=NP
        DO 9 NQ=1,NQM
          IPQ=IPQ+1
          IPQUT2=IPQUT2+1
          IPQUT3=IPQUT3+1
          IPQUT4=IPQUT4+1
*
*****************************************************************
*                                                               *
*         Read in a block of integrals  NPQ pq values           *
*                                                               *
*****************************************************************
*
          IF(LPQ.EQ.NPQ.and.(.not.NOFILE)) THEN
            CALL RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,
     &                 Buffer(ipB),LBUF,NPQ)
            Call GADSum(Buffer(ipB),LBuf)
            If (irc.ne.0) Then
               Write (6,*) 'TraMO_MCLR: error reading ORDINT!'
               Call Abend()
            End If
            IOPT=2
            LPQ=0
            IRSST=1-NBRS
          ENDIF
          LPQ=LPQ+1
          IRSST=IRSST+NBRS
*
*****************************************************************
*                                                               *
*        Start transformation of this pq pair -> (ij|kL) (ijKl) *
*                                                               *
*****************************************************************
*
          If ((nAS*nAR.ne.0).or.(nAS*(nAP+nAQ).ne.0)) Then
            If (iSR.eq.iSS) Then
              CALL SQUARE(Buffer(ipB+IRSST-1),X1,1,NBS,NBS)
              CALL DGEMM_('T','N',
     &                    NBR,NAS,NBS,
     &                    1.0d0,X1,NBS,
     &                    CMS,NBS,
     &                    0.0d0,X2,NBR)
            Else
              CALL DGEMM_('T','N',
     &                    NBR,NAS,NBS,
     &                    1.0d0,Buffer(ipB+IRSST-1),NBS,
     &                    CMS,NBS,
     &                    0.0d0,X2,NBR)
            End If
          End If
*
*****************************************************************
*
*     Sort and save integrals transformed in one index
*     that should be used for exchange type integrals
*
*****************************************************************
*
          If ((nAP+nAQ)*nAS.ne.0) Then
             If (IPQUT2.GT.iMax) Then
               iPQUT2=1
               iST2=1
               Do i=1,nAS*nBR
                  iBuf2=iBuf2+1
                  IF ( iBuf2.GT.LIOTAB ) THEN
                     Write (6,*) 'TraMO_MCLR: iBuf2.GT.LIOTAB'
                     Write (6,*) 'iBuf2,LIOTAB=',iBuf2,LIOTAB
                     Call Abend()
                  ENDIF
                  CALL PKR8(0,iMax,NBYTES,Buffer(ip2+ist2-1),
     &                                    Buffer(ip2+ist2-1))
                  LPKREC=(NBYTES+itob-1)/itob
*                 Save the length of this record
                  IRLHLF2(iBUF2)=LPKREC
*                 Save the address of this record
                  IDAHLF2(iBUF2)=iAD2
                  Call C_F_Pointer(C_Loc(Buffer(ip2+IST2-1)),iBuffer,
     &                             [LPKREC])
                  CALL iDAFILE(LUHLF2,1,iBuffer,LPKREC,IAD2)
                  Nullify(iBuffer)
                  iST2=iST2+iMax
               End Do
             End If
*
             call dcopy_(nAs*nBR,X2,1,Buffer(ip2+iPQUT2-1),iMax)
*
          End If
*****************************************************************
*
          IF(ISR.NE.ISS) THEN
            CALL DGEMM_('N','N',
     &                  NBS,NAR,NBR,
     &                  1.0d0,Buffer(ipB+IRSST-1),NBS,
     &                  CMR,NBR,
     &                  0.0d0,X3,NBS)
*
*****************************************************************
          If ((nAP+nAQ)*nAR.ne.0) Then
            If(IPQUT3.GT.iMax) Then
              IPQUT3=1
              IST3=1
              DO 34 I=1,nAR*nBS
                IBUF3=IBUF3+1
                IF ( IBUF3.GT.LIOTAB ) THEN
                   Write (6,*) 'TraMO_MCLR: iBuf3.GT.LIOTAB'
                   Write (6,*) 'iBuf3,LIOTAB=',iBuf3,LIOTAB
                   Call Abend()
                ENDIF
                CALL PKR8(0,IMAX,NBYTES,Buffer(ip3+IST3-1),
     &                                    Buffer(ip3+IST3-1))
                LPKREC=(NBYTES+itob-1)/itob
                IRLHLF3(IBUF3)=LPKREC
                IDAHLF3(IBUF3)=IAD3
                Call C_F_Pointer(C_Loc(Buffer(ip3+IST3-1)),iBuffer,
     &                           [LPKREC])
                CALL iDAFILE(LUHLF3,1,iBuffer,LPKREC,IAD3)
                Nullify(iBuffer)
                IST3=IST3+iMax
34             CONTINUE
            End If
          End If
*
          call dcopy_(nAR*nBS,X3,1,Buffer(ip3+iPQUT3-1),iMax)
*
          End If
*
*****************************************************************
*
*  Transform to Coulomb type integrals
*
*****************************************************************
*
          If (NAR*NAS.ne.0) Then
          IF(ISR.EQ.ISS) THEN
            CALL MXMT(X2,        NBR,1,
     *                CMR, 1,NBR,
     *                X4,
     *                NAR,NBR)
          ELSE
            CALL DGEMM_('T','N',
     &                  NAS,NAR,NBR,
     &                  1.0d0,X2,NBR,
     &                  CMR,NBR,
     &                  0.0d0,X4,NAS)
          End If
          End If
*
*****************************************************************
*
          If (nARS.ne.0) Then
            call dcopy_(nARS,X4,1,Buffer(ip4+iPQUT4-1),nBPQ)
          End If
*
*****************************************************************
*
9       CONTINUE
10    CONTINUE
*
*
*****************************************************************
*
*     Empty last buffers
*
*****************************************************************
*
      IF(INCORE.EQ.1) THEN
        IST=1
        DO 61 I=1,nBR*nAS
          IBUF2=IBUF2+1
          IF ( IBUF2.GT.LIOTAB ) THEN
             Write (6,*) 'TraMO_MCLR: iBuf2.GT.LIOTAB'
             Write (6,*) 'iBuf2,LIOTAB=',iBuf2,LIOTAB
             Call Abend()
          ENDIF
          CALL PKR8(0,iMax,NBYTES,Buffer(ip2+IST-1),Buffer(ip2+IST-1))
          LPKREC=(NBYTES+itob-1)/ItoB
          IDAHLF2(IBUF2)=IAD2
          IRLHLF2(IBUF2)=LPKREC
          Call C_F_Pointer(C_Loc(Buffer(ip2+IST-1)),iBuffer,[LPKREC])
          CALL iDAFILE(LUHLF2,1,iBuffer,LPKREC,IAD2)
          Nullify(iBuffer)
          IST=IST+IMAX
 61     CONTINUE
        If (iSS.ne.iSR) Then
         IST=1
         DO 71 I=1,nBS*nAR
          IBUF3=IBUF3+1
          IF ( IBUF3.GT.LIOTAB ) THEN
             Write (6,*) 'TraMO_MCLR: iBuf3.GT.LIOTAB'
             Write (6,*) 'iBuf3,LIOTAB=',iBuf3,LIOTAB
             Call Abend()
          ENDIF
          CALL PKR8(0,iMax,NBYTES,Buffer(ip3+IST-1),Buffer(ip3+IST-1))
          LPKREC=(NBYTES+itob-1)/itob
          IDAHLF3(IBUF3)=IAD3
          IRLHLF3(IBUF3)=LPKREC
          Call C_F_Pointer(C_Loc(Buffer(ip3+IST-1)),iBuffer,[LPKREC])
          CALL iDAFILE(LUHLF3,1,iBuffer,LPKREC,IAD3)
          Nullify(iBuffer)
          IST=IST+IMAX
71       CONTINUE
        End If
      End If
      End If ! nofile
      If (NOFILE)  ipqut4=NBPQ
      If (nars.ne.0)
     &  Call dDafile(LUTRI1,iOne,Buffer(ip4),nARS*NBPQ,iAD34)
      If (Buffer(ip3-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: Buffer(ip3-1).ne.-99999.0d0'
         Write (6,*) 'Buffer(ip3-1)=',Buffer(ip3-1)
         Call Abend()
      End If
      If (Buffer(ip4-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: Buffer(ip4-1).ne.-99999.0d0'
         Write (6,*) 'Buffer(ip4-1)=',Buffer(ip4-1)
         Call Abend()
      End If
      If (Buffer(ipB-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: Buffer(ipB-1).ne.-99999.0d0'
         Write (6,*) 'Buffer(ipB-1)=',Buffer(ipB-1)
         Call Abend()
      End If
      If (nAP+nAQ.eq.0) Goto 987
*****************************************************************
*
*     second transformation
*     =====================
*
*    The integrals are transformed and written to disk
*    this is done so that the integrals are written to
*    disk in batches of nB*nB*nA integrals. This is the
*    Least number of integrals that have to match in buffer
*    because we want to write the integrals sorted onto disk
*
*****************************************************************
*
      ipX=ip4
      ipY=ipX+nBR*NAS*nBP*nAQ+1   ! (Ij|kL)
      ipz=ipy+1
      If (isp.ne.isq) ipZ=ipY+nBR*NAS*nBQ*nAP+1
      If (ipz.gt.memx) Then
         Write (6,*) 'TraMO_MCLR: ipz.gt.memx'
         Write (6,*) 'ipz,memx=',ipz,memx
         Call Abend()
      End If
      Buffer(ipX-1)=-99999.0d0
      Buffer(ipY-1)=-99999.0d0
      Buffer(ipZ-1)=-99999.0d0

      IAD2=0
      IAD3=0
      IVX=0
      DO 125 NS=1,nAS
        nr2=(NS-1)*NBR*NBP*NAQ
        nr3=(NS-1)*nBR*NBQ*NAP
        If (.not.NoFile) Then
        DO 120 NR=1,nBR
         IPQST=1+NBPQ*IVX
         IVX=IVX+1
         IBUF=IVX
*
*****************************************************************
*
*      Read one buffer of integrals back into core if INCORE=1
*
*****************************************************************
*
         IF(INCORE.EQ.1) THEN
           ipq=0
           ip2=1
           ip5=ip2+((nBPQ-1)/iMax+1)*iMax
           ipX=ip5+iMax+1
           ipY=ipX+nBP*nAQ*nBR*NAS+1
           ipZ=ipY+1
           If (isp.ne.isq) ipZ=ipX+nAQ*NBP*nBR*NAS+1
           Buffer(ipX-1)=-99999.0d0
           Buffer(ipY-1)=-99999.0d0
           Buffer(ipZ-1)=-99999.0d0
           ipq=1
112        CONTINUE
           IAD2=IDAHLF2(IBUF)
           LPKREC=IRLHLF2(IBUF)
           Call C_F_Pointer(C_Loc(Buffer(ip5)),iBuffer,[LPKREC])
           CALL iDAFILE(LUHLF2,2,iBuffer,LPKREC,IAD2)
           Nullify(iBuffer)
           CALL UPKR8(0,iMax,NBYTES,Buffer(ip5),Buffer(ip2+IPQ-1))
           IPQ=IPQ+iMax
           iBuf=iBuf+nAS*nBR
           IF(IPQ.LE.NBPQ) GO TO 112
           IPQST=1
         ENDIF
*
*****************************************************************
*                                                               *
* Transform (ij|kL) -> (Ij|kL) & (iJ|kL)                        *
*                                                               *
*****************************************************************
*
         IF(ISP.EQ.ISQ) THEN
         If (NAQ.ne.0) Then
            CALL SQUARE(Buffer(ip2+IPQST-1), X1,1,NBP, NBQ)
            CALL DGEMM_('T','N',
     &                  NBP,NAQ,NBQ,
     &                  1.0d0,X1,NBQ,
     &                  CMQ,NBQ,
     &                  0.0d0,X2,NBP)
         ENDIF
         ELSE
         If (NAQ.ne.0) Then
            CALL DGEMM_('T','N',
     &                  NBP,NAQ,NBQ,
     &                  1.0d0,Buffer(ip2+IPQST-1),NBQ,
     &                  CMQ,NBQ,
     &                  0.0d0,X2,NBP)
         ENDIF
         If (NAP.ne.0) Then
            CALL DGEMM_('N','N',
     &                  NBQ,NAP,NBP,
     &                  1.0d0,Buffer(ip2+IPQST-1),NBQ,
     &                  CMP,NBP,
     &                  0.0d0,X3,NBQ)
         ENDIF
         ENDIF
*
*****************************************************************
*
*      Move integrals to output buffer
*
*****************************************************************
*

         Do i=0,nAQ-1
           call dcopy_(nBP,X2(i*nBP+1),1,
     &                Buffer(ipX+nr2+i*nBR*nBP),1)
         End Do
         nr2=nr2+nBP

         If (isp.ne.isq) Then
          Do i=0,nAP-1
           call dcopy_(nBQ,X3(i*nBQ+1),1,
     &                Buffer(ipY+nr3+i*nBR*nBQ),1)
          End Do
          nr3=nr3+nBQ
         End If
*
120    CONTINUE
*
*****************************************************************
*
*       Empty last buffers
*
*****************************************************************
*
        End If ! nofile
125    CONTINUE
        If (nAQ.ne.0) Then
*          Call GADSum(Buffer(ipX),nAQ*NAS*NBR*NBP)
           Call dDafile(LUTRI2,ione,Buffer(ipX),
     &                 nAQ*NAS*NBR*NBP,iAD24)
        End If
        If (iSP.ne.iSQ.and.nAP.ne.0) Then
*          Call GADSum(Buffer(ipY),NAP*NAS*NBR*NBQ)
           Call dDafile(LUTRI3,ione,Buffer(ipY),
     &              NAP*NAS*NBR*NBQ,iAD14)
        End If
*

*
*****************************************************************
*
*
*****************************************************************
      If (buffer(ipX-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipX-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipX-1)=',buffer(ipX-1)
         Call Abend()
      End If
      If (buffer(ipY-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipY-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipY-1)=',buffer(ipY-1)
         Call Abend()
      End If
      If (buffer(ipZ-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipZ-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipZ-1)=',buffer(ipZ-1)
         Call Abend()
      End If
      If (iSS.ne.iSR) Then
       ipX=ip4
       Buffer(ipX-1)=-99999.0d0
       ipY=ipX+nAR*NBS*nBP*nAQ+1   ! (Ij|kL)
       Buffer(ipY-1)=-99999.0d0
       ipZ=ipy+1
       if (isp.ne.isq) ipZ=ipY+nAR*NBS*nBQ*nAP+1
       Buffer(ipZ-1)=-99999.0d0
       IVX=0
       DO 225 NR=1,nAR
        nr2=NBP*NBS*(nR-1)*NAQ
        nr3=NBQ*NBS*(nR-1)*NAP
        If (.not.NoFile) Then
        DO 220 NS=1,nBS
         IPQST=1+NBPQ*IVX
         IVX=IVX+1
         IBUF=IVX
*
*****************************************************************
*
*      Read one buffer of integrals back into core if INCORE=1
*
*****************************************************************
*
         IF(INCORE.EQ.1) THEN
           ip2=1
           ip3=ip2
           ip5=ip3+((nBPQ-1)/iMax+1)*iMax
           ipX=ip5+iMax+1
           ipY=ipX+nBP*nAQ*nAR*NBS+1
           ipZ=ipY+nBQ*nAP*NAR*NBS+1
           Buffer(ipX-1)=-99999.0d0
           Buffer(ipY-1)=-99999.0d0
           Buffer(ipZ-1)=-99999.0d0
           IPQ=1
212        CONTINUE
           IAD3=IDAHLF3(IBUF)
           LPKREC=IRLHLF3(IBUF)
           Call C_F_Pointer(C_Loc(Buffer(ip5)),iBuffer,[LPKREC])
           CALL iDAFILE(LUHLF3,2,iBuffer,LPKREC,IAD3)
           Nullify(iBuffer)
           CALL UPKR8(0,iMax,NBYTES,Buffer(ip5),Buffer(ip3+IPQ-1))
           IPQ=IPQ+iMax
           iBuf=iBuf+nAR*nBS
           IF(IPQ.LE.NBPQ) GO TO 212
           IPQST=1
         ENDIF
*
*****************************************************************
*                                                               *
* Transform (ij|kL) -> (Ij|kL) & (iJ|kL)                        *
*                                                               *
*****************************************************************
*
         IF(ISP.EQ.ISQ) THEN
         If (nBQ*nAP.ne.0) Then

            CALL SQUARE(Buffer(ip3+IPQST-1), X1,1,NBQ, NBQ)
            CALL DGEMM_('T','N',
     &                  NBP,NAQ,NBQ,
     &                  1.0d0,X1,NBQ,
     &                  CMQ,NBQ,
     &                  0.0d0,X2,NBP)
*           CALL xxDGEMUL(X1,NBQ,'N',CMP,NBP,'N',
*    &                  X3,NBQ,NBQ,NBP,NAP)
*           Call RecPrt('PqRs',' ',X3,nBQ,nAP)
         ENDIF
         ELSE
         If (nBP*nAQ.ne.0)
     &      CALL DGEMM_('T','N',
     &                  NBP,NAQ,NBQ,
     &                  1.0d0,Buffer(ip3+IPQST-1),NBQ,
     &                  CMQ,NBQ,
     &                  0.0d0,X2,NBP)
         If (nBQ*nAP.ne.0)
     &      CALL DGEMM_('N','N',
     &                  NBQ,NAP,NBP,
     &                  1.0d0,Buffer(ip3+IPQST-1),NBQ,
     &                  CMP,NBP,
     &                  0.0d0,X3,NBQ)
         ENDIF
*
*****************************************************************
*
*      Move integrals to output buffer
*
*****************************************************************
*

         Do i=0,nAQ-1
           call dcopy_(nBP,X2(i*nBP+1),1,
     &                Buffer(ipX+i*nBS*nBP+nr2),1)
         End Do
         nr2=nr2+nBP
         If (iSP.ne.iSQ) Then
          Do i=0,nAP-1
           call dcopy_(nBQ,X3(i*nBQ+1),1,
     &                Buffer(ipY+i*NBS*nBQ+nR3),1)
          End Do
          nr3=nr3+nBQ
         End If
*
220    CONTINUE
      If (buffer(ipX-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipX-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipX-1)=',buffer(ipX-1)
         Call Abend()
      End If
      If (buffer(ipY-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipY-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipY-1)=',buffer(ipY-1)
         Call Abend()
      End If
      If (buffer(ipZ-1).ne.-99999.0d0) Then
         Write (6,*) 'TraMO_MCLR: buffer(ipZ-1).ne.-99999.0d0'
         Write (6,*) 'buffer(ipZ-1)=',buffer(ipZ-1)
         Call Abend()
      End If

*****************************************************************
*
*      Write buffer to disk (iJ|Kl) & (Ij|Kl)
*
*****************************************************************
*
        End If ! nofile
225     CONTINUE
        If (nAQ.ne.0) Then
*          Call GADSum(Buffer(ipX),nAR*nAQ*nBS*nBP)
           Call dDaFile(LuTri4,ione,Buffer(ipX),
     &                  nAR*nAQ*nBS*nBP,iad23)
        End If
        If (iSP.ne.iSQ.and.nap.ne.0) Then
*          Call GADSum(Buffer(ipY),nAR*nAP*nBS*nBQ)
           Call dDaFile(LuTri5,ione,Buffer(ipY),
     &                  nAR*nAP*nBS*nBQ,iad13)
        End If
       End If
*
 987  Continue
      IF ( INCORE.EQ.1 ) THEN
        CALL DACLOS(LUHLF2)
        CALL DACLOS(LUHLF3)
      END IF
*
      RETURN
      End Subroutine TRAMO_MCLR_INTERNAL
*
      END
