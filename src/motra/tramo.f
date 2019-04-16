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
      SUBROUTINE TRAMO(LBUF,OUTBUF,nOUTBUF,
     &                 X1,nX1,X2,nX2,X3,nX3,VXPQ,nVXPQ,CMO,iDsk,mOVX)
*
*     Purpose: two-electron transformation routine.
*              Transformation is made in core if all half-transformed
*              integrals fit. Otherwise sorted integrals
*              are written onto unit LUHALF.
*
      !> module dependencies
#ifdef _HDF5_QCM_
      use hdf5_utils
#endif

      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
#include "SysDef.fh"
*
      Real*8 OUTBUF(nOUTBUF),X1(nX1),X2(nX2),X3(nX3),VXPQ(nVXPQ)
      Real*8 CMO(*)
*
      Integer iDsk(3,mOVX)
#ifdef _HDF5_QCM_
      integer, save :: total_number_2ints = 0
      real*8, allocatable :: tmpbuf(:)
      integer             :: iout_total
#endif
*
      Call qEnter('Tramo')
*
*     Set some constants
*
      Call ICopy(mOVX,[-1],0,iDsk(1,1),3)
      Call ICopy(mOVX, [0],0,iDsk(2,1),3)
      Call ICopy(mOVX, [1],0,iDsk(3,1),3)

#ifdef _HDF5_QCM_
      if(ihdf5 == 1)then
!         allocate(tmpbuf(nx1));
! Stefan's dirty hack to increase the # of BF working in MOTRA with HDF5
        allocate(tmpbuf(((nop+1)*(noq+1)*(nor+1)*(nos+1))));tmpbuf = 0
        write(6,*) 'size of tmpbuf:',
     &  ((nop+1)*(noq+1)*(nor+1)*(nos+1))
        iout_total = 0
      end if
#endif
      KBUF1=1
      LTUVX=0
      IOUT=0
      IAD14=0
      IPQUT=0
      IPQMAX=NBPQ
      INCORE=0
      IF ( NBPQ*NOVX.GT.nVXPQ ) THEN
        INCORE=1
        IPQMAX=nVXPQ/NOVX
        IPQMAX=(nVXPQ-IPQMAX)/NOVX
        IPQMAX=(nVXPQ-IPQMAX)/NOVX
        IPQMAX=(nVXPQ-IPQMAX)/NOVX
        WRITE(6,'(6X,A)')'OUT OF CORE TRANSFORMATION'
        CALL DANAME_MF(LUHALF,FNHALF)
      ENDIF
*
*     Start loop over sorted AO-integrals: npq pq-pairs in each buffer
*
      IF(IPRINT.GE.20) WRITE(6,*) 'TRAMO 001'
      IRC=0
      IOPT=1
      IPQ=0
      LPQ=0
      NPQ=0
      IRSST=1-NBRS
      Do NP=1,NBP
        NQM=NBQ
        IF(ISP.EQ.ISQ) NQM=NP
        Do NQ=1,NQM
          IPQ=IPQ+1
          IPQUT=IPQUT+1
*
*         Read in a block of integrals  NPQ pq values
*
          IF(LPQ.EQ.NPQ) THEN
            CALL RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
            Call GADSum(X1,LBUF)
            IOPT=2
            LPQ=0
            IRSST=1-NBRS
          ENDIF
          LPQ=LPQ+1
          IRSST=IRSST+NBRS
*
*         Start transformation of this pq pair
*
          IF(ISR.EQ.ISS) THEN
            CALL SQUARE(X1(IRSST),X2,1,NBS,NBS)
            IF( IPRINT.GE.30 ) THEN
              WRITE(6,'(6X,A,2I4)')'AO Integrals for the pair',NP,NQ
              CALL SQPRT(X2,NBR)
            ENDIF
            IF ( NBR*NBS*NOS.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NBR,NOS,NBS,
     &                  1.0d0,X2,NBS,
     &                  CMO(LMOS),NBS,
     &                  0.0d0,X3,NBR)
            CALL MXMT(X3,        NBR,1,
     &                CMO(LMOR), 1,NBR,
     &                X2,
     &                NOR,NBR)
          ELSE
            IF ( NBR*NBS*NOS.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NBR,NOS,NBS,
     &                  1.0d0,X1(IRSST),NBS,
     &                  CMO(LMOS),NBS,
     &                  0.0d0,X3,NBR)
            IF ( NOS*NBR*NOR.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NOS,NOR,NBR,
     &                  1.0d0,X3,NBR,
     &                  CMO(LMOR),NBR,
     &                  0.0d0,X2,NOS)
          ENDIF
*
*         Sort the matrix X2 into VXPQ (sort after PQ instead of VX)
*         if INCORE=1 also write sorted integrals on unit LUHALF.
*
          IF(IPQUT.GT.IPQMAX) THEN
            IPQUT=1
            IST=1
            Do I=1,NOVX
              CALL PKR8(0,IPQMAX,NBYTES,VXPQ(IST),VXPQ(IST))
              LPKREC=(NBYTES+RtoB-1)/RtoB
*
              IAD14_=IAD14
              CALL dDAFILE(LUHALF,1,VXPQ(IST),LPKREC,IAD14)
              CALL iDAFILE(LUHALF,1,iDsk(1,I),2,IAD14)
              iDsk(1,I)=IAD14_
              iDsk(2,I)=LPKREC
              iDsk(3,I)=iDsk(3,I)+IPQMAX
*
              IST=IST+IPQMAX
            End Do
          END IF
*
          IPQX=IPQUT
          Do I=1,NOVX
            VXPQ(IPQX)=X2(I)
            IPQX=IPQX+IPQMAX
          End Do
        End Do !  NQ
      End Do   !  NP
*
      LVXPQ=NBPQ*NOVX
      IF(IPRINT.GE.20) WRITE(6,*) 'TRAMO 002'
      IF(IPRINT.GE.30) THEN
         WRITE(6,'(A,I6)') ' HALF TRANSFORMED INTEGRALS:',LVXPQ
         WRITE(6,'(1X,10F11.6)') (VXPQ(I),I=1,LVXPQ)
      END IF
*
*     Empty last buffers
*
      IF(INCORE.EQ.1) THEN
        IST=1
        DO I=1,NOVX
          CALL PKR8(0,IPQMAX,NBYTES,VXPQ(IST),VXPQ(IST))
          LPKREC=(NBYTES+RtoB-1)/RtoB
*
          IAD14_=IAD14
          CALL dDAFILE(LUHALF,1,VXPQ(IST),LPKREC,IAD14)
          CALL iDAFILE(LUHALF,1,iDsk(1,I),2,IAD14)
          iDsk(1,I)=IAD14_
          iDsk(2,I)=LPKREC
          iDsk(3,I)=iDsk(3,I)+IPQMAX
*
          IST=IST+IPQMAX
        End Do
        IAD14=0
      ENDIF
*
*     First half transformation is now done.VXPQ contains half trans-
*     formed integrals for this symmetry block, if INCORE=0
*     otherwise integrals or one VX pair at a time will be read
*     from unit LUHALF.
*
      IF(IPRINT.GE.20) WRITE(6,*)'TRAMO 004'
      IVX=0
      DO 20 NV=1,NOR
        NXM=NV
        IF(ISS.NE.ISR) NXM=NOS
        DO 25 NX=1,NXM
         IPQST=1+NBPQ*IVX
         IVX=IVX+1
*
*      Read one buffer of integrals back into core if INCORE=1
*
         IF(INCORE.EQ.1) THEN
*
*----------Back chain buffers
*
           inBuf=iDsk(3,iVX)
           IPQ = inBuf-IPQMAX
12         IAD14  = iDsk(1,iVX)  ! Disk Address of the next buffer
           If (IAD14.ge.0) Then
              LPKREC = iDsk(2,iVX)  ! Length of the  buffer
              CALL dDAFILE(LUHALF,2,VXPQ(inBuf),LPKREC,IAD14)
              CALL iDAFILE(LUHALF,2,iDsk(1,iVX),2,IAD14)
              CALL UPKR8(0,IPQMAX,NBYTES,VXPQ(inBuf),VXPQ(IPQ))
              IPQ=IPQ-IPQMAX
              GO TO 12
           End If
*
           IPQST=1
         END IF
         IF(ISP.NE.ISR) THEN
           IF(ISP.EQ.ISQ) THEN
             CALL SQUARE(VXPQ(IPQST), X2,1,NBQ, NBQ)
            IF( NBP*NBQ*NOQ.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NBP,NOQ,NBQ,
     &                  1.0d0,X2,NBQ,
     &                  CMO(LMOQ),NBQ,
     &                  0.0d0,X1,NBP)
             CALL MXMT(X1,        NBP,1,
     &                 CMO(LMOP), 1,NBP,
     &                 X2,
     &                 NOP,NBP)
             IX2=(NOP+NOP**2)/2
           ELSE
            IF ( NBP*NBQ*NOQ.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NBP,NOQ,NBQ,
     &                  1.0d0,VXPQ(IPQST),NBQ,
     &                  CMO(LMOQ),NBQ,
     &                  0.0d0,X1,NBP)
            IF ( NOQ*NBP*NOP.GT.0 )
     &      CALL DGEMM_('T','N',
     &                  NOQ,NOP,NBP,
     &                  1.0d0,X1,NBP,
     &                  CMO(LMOP),NBP,
     &                  0.0d0,X2,NOQ)
             IX2=NOP*NOQ
           ENDIF
         ELSE
           IF(ISP.EQ.ISQ) THEN
             CALL SQUARE(VXPQ(IPQST),X2,1,NBP,NBP)
             IF ( NBP*NBQ*NOQ.GT.0 )
     &       CALL DGEMM_('T','N',
     &                   NBP,NOQ,NBQ,
     &                   1.0d0,X2,NBQ,
     &                   CMO(LMOQ),NBQ,
     &                   0.0d0,X1,NBP)
           ELSE
             IF ( NBP*NBQ*NOQ.GT.0 )
     &       CALL DGEMM_('T','N',
     &                   NBP,NOQ,NBQ,
     &                   1.0d0,VXPQ(IPQST),NBQ,
     &                   CMO(LMOQ),NBQ,
     &                   0.0d0,X1,NBP)
           END IF
*
*          X1 now contains the matrix (pu/vx) for all p and all u
*          for a fixed pair vx.
*          If the first index t=v then the range of u is x,umax
*          Otherwise the range is 1,umax. a loop over t is necessary her
*
           ISTMOT=LMOP+NBP*(NV-1)
           IX2=1
           DO 14 NT=NV,NOP
             NUMAX=NOQ
             IF(ISP.EQ.ISQ) NUMAX=NT
             LOQ=NUMAX
             IF(NT.EQ.NV) LOQ=NUMAX-NX+1
             IX1=1
             IF(NT.EQ.NV) IX1=1+NBP*(NX-1)
             IF ( LOQ.GT.0 ) THEN
               IF ( NBP.EQ.0 ) THEN
                 CALL DCOPY_(LOQ,[0.0D0],0,X2(IX2),1)
               ELSE
                 CALL DGEMM_('T','N',
     &                       LOQ,1,NBP,
     &                       1.0d0,X1(IX1),NBP,
     &                       CMO(ISTMOT),NBP,
     &                       0.0d0,X2(IX2),LOQ)
               END IF
             END IF
             IX2=IX2+LOQ
             ISTMOT=ISTMOT+NBP
14         CONTINUE
           IX2=IX2-1
         ENDIF
*
*      Move integrals to output buffer and write them on LUTWOMO
*
         DO 16 NTUVX=1,IX2
           IF(IOUT.EQ.KBUF) THEN
             IOUT=0
             IF (IPRINT.GE.5) THEN
                WRITE(6,'(A)')
     &               'SAVE TRANSFORMED INTEGRALS:'
                WRITE(6,'(A,3I8)')
     &               'LU,KBUF,IDISK =',LUTWOMO,KBUF,IAD13
             END IF
             CALL dDAFILE(LUTWOMO,1,OUTBUF(KBUF1),KBUF,IAD13)
#ifdef _HDF5_QCM_
      !WRITE(6,*)
      !&'SAVE TRANSFORMED INTEGRALS:',
      !& KBUF,KBUF1,iout_total,iout_total+kbuf
             if(ihdf5 == 1)then
               call dcopy_(KBUF,OUTBUF(KBUF1),1,tmpbuf(iout_total+1),1)
               iout_total = iout_total + kbuf
             end if
#endif
             KBUF1=KBUF+2-KBUF1
           ENDIF
           IOUT=IOUT+1
           LTUVX=LTUVX+1
           OUTBUF(IOUT+KBUF1-1)=X2(NTUVX)
16       CONTINUE
25      CONTINUE
20    CONTINUE
*
*     Empty last buffer (which is never empty)
*
      IF (IPRINT.GE.5) THEN
         WRITE(6,'(A)')
     &        'SAVE TRANSFORMED INTEGRALS:'
         WRITE(6,'(A,3I8)')
     &        'LU,KBUF,IDISK =',LUTWOMO,KBUF,IAD13
      END IF
      CALL dDAFILE(LUTWOMO,1,OUTBUF(KBUF1),KBUF,IAD13)
      IF(IPRINT.GE.10) THEN
        WRITE(6,'(1X,A,4I2,A,I4)')
     *  'TRANSFORMED INTEGRALS FOR SYMMETRY BLOCK',ISP,ISQ,ISR,ISS,
     *  ' IOUT=',IOUT
        WRITE(6,'(1X,10F12.6)') (OUTBUF(I+KBUF1-1),I=1,IOUT)
      ENDIF

#ifdef _HDF5_QCM_
      if(ihdf5 == 1)then
        !> and the last buffer here as well...
        call dcopy_(iout,OUTBUF(KBUF1),1,tmpbuf(iout_total+1),1)
        iout_total = iout_total + iout
        !> put data to file
        write(datatag,'(a1,i3,a1,i3,a1,i3,a1,i3)')
     &  'p',isp,'q',isq,'r',isr,'s',iss
        tagx = 16
        datadim(1)        = 1; datadim_bound = 1
        call hdf5_put_data(file_id(1),"XXXXXX",datadim,iout_total)
        datadim(1)        = iout_total; datadim_bound = 1
        call hdf5_put_data(file_id(1),"XXXXXX",datadim,tmpbuf)

        IF(IPRINT.GE.5) THEN
          total_number_2ints = total_number_2ints + iout_total
          WRITE(6,'(1X,a,i10)') ' total number of integrals (sym)  ...',
     &    iout_total
          WRITE(6,'(1X,a,i10)') ' total number of integrals so far ...',
     &    total_number_2ints
        END IF

        deallocate(tmpbuf)
      end if
#endif
*
      IF ( INCORE.EQ.1 ) CALL DACLOS(LUHALF)
*
      Call qExit('Tramo')
      RETURN
      END
