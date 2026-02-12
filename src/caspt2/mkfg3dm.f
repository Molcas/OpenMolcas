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
C Procedure for computing 1-body, 2-body, and 3-body
C density elements with active indices only,
C and related matrices obtained from contractions with
C the diagonal one-electron Hamiltonian.
C
C Out: density matrices, denoted here G1, G2 and G3.
C Storage: G1 and G2 are simple two- and four-index arrays, and
C includes also such zeroes that are implied by symmetry.
C But G3 is quite large, and while it is stored with zeroes, it
C is made more compact by calculating only the minimum amount of unique values and storing the
C active indices in the array idxG3.  Later, the full matrix can be restored
C on the fly by using the full permutational symmetry.
C
C ----
C
C.NN15 NOTE:
C In DMRG-CASPT2, full G3 is computed from BLOCK code
C and F3 is approximated by cumulant reconstruction scheme from G1, G2, and G3 (cu4),
C to avoid computation of 4-particle density matrix.
C However, because this introduces instability of CASPT2 calculation
C (lots of negative denominators appear), relatively large IPEA and imaginary shifts
C are required to converge CASPT2 iteration.
C

#include "compiler_features.h"

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined _DMRG_
      SUBROUTINE MKFG3DM(IFF,G1,F1,G2,F2,G3,F3,idxG3,NLEV)
      use Symmetry_Info, only: Mul
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: DEBUG, VERBOSE
      use gugx, only: CIS, SGS, L2ACT
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      use Definitions, only: RtoB
      use caspt2_module, only: nActEl, nSym, STSym
#ifdef _DMRG_
      use caspt2_module, only: DMRG
#endif
      use gugx, only: MxLev
      use pt2_guga, only: MxCI, nG1, nG2, nG3
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IFF, NLEV
      REAL*8, INTENT(OUT) :: G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: G3(*), F3(*)
      INTEGER*1, INTENT(OUT) :: idxG3(6,*)

      LOGICAL RSV_TSK

      INTEGER, PARAMETER :: I1=KIND(idxG3)

      REAL*8 DG1,DG2,DG3,DF1,DF2,DF3

      INTEGER I,J,IDX,JDX
      INTEGER IB,IBMN,IBMX,IBUF,NB,NBTOT,IBUF1
      INTEGER IP1,IP2,IP3,IP1MN,IP1MX,IP1I,IP1STA,IP1END,IP3MX,IQ1
      INTEGER IG3,IG3OFF
      INTEGER ISTU,ISVX,ISYZ
      INTEGER IT,IU,IV,IX,IY,IZ
      INTEGER ITLEV,IULEV,IVLEV,IXLEV,IYLEV,IZLEV
      INTEGER NBUF1
      INTEGER IOFFSET
      INTEGER ISSG1,ISSG2,ISP1
      INTEGER ITASK,ISUBTASK,ID,NTASKS,NSUBTASKS,MXTASK,MYTASK,MYBUFFER
      INTEGER NTRI1,NTRI2
      INTEGER MEMMAX, MEMMAX_SAFE
      INTEGER NLEV2
#ifdef _ENABLE_BLOCK_DMRG_
      INTEGER NLEV4
      REAL*8, ALLOCATABLE:: G3Tmp(:)
#endif
      INTEGER NCI

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

      ! translation tables for levels i,j to and from pair indices idx
      INTEGER IJ2IDX(MXLEV,MXLEV)
      INTEGER IDX2IJ(2,MXLEV**2)
      INTEGER ICNJ(MXLEV**2)
      INTEGER IP1_BUF(MXLEV**2)

      REAL*8, ALLOCATABLE:: BUF1(:,:), BUF2(:), BUFT(:), BUFD(:)
      INTEGER, ALLOCATABLE:: TaskList(:,:)


      ! result buffer, maximum size is the largest possible ip1 range,
      ! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
*     REAL*8 BUFR(MXLEV**2)


C Put in zeroes. Recognize special cases:
      IF(nlev.EQ.0) GOTO 999

      CALL DCOPY_(NG1,[0.0D0],0,G1,1)
      CALL DCOPY_(NG2,[0.0D0],0,G2,1)
      CALL DCOPY_(NG3,[0.0D0],0,G3,1)
      IF(IFF.ne.0) THEN
        CALL DCOPY_(NG1,[0.0D0],0,F1,1)
        CALL DCOPY_(NG2,[0.0D0],0,F2,1)
        CALL DCOPY_(NG3,[0.0D0],0,F3,1)
      END IF

      IF(NACTEL.EQ.0) GOTO 999

      NCI=CIS%NCSF(STSYM)
* This should not happen, but...
      IF(NCI.EQ.0) GOTO 999

C Here, for regular CAS or RAS cases.

C Special pair index idx2ij allows true RAS cases to be handled:
      nlev2=nlev**2
      ntri1=(nlev2-nlev)/2
      ntri2=(nlev2+nlev)/2
      idx=0
      do i=1,nlev-1
        do j=i+1,nlev
          idx=idx+1
          ij2idx(i,j)=idx
          idx2ij(1,idx)=i
          idx2ij(2,idx)=j
          jdx=nlev2+1-idx
          ij2idx(j,i)=jdx
          idx2ij(1,jdx)=j
          idx2ij(2,jdx)=i
        end do
      end do
      do i=1,nlev
        idx=ntri1+i
        ij2idx(i,i)=idx
        idx2ij(1,idx)=i
        idx2ij(2,idx)=i
      end do
      do idx=1,nlev2
        i=idx2ij(1,idx)
        j=idx2ij(2,idx)
        jdx=ij2idx(j,i)
        icnj(idx)=jdx
      end do

      call mma_MaxDBLE(memmax)

* Use *almost* all remaining memory:
      memmax_safe=int(dble(memmax)*0.95D0)

* Buffers to compute CI expansion vectors into:
*
*
*
*
*
*
      nbuf1=max(1,min(nlev2,(memmax_safe-3*mxci)/mxci)) ! -> 1 w/ DMRG?
      CALL mma_allocate(BUF1,MXCI,NBUF1,LABEL='BUF1')
      CALL mma_allocate(BUF2,MXCI,LABEL='BUF2')
      CALL mma_allocate(BUFT,MXCI,LABEL='BUFT')
      CALL mma_allocate(BUFD,MXCI,LABEL='BUFD')

C-SVC20100301: calculate maximum number of tasks possible
      MXTASK=(NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
      CALL mma_allocate (TaskList,mxTask,4,LABEL='TaskList')

      IF(iPrGlb.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(2X,A)') 'Constructing G3/F3'
        WRITE(6,'(2X,A,F16.9,A)') ' memory avail: ',
     &    (memmax*RtoB)/1.0D9, ' GB'
        WRITE(6,'(2X,A,F16.9,A)') ' memory used:  ',
     &    (((nbuf1+3)*MXCI)*RtoB)/1.0D9, ' GB'
        call xFlush(6)
      ENDIF

      iG3OFF=0
* A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.
*
      DO issg1=1,nsym
       isp1=Mul(issg1,stsym)
*      nsgm1=CIS%ncsf(issg1)
*      CALL H0DIAG_CASPT2(ISSG1,BUFD,NOW1,IOW1,NMIDV)

C-SVC20100301: calculate number of larger tasks for this symmetry, this
C-is basically the number of buffers we fill with sigma1 vectors.
      iTask=1
      ibuf1=0
      DO ip1=1,nlev2
        itlev=idx2ij(1,ip1)
        iulev=idx2ij(2,ip1)
        istu=Mul(SGS%ism(itlev),SGS%ism(iulev))
        IF (istu.EQ.isp1) THEN
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1
          IF (ibuf1.EQ.1) TaskList(iTask,1)=ip1
        ENDIF
        IF (ibuf1.EQ.nbuf1.OR.(ibuf1.GT.0.AND.
     &         (ip1.EQ.ntri2.OR.ip1.EQ.nlev2))) THEN
            TaskList(iTask,2)=ip1_buf(ibuf1)
            TaskList(iTask,3)=ibuf1
            iTask=iTask+1
            ibuf1=0
        ENDIF
      ENDDO
      nTasks=iTask
      IF (ibuf1.EQ.0) nTasks=nTasks-1
C-SVC20100309: calculate number of inner loop iteration tasks.
      iOffSet=0
      DO iTask=1,nTasks
        TaskList(iTask,4)=iOffSet
        ip1sta=TaskList(iTask,1)
        ip1end=TaskList(iTask,2)
        ip3mx=ntri2
        if(ip1end.le.ntri2) ip3mx=ip1end
        if(ip1sta.gt.ntri2) ip3mx=ntri1
C-SVC20100309: Currently -we are going to limit this to the ip3-loop and
C-leave the ip2-loop intact.  This was based on the large overhead which
C-was observed for a very large number of small tasks.
C       iOffSet=iOffSet+ip3mx*ntri2-((ip3mx**2-ip3mx)/2)
        iOffSet=iOffSet+ip3mx
      ENDDO
      nSubTasks=iOffSet

      IF(iPrGlb.GE.VERBOSE) THEN
        WRITE(6,'(2X,A,I3,A,I6)') 'Sym: ',issg1,', #Tasks: ',nSubTasks
        call xFlush(6)
      ENDIF

      IF(iPrGlb.GE.DEBUG) THEN
        IF (nSubTasks .GT. 0) THEN
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "task ID ",
     &    " ip1 range  ",
     &    "ip3 ",
     &    "#elements"
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          call xFlush(6)
        END IF
      END IF

C-SVC20100301: initialize the series of subtasks
      Call Init_Tsk(ID, nSubTasks)

      myBuffer=0

 500  CONTINUE
C-SVC20100908: first check: can I actually do any task?
      IF ((NG3-iG3OFF).LT.nbuf1*ntri2) GOTO 501
C-SVC20100831: initialize counter for offset into G3
C-SVC20100302: BEGIN SEPARATE TASK EXECUTION
      If (.NOT.Rsv_Tsk(ID,iSubTask)) GOTO 501

      myTask=nTasks
      DO iTask=1,nTasks
        iBuf=iSubTask-TaskList(iTask,4)
        IF (iBuf.LE.0) THEN
          myTask=iTask-1
          goto 666
        ENDIF
      ENDDO
666   continue
      iTask=myTask

      iOffSet=TaskList(iTask,4)

C-SVC20100310: one task handles a range of ip1 values
C-that are in the buffer and one ip3 value, for which
C-a loop over ip2 values is then executed.
      ip1sta=TaskList(iTask,1)
      ip1end=TaskList(iTask,2)
      ip3=iSubTask-iOffSet

C-SVC20100301: fill the buffer with sigma vectors if they
C-have not been computed yet, else just get the number of
C-sigma vectors in the buffer.
      IF (myBuffer.NE.iTask) THEN
        ibuf1=0
        do ip1i=ip1sta,ip1end
         itlev=idx2ij(1,ip1i)
         iulev=idx2ij(2,ip1i)
         istu=Mul(SGS%ism(itlev),SGS%ism(iulev))
         it=L2ACT(itlev)
         iu=L2ACT(iulev)
         if(istu.eq.isp1) then
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1i
*         call dcopy_(nsgm1,[0.0D0],0,BUF1(:,ibuf1),1)
*         CALL SIGMA1(IULEV,ITLEV,1.0D00,STSYM,CI,BUF1(:,ibuf1))
         end if
        end do
        myBuffer=iTask
      ELSE
        ibuf1=TaskList(iTask,3)
      ENDIF
C-SVC20100301: necessary batch of sigma vectors is now in the buffer

*     ! The ip1 buffer could be the same on different processes
*     ! so only compute the G1 contribution when ip3 is 1, as
*     ! this will only be one task per buffer.
*     if (issg1.eq.stsym.AND.ip3.eq.1) then
*       do ib=1,ibuf1
*         idx=ip1_buf(ib)
*         itlev=idx2ij(1,idx)
*         iulev=idx2ij(2,idx)
*         it=L2ACT(itlev)
*         iu=L2ACT(iulev)
*         G1(it,iu)=DDOT_(nsgm1,ci,1,BUF1(:,ib),1)
*         IF(IFF.ne.0) then
*           F1sum=0.0D0
*           do i=1,nsgm1
*             F1sum=F1sum+CI(i)*BUF1(i,ib)*bufd(i)
*           end do
*           F1(it,iu)=F1sum-EPSA(iu)*G1(it,iu)
*         end if
*       end do
*     end if

C     ip3mx=ntri2
C     if(ip1end.le.ntri2) ip3mx=ip1end
C     if(ip1sta.gt.ntri2) ip3mx=ntri1
C-SVC20100309: loop over ip3, ip2
C     do ip3=1,ip3mx

C-SVC20100309: PAM's magic formula
*     iCnt=iSubTask-iOffSet
*     ip3=int(dble(ntri2)+1.5D0 -
*    &     sqrt((dble(ntri2)+0.5d0)**2-2*iCnt+0.000001D0))
*     ip2=iCnt-((ip3-1)*ntri2-((ip3-1)*(ip3-2))/2 )+ip3-1

C-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

C NN.14 TODO:
C To avoid storing full G3(tmp) in memory, need to store
C G3(:,:,it,iu,iy,iz) loaded from disk, for each process...

      iq1=icnj(ip3)
* The indices corresponding to pair index p3:
      iylev=idx2ij(1,ip3)
      izlev=idx2ij(2,ip3)
      isyz=Mul(SGS%ism(iylev),SGS%ism(izlev))
      issg2=Mul(isyz,stsym)
*     nsgm2=CIS%ncsf(issg2)
      iy=L2ACT(iylev)
      iz=L2ACT(izlev)
*     call dcopy_(nsgm2,0.0D0,0,BUF2,1)
*     CALL SIGMA1(IYLEV,IZLEV,1.0D00,STSYM,CI,BUF2)
*     if(issg2.eq.issg1) then
*       do ib=1,ibuf1
*         idx=ip1_buf(ib)
*         itlev=idx2ij(1,idx)
*         iulev=idx2ij(2,idx)
*         it=L2ACT(itlev)
*         iu=L2ACT(iulev)
*         G2(it,iu,iy,iz)=DDOT_(nsgm1,BUF2,1,BUF1(:,ib),1)
*         IF(IFF.ne.0) THEN
*           F2sum=0.0D0
*           do i=1,nlev
*             F2sum=F2sum+BUF2(i)*bufd(i)*BUF1(i,ib)
*           end do
*           F2(it,iu,iy,iz)=F2sum
*         END IF
*       end do
*     end if
      nbtot=0
      do ip2=ip3,ntri2
        ivlev=idx2ij(1,ip2)
        ixlev=idx2ij(2,ip2)
        isvx=Mul(SGS%ism(ivlev),SGS%ism(ixlev))
        iv=L2ACT(ivlev)
        ix=L2ACT(ixlev)
        if(isvx.ne.Mul(issg1,issg2)) goto 99
*       call dcopy_(nsgm1,[0.0D0],0,BUFT,1)
*       CALL SIGMA1(IVLEV,IXLEV,1.0D00,ISSG2,BUF2,BUFT)
*-----------
* Max and min values of index p1:
        ip1mx=ntri2
        if(ip3.le.ntri1) then
          ip1mx=nlev2
          if(ip2.gt.ntri1) ip1mx=iq1
        end if
        ip1mn=max(ip2,ip1sta)
        ip1mx=min(ip1mx,ip1end)
* The corresponding locations in the Sgm1 buffer:
        ibmn=999999
        ibmx=-999999
        do ib=ibuf1,1,-1
          ip1=ip1_buf(ib)
          if(ip1.ge.ip1mn)ibmn=ib
        end do
        do ib=1,ibuf1
          ip1=ip1_buf(ib)
          if(ip1.le.ip1mx)ibmx=ib
        end do
        nb=ibmx-ibmn+1
        if(nb.le.0) goto 99

*-----------
* Contract the Sgm1 wave functions with the Tau wave function.
*       call DGEMV_('T',nsgm1,nb,1.0D0,BUF1(:,ibmn),mxci,
*    &       buft,1,0.0D0,bufr,1)
* and distribute this result into G3:
*       call DCOPY_(nb,bufr,1,G3(iG3OFF+1),1)
* and copy the active indices into idxG3:
        do ib=1,nb
         iG3=iG3OFF+ib
         idx=ip1_buf(ibmn-1+ib)
         itlev=idx2ij(1,idx)
         iulev=idx2ij(2,idx)
         iT=l2act(itlev)
         iU=l2act(iulev)
         idxG3(1,iG3)=int(iT,I1)
         idxG3(2,iG3)=int(iU,I1)
         idxG3(3,iG3)=int(iV,I1)
         idxG3(4,iG3)=int(iX,I1)
         idxG3(5,iG3)=int(iY,I1)
         idxG3(6,iG3)=int(iZ,I1)
        end do
*       IF(IFF.ne.0) THEN
* Elementwise multiplication of Tau with H0 diagonal - EPSA(IV):
*         do icsf=1,nsgm1
*           buft(icsf)=
*    &           (bufd(icsf)-epsa(iv))*buft(icsf)
*         end do
* so Tau is now = Sum(eps(w)*E_vxww) Psi. Contract and distribute:
*         call DGEMV_('T',nsgm1,nb,1.0D0,BUF1(:ibmn),mxci,
*    &         buft,1,0.0D0,bufr,1)
*         call dcopy_(nb,bufr,1,F3(iG3OFF+1),1)
*       END IF
        iG3OFF=iG3OFF+nb
        nbtot=nbtot+nb
 99     continue
      end do
*     end do

      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)')
     &    iSubTask, ip1sta, ip1end, ip3, nbtot
        call xFlush(6)
      END IF

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.

C-SVC20100301: end of the task
      GOTO 500

 501  CONTINUE

C-SVC20100302: no more tasks, wait here for the others, then proceed
C with next symmetry
      CALL Free_Tsk(ID)

      IF(iPrGlb.GE.DEBUG) THEN
        IF (nSubTasks .GT. 0) THEN
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
        END IF
      END IF

* End of sectioning loop over symmetry of Sgm1 wave functions.
      END DO
C-SVC20100831: set correct number of elements in new G3
      NG3=iG3OFF

      CALL mma_deallocate(TASKLIST)
      ! free CI buffers
      CALL mma_deallocate(BUF1)
      CALL mma_deallocate(BUF2)
      CALL mma_deallocate(BUFT)
      CALL mma_deallocate(BUFD)

C-SVC20100302: Synchronized add into the densitry matrices
C  only for the G1 and G2 replicate arrays
      CALL GADSUM(G1,NG1)
      CALL GADSUM(G2,NG2)

      CALL GADSUM(F1,NG1)
      CALL GADSUM(F2,NG2)

#ifdef _ENABLE_BLOCK_DMRG_
      NLEV4=NLEV2**2
C
C allocate work space to store 3RDM
      Call mma_allocate(G3TMP,NLEV4,Label='G3Tmp')
C
C TODO: Here, several options to compute F3.
C Currently implemented only cu4, but cu34 and F3 from DMRG-sweep
C will be possible. They should be implemented at this section.
C
C MKFG3CU4 is located under block_dmrg_util/
      Call MKFG3CU4(IFF,G1,F1,G2,F2,G3,F3,idxG3,G3TMP)
C
      Call mma_deallocate(G3TMP)
#endif

! TODO: @kszenes: this should be wrapped in an if statement
#ifdef _ENABLE_CHEMPS2_DMRG_
      Call mkfg3chemps2(IFF,NLEV,G1,F1,G2,F2,G3,F3,idxG3)
#endif

#ifdef _DMRG_
      if (DMRG) then
        call mkfg3qcm(IFF,G1,F1,G2,F2,G3,F3,idxG3)
      endif
#endif

      IF(iPrGlb.GE.DEBUG) THEN
CSVC: if running parallel, G3/F3 are spread over processes,
C     so make sure that the _total_ fingerprint is computed
        dG1=DNRM2_(NG1,G1,1)
        dG2=DNRM2_(NG2,G2,1)
        dG3=DDOT_(NG3,G3,1,G3,1)
        CALL GADGOP_SCAL(dG3,'+')
        dG3=SQRT(dG3)
        dF1=DNRM2_(NG1,F1,1)
        dF2=DNRM2_(NG2,F2,1)
        dF3=DDOT_(NG3,F3,1,F3,1)
        CALL GADGOP_SCAL(dF3,'+')
        dF3=SQRT(dF3)

        WRITE(6,'("DEBUG> ",A)') "norms of the density matrices:"
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:", dG1
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G2:", dG2
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G3:", dG3
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "F1:", dF1
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "F2:", dF2
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "F3:", dF3
      ENDIF

 999  continue
      RETURN
      END

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#     include "macros.fh"
      subroutine empty_MKFG3DM()
      end subroutine empty_MKFG3DM

#endif
