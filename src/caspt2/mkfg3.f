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
* Copyright (C) Per Ake Malmqvist                                      *
*               Steven Vancoillie                                      *
************************************************************************
C> @brief
C>  Procedure for computing 1-body, 2-body, and 3-body
C>  density elements with active indices only,
C>  and related matrices obtained from contractions with
C>  the diagonal one-electron Hamiltonian.
C> @author Per &Aring;ke Malmqvist
C> @modified_by Steven Vancoillie
C>
C> @details
C> Computation of the 1-, 2-, and 3-body density matrices defined as
C> \f{align}{
C> G1(t,u)         &= \langle 0 \lvert E_{tu} \rvert 0 \rangle \\
C> G2(t,u,v,x)     &= \langle 0 \lvert E_{tuvx} \rvert 0 \rangle \\
C> G3(t,u,v,x,y,z) &= \langle 0 \lvert E_{tuvxyz} \rvert 0 \rangle \\
C> \f}
C> and the contractions with the diagonal 1-el Hamiltonian
C> \f{align}{
C> F1(t,u)         &= \sum_w \langle 0 \lvert E_{tuww} \rvert 0 \rangle e_w \\
C> F2(t,u,v,x)     &= \sum_w \langle 0 \lvert E_{tuvxww} \rvert 0 \rangle e_w \\
C> F3(t,u,v,x,y,z) &= \sum_w \langle 0 \lvert E_{tuvxyzww} \rvert 0 \rangle e_w \\
C> \f}
C> Storage: \p G1 and \p G2 are simple two- and four-index arrays, and
C> includes also such zeroes that are implied by symmetry.
C> But \p G3 is quite large, and while it is stored with zeroes, it
C> is made more compact by calculating only the minimum amount of
C> unique values and storing the active indices in the array \p idxG3.
C> Later, the full matrix can be restored on the fly by using the
C> full permutational symmetry (see ::mksmat and ::mkbmat). The
C> same storage applies to the \f$ F \f$ matrices.
C>
C> @param[in]  IFF   switch to activate computation of \f$ F \f$ matrices
C> @param[in]  CI    wave function CI coefficients, with symmetry \c STSYM
C> @param[out] G1    1-body active density matrix
C> @param[out] G2    2-body active density matrix
C> @param[out] G3    process-local part of 3-body active density matrix
C> @param[out] F1    1-body active density matrix contracted with
C>                   diagonal 1-el Hamiltonian
C> @param[out] F2    2-body active density matrix contracted with
C>                   diagonal 1-el Hamiltonian
C> @param[out] F3    process-local part of 3-body active density matrix
C>                   contracted with diagonal 1-el Hamiltonian
C> @param[out] idxG3 table to translate from process-local array index
C>                   to active indices

      SUBROUTINE MKFG3(IFF,CI,G1,F1,G2,F2,G3,F3,idxG3,NLEV)
      use caspt2_global, only: iPrGlb
      use fciqmc_interface, only: DoFCIQMC, mkfg3fciqmc
      use caspt2_global, only: do_grad, nbuf1_grad, nStpGrd,
     *                         iTasks_grad,nTasks_grad
      use PrintLevel, only: debug, verbose
      use gugx, only: CIS, SGS, L2ACT, EXS
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      IMPLICIT NONE
#include "caspt2.fh"
#include "SysDef.fh"
#include "pt2_guga.fh"


      INTEGER, INTENT(IN) :: IFF, NLEV
      REAL*8, INTENT(IN) :: CI(MXCI)
      REAL*8, INTENT(OUT) :: G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: G3(*), F3(*)
      INTEGER*1, INTENT(OUT) :: idxG3(6,*)

      INTEGER, PARAMETER :: I1=KIND(idxG3)
      LOGICAL RSV_TSK
      REAL*8 DG1,DG2,DG3,DF1,DF2,DF3
      REAL*8 F1SUM,F2SUM
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
      INTEGER NSGM1,NSGM2
      INTEGER NTRI1,NTRI2
      INTEGER MEMMAX, MEMMAX_SAFE
      INTEGER NLEV2
      INTEGER NCI,ICSF

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

      ! translation tables for levels i,j to and from pair indices idx
      INTEGER IJ2IDX(MXLEV,MXLEV)
      INTEGER IDX2IJ(2,MXLEV**2)
      INTEGER ICNJ(MXLEV**2)
      INTEGER IP1_BUF(MXLEV**2)

      ! result buffer, maximum size is the largest possible ip1 range,
      ! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
      REAL*8 BUFR(MXLEV**2)
      REAL*8, ALLOCATABLE:: BUF1(:,:), BUF2(:), BUFT(:), BUFD(:)
      INTEGER, ALLOCATABLE:: TASKLIST(:,:)

      Integer :: nMidV
      nMidV = CIS%nMidV

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

* Dummy values necessary for fooling syntax checkers:
      call mma_MaxDBLE(memmax)

* Use *almost* all remaining memory:
      memmax_safe=int(dble(memmax)*0.95D0)

* Buffers to compute CI expansion vectors into:
* <Psi0|E_ip1 | E_ip2 E_ip3|Psi0>
* buf1: bra buffer with E_ip1 excitations of Psi0
*       holds multiple CI vectors (allowed by memory)
* buf2: ket buffer for an E_ip3 excitation of Psi0
* buft: ket buffer for an E_ip2 excitation of E_ip3|Psi0>
* bufd: diagonal matrix elements to compute the F matrix
      nbuf1=max(1,min(nlev2,(memmax_safe-3*mxci)/mxci))
      !! if gradient, nbuf1 must be consistent here and in derfg3.f
      if (do_grad .or. nStpGrd==2) then
        !! Compute approximate available memory in derfg3.f
        !! allocated in DENS
        memmax_safe = memmax_safe - 16*NBAST**2 - 3*NASHT**2
        !! allocated (additionally) in CLagX
        memmax_safe = memmax_safe - 2*(NG1+NG2+NG3)
        nbuf1=max(1,min(nlev2,(memmax_safe-(6+nlev)*mxci)/mxci/3))
        nbuf1_grad = nbuf1
        nTasks_grad = 0
        iTasks_grad(:) = 0
      end if
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
* This also allows precomputing the Hamiltonian (H0) diagonal elements.
      DO issg1=1,nsym
        isp1=mul(issg1,stsym)
        if (.not. DoFCIQMC) then
          nsgm1=CIS%ncsf(issg1)
          CALL H0DIAG_CASPT2(ISSG1,BUFD,CIS%NOW,CIS%IOW,NMIDV)
        end if

C-SVC20100301: calculate number of larger tasks for this symmetry, this
C-is basically the number of buffers we fill with sigma1 vectors.
      iTask=1
      ibuf1=0
      DO ip1=1,nlev2
        itlev=idx2ij(1,ip1)
        iulev=idx2ij(2,ip1)
        istu=mul(SGS%ism(itlev),SGS%ism(iulev))
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
         istu=mul(SGS%ism(itlev),SGS%ism(iulev))
         it=L2ACT(itlev)
         iu=L2ACT(iulev)
         if(istu.eq.isp1) then
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1i
          if (.not. DoFCIQMC) then
              call dcopy_(nsgm1,[0.0D0],0,BUF1(:,ibuf1),1)
              CALL SIGMA1(SGS,CIS,EXS,
     &                    IULEV,ITLEV,1.0D00,STSYM,CI,BUF1(:,ibuf1))
          end if
         end if
        end do
        myBuffer=iTask
      ELSE
        ibuf1=TaskList(iTask,3)
      ENDIF
C-SVC20100301: necessary batch of sigma vectors is now in the buffer
      if (.not. DoFCIQMC) then
          ! The ip1 buffer could be the same on different processes
          ! so only compute the G1 contribution when ip3 is 1, as
          ! this will only be one task per buffer.
          if (issg1.eq.stsym.AND.ip3.eq.1) then
            do ib=1,ibuf1
              idx=ip1_buf(ib)
              itlev=idx2ij(1,idx)
              iulev=idx2ij(2,idx)
              it=L2ACT(itlev)
              iu=L2ACT(iulev)
              G1(it,iu)=DDOT_(nsgm1,ci,1,BUF1(:,ib),1)
              IF(IFF.ne.0) then
                F1sum=0.0D0
                do i=1,nsgm1
                  F1sum=F1sum+CI(i)*BUF1(i,ib)*bufd(i)
                end do
                F1(it,iu)=F1sum-EPSA(iu)*G1(it,iu)
              end if
            end do
          end if
      end if

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

      iq1=icnj(ip3)
* The indices corresponding to pair index p3:
      iylev=idx2ij(1,ip3)
      izlev=idx2ij(2,ip3)
      isyz=mul(SGS%ism(iylev),SGS%ism(izlev))
      issg2=mul(isyz,stsym)
      if (.not. DoFCIQMC) then
         nsgm2=CIS%ncsf(issg2)
      end if
      iy=L2ACT(iylev)
      iz=L2ACT(izlev)
      if (.not. DoFCIQMC) then
          call dcopy_(nsgm2,[0.0D0],0,BUF2,1)
          CALL SIGMA1(SGS,CIS,EXS,
     &                IYLEV,IZLEV,1.0D00,STSYM,CI,BUF2)
          if(issg2.eq.issg1) then
            do ib=1,ibuf1
              idx=ip1_buf(ib)
              itlev=idx2ij(1,idx)
              iulev=idx2ij(2,idx)
              it=L2ACT(itlev)
              iu=L2ACT(iulev)
              G2(it,iu,iy,iz)=DDOT_(nsgm1,BUF2,1,BUF1(:,ib),1)
              IF(IFF.ne.0) THEN
                F2sum=0.0D0
                do i=1,nsgm1
                  F2sum=F2sum+BUF2(i)*bufd(i)*buf1(i,ib)
                end do
                F2(it,iu,iy,iz)=F2sum
              END IF
            end do
          end if
      end if
      nbtot=0
      do ip2=ip3,ntri2
        ivlev=idx2ij(1,ip2)
        ixlev=idx2ij(2,ip2)
        isvx=mul(SGS%ism(ivlev),SGS%ism(ixlev))
        iv=L2ACT(ivlev)
        ix=L2ACT(ixlev)
        if(isvx.ne.mul(issg1,issg2)) goto 99
        if (.not. DoFCIQMC) then
            call dcopy_(nsgm1,[0.0D0],0,BUFT,1)
            CALL SIGMA1(SGS,CIS,EXS,
     &                  IVLEV,IXLEV,1.0D00,ISSG2,BUF2,BUFT)
        end if
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
        if (.not. DoFCIQMC) then
            call DGEMV_ ('T',nsgm1,nb,1.0D0,BUF1(:,ibmn),mxci,
     &           buft,1,0.0D0,bufr,1)
* and distribute this result into G3:
            call dcopy_(nb,bufr,1,G3(iG3OFF+1),1)
* and copy the active indices into idxG3:
        end if
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
        if (.not. DoFCIQMC) then
            IF(IFF.ne.0) THEN
* Elementwise multiplication of Tau with H0 diagonal - EPSA(IV):
                do icsf=1,nsgm1
                  buft(icsf)=
     &                 (bufd(icsf)-epsa(iv))*buft(icsf)
                end do
* so Tau is now = Sum(eps(w)*E_vxww) Psi. Contract and distribute:
                call DGEMV_ ('T',nsgm1,nb,1.0D0,BUF1(:,ibmn),mxci,
     &           buft,1,0.0D0,bufr,1)
                call dcopy_(nb,bufr,1,F3(iG3OFF+1),1)
            END IF
        end if
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
      !! consistent tasks must be executed here and in derfg3.f for grad
      if (do_grad .or. nStpGrd==2) then
        nTasks_grad = nTasks_grad + 1
        iTasks_grad(nTasks_grad) = iSubTask
      end if

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

      if (DoFCIQMC) then
          call mkfg3fciqmc(G1,G2,G3,F1,F2,F3,idxG3,nLev)
      else
          ! Correction to G2: It is now = <0| E_tu E_yz |0>
          do iu=1,nlev
           do iz=1,nlev
            do it=1,nlev
             G2(it,iu,iu,iz)=G2(it,iu,iu,iz)-G1(it,iz)
            end do
           end do
          end do
        ! SVC20100310: took some spurious mirroring of G2 values out
        ! of the loops and put them here, after the parallel section has
        ! finished, so that GAdSUM works correctly.
          do ip1=ntri2+1,nlev2
           itlev=idx2ij(1,ip1)
           iulev=idx2ij(2,ip1)
           it=L2ACT(itlev)
           iu=L2ACT(iulev)
           do ip3=ntri1+1,ip1
            iylev=idx2ij(1,ip3)
            izlev=idx2ij(2,ip3)
            iy=L2ACT(iylev)
            iz=L2ACT(izlev)
            G2(it,iu,iy,iz)=G2(iz,iy,iu,it)
           end do
          end do
       ! Correction to G2: Some values not computed follow from symmetry
          do ip1=1,nlev2-1
           itlev=idx2ij(1,ip1)
           iulev=idx2ij(2,ip1)
           it=L2ACT(itlev)
           iu=L2ACT(iulev)
           do ip3=ip1+1,nlev2
            iylev=idx2ij(1,ip3)
            izlev=idx2ij(2,ip3)
            iy=L2ACT(iylev)
            iz=L2ACT(izlev)
            G2(it,iu,iy,iz)=G2(iy,iz,it,iu)
           end do
          end do
          IF(IFF.ne.0) THEN
           ! Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
           do iz=1,nlev
            do iy=1,nlev
             do iu=1,nlev
              do it=1,nlev
               F2(it,iu,iy,iz)=F2(it,iu,iy,iz)-
     &               (EPSA(iu)+EPSA(iy))*G2(it,iu,iy,iz)
              end do
             end do
            end do
           end do
           do iz=1,nlev
            do iu=1,nlev
             do it=1,nlev
              F2(it,iu,iu,iz)=F2(it,iu,iu,iz)-
     &              (F1(it,iz)+EPSA(iu)*G1(it,iz))
             end do
            end do
           end do
        ! SVC20100310: took some spurious mirroring of F2 values out
        ! of the loops and put them here, after the parallel section has
        ! finished, so that GAdSUM works correctly.
           do ip1=ntri2+1,nlev2
            itlev=idx2ij(1,ip1)
            iulev=idx2ij(2,ip1)
            it=L2ACT(itlev)
            iu=L2ACT(iulev)
            do ip3=ntri1+1,ip1
             iylev=idx2ij(1,ip3)
             izlev=idx2ij(2,ip3)
             iy=L2ACT(iylev)
             iz=L2ACT(izlev)
             F2(it,iu,iy,iz)=F2(iz,iy,iu,it)
            end do
           end do
      ! Correction to F2: Some values not computed follow from symmetry
           do ip1=1,nlev2-1
            itlev=idx2ij(1,ip1)
            iulev=idx2ij(2,ip1)
            it=L2ACT(itlev)
            iu=L2ACT(iulev)
            do ip3=ip1+1,nlev2
             iylev=idx2ij(1,ip3)
             izlev=idx2ij(2,ip3)
             iy=L2ACT(iylev)
             iz=L2ACT(izlev)
             F2(it,iu,iy,iz)=F2(iy,iz,it,iu)
            end do
           end do
          END IF

      ! Correction to G3: It is now <0| E_tu E_vx E_yz |0>
      ! Similar for F3 values.
          DO iG3=1,NG3
           iT=idxG3(1,iG3)
           iU=idxG3(2,iG3)
           iV=idxG3(3,iG3)
           iX=idxG3(4,iG3)
           iY=idxG3(5,iG3)
           iZ=idxG3(6,iG3)
      ! Correction: From <0| E_tu E_vx E_yz |0>, form <0| E_tuvxyz |0>
           if(iY.eq.iX) then
            G3(iG3)=G3(iG3)-G2(iT,iU,iV,iZ)
            IF(IFF.ne.0) F3(iG3)=
     &               F3(iG3)-(F2(iT,iU,iV,iZ)+EPSA(iu)*G2(iT,iU,iV,iZ))
            if(iv.eq.iu) then
             G3(iG3)=G3(iG3)-G1(iT,iZ)
             IF(IFF.ne.0) F3(iG3)=F3(iG3)-F1(iT,iZ)
            end if
           end if
           if(iV.eq.iU) then
             G3(iG3)=G3(iG3)-G2(iT,iX,iY,iZ)
             IF(IFF.ne.0) F3(iG3)=F3(iG3)-
     &                (F2(iT,iX,iY,iZ)+EPSA(iY)*G2(iT,iX,iY,iZ))
           end if
           if(iY.eq.iU) then
             G3(iG3)=G3(iG3)-G2(iV,iX,iT,iZ)
             IF(IFF.ne.0) F3(iG3)=
     &                F3(iG3)-(F2(iV,iX,iT,iZ)+EPSA(iU)*G2(iV,iX,iT,iZ))
           end if
           IF(IFF.ne.0) F3(iG3)=F3(iG3)-(EPSA(iU)+EPSA(iY))*G3(iG3)
          END DO
      end if

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

        WRITE(6,'("DEBUG> ",A)') "MKFG3: norms of the density matrices:"
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
