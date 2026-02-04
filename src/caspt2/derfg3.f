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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      SUBROUTINE DERFG3(CI,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,
     *                  DEPSA,G1,G2,nLev)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug, verbose
#ifdef _MOLCAS_MPP_
      use caspt2_global, only: nbuf1_grad,iTasks_grad,nTasks_grad
#else
      use caspt2_global, only: nbuf1_grad
#endif
      use gugx, only: CIS, L2ACT, SGS, EXS
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      use definitions, only: iwp,wp,u6
      use SysDef, only: RtoB
      use caspt2_module, only: nConf, nActEl, nSym, STSym, EPSA, Mul
      use pt2_guga, only: MxCI, MxLev
      use Constants, only: Zero, One, Half

      IMPLICIT NONE

      logical(kind=iwp), external :: RSV_TSK

      real(kind=wp), intent(in) :: CI(MXCI), DG3(*), DF3(*),
     &  G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
      real(kind=wp), intent(inout) :: CLAG(NCONF), DG1(NLEV,NLEV),
     &  DG2(NLEV,NLEV,NLEV,NLEV), DF1(NLEV,NLEV),
     &  DF2(NLEV,NLEV,NLEV,NLEV), DEPSA(NLEV,NLEV)
      integer(kind=iwp), intent(in) :: nLev

      integer(kind=iwp) :: I,J,IDX,JDX
      integer(kind=iwp) :: IB,IBMN,IBMX,IBUF,NB,NBTOT,IBUF1
      integer(kind=iwp) :: IP1,IP2,IP3,IP1MN,IP1MX,IP1I,IP1STA,IP1END,
     &                     IP3MX,IQ1
      integer(kind=iwp) :: IG3,IG3OFF,IG3BK
      integer(kind=iwp) :: ISTU,ISVX,ISYZ
      integer(kind=iwp) :: IT,IU,IV,IX,IY,IZ
      integer(kind=iwp) :: ITLEV,IULEV,IVLEV,IXLEV,IYLEV,IZLEV,IXLEV0
      integer(kind=iwp) :: NBUF1,NBUFX,NDTU,NDAB
      integer(kind=iwp) :: IOFFSET
      integer(kind=iwp) :: ISSG1,ISSG2,ISP1
      integer(kind=iwp) :: ITASK,ISUBTASK,ID,NTASKS,NSUBTASKS,
     &                     MXTASK,MYTASK,MYBUFFER
      integer(kind=iwp) :: NSGM1,NSGM2
      integer(kind=iwp) :: NTRI1,NTRI2
      integer(kind=iwp) :: MEMMAX !, MEMMAX_SAFE
      integer(kind=iwp) :: NLEV2
      integer(kind=iwp) :: NCI,ICSF

      real(kind=wp), external :: DDOT_

      ! translation tables for levels i,j to and from pair indices idx
      integer(kind=iwp) :: IJ2IDX(MXLEV,MXLEV)
      integer(kind=iwp) :: IDX2IJ(2,MXLEV**2)
      integer(kind=iwp) :: ICNJ(MXLEV**2)
      integer(kind=iwp) :: IP1_BUF(MXLEV**2)

      REAL(KIND=WP),ALLOCATABLE:: BUF1(:,:),BUF2(:),BUFT(:),BUFD(:),
     *  DTU(:,:),DYZ(:),DAB(:,:),BUF3(:),BUF4(:),BUFX(:,:)
      INTEGER(KIND=IWP),ALLOCATABLE:: TASKLIST(:,:)

      ! result buffer, maximum size is the largest possible ip1 range,
      ! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
C     REAL*8 BUFR(MXLEV**2)
C
C     integer(kind=iwp) :: LFCDer1,LFCDer2
      integer(kind=iwp) :: iTask_loc
      real(kind=wp) :: SCAL,ScalG,ScalF
C     real(kind=wp) :: tmp,tmp2
      logical(kind=iwp) :: first
      integer(kind=iwp) :: nMidV
      nMidV = CIS%nMidV
C
C Put in zeroes. Recognize special cases:
      IF(nlev == 0) return
      IF(NACTEL == 0) return

      NCI=CIS%NCSF(STSYM)
* This should not happen, but...
      IF(NCI == 0) return

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
C
* Correction to G3: It is now <0| E_tu E_vx E_yz |0>
* Similar for F3 values.
* These corrections are already done in CLagDXA_FG3 and CLagDXC_FG3
C
* Correction to F2: Some values not computed follow from symmetry
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
        SCAL=DF2(iy,iz,it,iu)+DF2(it,iu,iy,iz)
        DF2(it,iu,iy,iz)=Zero
        DF2(iy,iz,it,iu)=Scal
       end do
      end do
C-SVC20100310: took some spurious mirroring of F2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
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
        Scal=DF2(iz,iy,iu,it)+DF2(it,iu,iy,iz)
        DF2(it,iu,iy,iz)=Zero
        DF2(iz,iy,iu,it)=Scal
       end do
      end do
* Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
       do iz=1,nlev
        do iy=1,nlev
         do iu=1,nlev
          do it=1,nlev
           DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)
     *      - DF2(it,iu,iy,iz)*(EPSA(iu)+EPSA(iy))
           !! DEPSA-related operations in the preparation step
           !! are done only on the master node
#ifdef _MOLCAS_MPP_
           if (king()) then
#endif
            do ix=1,nlev
             DEPSA(iu,ix) = DEPSA(iu,ix)
     *         - DF2(it,iu,iy,iz)*G2(it,ix,iy,iz)
             DEPSA(ix,iy) = DEPSA(ix,iy)
     *         - DF2(it,iu,iy,iz)*G2(it,iu,ix,iz)
            end do
#ifdef _MOLCAS_MPP_
           end if
#endif
          end do
         end do
        end do
       end do
       do iz=1,nlev
        do iu=1,nlev
         do it=1,nlev
          DG1(it,iz) = DG1(it,iz) - EPSA(iu)*DF2(it,iu,iu,iz)
          DF1(it,iz) = DF1(it,iz) - DF2(it,iu,iu,iz)
#ifdef _MOLCAS_MPP_
          if (king()) then
#endif
           do iy = 1, nlev
             DEPSA(iu,iy) = DEPSA(iu,iy) - G1(it,iz)*DF2(it,iu,iy,iz)
           end do
#ifdef _MOLCAS_MPP_
          end if
#endif
         end do
        end do
       end do
C
C
C
* Correction to G2: Some values not computed follow from symmetry
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
        SCAL=DG2(iy,iz,it,iu)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=Zero
        DG2(iy,iz,it,iu)=Scal
       end do
      end do
C-SVC20100310: took some spurious mirroring of G2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
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
        Scal=DG2(iz,iy,iu,it)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=Zero
        DG2(iz,iy,iu,it)=Scal
       end do
      end do
* Correction to G2: It is now = <0| E_tu E_yz |0>
      do iu=1,nlev
       do iz=1,nlev
        do it=1,nlev
         DG1(it,iz) = DG1(it,iz) - DG2(it,iu,iu,iz)
        end do
       end do
      end do
* Additional correction terms for F1
      Do iT = 1, NLEV
        Do iU = 1, NLEV
          !! With the improved algorithm, symmetrizing the DF1
          !! contribution to DF1 here is somehow required
          DG1(iT,iU) = DG1(iT,iU)
     *      - (DF1(iT,iU)+DF1(iU,iT))*EPSA(iU)*Half
        End Do
      End Do
#ifdef _MOLCAS_MPP_
      if (king()) then
#endif
        Do iT = 1, NLEV
          Do iU = 1, NLEV
            Do iV = 1, NLEV
              Do iX = 1, nLEV
                DEPSA(iV,iX) = DEPSA(iV,iX) + G2(iT,iU,iV,iX)*DF1(iT,iU)
              End Do
            End DO
          End DO
        End Do
#ifdef _MOLCAS_MPP_
      endif
#endif
C
      Call CLagSym(nLev,DG1,DG2,DF1,DF2,2)
C
      call mma_MaxDBLE(memmax)

* Use *almost* all remaining memory:
C     memmax_safe=int(dble(memmax)*0.95D0)

* Buffers to compute CI expansion vectors into:
* <Psi0|E_ip1 | E_ip2 E_ip3|Psi0>
* buf1: bra buffer with E_ip1 excitations of Psi0
*       holds multiple CI vectors (allowed by memory)
* buf2: ket buffer for an E_ip3 excitation of Psi0
* buft: ket buffer for an E_ip2 excitation of E_ip3|Psi0>
* bufd: diagonal matrix elements to compute the F matrix
      nbuf1= nbuf1_grad
      ndtu = nbuf1
      ndab = nbuf1
      nbufx= nlev
      call mma_allocate(BUF1,MXCI,NBUF1,Label='BUF1')
      call mma_allocate(BUF2,MXCI,Label='BUF2')
      call mma_allocate(BUFT,MXCI,Label='BUFT')
      call mma_allocate(BUFD,MXCI,Label='BUFD')
C
      call mma_allocate(DTU,MXCI,NDTU,Label='DTU')
      call mma_allocate(DYZ,MXCI,Label='DYZ')
      call mma_allocate(DAB,MXCI,NDAB,Label='DAB')
      call mma_allocate(BUF3,MXCI,Label='BUF3')
      call mma_allocate(BUF4,MXCI,Label='BUF4')
C
      call mma_allocate(BUFX,MXCI,NBUFX,Label='BUFX')

C-SVC20100301: calculate maximum number of tasks possible
      MXTASK=(NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
      CALL mma_allocate (TaskList,mxTask,4,LABEL='TaskList')

      IF(iPrGlb >= VERBOSE) THEN
        WRITE(u6,*)
        WRITE(u6,'(2X,A)') 'Constructing derivatives of G3/F3'
        WRITE(u6,'(2X,A,F16.9,A)') ' memory avail: ',
     &    (memmax*RtoB)/1.0e+09_wp, ' GB'
        WRITE(u6,'(2X,A,F16.9,A)') ' memory used:  ',
     &    (((3*nbuf1+6+nlev)*MXCI)*RtoB)/1.0e+09_wp, ' GB'
        call xFlush(u6)
      ENDIF
C     CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
C     CPUT =CPTF10-CPTF0
C     WALLT=TIOTF10-TIOTF0
C     write(6,*) "PREP    : CPU/WALL TIME=", cput,wallt

      iG3OFF=0
      iTask_loc = 1
      first = .true.
* A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.
* This also allows precomputing the Hamiltonian (H0) diagonal elements.
      DO issg1=1,nsym
       isp1=mul(issg1,STSYM)
       nsgm1=CIS%ncsf(issg1)
       !! Work(LBufD) = \sum_t <I|E_{tt}|I>*f_{tt}
       CALL H0DIAG_CASPT2(ISSG1,BUFD,CIS%NOW,CIS%IOW,nMidV)

C-SVC20100301: calculate number of larger tasks for this symmetry, this
C-is basically the number of buffers we fill with sigma1 vectors.
      iTask=1
      ibuf1=0
      DO ip1=1,nlev2
        itlev=idx2ij(1,ip1)
        iulev=idx2ij(2,ip1)
        istu=mul(SGS%ism(itlev),SGS%ism(iulev))
        IF (istu == isp1) THEN
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1
          IF (ibuf1 == 1) TaskList(iTask,1)=ip1
        ENDIF
        IF (ibuf1 == nbuf1 .OR. (ibuf1 > 0 .AND.
     &         (ip1 == ntri2 .OR. ip1 == nlev2))) THEN
            TaskList(iTask,2)=ip1_buf(ibuf1)
            TaskList(iTask,3)=ibuf1
            iTask=iTask+1
            ibuf1=0
        ENDIF
      ENDDO
      nTasks=iTask
C     write(6,*) "nTasks = ", nTasks
      IF (ibuf1 == 0) nTasks=nTasks-1
C-SVC20100309: calculate number of inner loop iteration tasks.
      iOffSet=0
      DO iTask=1,nTasks
        TaskList(iTask,4)=iOffSet
        ip1sta=TaskList(iTask,1)
        ip1end=TaskList(iTask,2)
        ip3mx=ntri2
        if(ip1end <= ntri2) ip3mx=ip1end
        if(ip1sta > ntri2) ip3mx=ntri1
C       write(6,*) "iTask = ", iTask
C       write(6,*) "start,end=",ip1sta,ip1end
C       write(6,*) "ip3mx = ", ip3mx
C-SVC20100309: Currently -we are going to limit this to the ip3-loop and
C-leave the ip2-loop intact.  This was based on the large overhead which
C-was observed for a very large number of small tasks.
C       iOffSet=iOffSet+ip3mx*ntri2-((ip3mx**2-ip3mx)/2)
        iOffSet=iOffSet+ip3mx
      ENDDO
      nSubTasks=iOffSet
C     write(6,*) "nSubTasks = ", nSubTasks

      IF(iPrGlb >= VERBOSE) THEN
        WRITE(u6,'(2X,A,I3,A,I6)') 'Sym: ',issg1,', #Tasks: ',nSubTasks
        call xFlush(u6)
      ENDIF

      IF(iPrGlb >= DEBUG) THEN
        IF (nSubTasks > 0) THEN
          WRITE(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          WRITE(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "task ID ",
     &    " ip1 range  ",
     &    "ip3 ",
     &    "#elements"
          WRITE(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          call xFlush(u6)
        END IF
      END IF

C-SVC20100301: initialize the series of subtasks
      Call Init_Tsk(ID, nSubTasks)

      myBuffer=0

      !! loop start
      do
C-SVC20100908: first check: can I actually do any task?
C     IF ((NG3-iG3OFF) < nbuf1*ntri2) GOTO 501
C-SVC20100831: initialize counter for offset into G3
C-SVC20100302: BEGIN SEPARATE TASK EXECUTION
C     write(6,*) rsv_tsk(id,isubtask)
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          !! do the same tasks here and in mkfg3.f
          iSubTask = iTasks_grad(iTask_loc)
          if (iSubTask == 0) exit
        else
#endif
          If (.NOT.Rsv_Tsk(ID,iSubTask)) exit
#ifdef _MOLCAS_MPP_
        end if
#endif

        myTask=nTasks
        DO iTask=1,nTasks
          iBuf=iSubTask-TaskList(iTask,4)
          IF (iBuf <= 0) THEN
            myTask=iTask-1
            exit
          ENDIF
        ENDDO
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
C     write(6,*) "myBuffer,iTask = ", myBuffer,iTask
        IF (myBuffer /= iTask) THEN
          if (.not.first) then
            !! Compute left derivative and DEPSA contributions before
            !! the TASK is completely switched
            CALL LEFT_DEPSA()
          end if
          ibuf1=0
          do ip1i=ip1sta,ip1end
           itlev=idx2ij(1,ip1i)
           iulev=idx2ij(2,ip1i)
           istu=mul(SGS%ism(itlev),SGS%ism(iulev))
           it=L2ACT(itlev)
           iu=L2ACT(iulev)
           if(istu == isp1) then
            ibuf1=ibuf1+1
            ip1_buf(ibuf1)=ip1i
            BUF1(1:nsgm1,ibuf1) = Zero
            CALL SIGMA1(SGS,CIS,EXS,
     &                  IULEV,ITLEV,One,STSYM,CI,BUF1(:,ibuf1))
           end if
          end do
          myBuffer=iTask
          DTU(1:MXCI,1:ibuf1) = Zero
          DAB(1:MXCI,1:ibuf1) = Zero
          if (first) first = .false.
        ELSE
          ibuf1=TaskList(iTask,3)
        ENDIF
C-SVC20100301: necessary batch of sigma vectors is now in the buffer

        ! The ip1 buffer could be the same on different processes
        ! so only compute the G1 contribution when ip3 is 1, as
        ! this will only be one task per buffer.
        if (issg1 == STSYM .AND. ip3 == 1) then
          !! buf1 = <Psi0|E_ip1|I>
          !! <0|E_{tu}I> = <I|E_{ut}|0>
C         write(6,*) "ib loop"
          do ib=1,ibuf1
            idx=ip1_buf(ib)
            itlev=idx2ij(1,idx)
            iulev=idx2ij(2,idx)
            it=L2ACT(itlev)
            iu=L2ACT(iulev)
C           write(6,'("itlev,iulev,it,iu = ",4i3)') itlev,iulev,it,iu
            !! DG1 contribution
            SCAL = DG1(iT,iU) + DG1(iT,iU)
            CLag(1:nsgm1) = CLag(1:nsgm1) + SCAL*BUF1(1:nsgm1,ib)
C
            !! left derivative of DF1
            Do icsf = 1, nsgm1
              DTU(icsf,ib) = DTU(icsf,ib)
     &          + DF1(it,iu)*BUFD(icsf)*CI(icsf)
            End Do
            !! right derivative of DF1
            Do icsf = 1, nsgm1
              CLag(icsf) = CLag(icsf)
     *          + DF1(it,iu)*BUF1(icsf,ib)*BUFD(icsf)
            End Do
C           G1(it,iu)=DDOT_(nsgm1,ci,1,work(lto),1)
C           IF(IFF.ne.0) then
C             F1sum=0.0D0
C             do i=1,nsgm1
C               F1sum=F1sum+CI(i)*work(lto-1+i)*bufd(i)
C             end do
C             F1(it,iu)=F1sum-EPSA(iu)*G1(it,iu)
C           end if
          end do
        end if

C       ip3mx=ntri2
C       if(ip1end <= ntri2) ip3mx=ip1end
C       if(ip1sta > ntri2) ip3mx=ntri1
C-SVC20100309: loop over ip3, ip2
C       do ip3=1,ip3mx

C-SVC20100309: PAM's magic formula
*       iCnt=iSubTask-iOffSet
*       ip3=int(dble(ntri2)+1.5D0 -
*    &       sqrt((dble(ntri2)+0.5d0)**2-2*iCnt+0.000001D0))
*       ip2=iCnt-((ip3-1)*ntri2-((ip3-1)*(ip3-2))/2 )+ip3-1

C-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

C       write(6,*) "ip3 = ", ip3
C       CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        iq1=icnj(ip3)
* The indices corresponding to pair index p3:
        iylev=idx2ij(1,ip3)
        izlev=idx2ij(2,ip3)
        isyz=mul(SGS%ism(iylev),SGS%ism(izlev))
        issg2=mul(isyz,STSYM)
        nsgm2=CIS%ncsf(issg2)
        iy=L2ACT(iylev)
        iz=L2ACT(izlev)
        buf2(1:nsgm2) = Zero
        CALL SIGMA1(SGS,CIS,EXS,IYLEV,IZLEV,One,STSYM,CI,BUF2)
        DYZ(1:nsgm1) = Zero
        if(issg2 == issg1) then
          buf3(1:nsgm2) = Zero
          do ib=1,ibuf1
            idx=ip1_buf(ib)
            itlev=idx2ij(1,idx)
            iulev=idx2ij(2,idx)
            it=L2ACT(itlev)
            iu=L2ACT(iulev)
C
            ScalG=DG2(iT,iU,iY,iZ)
            ScalF=DF2(iT,iU,iY,iZ)
            If (ScalG == Zero .and. ScalF == Zero) cycle
C
            !! left derivative
            Do icsf = 1, nsgm1
             BUFT(icsf) = ScalG*BUF2(icsf) + ScalF*BUF2(icsf)*BUFD(icsf)
            End Do
            DTU(1:nsgm1,ib) = DTU(1:nsgm1,ib) + BUFT(1:nsgm1)
C
            !! right derivative
            buft(1:mxci) = buf1(1:mxci,ib)
            Do icsf = 1, nsgm1
              BUF3(icsf) = BUF3(icsf)
     *          + ScalG*BUFT(icsf) + ScalF*BUFT(icsf)*BUFD(icsf)
            End Do
C
            !! For DEPSA
            DAB(1:nsgm1,ib) = DAB(1:nsgm1,ib) + ScalF*BUF2(1:nsgm1)
          end do
          !! Save for Eyz
          DYZ(1:nsgm2) = DYZ(1:nsgm2) + BUF3(1:nsgm2)
        end if
        nbtot=0
C
C
C
        !! Prepare for DEPSA for the -epsa(iv) term with square
        iG3bk = iG3OFF
        Do ixlev0 = 1, nlev
          Do ivlev = 1, nlev
            BUFX(1:nsgm1,ivlev) = Zero
            CALL SIGMA1(SGS,CIS,EXS,
     &                  IVLEV,IXLEV0,One,STSYM,BUF2,BUFX(1,ivlev))
          End Do
          iG3OFF = iG3bk
          do ip2=ip3,ntri2
            ivlev=idx2ij(1,ip2)
            ixlev=idx2ij(2,ip2)
            isvx=mul(SGS%ism(ivlev),SGS%ism(ixlev))
            iv=L2ACT(ivlev)
            ix=L2ACT(ixlev)
            if(isvx /= mul(issg1,issg2)) cycle
            !! <I|EvxEyz|0>
            If (IXLEV == IXLEV0) BUFT(1:nsgm1) = BUFX(1:nsgm1,IVLEV)
*-----------
* Max and min values of index p1:
            ip1mx=ntri2
            if(ip3 <= ntri1) then
              ip1mx=nlev2
              if(ip2 > ntri1) ip1mx=iq1
            end if
            ip1mn=max(ip2,ip1sta)
            ip1mx=min(ip1mx,ip1end)
* The corresponding locations in the Sgm1 buffer:
            ibmn=999999
            ibmx=-999999
            do ib=ibuf1,1,-1
              ip1=ip1_buf(ib)
              if(ip1 >= ip1mn)ibmn=ib
            end do
            do ib=1,ibuf1
              ip1=ip1_buf(ib)
              if(ip1 <= ip1mx)ibmx=ib
            end do
            nb=ibmx-ibmn+1
            if(nb <= 0) cycle
            if (ixlev /= ixlev0) then
              iG3OFF = iG3OFF + nb
              cycle
            end if
C
C         ----- left derivative
C
            do icsf = 1, nsgm1
              ! BUF3 = (<I|Ett|I>-EPSA(V))*<I|EvxEyz|0> = <I|fEvxEyz|0>
              buf3(icsf) = (bufd(icsf)-epsa(iv))*buft(icsf)
            end do
            do ib=1,nb
              iG3=iG3OFF+ib
              idx=ip1_buf(ibmn-1+ib)
C
              !! <I|EvxEyz|0>*Dtuvxyz + <I|fEvxEyz|0>*Ftuvxyz
              DTU(1:nsgm1,ibmn+ib-1) = DTU(1:nsgm1,ibmn+ib-1)
     &          + DG3(iG3)*BUFT(1:nsgm1) + DF3(iG3)*BUF3(1:nsgm1)
              !! DEPSA of the BUFD term
              DAB(1:nsgm1,ibmn+ib-1) = DAB(1:nsgm1,ibmn+ib-1)
     &          + DF3(iG3)*BUFT(1:nsgm1)
            end do
C
C         ----- right derivative
C
            BUF3(1:nsgm1) = Zero
            BUF4(1:nsgm1) = Zero
            !! right derivative (1):
            !! <0|Etu|I>*Dtuvxyz and <0|Etu|I>*Ftuvxyz
            !! <0|EtuEvxEyz|I> -> <I|EzyExvEut|0>
            do ib=1,nb
              iG3=iG3OFF+ib
              idx=ip1_buf(ibmn-1+ib)
C
              !! BUF3 = <0|Etu|I>*Dtuvxyz
              BUF3(1:nsgm1) = BUF3(1:nsgm1)
     &          + DG3(iG3)*BUF1(1:nsgm1,ibmn+ib-1)
              !! BUFC = <0|Etu|I>*Ftuvxyz
              BUF4(1:nsgm1) = BUF4(1:nsgm1)
     &          + DF3(iG3)*BUF1(1:nsgm1,ibmn+ib-1)
            end do
C
            !! DEPSA of the -EPSA(iv) term
            Call DGEMV_('T',nsgm1,NLEV,
     *                 -One,BUFX,mxci,buf4,1,
     *                  One,DEPSA(1,IVLEV),1)
C
            !! Scale the DF3 contribution with the diagonal Fock
            !! and add to the DG3 contribution
            do icsf = 1, nsgm1
              buf3(icsf) = buf3(icsf) + buf4(icsf)*(bufd(icsf)-epsa(iv))
            end do
            !! right derivative (2): <0|EtuEvx|I>*Dtuvxyz
            CALL SIGMA1(SGS,CIS,EXS,IXLEV,IVLEV,One,STSYM,BUF3,DYZ)
C
            iG3OFF=iG3OFF+nb
            nbtot=nbtot+nb
          end do !! end of ip2 loop
        End Do !! end of ixlev0 loop
C
        !! Complete the right derivative contribution:
        !! <0|EtuEyz|I> and <0|EtuEvxEyz|I>
        CALL SIGMA1(SGS,CIS,EXS,IZLEV,IYLEV,One,STSYM,DYZ,CLAG)
C
        IF(iPrGlb >= DEBUG) THEN
          WRITE(u6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)')
     &      iSubTask, ip1sta, ip1end, ip3, nbtot
          call xFlush(u6)
        END IF
        iTask_loc = iTask_loc + 1

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.

      end do
C-SVC20100301: end of the task

      !! Final (the last task) left derivative and DEPSA contributions
      CALL LEFT_DEPSA()

C-SVC20100302: no more tasks, wait here for the others, then proceed
C with next symmetry
      CALL Free_Tsk(ID)

      IF(iPrGlb >= DEBUG) THEN
        IF (nSubTasks > 0) THEN
          WRITE(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
        END IF
      END IF

* End of sectioning loop over symmetry of Sgm1 wave functions.
      END DO
C
#ifdef _MOLCAS_MPP_
      if (is_real_par() .and. (iTask_loc-1 /= nTasks_grad)) then
        write (u6,*)
        write (u6,*) "Somehow, the number of tasks in mkfg3.f and ",
     *              "derfg3.f is not consistent..."
        write (u6,*) "probably, bug"
        write (u6,*) "# of tasks in  mkfg3.f = ", nTasks_grad
        write (u6,*) "# of tasks in derfg3.f = ", iTask_loc-1
        call abend()
      end if
#endif
C
      CALL mma_deallocate(TASKLIST)
      ! free CI buffers
      call mma_deallocate(BUF1)
      call mma_deallocate(BUF2)
      call mma_deallocate(BUFT)
      call mma_deallocate(BUFD)
C
      call mma_deallocate(DTU)
      call mma_deallocate(DYZ)
      call mma_deallocate(DAB)
      call mma_deallocate(BUF3)
      call mma_deallocate(BUF4)
C
      call mma_deallocate(BUFX)
C
      contains

      SUBROUTINE LEFT_DEPSA()

      IMPLICIT NONE

      integer(kind=iwp) :: IALEVloc, IBLEVloc, ibloc

      do ibloc=1,ibuf1
        idx=ip1_buf(ibloc)
        itlev=idx2ij(1,idx)
        iulev=idx2ij(2,idx)
        !! left derivative
        CALL SIGMA1(SGS,CIS,EXS,ITLEV,IULEV,One,STSYM,DTU(1,ibloc),CLAG)
        !! the rest is DEPSA contribution
        Do IALEVloc = 1, NLEV
          Do IBLEVloc = 1, NLEV
            BUF2(:) = Zero
            CALL SIGMA1(SGS,CIS,EXS,
     &                 IALEVloc,IBLEVloc,One,STSYM,DAB(1,ibloc),BUF2)
            DEPSA(IALEVloc,IBLEVloc) = DEPSA(IALEVloc,IBLEVloc)
     *        + DDot_(nsgm1,BUF1(1,IBloc),1,BUF2,1)
          End Do
        End Do
      end do

      END SUBROUTINE LEFT_DEPSA

      END SUBROUTINE DERFG3
C
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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE DERSPE(DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use gugx, only: SGS, LEVEL
      use caspt2_module, only: NACTEL, NASHT, ISCF
      use pt2_guga, only: ETA, NG3
      use Constants, only: Zero, One, Two
      use definitions, only: wp, iwp, byte, u6

      implicit none

      real(kind=wp), intent(in) :: DF1(NASHT,NASHT),
     &  DF2(NASHT,NASHT,NASHT,NASHT), DF3(*), G1(NASHT,NASHT),
     &  G2(NASHT,NASHT,NASHT,NASHT), G3(*)
      integer(kind=byte), intent(inout) :: idxG3(6,*)
      real(kind=wp), intent(inout) :: DEPSA(NASHT,NASHT)

      integer(kind=iwp), parameter :: I1=KIND(idxG3)
C SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
C OR CLOSED-SHELL SCF CASE.

      logical(kind=iwp), external ::  RSV_TSK
      integer(kind=iwp) :: nLev, I, NLEV2, NLEV4, iG3, nTask, ID, iTask,
     &  IND1, IND2, IT1, IU1, LU1, IT2, IU2, LU2, IT3, IU3, IND3, LU3,
     &  LT, IT, IU, LU, IV, LV
      real(kind=wp) :: ESUM, DESUM, OCC

      nLev = SGS%nLev

      ESUM=Zero
      DESUM=Zero
      DO I=1,NLEV
        ESUM=ESUM+ETA(I)
      END DO
C ISCF=1 for closed-shell, =2 for hispin
      OCC=Two
      IF(ISCF == 2) OCC=One

      if (NACTEL == 1 .or. NACTEL == 2) NG3 = 0
      if (NACTEL /= 1) then
        if (NACTEL /= 2) then
          write(u6,*) "I have not implemented for non-standard Psi0, ",
     *      "when A and C subspaces contribute to the energy, ",
     *      "in particular"
          write(u6,*) "I cannot debug, ",
     &      "because I do not know when it happens"
C         call abend
C
          NLEV2=NLEV**2
          NLEV4=NLEV**4

          iG3=0
          nTask=NLEV4
C SVC20100908 initialize the series of tasks
          Call Init_Tsk(ID, nTask)

          do
#ifdef _MOLCAS_MPP_
            IF ((NG3-iG3) < NLEV2) exit
#endif
            IF (.NOT.Rsv_Tsk(ID,iTask)) exit

            IND1=MOD(iTask-1,NLEV2)+1
            IND2=((iTask-IND1)/(NLEV2))+1
            IF(IND2 > IND1) cycle

            IT1=MOD(IND1-1,NASHT)+1
            IU1=(IND1-IT1)/NASHT+1
            LU1=LEVEL(IU1)
            IT2=MOD(IND2-1,NASHT)+1
            IU2=(IND2-IT2)/NASHT+1
            LU2=LEVEL(IU2)

            DO IT3=1,NLEV
             DO IU3=1,NLEV
              IND3=IT3+NASHT*(IU3-1)
              IF(IND3 > IND2) cycle
              LU3=LEVEL(IU3)
C             VAL=G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

C Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Add here the necessary Kronecker deltas times 2-body matrix
C elements and lower, so we get a true normal-ordered density matrix
C element.

C <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
C = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
C -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
C -D(T2,U1)*G2(T1,U2,T3,U3)
C -D(T3,U1)*G2(T2,U2,T1,U3)

C             IF(IT3 == IU2) THEN
C               VAL=VAL-G2(IT1,IU1,IT2,IU3)
C               IF(IT2 == IU1) THEN
C                 VAL=VAL-G1(IT1,IU3)
C               END IF
C             END IF
C             IF(IT2 == IU1) THEN
C               VAL=VAL-G2(IT1,IU2,IT3,IU3)
C             END IF
C             IF(IT3 == IU1) THEN
C               VAL=VAL-G2(IT2,IU2,IT1,IU3)
C             END IF

C VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
              iG3=iG3+1
              idxG3(1,iG3)=INT(iT1,I1)
              idxG3(2,iG3)=INT(iU1,I1)
              idxG3(3,iG3)=INT(iT2,I1)
              idxG3(4,iG3)=INT(iU2,I1)
              idxG3(5,iG3)=INT(iT3,I1)
              idxG3(6,iG3)=INT(iU3,I1)
C             G3(iG3)=VAL
C             F3(iG3)=(ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL
              DESUM = DESUM + OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU1,LU1) = DEPSA(LU1,LU1) - OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU2,LU2) = DEPSA(LU2,LU2) - OCC*G3(iG3)*DF3(iG3)
              DEPSA(LU3,LU3) = DEPSA(LU3,LU3) - OCC*G3(iG3)*DF3(iG3)
              END DO
            END DO

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
          end do

C SVC2010: no more tasks, wait here for the others.
          CALL Free_Tsk(ID)

          NG3=iG3
        end if
        DO IT=1,NASHT
         LT=LEVEL(IT)
         DO IU=1,NASHT
          LU=LEVEL(IU)
C         G2(IT,IT,IU,IU)=G1(IT,IT)*G1(IU,IU)
C         IF(IU == IT) THEN
C          G2(IT,IT,IU,IU)=G2(IT,IT,IU,IU)-G1(IT,IU)
C         ELSE
C          G2(IT,IU,IU,IT)=-G1(IT,IT)
C         END IF
C         F2(IT,IT,IU,IU)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
C         F2(IT,IU,IU,IT)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
          DESUM = DESUM + OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IU)
          DESUM = DESUM + OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IT)
          DO IV=1,NASHT
          LV=LEVEL(IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IV,IU,IU)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IV)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IV,IT)
          END DO
         END DO
        END DO
      end if

      DO IT=1,NASHT
C       G1(IT,IT)=OCC
        LT=LEVEL(IT)
C       F1(IT,IT)=(ESUM*OCC-ETA(LT))*G1(IT,IT)
        DESUM = DESUM + OCC*G1(IT,IT)*DF1(IT,IT)
        Do IU=1, NASHT
          LU=LEVEL(IU)
          DEPSA(LT,LU) = DEPSA(LT,LU) - OCC*G1(IT,IT)*DF1(IT,IU)
        End Do
      END DO

      RETURN
      END
