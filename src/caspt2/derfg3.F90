!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine DERFG3(CI,NCONF,NLEV,NG3,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G1,G2)

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG, VERBOSE
use caspt2_global, only: nbuf1_grad
#ifdef _MOLCAS_MPP_
use caspt2_global, only: iTasks_grad, nTasks_grad
#endif
use Symmetry_Info, only: Mul
use sguga, only: CIS, L2ACT, SGS, EXS
use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
use definitions, only: iwp, wp, u6, RtoB
use caspt2_module, only: nActEl, nSym, STSym, EPSA
use Molcas, only: MxLev
use caspt2_module, only: MxCI
use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk
use Constants, only: Zero, One, Half

implicit none
integer(kind=iwp), intent(in) :: nCONF, NLEV, NG3
real(kind=wp), intent(in) :: CI(nCONF), DG3(NG3), DF3(NG3), G1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV)
real(kind=wp), intent(inout) :: CLAG(NCONF), DG1(NLEV,NLEV), DG2(NLEV,NLEV,NLEV,NLEV), DF1(NLEV,NLEV), DF2(NLEV,NLEV,NLEV,NLEV), &
                                DEPSA(NLEV,NLEV)
integer(kind=iwp) :: I, J, IDX, JDX
integer(kind=iwp) :: IB, IBMN, IBMX, IBUF, NB, NBTOT, IBUF1
integer(kind=iwp) :: IP1, IP2, IP3, IP1MN, IP1MX, IP1I, IP1STA, IP1END, IP3MX, IQ1
integer(kind=iwp) :: IG3, IG3OFF, IG3BK
integer(kind=iwp) :: ISTU, ISVX, ISYZ
integer(kind=iwp) :: IT, IU, IV, IX, IY, IZ
integer(kind=iwp) :: ITLEV, IULEV, IVLEV, IXLEV, IYLEV, IZLEV, IXLEV0
integer(kind=iwp) :: NBUF1, NBUFX, NDTU, NDAB
integer(kind=iwp) :: IOFFSET
integer(kind=iwp) :: ISSG1, ISSG2, ISP1
integer(kind=iwp) :: ITASK, ISUBTASK, ID, NTASKS, NSUBTASKS, MXTASK, MYTASK, MYBUFFER
integer(kind=iwp) :: NSGM1, NSGM2
integer(kind=iwp) :: NTRI1, NTRI2
integer(kind=iwp) :: MEMMAX !, MEMMAX_SAFE
integer(kind=iwp) :: NLEV2
integer(kind=iwp) :: ICSF
real(kind=wp), external :: DDOT_
! translation tables for levels i,j to and from pair indices idx
integer(kind=iwp) :: IJ2IDX(MXLEV,MXLEV)
integer(kind=iwp) :: IDX2IJ(2,MXLEV**2)
integer(kind=iwp) :: ICNJ(MXLEV**2)
integer(kind=iwp) :: IP1_BUF(MXLEV**2)
real(kind=wp), allocatable :: BUF1(:,:), BUF2(:), BUFT(:), BUFD(:), DTU(:,:), DYZ(:), DAB(:,:), BUF3(:), BUF4(:), BUFX(:,:)
integer(kind=iwp), allocatable :: TASKLIST(:,:)
! result buffer, maximum size is the largest possible ip1 range,
! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
!real(kind=wp) BUFR(MXLEV**2)
!integer(kind=iwp) :: LFCDer1,LFCDer2
integer(kind=iwp) :: iTask_loc
real(kind=wp) :: SCAL, ScalG, ScalF
!real(kind=wp) :: tmp,tmp2
logical(kind=iwp) :: first
integer(kind=iwp) :: nMidV

nMidV = CIS%nMidV

! Put in zeroes. Recognize special cases:
if (nlev == 0) return
if (NACTEL == 0) return

! This should not happen, but...
if (NCONF == 0) return

! Here, for regular CAS or RAS cases.

! Special pair index idx2ij allows true RAS cases to be handled:
nlev2 = nlev**2
ntri1 = (nlev2-nlev)/2
ntri2 = (nlev2+nlev)/2
idx = 0
do i=1,nlev-1
  do j=i+1,nlev
    idx = idx+1
    ij2idx(i,j) = idx
    idx2ij(1,idx) = i
    idx2ij(2,idx) = j
    jdx = nlev2+1-idx
    ij2idx(j,i) = jdx
    idx2ij(1,jdx) = j
    idx2ij(2,jdx) = i
  end do
end do
do i=1,nlev
  idx = ntri1+i
  ij2idx(i,i) = idx
  idx2ij(1,idx) = i
  idx2ij(2,idx) = i
end do
do idx=1,nlev2
  i = idx2ij(1,idx)
  j = idx2ij(2,idx)
  jdx = ij2idx(j,i)
  icnj(idx) = jdx
end do

! Correction to G3: It is now <0| E_tu E_vx E_yz |0>
! Similar for F3 values.
! These corrections are already done in CLagDXA_FG3 and CLagDXC_FG3

! Correction to F2: Some values not computed follow from symmetry
do ip1=1,nlev2-1
  itlev = idx2ij(1,ip1)
  iulev = idx2ij(2,ip1)
  it = L2ACT(itlev)
  iu = L2ACT(iulev)
  do ip3=ip1+1,nlev2
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    SCAL = DF2(iy,iz,it,iu)+DF2(it,iu,iy,iz)
    DF2(it,iu,iy,iz) = Zero
    DF2(iy,iz,it,iu) = Scal
  end do
end do
!-SVC20100310: took some spurious mirroring of F2 values out
!-of the loops and put them here, after the parallel section has
!-finished, so that GAdGOP works correctly.
do ip1=ntri2+1,nlev2
  itlev = idx2ij(1,ip1)
  iulev = idx2ij(2,ip1)
  it = L2ACT(itlev)
  iu = L2ACT(iulev)
  do ip3=ntri1+1,ip1
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    Scal = DF2(iz,iy,iu,it)+DF2(it,iu,iy,iz)
    DF2(it,iu,iy,iz) = Zero
    DF2(iz,iy,iu,it) = Scal
  end do
end do
! Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
do iz=1,nlev
  do iy=1,nlev
    do iu=1,nlev
      do it=1,nlev
        DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)-DF2(it,iu,iy,iz)*(EPSA(iu)+EPSA(iy))
        !! DEPSA-related operations in the preparation step
        !! are done only on the master node
#       ifdef _MOLCAS_MPP_
        if (king()) then
#       endif
          do ix=1,nlev
            DEPSA(iu,ix) = DEPSA(iu,ix)-DF2(it,iu,iy,iz)*G2(it,ix,iy,iz)
            DEPSA(ix,iy) = DEPSA(ix,iy)-DF2(it,iu,iy,iz)*G2(it,iu,ix,iz)
          end do
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end do
    end do
  end do
end do
do iz=1,nlev
  do iu=1,nlev
    do it=1,nlev
      DG1(it,iz) = DG1(it,iz)-EPSA(iu)*DF2(it,iu,iu,iz)
      DF1(it,iz) = DF1(it,iz)-DF2(it,iu,iu,iz)
#     ifdef _MOLCAS_MPP_
      if (king()) then
#     endif
        do iy=1,nlev
          DEPSA(iu,iy) = DEPSA(iu,iy)-G1(it,iz)*DF2(it,iu,iy,iz)
        end do
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end do
  end do
end do

! Correction to G2: Some values not computed follow from symmetry
do ip1=1,nlev2-1
  itlev = idx2ij(1,ip1)
  iulev = idx2ij(2,ip1)
  it = L2ACT(itlev)
  iu = L2ACT(iulev)
  do ip3=ip1+1,nlev2
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    SCAL = DG2(iy,iz,it,iu)+DG2(it,iu,iy,iz)
    DG2(it,iu,iy,iz) = Zero
    DG2(iy,iz,it,iu) = Scal
  end do
end do
!-SVC20100310: took some spurious mirroring of G2 values out
!-of the loops and put them here, after the parallel section has
!-finished, so that GAdSUM works correctly.
do ip1=ntri2+1,nlev2
  itlev = idx2ij(1,ip1)
  iulev = idx2ij(2,ip1)
  it = L2ACT(itlev)
  iu = L2ACT(iulev)
  do ip3=ntri1+1,ip1
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    Scal = DG2(iz,iy,iu,it)+DG2(it,iu,iy,iz)
    DG2(it,iu,iy,iz) = Zero
    DG2(iz,iy,iu,it) = Scal
  end do
end do
! Correction to G2: It is now = <0| E_tu E_yz |0>
do iu=1,nlev
  do iz=1,nlev
    do it=1,nlev
      DG1(it,iz) = DG1(it,iz)-DG2(it,iu,iu,iz)
    end do
  end do
end do
! Additional correction terms for F1
do iT=1,NLEV
  do iU=1,NLEV
    !! With the improved algorithm, symmetrizing the DF1
    !! contribution to DF1 here is somehow required
    DG1(iT,iU) = DG1(iT,iU)-(DF1(iT,iU)+DF1(iU,iT))*EPSA(iU)*Half
  end do
end do
#ifdef _MOLCAS_MPP_
if (king()) then
#endif
  do iT=1,NLEV
    do iU=1,NLEV
      do iV=1,NLEV
        do iX=1,nLEV
          DEPSA(iV,iX) = DEPSA(iV,iX)+G2(iT,iU,iV,iX)*DF1(iT,iU)
        end do
      end do
    end do
  end do
#ifdef _MOLCAS_MPP_
end if
#endif

call CLagSym(nLev,DG1,DG2,DF1,DF2,2)

call mma_MaxDBLE(memmax)

! Use *almost* all remaining memory:
!memmax_safe = int(dble(memmax)*0.95D0)

! Buffers to compute CI expansion vectors into:
! <Psi0|E_ip1 | E_ip2 E_ip3|Psi0>
! buf1: bra buffer with E_ip1 excitations of Psi0
!       holds multiple CI vectors (allowed by memory)
! buf2: ket buffer for an E_ip3 excitation of Psi0
! buft: ket buffer for an E_ip2 excitation of E_ip3|Psi0>
! bufd: diagonal matrix elements to compute the F matrix
nbuf1 = nbuf1_grad
ndtu = nbuf1
ndab = nbuf1
nbufx = nlev
call mma_allocate(BUF1,MXCI,NBUF1,Label='BUF1')
call mma_allocate(BUF2,MXCI,Label='BUF2')
call mma_allocate(BUFT,MXCI,Label='BUFT')
call mma_allocate(BUFD,MXCI,Label='BUFD')

call mma_allocate(DTU,MXCI,NDTU,Label='DTU')
call mma_allocate(DYZ,MXCI,Label='DYZ')
call mma_allocate(DAB,MXCI,NDAB,Label='DAB')
call mma_allocate(BUF3,MXCI,Label='BUF3')
call mma_allocate(BUF4,MXCI,Label='BUF4')

call mma_allocate(BUFX,MXCI,NBUFX,Label='BUFX')

!-SVC20100301: calculate maximum number of tasks possible
MXTASK = (NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
call mma_allocate(TaskList,mxTask,4,LABEL='TaskList')

if (iPrGlb >= VERBOSE) then
  write(u6,*)
  write(u6,'(2X,A)') 'Constructing derivatives of G3/F3'
  write(u6,'(2X,A,F16.9,A)') ' memory avail: ',(memmax*RtoB)/1.0e+09_wp,' GB'
  write(u6,'(2X,A,F16.9,A)') ' memory used:  ',(((3*nbuf1+6+nlev)*MXCI)*RtoB)/1.0e+09_wp,' GB'
  call xFlush(u6)
end if
!call TIMING(CPTF10,CPE,TIOTF10,TIOE)
!CPUT = CPTF10-CPTF0
!WALLT = TIOTF10-TIOTF0
!write(u6,*) 'PREP    : CPU/WALL TIME=',cput,wallt

iG3OFF = 0
iTask_loc = 1
first = .true.
! A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.
! This also allows precomputing the Hamiltonian (H0) diagonal elements.
do issg1=1,nsym
  isp1 = Mul(issg1,STSYM)
  nsgm1 = CIS%ncsf(issg1)
  !!BufD_I = \sum_t <I|E_{tt}|I>*f_{tt}
  call H0DIAG_CASPT2(ISSG1,BUFD,nsgm1,CIS%NOW,CIS%IOW,nMidV)

  !-SVC20100301: calculate number of larger tasks for this symmetry, this
  !-is basically the number of buffers we fill with SG_Epq_Psi vectors.
  iTask = 1
  ibuf1 = 0
  do ip1=1,nlev2
    itlev = idx2ij(1,ip1)
    iulev = idx2ij(2,ip1)
    istu = Mul(SGS%ism(itlev),SGS%ism(iulev))
    if (istu == isp1) then
      ibuf1 = ibuf1+1
      ip1_buf(ibuf1) = ip1
      if (ibuf1 == 1) TaskList(iTask,1) = ip1
    end if
    if ((ibuf1 == nbuf1) .or. ((ibuf1 > 0) .and. ((ip1 == ntri2) .or. (ip1 == nlev2)))) then
      TaskList(iTask,2) = ip1_buf(ibuf1)
      TaskList(iTask,3) = ibuf1
      iTask = iTask+1
      ibuf1 = 0
    end if
  end do
  nTasks = iTask
  !write(u6,*) 'nTasks = ',nTasks
  if (ibuf1 == 0) nTasks = nTasks-1
  !-SVC20100309: calculate number of inner loop iteration tasks.
  iOffSet = 0
  do iTask=1,nTasks
    TaskList(iTask,4) = iOffSet
    ip1sta = TaskList(iTask,1)
    ip1end = TaskList(iTask,2)
    ip3mx = ntri2
    if (ip1end <= ntri2) ip3mx = ip1end
    if (ip1sta > ntri2) ip3mx = ntri1
    !write(u6,*) 'iTask = ',iTask
    !write(u6,*) 'start,end=',ip1sta,ip1end
    !write(u6,*) 'ip3mx = ',ip3mx
    !-SVC20100309: Currently -we are going to limit this to the ip3-loop and
    !-leave the ip2-loop intact.  This was based on the large overhead which
    !-was observed for a very large number of small tasks.
    !iOffSet = iOffSet+ip3mx*ntri2-((ip3mx**2-ip3mx)/2)
    iOffSet = iOffSet+ip3mx
  end do
  nSubTasks = iOffSet
  !write(u6,*) 'nSubTasks = ',nSubTasks

  if (iPrGlb >= VERBOSE) then
    write(u6,'(2X,A,I3,A,I6)') 'Sym: ',issg1,', #Tasks: ',nSubTasks
    call xFlush(u6)
  end if

  if (iPrGlb >= DEBUG) then
    if (nSubTasks > 0) then
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') 'task ID ',' ip1 range  ','ip3 ','#elements'
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
      call xFlush(u6)
    end if
  end if

  !-SVC20100301: initialize the series of subtasks
  call Init_Tsk(ID,nSubTasks)

  myBuffer = 0

  !! loop start
  do
    !-SVC20100908: first check: can I actually do any task?
    !if (NG3-iG3OFF < nbuf1*ntri2) goto 501
    !-SVC20100831: initialize counter for offset into G3
    !-SVC20100302: BEGIN SEPARATE TASK EXECUTION
    !write(u6,*) rsv_tsk(id,isubtask)
#   ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      !! do the same tasks here and in mkfg3
      iSubTask = iTasks_grad(iTask_loc)
      if (iSubTask == 0) exit
    else
#   endif
      if (.not. Rsv_Tsk(ID,iSubTask)) exit
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    myTask = nTasks
    do iTask=1,nTasks
      iBuf = iSubTask-TaskList(iTask,4)
      if (iBuf <= 0) then
        myTask = iTask-1
        exit
      end if
    end do
    iTask = myTask

    iOffSet = TaskList(iTask,4)

    !-SVC20100310: one task handles a range of ip1 values
    !-that are in the buffer and one ip3 value, for which
    !-a loop over ip2 values is then executed.
    ip1sta = TaskList(iTask,1)
    ip1end = TaskList(iTask,2)
    ip3 = iSubTask-iOffSet

    !-SVC20100301: fill the buffer with sigma vectors if they
    !-have not been computed yet, else just get the number of
    !-sigma vectors in the buffer.
    !write(u6,*) 'myBuffer,iTask = ',myBuffer,iTask
    if (myBuffer /= iTask) then
      !! Compute left derivative and DEPSA contributions before
      !! the TASK is completely switched
      if (.not. first) call LEFT_DEPSA()
      ibuf1 = 0
      do ip1i=ip1sta,ip1end
        itlev = idx2ij(1,ip1i)
        iulev = idx2ij(2,ip1i)
        istu = Mul(SGS%ism(itlev),SGS%ism(iulev))
        it = L2ACT(itlev)
        iu = L2ACT(iulev)
        if (istu == isp1) then
          ibuf1 = ibuf1+1
          ip1_buf(ibuf1) = ip1i
          BUF1(1:nsgm1,ibuf1) = Zero
          call SG_Epq_Psi(SGS,CIS,EXS,IULEV,ITLEV,One,STSYM,CI,BUF1(:,ibuf1))
        end if
      end do
      myBuffer = iTask
      DTU(1:MXCI,1:ibuf1) = Zero
      DAB(1:MXCI,1:ibuf1) = Zero
      if (first) first = .false.
    else
      ibuf1 = TaskList(iTask,3)
    end if
    !-SVC20100301: necessary batch of sigma vectors is now in the buffer

    ! The ip1 buffer could be the same on different processes
    ! so only compute the G1 contribution when ip3 is 1, as
    ! this will only be one task per buffer.
    if ((issg1 == STSYM) .and. (ip3 == 1)) then
      !! buf1 = <Psi0|E_ip1|I>
      !! <0|E_{tu}I> = <I|E_{ut}|0>
      !write(u6,*) 'ib loop'
      do ib=1,ibuf1
        idx = ip1_buf(ib)
        itlev = idx2ij(1,idx)
        iulev = idx2ij(2,idx)
        it = L2ACT(itlev)
        iu = L2ACT(iulev)
        !write(u6,'('itlev,iulev,it,iu = ',4i3)') itlev,iulev,it,iu
        !! DG1 contribution
        SCAL = DG1(iT,iU)+DG1(iT,iU)
        CLag(1:nsgm1) = CLag(1:nsgm1)+SCAL*BUF1(1:nsgm1,ib)

        !! left derivative of DF1
        do icsf=1,nsgm1
          DTU(icsf,ib) = DTU(icsf,ib)+DF1(it,iu)*BUFD(icsf)*CI(icsf)
        end do
        !! right derivative of DF1
        do icsf=1,nsgm1
          CLag(icsf) = CLag(icsf)+DF1(it,iu)*BUF1(icsf,ib)*BUFD(icsf)
        end do
        !G1(it,iu) = DDOT_(nsgm1,ci,1,work(lto),1)
        !if (IFF /= 0) then
        !  F1sum = Zero
        !  do i=1,nsgm1
        !    F1sum = F1sum+CI(i)*work(lto-1+i)*bufd(i)
        !  end do
        !  F1(it,iu) = F1sum-EPSA(iu)*G1(it,iu)
        !end if
      end do
    end if

    !ip3mx = ntri2
    !if (ip1end <= ntri2) ip3mx = ip1end
    !if (ip1sta > ntri2) ip3mx = ntri1
    !-SVC20100309: loop over ip3, ip2
    !do ip3=1,ip3mx

    !-SVC20100309: PAM's magic formula
    !iCnt = iSubTask-iOffSet
    !ip3 = int(dble(ntri2)+1.5D0-sqrt((dble(ntri2)+0.5d0)**2-2*iCnt+0.000001D0))
    !ip2 = iCnt-((ip3-1)*ntri2-((ip3-1)*(ip3-2))/2 )+ip3-1

    !-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

    !write(u6,*) 'ip3 = ',ip3
    !CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
    iq1 = icnj(ip3)
    ! The indices corresponding to pair index p3:
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    isyz = Mul(SGS%ism(iylev),SGS%ism(izlev))
    issg2 = Mul(isyz,STSYM)
    nsgm2 = CIS%ncsf(issg2)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    buf2(1:nsgm2) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,IYLEV,IZLEV,One,STSYM,CI,BUF2)
    DYZ(1:nsgm1) = Zero
    if (issg2 == issg1) then
      buf3(1:nsgm2) = Zero
      do ib=1,ibuf1
        idx = ip1_buf(ib)
        itlev = idx2ij(1,idx)
        iulev = idx2ij(2,idx)
        it = L2ACT(itlev)
        iu = L2ACT(iulev)

        ScalG = DG2(iT,iU,iY,iZ)
        ScalF = DF2(iT,iU,iY,iZ)
        if ((ScalG == Zero) .and. (ScalF == Zero)) cycle

        !! left derivative
        do icsf=1,nsgm1
          BUFT(icsf) = ScalG*BUF2(icsf)+ScalF*BUF2(icsf)*BUFD(icsf)
        end do
        DTU(1:nsgm1,ib) = DTU(1:nsgm1,ib)+BUFT(1:nsgm1)

        !! right derivative
        buft(1:mxci) = buf1(1:mxci,ib)
        do icsf=1,nsgm1
          BUF3(icsf) = BUF3(icsf)+ScalG*BUFT(icsf)+ScalF*BUFT(icsf)*BUFD(icsf)
        end do

        !! For DEPSA
        DAB(1:nsgm1,ib) = DAB(1:nsgm1,ib)+ScalF*BUF2(1:nsgm1)
      end do
      !! Save for Eyz
      DYZ(1:nsgm2) = DYZ(1:nsgm2)+BUF3(1:nsgm2)
    end if
    nbtot = 0

    !! Prepare for DEPSA for the -epsa(iv) term with square
    iG3bk = iG3OFF
    do ixlev0=1,nlev
      do ivlev=1,nlev
        BUFX(1:nsgm1,ivlev) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,IVLEV,IXLEV0,One,STSYM,BUF2,BUFX(1,ivlev))
      end do
      iG3OFF = iG3bk
      do ip2=ip3,ntri2
        ivlev = idx2ij(1,ip2)
        ixlev = idx2ij(2,ip2)
        isvx = Mul(SGS%ism(ivlev),SGS%ism(ixlev))
        iv = L2ACT(ivlev)
        ix = L2ACT(ixlev)
        if (isvx /= Mul(issg1,issg2)) cycle
        !! <I|EvxEyz|0>
        if (IXLEV == IXLEV0) BUFT(1:nsgm1) = BUFX(1:nsgm1,IVLEV)
        !-----------
        ! Max and min values of index p1:
        ip1mx = ntri2
        if (ip3 <= ntri1) then
          ip1mx = nlev2
          if (ip2 > ntri1) ip1mx = iq1
        end if
        ip1mn = max(ip2,ip1sta)
        ip1mx = min(ip1mx,ip1end)
        ! The corresponding locations in the Sgm1 buffer:
        ibmn = 999999
        ibmx = -999999
        do ib=ibuf1,1,-1
          ip1 = ip1_buf(ib)
          if (ip1 >= ip1mn) ibmn = ib
        end do
        do ib=1,ibuf1
          ip1 = ip1_buf(ib)
          if (ip1 <= ip1mx) ibmx = ib
        end do
        nb = ibmx-ibmn+1
        if (nb <= 0) cycle
        if (ixlev /= ixlev0) then
          iG3OFF = iG3OFF+nb
          cycle
        end if

        ! ----- left derivative

        do icsf=1,nsgm1
          ! BUF3 = (<I|Ett|I>-EPSA(V))*<I|EvxEyz|0> = <I|fEvxEyz|0>
          buf3(icsf) = (bufd(icsf)-epsa(iv))*buft(icsf)
        end do
        do ib=1,nb
          iG3 = iG3OFF+ib
          idx = ip1_buf(ibmn-1+ib)

          !! <I|EvxEyz|0>*Dtuvxyz + <I|fEvxEyz|0>*Ftuvxyz
          DTU(1:nsgm1,ibmn+ib-1) = DTU(1:nsgm1,ibmn+ib-1)+DG3(iG3)*BUFT(1:nsgm1)+DF3(iG3)*BUF3(1:nsgm1)
          !! DEPSA of the BUFD term
          DAB(1:nsgm1,ibmn+ib-1) = DAB(1:nsgm1,ibmn+ib-1)+DF3(iG3)*BUFT(1:nsgm1)
        end do

        ! ----- right derivative

        BUF3(1:nsgm1) = Zero
        BUF4(1:nsgm1) = Zero
        !! right derivative (1):
        !! <0|Etu|I>*Dtuvxyz and <0|Etu|I>*Ftuvxyz
        !! <0|EtuEvxEyz|I> -> <I|EzyExvEut|0>
        do ib=1,nb
          iG3 = iG3OFF+ib
          idx = ip1_buf(ibmn-1+ib)

          !! BUF3 = <0|Etu|I>*Dtuvxyz
          BUF3(1:nsgm1) = BUF3(1:nsgm1)+DG3(iG3)*BUF1(1:nsgm1,ibmn+ib-1)
          !! BUFC = <0|Etu|I>*Ftuvxyz
          BUF4(1:nsgm1) = BUF4(1:nsgm1)+DF3(iG3)*BUF1(1:nsgm1,ibmn+ib-1)
        end do

        !! DEPSA of the -EPSA(iv) term
        call DGEMV_('T',nsgm1,NLEV,-One,BUFX,mxci,buf4,1,One,DEPSA(1,IVLEV),1)

        !! Scale the DF3 contribution with the diagonal Fock
        !! and add to the DG3 contribution
        do icsf=1,nsgm1
          buf3(icsf) = buf3(icsf)+buf4(icsf)*(bufd(icsf)-epsa(iv))
        end do
        !! right derivative (2): <0|EtuEvx|I>*Dtuvxyz
        call SG_Epq_Psi(SGS,CIS,EXS,IXLEV,IVLEV,One,STSYM,BUF3,DYZ)

        iG3OFF = iG3OFF+nb
        nbtot = nbtot+nb
      end do !! end of ip2 loop
    end do !! end of ixlev0 loop

    !! Complete the right derivative contribution:
    !! <0|EtuEyz|I> and <0|EtuEvxEyz|I>
    call SG_Epq_Psi(SGS,CIS,EXS,IZLEV,IYLEV,One,STSYM,DYZ,CLAG)

    if (iPrGlb >= DEBUG) then
      write(u6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)') iSubTask,ip1sta,ip1end,ip3,nbtot
      call xFlush(u6)
    end if
    iTask_loc = iTask_loc+1

    !SVC: The master node now continues to only handle task scheduling,
    !     needed to achieve better load balancing. So it exits from the task
    !     list.  It has to do it here since each process gets at least one
    !     task.

  end do
  !-SVC20100301: end of the task

  !! Final (the last task) left derivative and DEPSA contributions
  call LEFT_DEPSA()

  !-SVC20100302: no more tasks, wait here for the others, then proceed
  ! with next symmetry
  call Free_Tsk(ID)

  if (iPrGlb >= DEBUG) then
    if (nSubTasks > 0) write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
  end if

! End of sectioning loop over symmetry of Sgm1 wave functions.
end do

#ifdef _MOLCAS_MPP_
if (is_real_par() .and. (iTask_loc-1 /= nTasks_grad)) then
  write(u6,*)
  write(u6,*) 'Somehow, the number of tasks in mkfg3 and derfg3 is not consistent...'
  write(u6,*) 'probably, bug'
  write(u6,*) '# of tasks in  mkfg3 = ',nTasks_grad
  write(u6,*) '# of tasks in derfg3 = ',iTask_loc-1
  call abend()
end if
#endif

call mma_deallocate(TASKLIST)
! free CI buffers
call mma_deallocate(BUF1)
call mma_deallocate(BUF2)
call mma_deallocate(BUFT)
call mma_deallocate(BUFD)

call mma_deallocate(DTU)
call mma_deallocate(DYZ)
call mma_deallocate(DAB)
call mma_deallocate(BUF3)
call mma_deallocate(BUF4)

call mma_deallocate(BUFX)

contains

subroutine LEFT_DEPSA()

  integer(kind=iwp) :: IALEVloc, IBLEVloc, ibloc

  do ibloc=1,ibuf1
    idx = ip1_buf(ibloc)
    itlev = idx2ij(1,idx)
    iulev = idx2ij(2,idx)
    !! left derivative
    call SG_Epq_Psi(SGS,CIS,EXS,ITLEV,IULEV,One,STSYM,DTU(1,ibloc),CLAG)
    !! the rest is DEPSA contribution
    do IALEVloc=1,NLEV
      do IBLEVloc=1,NLEV
        BUF2(:) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,IALEVloc,IBLEVloc,One,STSYM,DAB(1,ibloc),BUF2)
        DEPSA(IALEVloc,IBLEVloc) = DEPSA(IALEVloc,IBLEVloc)+DDot_(nsgm1,BUF1(1,IBloc),1,BUF2,1)
      end do
    end do
  end do

end subroutine LEFT_DEPSA

end subroutine DERFG3
