!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! Procedure for computing 1-body, 2-body, and 3-body
! density elements with active indices only,
! and related matrices obtained from contractions with
! the diagonal one-electron Hamiltonian.
!
! Out: density matrices, denoted here G1, G2 and G3.
! Storage: G1 and G2 are simple two- and four-index arrays, and
! includes also such zeroes that are implied by symmetry.
! But G3 is quite large, and while it is stored with zeroes, it
! is made more compact by calculating only the minimum amount of unique values and storing the
! active indices in the array idxG3.  Later, the full matrix can be restored
! on the fly by using the full permutational symmetry.
!
! ----
!
!.NN15 NOTE:
! In DMRG-CASPT2, full G3 is computed from BLOCK code
! and F3 is approximated by cumulant reconstruction scheme from G1, G2, and G3 (cu4),
! to avoid computation of 4-particle density matrix.
! However, because this introduces instability of CASPT2 calculation
! (lots of negative denominators appear), relatively large IPEA and imaginary shifts
! are required to converge CASPT2 iteration.
!

#include "compiler_features.h"
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined _DMRG_

subroutine MKFG3DM(mkF,G1,F1,G2,F2,G3,F3,idxG3,NLEV,mG3)

use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use Symmetry_Info, only: Mul
use caspt2_global, only: iPrGlb
use PrintLevel, only: DEBUG, VERBOSE
use sguga, only: CIS, L2ACT, SGS
use Molcas, only: MxLev
use caspt2_module, only: MxCI, nActEl, nG1, nG2, nG3, nSym, STSym
#ifdef _DMRG_
use caspt2_module, only: DMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, byte, RtoB

implicit none
logical(kind=iwp), intent(in) :: mkF
integer(kind=iwp), intent(in) :: NLEV, mG3
real(kind=wp), intent(out) :: G1(NLEV,NLEV), F1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV), F2(NLEV,NLEV,NLEV,NLEV), G3(mG3), F3(mG3)
integer(kind=byte), intent(out) :: idxG3(6,mG3)
integer(kind=iwp) :: I, IB, IBMN, IBMX, IBUF, IBUF1, ID, IDX, IG3, IG3OFF, IOFFSET, IP1, IP1END, IP1I, IP1MN, IP1MX, IP1STA, IP2, &
                     IP3, IP3MX, IQ1, ISP1, ISSG1, ISSG2, ISTU, ISUBTASK, ISVX, ISYZ, IT, ITASK, ITLEV, IU, IULEV, IV, IVLEV, IX, &
                     IXLEV, IY, IYLEV, IZ, IZLEV, J, JDX, MEMMAX, MEMMAX_SAFE, MXTASK, MYBUFFER, MYTASK, NB, NBTOT, NBUF1, NCI, &
                     NLEV2, NSUBTASKS, NTASKS, NTRI1, NTRI2
real(kind=wp) :: DF1, DF2, DF3, DG1, DG2, DG3
#ifdef _ENABLE_BLOCK_DMRG_
integer(kind=iwp) :: NLEV4
real(kind=wp), allocatable :: G3Tmp(:)
#endif
integer(kind=iwp), allocatable :: ICNJ(:), IDX2IJ(:,:), IJ2IDX(:,:), IP1_BUF(:), TaskList(:,:)
real(kind=wp), allocatable :: BUF1(:,:), BUF2(:), BUFD(:), BUFT(:)
real(kind=wp), external :: DDOT_, DNRM2_

! IJ2IDX, IDX2IJ, ICNJ, IP1_BUF: translation tables for levels i,j to and from pair indices idx

! Put in zeroes. Recognize special cases:
if (nlev == 0) return

G1(:,:) = Zero
G2(:,:,:,:) = Zero
call DCOPY_(NG3,[Zero],0,G3,1)
if (mkF) then
  F1(:,:) = Zero
  F2(:,:,:,:) = Zero
  call DCOPY_(NG3,[Zero],0,F3,1)
end if

if (NACTEL == 0) return

NCI = CIS%NCSF(STSYM)
! This should not happen, but...
if (NCI == 0) return

! Here, for regular CAS or RAS cases.

! Special pair index idx2ij allows true RAS cases to be handled:
nlev2 = nlev**2
ntri1 = (nlev2-nlev)/2
ntri2 = (nlev2+nlev)/2
call mma_allocate(ij2idx,nlev,nlev,Label='ij2idx')
call mma_allocate(idx2ij,2,nlev2,Label='idx2ij')
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
call mma_allocate(icnj,nlev2,Label='icnj')
do idx=1,nlev2
  i = idx2ij(1,idx)
  j = idx2ij(2,idx)
  jdx = ij2idx(j,i)
  icnj(idx) = jdx
end do
call mma_deallocate(ij2idx)

call mma_MaxDBLE(memmax)

! Use *almost* all remaining memory:
memmax_safe = int(real(memmax,kind=wp)*0.95_wp)

! Buffers to compute CI expansion vectors into:

nbuf1 = max(1,min(nlev2,(memmax_safe-3*mxci)/mxci)) ! -> 1 w/ DMRG?
call mma_allocate(BUF1,MXCI,NBUF1,LABEL='BUF1')
call mma_allocate(BUF2,MXCI,LABEL='BUF2')
call mma_allocate(BUFT,MXCI,LABEL='BUFT')
call mma_allocate(BUFD,MXCI,LABEL='BUFD')

!-SVC20100301: calculate maximum number of tasks possible
MXTASK = (NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
call mma_allocate(TaskList,mxTask,4,LABEL='TaskList')

if (iPrGlb >= VERBOSE) then
  write(u6,*)
  write(u6,'(2X,A)') 'Constructing G3/F3'
  write(u6,'(2X,A,F16.9,A)') ' memory avail: ',(memmax*RtoB)*1.0e-9_wp,' GB'
  write(u6,'(2X,A,F16.9,A)') ' memory used:  ',(((nbuf1+3)*MXCI)*RtoB)*1.0e-9_wp,' GB'
end if

call mma_allocate(ip1_buf,nlev2,Label='ip1_buf')

iG3OFF = 0
! A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.

do issg1=1,nsym
  isp1 = Mul(issg1,stsym)
  !nsgm1 = CIS%ncsf(issg1)
  !call H0DIAG_CASPT2(ISSG1,BUFD,nsgm1,NOW1,IOW1,NMIDV)

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
    !-SVC20100309: Currently -we are going to limit this to the ip3-loop and
    !-leave the ip2-loop intact.  This was based on the large overhead which
    !-was observed for a very large number of small tasks.
    !iOffSet = iOffSet+ip3mx*ntri2-((ip3mx**2-ip3mx)/2)
    iOffSet = iOffSet+ip3mx
  end do
  nSubTasks = iOffSet

  if (iPrGlb >= VERBOSE) write(u6,'(2X,A,I3,A,I6)') 'Sym: ',issg1,', #Tasks: ',nSubTasks

  if (iPrGlb >= DEBUG) then
    if (nSubTasks > 0) then
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') 'task ID ',' ip1 range  ','ip3 ','#elements'
      write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
    end if
  end if

  !-SVC20100301: initialize the series of subtasks
  call Init_Tsk(ID,nSubTasks)

  myBuffer = 0

  do
    !-SVC20100908: first check: can I actually do any task?
    if ((NG3-iG3OFF) < nbuf1*ntri2) exit
    !-SVC20100831: initialize counter for offset into G3
    !-SVC20100302: BEGIN SEPARATE TASK EXECUTION
    if (.not. Rsv_Tsk(ID,iSubTask)) exit

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
    if (myBuffer /= iTask) then
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
          !call dcopy_(nsgm1,[Zero],0,BUF1(:,ibuf1),1)
          !call SG_Epq_Psi(IULEV,ITLEV,One,STSYM,CI,BUF1(:,ibuf1))
        end if
      end do
      myBuffer = iTask
    else
      ibuf1 = TaskList(iTask,3)
    end if
    !-SVC20100301: necessary batch of sigma vectors is now in the buffer

    ! The ip1 buffer could be the same on different processes
    ! so only compute the G1 contribution when ip3 is 1, as
    ! this will only be one task per buffer.
    !if ((issg1 == stsym) .and. (ip3 == 1)) then
    !  do ib=1,ibuf1
    !    idx = ip1_buf(ib)
    !    itlev = idx2ij(1,idx)
    !    iulev = idx2ij(2,idx)
    !    it = L2ACT(itlev)
    !    iu = L2ACT(iulev)
    !    G1(it,iu) = DDOT_(nsgm1,ci,1,BUF1(:,ib),1)
    !    if (mkF) then
    !      F1sum = Zero
    !      do i=1,nsgm1
    !        F1sum = F1sum+CI(i)*BUF1(i,ib)*bufd(i)
    !      end do
    !      F1(it,iu) = F1sum-EPSA(iu)*G1(it,iu)
    !    end if
    !  end do
    !end if

    !ip3mx = ntri2
    !if (ip1end <= ntri2) ip3mx = ip1end
    !if (ip1sta > ntri2) ip3mx = ntri1
    !-SVC20100309: loop over ip3, ip2
    !do ip3=1,ip3mx

    !-SVC20100309: PAM's magic formula
    !iCnt = iSubTask-iOffSet
    !ip3 = int(real(ntri2,kind=wp)+OneHalf-sqrt((real(ntri2,kind=wp)+Half)**2-2*iCnt+1.0e-6_wp))
    !ip2 = iCnt-((ip3-1)*ntri2-((ip3-1)*(ip3-2))/2 )+ip3-1

    !-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

    ! NN.14 TODO:
    ! To avoid storing full G3(tmp) in memory, need to store
    ! G3(:,:,it,iu,iy,iz) loaded from disk, for each process...

    iq1 = icnj(ip3)
    ! The indices corresponding to pair index p3:
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    isyz = Mul(SGS%ism(iylev),SGS%ism(izlev))
    issg2 = Mul(isyz,stsym)
    !nsgm2 = CIS%ncsf(issg2)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    !call dcopy_(nsgm2,Zero,0,BUF2,1)
    !call SG_Epq_Psi(IYLEV,IZLEV,One,STSYM,CI,BUF2)
    !if (issg2 == issg1) then
    !  do ib=1,ibuf1
    !    idx = ip1_buf(ib)
    !    itlev = idx2ij(1,idx)
    !    iulev = idx2ij(2,idx)
    !    it = L2ACT(itlev)
    !    iu = L2ACT(iulev)
    !    G2(it,iu,iy,iz) = DDOT_(nsgm1,BUF2,1,BUF1(:,ib),1)
    !    if (mkF) then
    !      F2sum = Zero
    !      do i=1,nlev
    !        F2sum = F2sum+BUF2(i)*bufd(i)*BUF1(i,ib)
    !      end do
    !      F2(it,iu,iy,iz) = F2sum
    !    end if
    !  end do
    !end if
    nbtot = 0
    do ip2=ip3,ntri2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx = Mul(SGS%ism(ivlev),SGS%ism(ixlev))
      iv = L2ACT(ivlev)
      ix = L2ACT(ixlev)
      if (isvx == Mul(issg1,issg2)) then
        !call dcopy_(nsgm1,[Zero],0,BUFT,1)
        !call SG_Epq_Psi(IVLEV,IXLEV,One,ISSG2,BUF2,BUFT)
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
        if (nb > 0) then

          !-----------
          ! Contract the Sgm1 wave functions with the Tau wave function.
          !call DGEMV_('T',nsgm1,nb,One,BUF1(:,ibmn),mxci,buft,1,Zero,bufr,1)
          ! and distribute this result into G3:
          !call DCOPY_(nb,bufr,1,G3(iG3OFF+1),1)
          ! and copy the active indices into idxG3:
          do ib=1,nb
            iG3 = iG3OFF+ib
            idx = ip1_buf(ibmn-1+ib)
            itlev = idx2ij(1,idx)
            iulev = idx2ij(2,idx)
            iT = l2act(itlev)
            iU = l2act(iulev)
            idxG3(1,iG3) = int(iT,kind=byte)
            idxG3(2,iG3) = int(iU,kind=byte)
            idxG3(3,iG3) = int(iV,kind=byte)
            idxG3(4,iG3) = int(iX,kind=byte)
            idxG3(5,iG3) = int(iY,kind=byte)
            idxG3(6,iG3) = int(iZ,kind=byte)
          end do
          !if (mkF) then
          !  ! Elementwise multiplication of Tau with H0 diagonal - EPSA(IV):
          !  do icsf=1,nsgm1
          !    buft(icsf) = (bufd(icsf)-epsa(iv))*buft(icsf)
          !  end do
          !  ! so Tau is now = Sum(eps(w)*E_vxww) Psi. Contract and distribute:
          ! call DGEMV_('T',nsgm1,nb,One,BUF1(:ibmn),mxci,buft,1,Zero,bufr,1)
          ! call dcopy_(nb,bufr,1,F3(iG3OFF+1),1)
          !end if
          iG3OFF = iG3OFF+nb
          nbtot = nbtot+nb
        end if
      end if
    end do
    !end do

    if (iPrGlb >= DEBUG) write(u6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)') iSubTask,ip1sta,ip1end,ip3,nbtot

    !SVC: The master node now continues to only handle task scheduling,
    !     needed to achieve better load balancing. So it exits from the task
    !     list.  It has to do it here since each process gets at least one
    !     task.

    !-SVC20100301: end of the task
  end do

  !-SVC20100302: no more tasks, wait here for the others, then proceed
  ! with next symmetry
  call Free_Tsk(ID)

  if (iPrGlb >= DEBUG) then
    if (nSubTasks > 0) write(u6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)') '--------','------------','----','---------'
  end if

! End of sectioning loop over symmetry of Sgm1 wave functions.
end do
!-SVC20100831: set correct number of elements in new G3
NG3 = iG3OFF

call mma_deallocate(idx2ij)
call mma_deallocate(icnj)
call mma_deallocate(ip1_buf)
call mma_deallocate(TASKLIST)
! free CI buffers
call mma_deallocate(BUF1)
call mma_deallocate(BUF2)
call mma_deallocate(BUFT)
call mma_deallocate(BUFD)

!-SVC20100302: Synchronized add into the densitry matrices
!  only for the G1 and G2 replicate arrays
call GADGOP(G1,NG1,'+')
call GADGOP(G2,NG2,'+')

call GADGOP(F1,NG1,'+')
call GADGOP(F2,NG2,'+')

#ifdef _ENABLE_BLOCK_DMRG_
NLEV4 = NLEV2**2

! allocate work space to store 3RDM
call mma_allocate(G3TMP,NLEV4,Label='G3Tmp')

! TODO: Here, several options to compute F3.
! Currently implemented only cu4, but cu34 and F3 from DMRG-sweep
! will be possible. They should be implemented at this section.

! MKFG3CU4 is located under block_dmrg_util/
call MKFG3CU4(mkF,nLEV,G1,F1,G2,F2,G3,F3,idxG3,nG3,G3TMP)

call mma_deallocate(G3TMP)
#endif

! TODO: @kszenes: this should be wrapped in an if statement
#ifdef _ENABLE_CHEMPS2_DMRG_
call mkfg3chemps2(mkF,NLEV,G1,F1,G2,F2,G3,F3,idxG3,nG3)
#endif

#ifdef _DMRG_
if (DMRG) call mkfg3qcm(mkF,nLEV,G1,F1,G2,F2,G3,F3,idxG3,nG3)
#endif

if (iPrGlb >= DEBUG) then
  !SVC: if running parallel, G3/F3 are spread over processes,
  !     so make sure that the _total_ fingerprint is computed
  dG1 = DNRM2_(NG1,G1,1)
  dG2 = DNRM2_(NG2,G2,1)
  dG3 = DDOT_(NG3,G3,1,G3,1)
  call GADGOP_SCAL(dG3,'+')
  dG3 = sqrt(dG3)
  dF1 = DNRM2_(NG1,F1,1)
  dF2 = DNRM2_(NG2,F2,1)
  dF3 = DDOT_(NG3,F3,1,F3,1)
  call GADGOP_SCAL(dF3,'+')
  dF3 = sqrt(dF3)

  write(u6,'("DEBUG> ",A)') 'norms of the density matrices:'
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G1:',dG1
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G2:',dG2
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G3:',dG3
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F1:',dF1
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F2:',dF2
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F3:',dF3
end if

end subroutine MKFG3DM

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(MKFG3DM)

#endif
