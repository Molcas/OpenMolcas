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
! Copyright (C) Per Ake Malmqvist                                      *
!               Steven Vancoillie                                      *
!***********************************************************************
!> @brief
!>  Procedure for computing 1-body, 2-body, and 3-body
!>  density elements with active indices only,
!>  and related matrices obtained from contractions with
!>  the diagonal one-electron Hamiltonian.
!> @author Per &Aring;ke Malmqvist
!> @modified_by Steven Vancoillie
!>
!> @details
!> Computation of the 1-, 2-, and 3-body density matrices defined as
!> \f{align}{
!> G1(t,u)         &= \langle 0 \lvert E_{tu} \rvert 0 \rangle \\
!> G2(t,u,v,x)     &= \langle 0 \lvert E_{tuvx} \rvert 0 \rangle \\
!> G3(t,u,v,x,y,z) &= \langle 0 \lvert E_{tuvxyz} \rvert 0 \rangle \\
!> \f}
!> and the contractions with the diagonal 1-el Hamiltonian
!> \f{align}{
!> F1(t,u)         &= \sum_w \langle 0 \lvert E_{tuww} \rvert 0 \rangle e_w \\
!> F2(t,u,v,x)     &= \sum_w \langle 0 \lvert E_{tuvxww} \rvert 0 \rangle e_w \\
!> F3(t,u,v,x,y,z) &= \sum_w \langle 0 \lvert E_{tuvxyzww} \rvert 0 \rangle e_w \\
!> \f}
!> Storage: \p G1 and \p G2 are simple two- and four-index arrays, and
!> includes also such zeroes that are implied by symmetry.
!> But \p G3 is quite large, and while it is stored with zeroes, it
!> is made more compact by calculating only the minimum amount of
!> unique values and storing the active indices in the array \p idxG3.
!> Later, the full matrix can be restored on the fly by using the
!> full permutational symmetry (see ::mksmat and ::mkbmat). The
!> same storage applies to the \f$ F \f$ matrices.
!>
!> @param[in]  mkF   switch to activate computation of \f$ F \f$ matrices
!> @param[in]  CI    wave function CI coefficients, with symmetry \c STSYM
!> @param[in]  nCI   number of CI coefficients, with symmetry \c STSYM
!> @param[out] G1    1-body active density matrix
!> @param[out] G2    2-body active density matrix
!> @param[out] G3    process-local part of 3-body active density matrix
!> @param[out] F1    1-body active density matrix contracted with
!>                   diagonal 1-el Hamiltonian
!> @param[out] F2    2-body active density matrix contracted with
!>                   diagonal 1-el Hamiltonian
!> @param[out] F3    process-local part of 3-body active density matrix
!>                   contracted with diagonal 1-el Hamiltonian
!> @param[out] idxG3 table to translate from process-local array index
!>                   to active indices

subroutine MKFG3(mkF,CI,nCI,G1,F1,G2,F2,G3,F3,idxG3,NLEV,nG1,nG2,nG3)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use fciqmc_interface, only: DoFCIQMC, mkfg3fciqmc
use PrintLevel, only: DEBUG, VERBOSE
use sguga, only: CIS, EXS, L2ACT, SGS
use caspt2_global, only: do_grad, iPrGlb, iTasks_grad, nbuf1_grad, nStpGrd, nTasks_grad
use caspt2_module, only: EPSA, MxCI, nActEl, nAshT, nBasT, nSym, STSym
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, byte, RtoB

implicit none
logical(kind=iwp), intent(in) :: mkF
integer(kind=iwp), intent(in) :: nCI, NLEV, nG1, nG2
real(kind=wp), intent(in) :: CI(nCI)
integer(kind=iwp), intent(inout) :: nG3
real(kind=wp), intent(out) :: G1(NLEV,NLEV), F1(NLEV,NLEV), G2(NLEV,NLEV,NLEV,NLEV), F2(NLEV,NLEV,NLEV,NLEV), G3(nG3), F3(nG3)
integer(kind=byte), intent(out) :: idxG3(6,nG3)
integer(kind=iwp) :: I, IB, IBMN, IBMX, IBUF, IBUF1, ID, IDX, IG3, IG3OFF, IOFFSET, IP1, IP1END, IP1I, IP1MN, IP1MX, IP1STA, IP2, &
                     IP3, IP3MX, IQ1, ISP1, ISSG1, ISSG2, ISTU, ISUBTASK, ISVX, ISYZ, IT, ITASK, ITLEV, IU, IULEV, IV, IVLEV, IX, &
                     IXLEV, IY, IYLEV, IZ, IZLEV, J, JDX, MEMMAX, MEMMAX_SAFE, MXTASK, MYBUFFER, MYTASK, NB, NBTOT, NBUF1, &
                     NLEV2, NSGM1, NSGM2, NSUBTASKS, NTASKS, NTRI1, NTRI2
real(kind=wp) :: DF1, DF2, DF3, DG1, DG2, DG3
integer(kind=iwp), allocatable :: ICNJ(:), IDX2IJ(:,:), IJ2IDX(:,:), IP1_BUF(:), TASKLIST(:,:)
real(kind=wp), allocatable :: BUF1(:,:), BUF2(:), BUFD(:), BUFR(:), BUFT(:)
real(kind=wp), external :: DDOT_, DNRM2_

! IJ2IDX, IDX2IJ, ICNJ, IP1_BUF: translation tables for levels i,j to and from pair indices idx
! BUFR: result buffer, maximum size is the largest possible ip1 range,
! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2

! Put in zeroes. Recognize special cases:
if (nlev == 0) return

G1(:,:) = Zero
G2(:,:,:,:) = Zero
G3(:) = Zero
if (mkF) then
  F1(:,:) = Zero
  F2(:,:,:,:) = Zero
  F3(:) = Zero
end if

if (NACTEL == 0) return

! This should not happen, but...
if (NCI == 0) return

! Here, for regular CAS or RAS cases.

call mma_allocate(IJ2IDX,nLev,nLev,Label='IJ2IDX')
call mma_allocate(IDX2IJ,2,nLev**2,Label='IDX2IJ')
call mma_allocate(ICNJ,nLev**2,Label='ICNJ')
IJ2IDX(:,:) = -1
IDX2IJ(:,:) = -1
ICNJ(:) = -1

! Special pair index idx2ij allows true RAS cases to be handled:
nlev2 = nlev**2
ntri1 = nTri_Elem(nlev-1)
ntri2 = nTri_Elem(nlev)
idx = 0
do i=1,nlev-1
  do j=i+1,nlev
    ! i<j
    idx = idx+1
    ij2idx(i,j) = idx
    idx2ij(1,idx) = i
    idx2ij(2,idx) = j
    ! i>j
    jdx = nlev2+1-idx ! note the reverse indexation, ...-idx!
    ij2idx(j,i) = jdx
    idx2ij(1,jdx) = j
    idx2ij(2,jdx) = i
  end do
end do
! i=j
do i=1,nlev
  idx = ntri1+i
  ij2idx(i,i) = idx
  idx2ij(1,idx) = i
  idx2ij(2,idx) = i
end do
! Loop over the compond index, idx, corresponding to (i,j) and
! tabulate in icnj(idx) the compound index jdx that corresponds
! to the pair (j,i)
do idx=1,nlev2
  i = idx2ij(1,idx)
  j = idx2ij(2,idx)
  icnj(idx) = ij2idx(j,i)
end do

call mma_deallocate(IJ2IDX)

call mma_MaxDBLE(memmax)

! Use *almost* all remaining memory:
memmax_safe = int(real(memmax,kind=wp)*0.95_wp)

! Buffers to compute CI expansion vectors into:
! <Psi0|E_ip1 | E_ip2 E_ip3|Psi0>
! buf1: bra buffer with E_ip1 excitations of Psi0
!       holds multiple CI vectors (allowed by memory)
! buf2: ket buffer for an E_ip3 excitation of Psi0
! buft: ket buffer for an E_ip2 excitation of E_ip3|Psi0>
! bufd: diagonal matrix elements to compute the F matrix
nbuf1 = max(1,min(nlev2,(memmax_safe-3*mxci)/mxci))
!! if gradient, nbuf1 must be consistent here and in derfg3
if (do_grad .or. (nStpGrd == 2)) then
  !! Compute approximate available memory in derfg3
  !! allocated in DENS
  memmax_safe = memmax_safe-16*NBAST**2-3*NASHT**2
  !! allocated (additionally) in CLagX
  memmax_safe = memmax_safe-2*(NG1+NG2+NG3)
  nbuf1 = max(1,min(nlev2,(memmax_safe-(6+nlev)*mxci)/mxci/3))
  nbuf1_grad = nbuf1
end if
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
  write(u6,'(2X,A,F16.9,A)') ' memory avail: ',real(memmax*RtoB,kind=wp)*1.0e-9_wp,' GB'
  write(u6,'(2X,A,F16.9,A)') ' memory used:  ',real(((nbuf1+3)*MXCI)*RtoB,kind=wp)*1.0e-9_wp,' GB'
end if

call mma_allocate(ip1_buf,nlev2,Label='ip1_buf')
call mma_allocate(bufr,nlev2,Label='bufr')

!! For gradients, iTasks_grad has the sequence of subtasks executed below so that derfg3 processes exactly the same batching
if (do_grad .or. (nStpGrd == 2)) then
  nTasks_grad = 0
  !! This nTasks_grad is the largest number of tasks; the actual number is smaller if parallel
  do issg1 = 1, nsym
    call build_TaskList(issg1,nTasks,nSubTasks)
    nTasks_grad = nTasks_grad + nSubTasks
  end do
  if (allocated(iTasks_grad)) call mma_deallocate(iTasks_grad)
  call mma_allocate(iTasks_grad,max(1,nTasks_grad),Label='Tasks_grad')
  iTasks_grad(:) = 0
  nTasks_grad = 0
end if

!***********************************************************************
!                                                                      *
! A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.
! This also allows precomputing the Hamiltonian (H0) diagonal elements.

iG3OFF = 0
Symmetry_Loop: do issg1=1,nsym   ! Symmetry index of E_ut/0>
!                                                                      *
!***********************************************************************

  isp1 = Mul(issg1,stsym)       ! Symmetry index of E_ut

  if (.not. DoFCIQMC) then

    ! form: H0_I = \sum_t e_{t}  |I><I|E_{tt}|I><I|

    ! Note: this only works if the Fock matrix is presented in the
    ! block diagonal form. i.e. in the pseudo canonical basis.

    nsgm1 = CIS%ncsf(issg1)
    if (mkF) call H0DIAG_CASPT2(ISSG1,BUFD,nsgm1,CIS%NOW,CIS%IOW,CIS%nMidV)

  end if

  ! build the task list for this symmetry (TaskList and ip1_buf) and count the number of larger tasks/inner-loop subtasks.
  call build_TaskList(issg1,nTasks,nSubTasks)

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

      if (.not. DoFCIQMC) then
        do ip1i=ip1sta,ip1end
          itlev = idx2ij(1,ip1i)
          iulev = idx2ij(2,ip1i)

          istu = Mul(SGS%ism(itlev),SGS%ism(iulev))

          it = L2ACT(itlev)
          iu = L2ACT(iulev)

          if (istu == isp1) then
            ibuf1 = ibuf1+1
            ip1_buf(ibuf1) = ip1i
            ! form E_ut |0>
            BUF1(1:nSgm1,iBuf1) = Zero
            call SG_Epq_Psi(SGS,CIS,EXS,IULEV,ITLEV,One,STSYM,CI,BUF1(:,ibuf1))
          end if
        end do

      else

        do ip1i=ip1sta,ip1end
          itlev = idx2ij(1,ip1i)
          iulev = idx2ij(2,ip1i)
          istu = Mul(SGS%ism(itlev),SGS%ism(iulev))
          it = L2ACT(itlev)
          iu = L2ACT(iulev)
          if (istu == isp1) then
            ibuf1 = ibuf1+1
            ip1_buf(ibuf1) = ip1i
          end if
        end do
      end if

      myBuffer = iTask
    else
      ibuf1 = TaskList(iTask,3)
    end if

    !-SVC20100301: necessary batch of sigma vectors is now in the buffer
    if (.not. DoFCIQMC) then
      ! The ip1 buffer could be the same on different processes
      ! so only compute the G1 contribution when ip3 is 1, as
      ! this will only be one task per buffer.

      ! form <0| *  E_ut|0> = G_tu and
      !      <0| * H0 * E_ut|0> - e_t G_tu = F1_tu

      if ((issg1 == stsym) .and. (ip3 == 1)) then
        if (mkF) then
          do ib=1,ibuf1
            idx = ip1_buf(ib)
            itlev = idx2ij(1,idx)
            iulev = idx2ij(2,idx)
            it = L2ACT(itlev)
            iu = L2ACT(iulev)
            G1(it,iu) = dot_product(CI(1:nsgm1),BUF1(1:nsgm1,ib))
            F1(it,iu) = sum(CI(1:nsgm1)*BUF1(1:nsgm1,ib)*bufd(1:nsgm1))-EPSA(iu)*G1(it,iu)
          end do
        else
          do ib=1,ibuf1
            idx = ip1_buf(ib)
            itlev = idx2ij(1,idx)
            iulev = idx2ij(2,idx)
            it = L2ACT(itlev)
            iu = L2ACT(iulev)
            G1(it,iu) = dot_product(CI(1:nsgm1),BUF1(1:nsgm1,ib))
          end do
        end if
      end if
    end if

    !ip3mx = ntri2
    !if (ip1end <= ntri2) ip3mx =ip1end
    !if (ip1sta > ntri2) ip3mx =ntri1
    !-SVC20100309: loop over ip3, ip2
    !do ip3=1,ip3mx

    !-SVC20100309: PAM's magic formula
    !iCnt = iSubTask-iOffSet
    !ip3 = int(real(ntri2,kind=wp)+OneHalf-sqrt((real(ntri2,kind=wp)+Half)**2-2*iCnt+1.0e-6_wp))
    !ip2 = iCnt-((ip3-1)*ntri2-nTri_Elem(ip3-2)+ip3-1

    !-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

    iq1 = icnj(ip3)
    ! The indices corresponding to pair index p3:
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    isyz = Mul(SGS%ism(iylev),SGS%ism(izlev))
    issg2 = Mul(isyz,stsym)
    iy = L2ACT(iylev)
    iz = L2ACT(izlev)
    if (.not. DoFCIQMC) then
      nsgm2 = CIS%ncsf(issg2)
      BUF2(1:nSgm2) = Zero
      ! form <0| E_zy|
      call SG_Epq_Psi(SGS,CIS,EXS,IYLEV,IZLEV,One,STSYM,CI,BUF2)

      if (issg2 == issg1) then

        if (mkF) then
          do ib=1,ibuf1
            idx = ip1_buf(ib)
            itlev = idx2ij(1,idx)
            iulev = idx2ij(2,idx)
            it = L2ACT(itlev)
            iu = L2ACT(iulev)
            ! form <0| E_zy E_ut |0> = G_tu,yz
            !      <0| E_zy * H0 * E_ut |0> = F_tu,yz
            G2(it,iu,iy,iz) = dot_product(BUF2(1:nsgm1),BUF1(1:nsgm1,ib))
            F2(it,iu,iy,iz) = sum(BUF2(1:nsgm1)*bufd(1:nsgm1)*buf1(1:nsgm1,ib))
          end do
        else
          do ib=1,ibuf1
            idx = ip1_buf(ib)
            itlev = idx2ij(1,idx)
            iulev = idx2ij(2,idx)
            it = L2ACT(itlev)
            iu = L2ACT(iulev)
            G2(it,iu,iy,iz) = dot_product(BUF2(1:nsgm1),BUF1(1:nsgm1,ib))
          end do
        end if

      end if

    end if

    nbtot = 0
    do ip2=ip3,ntri2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx = Mul(SGS%ism(ivlev),SGS%ism(ixlev))
      iv = L2ACT(ivlev)
      ix = L2ACT(ixlev)
      if (isvx == Mul(issg1,issg2)) then
        if (.not. DoFCIQMC) then
          BUFT(1:nSgm1) = Zero
          ! form <0| E_zy E_xv
          call SG_Epq_Psi(SGS,CIS,EXS,IVLEV,IXLEV,One,ISSG2,BUF2,BUFT)
        end if
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
          if (.not. DoFCIQMC) then
            ! form <0| E_zy E_xv * E_ut|0>
            call DGEMV_('T',nsgm1,nb,One,BUF1(:,ibmn),mxci,buft,1,Zero,bufr,1)
            ! and distribute this result into G3:
            G3(iG3OFF+1:iG3OFF+nb) = Bufr(1:nb)
            ! and copy the active indices into idxG3:
          end if
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

          if (.not. DoFCIQMC) then
            if (mkF) then
              ! Elementwise multiplication of Tau with H0 diagonal - EPSA(IV):
              ! form <0| E_zy E_xv (H0 -e_v)
              BufT(1:nSgm1) = (BufD(1:nSgm1)-EpsA(iV))*BufT(1:nSgm1)
              !do icsf=1,nsgm1
              !  buft(icsf) = (bufd(icsf)-epsa(iv))*buft(icsf)
              !end do
              ! so Tau is now = Sum(eps(w)*E_vxww) Psi. Contract and distribute:

              ! form <0| E_zy E_xv (H0 -e_v) E_ut|0> = F3(iuv,xyz)
              call DGEMV_('T',nsgm1,nb,One,BUF1(:,ibmn),mxci,buft,1,Zero,bufr,1)
              F3(iG3OFF+1:iG3OFF+nb) = Bufr(1:nb)
            end if
          end if

          iG3OFF = iG3OFF+nb
          nbtot = nbtot+nb
        end if
      end if
    end do

    if (iPrGlb >= DEBUG) write(u6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)') iSubTask,ip1sta,ip1end,ip3,nbtot
    !! consistent tasks must be executed here and in derfg3 for grad
    if (do_grad .or. (nStpGrd == 2)) then
      nTasks_grad = nTasks_grad+1
      iTasks_grad(nTasks_grad) = iSubTask
    end if

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

!***********************************************************************
!                                                                      *
! End of sectioning loop over symmetry of Sgm1 wave functions.
end do Symmetry_Loop
!-SVC20100831: set correct number of elements in new G3
NG3 = iG3OFF
!                                                                      *
!***********************************************************************

call mma_deallocate(ip1_buf)
call mma_deallocate(bufr)
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

if (mkF) then
  call GADGOP(F1,NG1,'+')
  call GADGOP(F2,NG2,'+')
end if

if (DoFCIQMC) then
  call mkfg3fciqmc(G1,G2,G3,F1,F2,F3,idxG3,nLev)
else
  ! Correction to G2: It is now = <0| E_tu E_yz |0>
  do iu=1,nlev
    G2(:,iu,iu,:) = G2(:,iu,iu,:)-G1(:,:)
  end do
  ! SVC20100310: took some spurious mirroring of G2 values out
  ! of the loops and put them here, after the parallel section has
  ! finished, so that GAdGOP works correctly.
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
      G2(it,iu,iy,iz) = G2(iz,iy,iu,it)
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
      G2(it,iu,iy,iz) = G2(iy,iz,it,iu)
    end do
  end do
  if (mkF) then
    ! Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
    do iy=1,nlev
      do iu=1,nlev
        F2(:,iu,iy,:) = F2(:,iu,iy,:)-(EPSA(iu)+EPSA(iy))*G2(:,iu,iy,:)
      end do
    end do
    do iu=1,nlev
      F2(:,iu,iu,:) = F2(:,iu,iu,:)-(F1(:,:)+EPSA(iu)*G1(:,:))
    end do
    ! SVC20100310: took some spurious mirroring of F2 values out
    ! of the loops and put them here, after the parallel section has
    ! finished, so that GAdGOP works correctly.
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
        F2(it,iu,iy,iz) = F2(iz,iy,iu,it)
      end do
    end do
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
        F2(it,iu,iy,iz) = F2(iy,iz,it,iu)
      end do
    end do
  end if

  ! Correction to G3: It is now <0| E_tu E_vx E_yz |0>
  ! Similar for F3 values.
  if (mkF) then
    do iG3=1,NG3
      iT = idxG3(1,iG3)
      iU = idxG3(2,iG3)
      iV = idxG3(3,iG3)
      iX = idxG3(4,iG3)
      iY = idxG3(5,iG3)
      iZ = idxG3(6,iG3)
      ! Correction: From <0| E_tu E_vx E_yz |0>, form <0| E_tuvxyz |0>
      if (iY == iX) then
        G3(iG3) = G3(iG3)-G2(iT,iU,iV,iZ)
        F3(iG3) = F3(iG3)-(F2(iT,iU,iV,iZ)+EPSA(iu)*G2(iT,iU,iV,iZ))
        if (iv == iu) then
          G3(iG3) = G3(iG3)-G1(iT,iZ)
          F3(iG3) = F3(iG3)-F1(iT,iZ)
        end if
      end if
      if (iV == iU) then
        G3(iG3) = G3(iG3)-G2(iT,iX,iY,iZ)
        F3(iG3) = F3(iG3)-(F2(iT,iX,iY,iZ)+EPSA(iY)*G2(iT,iX,iY,iZ))
      end if
      if (iY == iU) then
        G3(iG3) = G3(iG3)-G2(iV,iX,iT,iZ)
        F3(iG3) = F3(iG3)-(F2(iV,iX,iT,iZ)+EPSA(iU)*G2(iV,iX,iT,iZ))
      end if
      F3(iG3) = F3(iG3)-(EPSA(iU)+EPSA(iY))*G3(iG3)
    end do
  else
    do iG3=1,NG3
      iT = idxG3(1,iG3)
      iU = idxG3(2,iG3)
      iV = idxG3(3,iG3)
      iX = idxG3(4,iG3)
      iY = idxG3(5,iG3)
      iZ = idxG3(6,iG3)
      ! Correction: From <0| E_tu E_vx E_yz |0>, form <0| E_tuvxyz |0>
      if (iY == iX) then
        G3(iG3) = G3(iG3)-G2(iT,iU,iV,iZ)
        if (iv == iu) G3(iG3) = G3(iG3)-G1(iT,iZ)
      end if
      if (iV == iU) G3(iG3) = G3(iG3)-G2(iT,iX,iY,iZ)
      if (iY == iU) G3(iG3) = G3(iG3)-G2(iV,iX,iT,iZ)
    end do
  end if
end if

call mma_deallocate(idx2ij)
call mma_deallocate(icnj)

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

  write(u6,'("DEBUG> ",A)') 'MKFG3: norms of the density matrices:'
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G1:',dG1
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G2:',dG2
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'G3:',dG3
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F1:',dF1
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F2:',dF2
  write(u6,'("DEBUG> ",A,1X,ES21.14)') 'F3:',dF3
end if

contains

! Build the task list for symmetry issg
subroutine build_TaskList(issg,nTasks_,nSubTasks_)

  integer(kind=iwp), intent(in) :: issg
  integer(kind=iwp), intent(out) :: nTasks_, nSubTasks_
  integer(kind=iwp) :: isp1, ibuf1, iTask, ip1, itlev, iulev, istu, iOffSet, ip1sta, ip1end, ip3mx

  isp1 = Mul(issg,stsym)  ! Symmetry index of E_ut

  !-SVC20100301: calculate number of larger tasks for this symmetry, this
  !-is basically the number of buffers we fill with SG_Epq_Psi vectors.

  iTask = 1
  ibuf1 = 0
  do ip1=1,nlev2
    itlev = idx2ij(1,ip1)
    iulev = idx2ij(2,ip1)
    istu = Mul(SGS%ism(itlev),SGS%ism(iulev)) ! Symmetry of E_tu

    if (istu == isp1) then
      ! Add pair index to ip1_buf
      ibuf1 = ibuf1+1
      ip1_buf(ibuf1) = ip1
      ! Add pair index to Tasklist(iTask,1) if it is the first entry in ip1_buf
      if (ibuf1 == 1) TaskList(iTask,1) = ip1  ! Start ip1 index
    end if

    ! Terminate if buffer is full, or, if
    if ((ibuf1 == nbuf1) .or. ((ibuf1 > 0) .and. ((ip1 == ntri2) .or. (ip1 == nlev2)))) then
      TaskList(iTask,2) = ip1_buf(ibuf1)  ! End ip1 index
      TaskList(iTask,3) = ibuf1           ! End ibuf value
      iTask = iTask+1
      ibuf1 = 0
    end if

  end do

  nTasks_ = iTask
  if (ibuf1 == 0) nTasks_ = nTasks_-1

  !-SVC20100309: calculate number of inner loop iteration tasks.
  iOffSet = 0
  do iTask=1,nTasks_
    TaskList(iTask,4) = iOffSet
    ip1sta = TaskList(iTask,1)
    ip1end = TaskList(iTask,2)
    ip3mx = ntri2
    if (ip1end <= ntri2) ip3mx = ip1end
    if (ip1sta > ntri2) ip3mx = ntri1
    !-SVC20100309: Currently -we are going to limit this to the ip3-loop and
    !-leave the ip2-loop intact.  This was based on the large overhead which
    !-was observed for a very large number of small tasks.
    !iOffSet = iOffSet+ip3mx*ntri2-nTri_Elem(ip3mx-1)
    iOffSet = iOffSet+ip3mx
  end do
  nSubTasks_ = iOffSet

end subroutine build_tasklist

end subroutine MKFG3
