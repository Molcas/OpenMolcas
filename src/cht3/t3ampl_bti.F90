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

subroutine T3AMPL_BTI(OEH,OEP)
!mp subroutine T3AMPL_BTI(W,OEH,OEP)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    *** PROGRAM TRIPLY-UHF/RHF - REDUCED DIMENSION and timing - DRIVER
!
!    CALCULATION OF TRIPLE EXCITATION CONTRIBUTION TO
!    THE 4TH ORDER CORRELATION ENERGY OF MB RSPT AND CC-TRIPLES
!
!    Provides T3AMPL_BLOCKEDIMPROVED:  Accelerated algorithm.
!
!     Combines triply w/ trick by JN and limited I/Os by buffered
!     blocking scheme by PV. Blocks the virtual space to pieces
!     optimal either for parallelization or in using the available
!     memory
!     Involves also more robust and more reliable memory allocator.
!
! History
!     Buffered blocking scheme: PV, LAOG Grenoble, 16 april 2003.
!     First version w/ trick: JN, June 12, 2003
!     Parallel version: PV, 15 oct 2003.
!     Implemented integer offsets, PV, 14 may 2004.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

use Para_Info, only: MyRank, nProcs

implicit none
integer IUHF, LU(6)
integer i, nuga, nugc, nga, ngb, ngc, vblock, it1
integer isp, krem
real*8 OEH(*), OEP(*), ddot_, ccsdt, ccsdt4, energ(4), tccsd, times(10), times_parr(10), totcpu, totwal, timerel
!real*8 cpu0,cpu1,wall0,wall1
character FN*6
logical ifvo
integer la, t1a, t1b
logical lastcall, scored, skip
!integer maxspace
integer nla, nlb
real*8 enx1
#include "param_cht3.fh"
#include "ioind.fh"
#include "uhf.fh"
#include "cht3_casy.fh"
#include "WrkSpc.fh"
#include "ndisk.fh"
#include "dupfiles.fh"
!mp !#include 'task_info_inc'
!mp !#include 'ws_conn_inc'
#include "cht3_ccsd1.fh"
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
real*8 :: real_buffer(1)
#endif
!??? #ifdef _MOLCAS_MPP_
!??? integer mytype,krems(nhmx)
!??? #endif
!mp !!integer me,err,len_trim
!? integer nprocs0
!?cmp common /my_mpi_world_com/ me, nprocs
integer itmp
integer imy_tsk
integer id, j
logical rsv_tsk
external rsv_tsk
real*8 e_ccsd, e_scf
real*8 cpu0_aaa, wall0_aaa, cpu0_aab, wall0_aab

!? nprocs0 = nprocs
! Uncomment the following to force sequential mode
! nprocs = 1

if (nprocs > 1) then
  write(6,'(A,i4,A)') ' Parallel run on ',nprocs,' nodes'
else
  write(6,'(A,i4,A)') ' Serial run'
end if
call xflush(6)

! t^ijk_abc . D^abc_ijk =
!
! alpha-alpha-alpha (K_ab^ir+K_ab^jr-K_ac^kr)*(L_rc^jk+L_rc^ki+L_rc^ij)
! or beta-beta-beta   + permutations bca, acb
!        - M_abc^iki -M_abc^iij -M_abc^jjk+M_abc^jij-M_abc^kki+M_abc^kjk
!         + permutations bca cab
!
!        where    M_abc^ijk=K_ab^ir*L_ra^jk
!
!  strategy:   blocks the number of virtuals to pieces by ng
!              then everything in the memory!
!
! free everything from W(1) including Hermit space
! (PV) done for slaves in calling routine slave_T3AMPL_BTI
!
! parallelization: distribute common UHF
! (PV) done for slaves in calling routine slave_T3AMPL_BTI

!mp !!call gettim(cpu0,wall0)

do I=1,6
  LU(I) = 90+I
end do

call zeroma(times,1,10)
call zeroma(energ,1,4)

! calculate the T2-T3 contraction?  ifvo is the answer ! open-shell stuff !
!mp call get3dm('FOC-VO',w,noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
!mp ifvo = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w,1,w,1) > 1d-14
ifvo = .false.

! parallelization:  do this on each host >>>>>>> start

! allocates for two arrays of the size T1

it1 = NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)

call GetMem('t3_t1a','Allo','Real',t1a,it1)
call GetMem('t3_t1b','Allo','Real',t1b,it1)

call zeroma(Work(t1a),1,it1)
call zeroma(Work(t1b),1,it1)

call GetMem('t3_t1','Allo','Real',itmp,noab(1)*nuab(1))
call GetMem('t3_la','Allo','Real',la,it1)

!mp ! read t1 amplitudes
if (printkey >= 10) then
  write(6,*) 'Reading ',it1,' t1_alpha, t1_beta amplitudes'
end if
call GetRest_t3(Work(la),Work(itmp),e_ccsd)

!mp !write(6,*) 'Ze t1 = ',ddot_(it1,Work(la),1,Work(la),1)

! t1a and t1b always remain in the memory -- small fields

! parallelization:   on each host <<<<<<< end  - block done

! Remaining space can be used for the blocking in-core algorithms

call GetMem('t3_la','Max','Real',krem,krem)

if (printkey >= 10) then
  write(6,*)
  write(6,*) 'Available Memory before v_size_t3 = ',krem
  call xflush(6)
end if

! determines the virtual block size
! vblock - the virtual orbitals block size
! look for optimal number with respect nprocs to v_size_t3

write(6,*)
write(6,'(2x,A)') 'Starting triply with blocked virtuals algorithm'
if (nprocs > 1) then
  write(6,*) ' Node Number ',MyRank
end if

! Compute memory pattern.

call v_size_t3(vblock,nprocs,krem,printkey)

IUHF = 1 ! open-shell stuff
if (IOPT(76) == 0) IUHF = 0

if (printkey >= 10) then
  write(6,'(A,i1)') ' Closed-Shell calculation, IUHF = ',iuhf
end if

! create K(beta-alpha,beta-alpha) and L(alpha-beta,beta-alpha)
ndup = 0 !?????

!mp - print information on number of steps in loopa and loopb

call check_loops(nuab(1),vblock,nla,nlb)
write(6,*)
write(6,'(A,i6)') 'Number of steps in loopa : ',nla
write(6,'(A,i6)') 'Number of steps in loopb : ',nlb

! checking :

if (t3_stopa > nla) then
  write(6,*) 'Too many steps in t3_loopa requested in input'
  write(6,*) 'T3_LOOPA from input : ',t3_starta,t3_stopa,t3_stopa-t3_startb+1
  write(6,*) 'total T3_LOOPA steps : ',nla
  call abend()
end if

if (t3_stopb > nlb) then
  write(6,*) 'Too many steps in t3_loopb requested in input'
  write(6,*) 'T3_LOOPB from input : ',t3_startb,t3_stopb,t3_stopb-t3_startb+1
  write(6,*) 'total T3_LOOPB steps : ',nlb
  call abend()
end if

!mp
if (gen_files) then

  write(6,*)
  write(6,*) 'Creating KMAT, LMAT Scratch Integral Files-----------------------------'
  write(6,*)
  !mp
  call create_klvab_t3(vblock)
  !mp
  call CWTime(TCpu,TWall)

  if (printkey > 1) then
    write(6,*)
    write(6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(6,*)
    write(6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0d0*TCpu/(TWall-TWall0)
    write(6,*)
  end if

  TCpu_l = TCpu
  TWall_l = TWall

  !mp

  write(6,*) 'Create of Integrals done-------------------------------------------------'
  write(6,*)
  !mp

else
  write(6,*)
  write(6,*) 'Skipping KMAT, LMAT scratch integral files generation as requested'
end if
!mp
cpu0_aaa = TCpu
wall0_aaa = TWall
!
if (.not. run_triples) then
  write(6,*)
  write(6,*) 'Exiting triples after scratch integral files generation as requested'
  return
  !? call abend()
end if
!mp

! creates K(alpha-alpha,alpha-alpha) (dummy for closed shell)
! (master only)
!mp
!mp ! open-shell stuff
!mp !
!mp !if (iuhf /= 0) then
!mp !  do isp=1,1+iuhf
!mp !    call create_klvaa_t3(w(la),vblock,isp)
!mp !  end do
!mp !end if

! parallelization: files KMATxy LMATxy x=A,B  y=A,B (BA - closed shell)
!                  generated on each host

! Sorting and initialization timing

!mp !call gettim(cpu1,wall1)  !  prerob na poriadne timingy
call CWTime(TCpu,TWall)
times(9) = times(9)+TCpu-TCpu0
times(10) = times(10)+TWall-TWall0

do isp=1,1+iuhf

  nuga = nuab(isp)/vblock
  if ((nuga*vblock) < nuab(isp)) nuga = nuga+1
  nugc = nuab(3-isp)/vblock
  if ((nugc*vblock) < nuab(3-isp)) nugc = nugc+1
  ! creates K(alpha-alpha,alpha-alpha) (dummy for closed shell)
  !!if (iuhf /= 0) then   ! open-shell stuff.
  !!  call create_klvaa_t3(w(la),vblock,isp)
  !!end if
  FN(5:5) = ich(isp)
  FN(6:6) = ich(isp)
  FN(1:4) = 'KMAT'
  call multi_opendir(FN,LU(1))
  FN(1:4) = 'LMAT'
  call multi_opendir(FN,LU(2))

  ! parallelization: files KMATxx LMATxx are assumed to be available
  ! to all tasks on each host

  !                  not for closed shell
  !     alpha-alpha-alpha   or closed shell

  ! - skip t3loopa if needed

  if ((t3_starta < 0) .and. (t3_startb > 0)) then
    write(6,*)
    write(6,*) 'Skipping t3loopa on user request '
    write(6,*)
  else

    i = 0
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,ngb
          i = i+1
        end do
      end do
    end do

    if (t3_starta < 0) then
      call GetMem('Imy_tsk','Allo','Inte',imy_tsk,3*i)
    else
      call GetMem('Imy_tsk','Allo','Inte',imy_tsk,(t3_stopa-t3_starta+1)*3)
    end if

    if (t3_starta < 0) then
      i = 0
      do nga=1,nuga
        do ngb=1,nga
          do ngc=1,ngb
            i = i+1
            iWork(imy_tsk-3+3*i) = nga
            iWork(imy_tsk-3+3*i+1) = ngb
            iWork(imy_tsk-3+3*i+2) = ngc
          end do
        end do
      end do
    else
      i = 0
      do nga=1,nuga
        do ngb=1,nga
          do ngc=1,ngb
            i = i+1

            if ((i >= t3_starta) .and. (i <= t3_stopa)) then
              iWork(imy_tsk-3+3*(i-t3_starta+1)) = nga
              iWork(imy_tsk-3+3*(i-t3_starta+1)+1) = ngb
              iWork(imy_tsk-3+3*(i-t3_starta+1)+2) = ngc
            end if

          end do
        end do
      end do
    end if

    if (t3_starta > 0) then ! for correct deallocation
      i = t3_stopa-t3_starta+1
    end if

    write(6,*)
    write(6,*) '# of tasks to be parallelized in t3loop a = ',i
    id = 666
    lastcall = .false.
    scored = .false.
    call init_tsk(id,i)
    do while (rsv_tsk(id,j))

      nga = iWork(imy_tsk-3+3*j)
      ngb = iWork(imy_tsk-3+3*j+1)
      ngc = iWork(imy_tsk-3+3*j+2)

      !mp
      !mp call GetMem('(T)','Max','Real',maxspace,maxspace)
      !mp write(6,*) 'maxspace before ',maxspace
      !mp
      call t3loopa(oeh(noab(1)*(isp-1)+1),oep(nuab(1)*(isp-1)+1),Work(t1a+noab(1)*nuab(1)*(isp-1)), &
                   Work(t1b+noab(1)*nuab(1)*(isp-1)),nga,ngb,ngc,vblock,energ,isp,LU,ifvo,lastcall,scored,j,enx1)
      !mp
      ! update 5th order terms

      !mp call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1,noab(1)*nuab(1))

      call daxpy_((noab(1)*nuab(1)),1.0d0,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1)
      ccsdt = 2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
      call daxpy_((noab(1)*nuab(1)),-1.0d0,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1)

      !mp
      call CWTime(TCpu,TWall)
      !mp
      if (t3_starta < 0) then
        write(6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j,nga,ngb,ngc,2.0d0*enx1,ccsdt, &
                                                                TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
      else
        write(6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j+t3_starta-1,nga,ngb,ngc,2.0d0*enx1, &
                                                                ccsdt,TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
      end if

      TCpu_l = TCpu
      TWall_l = TWall
      !mp
      !mp call GetMem('(T)','Max','Real',maxspace,maxspace)
      !mp write(6,*) 'maxspace after ',maxspace
      !mp

    end do
    write(6,*) 't3loopa finished'
    write(6,*)
    call Free_tsk(id)
    call GetMem('Imy_tsk','Free','Inte',imy_tsk,3*i)

    !mp !call gettim(cpu1,wall1)            ! dorob timingy !!!!!!
    call CWTime(TCpu,TWall)
    times(isp) = times(isp)+TCpu-cpu0_aaa
    times(isp+4) = times(isp+4)+TWall-wall0_aaa

    write(6,'(1X,5A,D12.4,A,D12.4,A)') 'Spin case ',ich(isp),ich(isp),ich(isp),' done:',TCpu-cpu0_aaa,' CPU [s]',TWall-wall0_aaa, &
                                       ' Wall [s]'
    call xflush(6)
    !mp
    call CWTime(TCpu,TWall)

    if (printkey > 1) then
      write(6,*)
      write(6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
      write(6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
      write(6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0d0*TCpu/(TWall-TWall0)
      write(6,*)
    end if

    TCpu_l = TCpu
    TWall_l = TWall
    !mp

  end if

  !mp
  cpu0_aab = TCpu
  wall0_aab = TWall
  !mp

  if ((t3_starta > 0) .and. (t3_startb < 0)) then
    write(6,*)
    write(6,*) 'Skipping t3loopb on user request '
    write(6,*)
    skip = .true.
    exit
  else
    skip = .false.
  end if

  ! alpha-alpha-beta or beta-beta-alpha only UHF
  !!if (IUHF /= 0)then
  !mp!!call gettim(cpu0,wall0)
  FN(5:5) = ich(3-isp)
  FN(6:6) = ich(isp)
  FN(1:4) = 'KMAT'
  call multi_opendir(FN,LU(3))
  FN(1:4) = 'LMAT'
  call multi_opendir(FN,LU(5))

  if (IUHF /= 0) then ! open-shell stuff
    FN(5:5) = ich(isp)
    FN(6:6) = ich(3-isp)
    FN(1:4) = 'KMAT'
    call multi_opendir(FN,LU(4))
    FN(1:4) = 'LMAT'
    call multi_opendir(FN,LU(6))
  else
    LU(4) = LU(3)
    LU(6) = LU(5)
  end if

  !mp
  i = 0
  do nga=1,nuga
    do ngb=1,nga
      do ngc=1,nugc
        i = i+1
      end do
    end do
  end do

  if (t3_startb < 0) then
    call GetMem('Imy_tsk','Allo','Inte',imy_tsk,3*i)
  else
    call GetMem('Imy_tsk','Allo','Inte',imy_tsk,(t3_stopb-t3_startb+1)*3)
  end if

  i = 0
  if (t3_startb < 0) then
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,nugc
          i = i+1
          iWork(imy_tsk-3+3*i) = nga
          iWork(imy_tsk-3+3*i+1) = ngb
          iWork(imy_tsk-3+3*i+2) = ngc
        end do
      end do
    end do
  else
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,nugc
          i = i+1
          if ((i >= t3_startb) .and. (i <= t3_stopb)) then
            iWork(imy_tsk-3+3*(i-t3_startb+1)) = nga
            iWork(imy_tsk-3+3*(i-t3_startb+1)+1) = ngb
            iWork(imy_tsk-3+3*(i-t3_startb+1)+2) = ngc
          end if
        end do
      end do
    end do
  end if

  if (t3_startb > 0) then
    i = t3_stopb-t3_startb+1
  end if

  write(6,*)
  write(6,*) '# of tasks to be parallelized in t3loopb = ',i
  id = 667
  lastcall = .false.
  scored = .false.
  call init_tsk(id,i)
  do while (rsv_tsk(id,j))

    nga = iWork(imy_tsk-3+3*j)
    ngb = iWork(imy_tsk-3+3*j+1)
    ngc = iWork(imy_tsk-3+3*j+2)
    !mp write(6,'(A,4(i5,2x))') 'Tsk, nga, ngb, ngc = ',j,nga,ngb,ngc

    !mp
    !mp call GetMem('(T)','Max','Real',maxspace,maxspace)
    !mp write(6,*) 'maxspace before ',maxspace
    !mp
    call t3loopb(oeh,oep,Work(t1a),Work(t1b),nga,ngb,ngc,vblock,energ(3),isp,LU,ifvo,lastcall,scored,j,enx1)

    !mp??? call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1,noab(1)*nuab(1))
    call daxpy_((noab(1)*nuab(1)),1.0d0,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1)
    ccsdt = 2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
    call daxpy_((noab(1)*nuab(1)),-1.0d0,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1)

    !mp
    call CWTime(TCpu,TWall)
    !mp
    if (t3_startb < 0) then
      write(6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j,nga,ngb,ngc,2.0d0*enx1,ccsdt, &
                                                              TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
    else
      write(6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j+t3_startb-1,nga,ngb,ngc,2.0d0*enx1, &
                                                              ccsdt,TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
    end if

    TCpu_l = TCpu
    TWall_l = TWall

    !mp
    !mp call GetMem('(T)','Max','Real',maxspace,maxspace)
    !mp write(6,*) 'maxspace before ',maxspace
    !mp

  end do
  write(6,*) 't3loopb finished'
  write(6,*)
  call Free_tsk(id)
  call GetMem('Imy_tsk','Free','Inte',imy_tsk,3*i)
  !mp write(6,*) ' energ = ',energ
  !mp
  ! - deallocate arrays in t3loopb
  !
  !mp !call gettim(cpu1,wall1)  ! urob poriadne timingy !!!!
  call CWTime(TCpu,TWall)
  times(2+isp) = times(2+isp)+TCpu-cpu0_aab
  times(6+isp) = times(6+isp)+TWall-wall0_aab
  write(6,'(1X,5A,D12.4,A,D12.4,A)') 'Spin case ',ich(isp),ich(isp),ich(3-isp),' done:',TCpu-cpu0_aab,' CPU [s]',TWall-wall0_aab, &
                                     ' Wall [s]'
  call xflush(6)
  do i=3,6
    close(LU(i))
  end do
  !!end if   !! IUHF
  close(LU(1))
  close(LU(2))
end do        ! isp
if (.not. skip) then
  !mp !!call gettim(cpu0,wall0) ! urob timingy!
  !mp
  call CWTime(TCpu,TWall)

  if (printkey > 1) then
    write(6,*)
    write(6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0d0*TCpu/(TWall-TWall0)
    write(6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
end if

do i=1,10
  times_parr(i) = times(i)
end do

#ifdef _MOLCAS_MPP_

it1 = NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)
real_buffer(1) = ccsdt
call gadgop(real_buffer,size(real_buffer),'+')
ccsdt = real_buffer(1)
call gadgop(energ,4,'+')
call gadgop(times_parr,10,'+')

#endif

call xflush(6)

! Add MPI_Reduce contribution to sorting timing
!mp !cmp call gettim(cpu1,wall1)
!mp !call CWTime(TCpu,TWall)
!mp !times(9) = times(9)+TCpu-TCpu0
!mp !times(10) = times(10)+TWall-TWall0

! T(CCSD)

if (IUHF == 0) then
  energ(2) = energ(1)
  energ(4) = energ(3)
end if
tccsd = energ(1)+energ(2)+energ(3)+energ(4)

!mp !if (IUHF == 0) then ! open-shell stuff
!mp !  call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,Work(t1a),1,noab(1)*nuab(1))
!mp !  ccsdt = 2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
!mp !else
!mp !  ccsdt = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),Work(la),1,Work(t1a),1)
!mp !  write(6,*) 'ze co do ... ?'
!mp !  stop
!mp !end if

!mp for MOLCAS verify
call Get_dScalar('SCF energy',e_scf)
call Add_Info('CHT3ene',tccsd+ccsdt,1,6)
call Add_Info('E_CHT3',tccsd+ccsdt,1,6)
call Add_Info('E_HYPE',e_scf+e_ccsd+tccsd+ccsdt,1,6)
! for NUMERICAL_GRADIENTS
call Put_cArray('Relax Method','CHT3    ',8)
call Store_Energies(1,tccsd+ccsdt+e_ccsd+e_scf,1)
!mp

if (printkey > 1) then
  write(6,*)
  write(6,*) '--------------------------------------------------'
  write(6,*)
  write(6,*) 'GADGOPed energ & ccsdt values: '
  write(6,*)
  write(6,'(A,4(f15.12,1x))') 'energ (t2-w-t3 = 2*e1 + 2*e3) ',energ
  write(6,*)
  write(6,'(A,f15.12)') 'Sum                           ',2.0d0*energ(1)+2.0d0*energ(3)
  write(6,'(A,f15.12)') 'ccsdt (e5th ord. ST)          ',ccsdt
  write(6,*)
  write(6,*) '--------------------------------------------------'
end if

! Contents of times array:
!
!              CPU      WALL
! sorting
! aaa           1        4
! bbb           2        6
! aab           3        7
! bba           4        8

! Master has only to sum up the results
do i=1,10
  times(i) = times(i)/60.d0 ! now in minutes
end do
if (nprocs > 1) then
  do i=1,10
    times_parr(i) = times_parr(i)/60.d0
  end do
end if

ccsdt4 = 0.d0
if (ifvo) then
  write(6,*) 'ifvo correspond to open-shell system. Exiting'
  call abend()

  !mp !call get3dm('FOC-VO',w(la),noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
  !mp !call transm(w(la),w(t1a),nuab(1),noab(1))
  !mp !call transm(w(noab(1)*nuab(1)+la),w(noab(1)*nuab(1)+t1a),nuab(2),noab(2))
  ! E4 T2FT3
  !mp !if (IUHF == 0) then
  !mp !  call vadd(w(t1b),1,w(t1b+noab(1)*nuab(1)),1,w(t1b),1,noab(1)*nuab(1))
  !mp !  ccsdt4 = 2.0*ddot_(noab(1)*nuab(1),w(t1a),1,w(t1b),1)
  !mp !else
  !mp !  ccsdt4 = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w(t1a),1,w(t1b),1)
  !mp !end if
end if

!RESULT(IT+3,5) = ccsdt4+ccsdt
if (IT == 0) then
  !RESULT(IT+1,5) = tccsd+ccsdt4+ccsdt
  if (ifvo) then
    write(6,9993) TCCSD
    write(6,9991) ccsdt4
    write(6,9990) ccsdt
  end if
  write(6,9994) TCCSD+ccsdt+ccsdt4
else
  !RESULT(IT+2,5) = TCCSD
  write(6,9993) TCCSD
  if (ifvo) then
    write(6,9991) ccsdt4
    write(6,9990) ccsdt
    write(6,9995) CCSDT+ccsdt4
  else
    write(6,9995) CCSDT
  end if
end if

write(6,*) '*************************************************'
write(6,*)
write(6,*) 'Final Results : '
write(6,*)
write(6,'(A,f15.12)') '   (T)     corr. = ',tccsd+ccsdt
write(6,'(A,f15.12)') '   CCSD(T) corr. = ',e_ccsd+tccsd+ccsdt
write(6,*)
write(6,*) '*************************************************'
write(6,*)

! Timing
if (nprocs > 1) write(6,*) 'Master timings:'
if (times(10) <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = times(9)/times(10)
end if
write(6,'(1x,a,2f13.3,a,f6.3)') 'sorting cpu & wall',times(9),times(10),' (min);   cpu/wall=',timerel
if (times(5) <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = times(1)/times(5)
end if
write(6,'(1x,a,2f13.3,a,f6.3)') 'aaa     cpu & wall',times(1),times(5),' (min);   cpu/wall=',timerel
if (times(6) <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = times(2)/times(6)
end if
if (IUHF /= 0) write(6,'(1x,a,2f13.3,a,f6.3)') 'bbb     cpu & wall',times(2),times(6),' (min);   cpu/wall=',timerel
if (times(7) <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = times(3)/times(7)
end if
write(6,'(1x,a,2f13.3,a,f6.3)') 'aab     cpu & wall',times(3),times(7),' (min);   cpu/wall=',timerel
if (times(8) <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = times(4)/times(8)
end if
if (IUHF /= 0) write(6,'(1x,a,2f13.3,a,f6.3)') 'bba     cpu & wall',times(4),times(8),' (min);   cpu/wall=',timerel
totcpu = times(9)
totwal = times(10)
do i=1,4
  totcpu = totcpu+times(i)
  totwal = totwal+times(4+i)
end do
if (totwal <= 0.0d0) then
  timerel = 1.0d0
else
  timerel = totcpu/totwal
end if
write(6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',totcpu,totwal,' (min);   cpu/wall=',timerel

if (nprocs > 1) then
  ! Parallel timing
  totcpu = times_parr(9)
  totwal = times_parr(10)
  times(9) = times_parr(9)
  do i=1,4
    totcpu = totcpu+times_parr(i)
    totwal = totwal+times_parr(4+i)
  end do
  if (totwal <= 0.0d0) then
    timerel = 1.0d0
  else
    timerel = totcpu/totwal
  end if
  write(6,*)
  write(6,*) 'Aggregate parallel timings:'
  write(6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',totcpu,totwal,' (min);   cpu/wall=',timerel
end if

! Return
!mp call w_debug(.false.,.false.,'Triply done')
!? nprocs = nprocs0

call GetMem('t3_la','Free','Real',la,NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2))
call GetMem('t3_t1','Free','Real',itmp,noab(1)*nuab(1))
call GetMem('t3_t1a','Free','Real',t1a,it1)
call GetMem('t3_t1b','Free','Real',t1b,it1)
return

9993 format(/1X,'T2-W-T3 contribution from current amplitudes ',D18.10)
9991 format(1X,'T2-F-T3 contribution from current amplitudes ',D18.10)
9990 format(1X,'T1-W-T3 contribution from current amplitudes ',D18.10)
9994 format(/1X,'4th order MBPT tripleexcitation contribution ',D18.10)
9995 format(/1X,'5th order noniterative E[5] ST correction ',D18.10/)

end subroutine T3AMPL_BTI
