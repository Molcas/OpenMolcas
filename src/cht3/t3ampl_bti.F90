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

use ChT3_global, only: gen_files, ICH, IOPT, IT, NOAB, NUAB, printkey, run_triples, t3_starta, t3_startb, t3_stopa, t3_stopb, &
                       TCpu, TCpu_l, TCpu0, TWall, TWall_l, TWall0
use Para_Info, only: MyRank, nProcs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: OEH(*), OEP(*)
integer(kind=iwp) :: i, id, isp, it1, IUHF, j, krem, LU(6), n, nga, ngb, ngc, nla, nlb, nuga, nugc, vblock
real(kind=wp) :: ccsdt, ccsdt4, cpu0_aaa, cpu0_aab, e_ccsd, e_scf, energ(4), enx1, tccsd, timerel, times(10), times_parr(10), &
                 totcpu, totwal, wall0_aaa, wall0_aab
#ifdef _MOLCAS_MPP_
real(kind=wp) :: real_buffer(1)
#endif
logical(kind=iwp) :: ifvo, scored, skip
character(len=6) :: FN
integer(kind=iwp), allocatable :: my_tsk(:,:)
real(kind=wp), allocatable :: la(:), t1a(:), t1b(:), tmp(:)
real(kind=wp), external :: ddot_
logical(kind=iwp), external :: rsv_tsk

!? nprocs0 = nprocs
! Uncomment the following to force sequential mode
! nprocs = 1

if (nprocs > 1) then
  write(u6,'(A,i4,A)') ' Parallel run on ',nprocs,' nodes'
else
  write(u6,'(A,i4,A)') ' Serial run'
end if
call xflush(u6)

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

times(:) = Zero
energ(:) = Zero

! calculate the T2-T3 contraction?  ifvo is the answer ! open-shell stuff !
!mp call get3dm('FOC-VO',w,noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
!mp ifvo = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w,1,w,1) > 1.0e-14_wp
ifvo = .false.

! parallelization:  do this on each host >>>>>>> start

! allocates for two arrays of the size T1

it1 = NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)

call mma_allocate(t1a,it1,label='t3_t1a')
call mma_allocate(t1b,it1,label='t3_t1b')

t1a(:) = Zero
t1b(:) = Zero

call mma_allocate(tmp,noab(1)*nuab(1),label='t3_t1')
call mma_allocate(la,it1,label='t3_la')

!mp ! read t1 amplitudes
if (printkey >= 10) write(u6,*) 'Reading ',it1,' t1_alpha, t1_beta amplitudes'
call GetRest_t3(la,tmp,e_ccsd)

!mp !write(u6,*) 'Ze t1 = ',ddot_(it1,la,1,la,1)

! t1a and t1b always remain in the memory -- small fields

! parallelization:   on each host <<<<<<< end  - block done

! Remaining space can be used for the blocking in-core algorithms

call mma_maxDBLE(krem)

if (printkey >= 10) then
  write(u6,*)
  write(u6,*) 'Available Memory before v_size_t3 = ',krem
  call xflush(u6)
end if

! determines the virtual block size
! vblock - the virtual orbitals block size
! look for optimal number with respect nprocs to v_size_t3

write(u6,*)
write(u6,'(2x,A)') 'Starting triply with blocked virtuals algorithm'
if (nprocs > 1) write(u6,*) ' Node Number ',MyRank

! Compute memory pattern.

call v_size_t3(vblock,nprocs,krem,printkey)

IUHF = 1 ! open-shell stuff
if (IOPT(1) == 0) IUHF = 0

if (printkey >= 10) write(u6,'(A,i1)') ' Closed-Shell calculation, IUHF = ',iuhf

! create K(beta-alpha,beta-alpha) and L(alpha-beta,beta-alpha)
!ndup = 0 !?????

!mp - print information on number of steps in loopa and loopb

call check_loops(nuab(1),vblock,nla,nlb)
write(u6,*)
write(u6,'(A,i6)') 'Number of steps in loopa : ',nla
write(u6,'(A,i6)') 'Number of steps in loopb : ',nlb

! checking :

if (t3_stopa > nla) then
  write(u6,*) 'Too many steps in t3_loopa requested in input'
  write(u6,*) 'T3_LOOPA from input : ',t3_starta,t3_stopa,t3_stopa-t3_startb+1
  write(u6,*) 'total T3_LOOPA steps : ',nla
  call abend()
end if

if (t3_stopb > nlb) then
  write(u6,*) 'Too many steps in t3_loopb requested in input'
  write(u6,*) 'T3_LOOPB from input : ',t3_startb,t3_stopb,t3_stopb-t3_startb+1
  write(u6,*) 'total T3_LOOPB steps : ',nlb
  call abend()
end if

!mp
if (gen_files) then

  write(u6,*)
  write(u6,*) 'Creating KMAT, LMAT Scratch Integral Files-----------------------------'
  write(u6,*)
  !mp
  call create_klvab_t3(vblock)
  !mp
  call CWTime(TCpu,TWall)

  if (printkey > 1) then
    write(u6,*)
    write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
    write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if

  TCpu_l = TCpu
  TWall_l = TWall

  !mp

  write(u6,*) 'Create of Integrals done-------------------------------------------------'
  write(u6,*)
  !mp

else
  write(u6,*)
  write(u6,*) 'Skipping KMAT, LMAT scratch integral files generation as requested'
end if
!mp
cpu0_aaa = TCpu
wall0_aaa = TWall
!
if (.not. run_triples) then
  write(u6,*)
  write(u6,*) 'Exiting triples after scratch integral files generation as requested'
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

skip = .false.
do isp=1,1+iuhf

  nuga = nuab(isp)/vblock
  if ((nuga*vblock) < nuab(isp)) nuga = nuga+1
  nugc = nuab(3-isp)/vblock
  if ((nugc*vblock) < nuab(3-isp)) nugc = nugc+1
  ! creates K(alpha-alpha,alpha-alpha) (dummy for closed shell)
  !!if (iuhf /= 0) then   ! open-shell stuff.
  !!  call create_klvaa_t3(w(la),vblock,isp)
  !!end if
  FN = 'KMAT'//ich(isp)//ich(isp)
  call multi_opendir(FN,LU(1))
  FN = 'LMAT'//ich(isp)//ich(isp)
  call multi_opendir(FN,LU(2))

  ! parallelization: files KMATxx LMATxx are assumed to be available
  ! to all tasks on each host

  !                  not for closed shell
  !     alpha-alpha-alpha   or closed shell

  ! - skip t3loopa if needed

  if ((t3_starta < 0) .and. (t3_startb > 0)) then
    write(u6,*)
    write(u6,*) 'Skipping t3loopa on user request '
    write(u6,*)
  else

    if (t3_starta < 0) then
      i = 0
      do nga=1,nuga
        do ngb=1,nga
          do ngc=1,ngb
            i = i+1
          end do
        end do
      end do
    else
      i = t3_stopa-t3_starta+1
    end if
    call mma_allocate(my_tsk,3,i,label='my_tsk')

    i = 0
    if (t3_starta < 0) then
      do nga=1,nuga
        do ngb=1,nga
          do ngc=1,ngb
            i = i+1
            my_tsk(1,i) = nga
            my_tsk(2,i) = ngb
            my_tsk(3,i) = ngc
          end do
        end do
      end do
    else
      do nga=1,nuga
        do ngb=1,nga
          do ngc=1,ngb
            i = i+1
            if ((i >= t3_starta) .and. (i <= t3_stopa)) then
              my_tsk(1,i-t3_starta) = nga
              my_tsk(2,i-t3_starta) = ngb
              my_tsk(3,i-t3_starta) = ngc
            end if
          end do
        end do
      end do
    end if

    if (t3_starta > 0) then ! for correct deallocation
      i = t3_stopa-t3_starta+1
    end if

    write(u6,*)
    write(u6,*) '# of tasks to be parallelized in t3loop a = ',i
    id = 666
    scored = .false.
    call init_tsk(id,i)
    do while (rsv_tsk(id,j))

      nga = my_tsk(1,j)
      ngb = my_tsk(2,j)
      ngc = my_tsk(3,j)

      !mp
      !mp call mma_maxDBLE(maxspace)
      !mp write(u6,*) 'maxspace before ',maxspace
      !mp
      call t3loopa(oeh(noab(1)*(isp-1)+1),oep(nuab(1)*(isp-1)+1),t1a(noab(1)*nuab(1)*(isp-1)+1), &
                   t1b(noab(1)*nuab(1)*(isp-1)+1),nga,ngb,ngc,vblock,energ,isp,LU,ifvo,scored,enx1)
      !mp
      ! update 5th order terms

      n = noab(1)*nuab(1)
      !mp t1a(1:n) = t1a(1:n)+t1a(n+1:2*n)
      ccsdt = Two*(ddot_(n,la,1,t1a,1)+ddot_(n,la,1,t1a(n+1),1))

      !mp
      call CWTime(TCpu,TWall)
      !mp
      if (t3_starta < 0) then
        write(u6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j,nga,ngb,ngc,Two*enx1,ccsdt, &
                                                                 TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
      else
        write(u6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j+t3_starta-1,nga,ngb,ngc,Two*enx1, &
                                                                 ccsdt,TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
      end if

      TCpu_l = TCpu
      TWall_l = TWall
      !mp
      !mp call mma_maxDBLE(maxspace)
      !mp write(u6,*) 'maxspace after ',maxspace
      !mp

    end do
    write(u6,*) 't3loopa finished'
    write(u6,*)
    call Free_tsk(id)
    call mma_deallocate(my_tsk)

    !mp !call gettim(cpu1,wall1)            ! dorob timingy !!!!!!
    call CWTime(TCpu,TWall)
    times(isp) = times(isp)+TCpu-cpu0_aaa
    times(isp+4) = times(isp+4)+TWall-wall0_aaa

    write(u6,'(1X,5A,ES12.4,A,ES12.4,A)') 'Spin case ',ich(isp),ich(isp),ich(isp),' done:',TCpu-cpu0_aaa,' CPU [s]', &
                                          TWall-wall0_aaa,' Wall [s]'
    call xflush(u6)
    !mp
    call CWTime(TCpu,TWall)

    if (printkey > 1) then
      write(u6,*)
      write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
      write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
      write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
      write(u6,*)
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
    write(u6,*)
    write(u6,*) 'Skipping t3loopb on user request '
    write(u6,*)
    skip = .true.
    exit
  end if

  ! alpha-alpha-beta or beta-beta-alpha only UHF
  !!if (IUHF /= 0)then
  !mp!!call gettim(cpu0,wall0)
  FN = 'KMAT'//ich(3-isp)//ich(isp)
  call multi_opendir(FN,LU(3))
  FN = 'LMAT'//ich(3-isp)//ich(isp)
  call multi_opendir(FN,LU(5))

  if (IUHF /= 0) then ! open-shell stuff
    FN = 'KMAT'//ich(isp)//ich(3-isp)
    call multi_opendir(FN,LU(4))
    FN = 'LMAT'//ich(isp)//ich(3-isp)
    call multi_opendir(FN,LU(6))
  else
    LU(4) = LU(3)
    LU(6) = LU(5)
  end if

  !mp
  i = 0
  if (t3_startb < 0) then
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,nugc
          i = i+1
        end do
      end do
    end do
  else
    i = t3_stopb-t3_startb+1
  end if
  call mma_allocate(my_tsk,3,i,label='my_tsk')

  i = 0
  if (t3_startb < 0) then
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,nugc
          i = i+1
          my_tsk(1,i) = nga
          my_tsk(2,i) = ngb
          my_tsk(3,i) = ngc
        end do
      end do
    end do
  else
    do nga=1,nuga
      do ngb=1,nga
        do ngc=1,nugc
          i = i+1
          if ((i >= t3_startb) .and. (i <= t3_stopb)) then
            my_tsk(1,i-t3_startb) = nga
            my_tsk(2,i-t3_startb) = ngb
            my_tsk(3,i-t3_startb) = ngc
          end if
        end do
      end do
    end do
  end if

  if (t3_startb > 0) then
    i = t3_stopb-t3_startb+1
  end if

  write(u6,*)
  write(u6,*) '# of tasks to be parallelized in t3loopb = ',i
  id = 667
  scored = .false.
  call init_tsk(id,i)
  do while (rsv_tsk(id,j))

    nga = my_tsk(1,j)
    ngb = my_tsk(2,j)
    ngc = my_tsk(3,j)
    !mp write(u6,'(A,4(i5,2x))') 'Tsk, nga, ngb, ngc = ',j,nga,ngb,ngc

    !mp
    !mp call mma_maxDBLE(maxspace)
    !mp write(u6,*) 'maxspace before ',maxspace
    !mp
    call t3loopb(oeh,oep,t1a,t1b,nga,ngb,ngc,vblock,energ(3),isp,LU,ifvo,scored,enx1)

    n = noab(1)*nuab(1)
    !mp??? t1a(1:n) = t1a(1:n)+t1a(n+1:2*n)
    ccsdt = Two*(ddot_(n,la,1,t1a,1)+ddot_(n,la,1,t1a(n+1),1))

    !mp
    call CWTime(TCpu,TWall)
    !mp
    if (t3_startb < 0) then
      write(u6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j,nga,ngb,ngc,Two*enx1,ccsdt, &
                                                               TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
    else
      write(u6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))') 'Tsk, nga, ngb, ngc, inc = ',j+t3_startb-1,nga,ngb,ngc,Two*enx1, &
                                                               ccsdt,TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
    end if

    TCpu_l = TCpu
    TWall_l = TWall

    !mp
    !mp call mma_maxDBLE(maxspace)
    !mp write(u6,*) 'maxspace before ',maxspace
    !mp

  end do
  write(u6,*) 't3loopb finished'
  write(u6,*)
  call Free_tsk(id)
  call mma_deallocate(my_tsk)
  !mp write(u6,*) ' energ = ',energ
  !mp
  ! - deallocate arrays in t3loopb
  !
  !mp !call gettim(cpu1,wall1)  ! urob poriadne timingy !!!!
  call CWTime(TCpu,TWall)
  times(2+isp) = times(2+isp)+TCpu-cpu0_aab
  times(6+isp) = times(6+isp)+TWall-wall0_aab
  write(u6,'(1X,5A,ES12.4,A,ES12.4,A)') 'Spin case ',ich(isp),ich(isp),ich(3-isp),' done:',TCpu-cpu0_aab,' CPU [s]', &
                                        TWall-wall0_aab,' Wall [s]'
  call xflush(u6)
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
    write(u6,*)
    write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
    write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
    write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0_wp*TCpu/(TWall-TWall0)
    write(u6,*)
  end if
  TCpu_l = TCpu
  TWall_l = TWall
  !mp
end if

times_parr(:) = times

#ifdef _MOLCAS_MPP_

it1 = NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)
real_buffer(1) = ccsdt
call gadgop(real_buffer,size(real_buffer),'+')
ccsdt = real_buffer(1)
call gadgop(energ,4,'+')
call gadgop(times_parr,10,'+')

#endif

call xflush(u6)

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
!mp !  call n = noab(1)*nuab(1)
!mp !  t1a(1:n) = t1a(1:n)+t1a(n+1:2*n)
!mp !  ccsdt = Two*ddot_(n,la,1,t1a,1)
!mp !else
!mp !  ccsdt = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),la,1,t1a,1)
!mp !  write(u6,*) 'ze co do ... ?'
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
  write(u6,*)
  write(u6,*) '--------------------------------------------------'
  write(u6,*)
  write(u6,*) 'GADGOPed energ & ccsdt values: '
  write(u6,*)
  write(u6,'(A,4(f15.12,1x))') 'energ (t2-w-t3 = 2*e1 + 2*e3) ',energ
  write(u6,*)
  write(u6,'(A,f15.12)') 'Sum                           ',Two*energ(1)+Two*energ(3)
  write(u6,'(A,f15.12)') 'ccsdt (e5th ord. ST)          ',ccsdt
  write(u6,*)
  write(u6,*) '--------------------------------------------------'
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
times(:) = times/60.0_wp ! now in minutes
if (nprocs > 1) times_parr(:) = times_parr/60.0_wp

ccsdt4 = Zero
if (ifvo) then
  write(u6,*) 'ifvo correspond to open-shell system. Exiting'
  call abend()

  !mp !call get3dm('FOC-VO',w(la),noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
  !mp !call map2_21_t3(w(la),w(t1a),nuab(1),noab(1))
  !mp !call map2_21_t3(w(noab(1)*nuab(1)+la),w(t1a+noab(1)*nuab(1)),nuab(2),noab(2))
  ! E4 T2FT3
  !mp !if (IUHF == 0) then
  !mp !  n = noab(1)*nuab(1)
  !mp !  w(t1b:t1b+n-1) = w(t1b:t1b+n-1)+w(t1b+n:n1b+2*n-1)
  !mp !  ccsdt4 = Two*ddot_(noab(1)*nuab(1),w(t1a),1,w(t1b),1)
  !mp !else
  !mp !  ccsdt4 = ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w(t1a),1,w(t1b),1)
  !mp !end if
end if

!RESULT(IT+3,5) = ccsdt4+ccsdt
if (IT == 0) then
  !RESULT(IT+1,5) = tccsd+ccsdt4+ccsdt
  if (ifvo) then
    write(u6,9993) TCCSD
    write(u6,9991) ccsdt4
    write(u6,9990) ccsdt
  end if
  write(u6,9994) TCCSD+ccsdt+ccsdt4
else
  !RESULT(IT+2,5) = TCCSD
  write(u6,9993) TCCSD
  if (ifvo) then
    write(u6,9991) ccsdt4
    write(u6,9990) ccsdt
    write(u6,9995) CCSDT+ccsdt4
  else
    write(u6,9995) CCSDT
  end if
end if

write(u6,*) '*************************************************'
write(u6,*)
write(u6,*) 'Final Results : '
write(u6,*)
write(u6,'(A,f15.12)') '   (T)     corr. = ',tccsd+ccsdt
write(u6,'(A,f15.12)') '   CCSD(T) corr. = ',e_ccsd+tccsd+ccsdt
write(u6,*)
write(u6,*) '*************************************************'
write(u6,*)

! Timing
if (nprocs > 1) write(u6,*) 'Master timings:'
if (times(10) <= Zero) then
  timerel = One
else
  timerel = times(9)/times(10)
end if
write(u6,'(1x,a,2f13.3,a,f6.3)') 'sorting cpu & wall',times(9),times(10),' (min);   cpu/wall=',timerel
if (times(5) <= Zero) then
  timerel = One
else
  timerel = times(1)/times(5)
end if
write(u6,'(1x,a,2f13.3,a,f6.3)') 'aaa     cpu & wall',times(1),times(5),' (min);   cpu/wall=',timerel
if (times(6) <= Zero) then
  timerel = One
else
  timerel = times(2)/times(6)
end if
if (IUHF /= 0) write(u6,'(1x,a,2f13.3,a,f6.3)') 'bbb     cpu & wall',times(2),times(6),' (min);   cpu/wall=',timerel
if (times(7) <= Zero) then
  timerel = One
else
  timerel = times(3)/times(7)
end if
write(u6,'(1x,a,2f13.3,a,f6.3)') 'aab     cpu & wall',times(3),times(7),' (min);   cpu/wall=',timerel
if (times(8) <= Zero) then
  timerel = One
else
  timerel = times(4)/times(8)
end if
if (IUHF /= 0) write(u6,'(1x,a,2f13.3,a,f6.3)') 'bba     cpu & wall',times(4),times(8),' (min);   cpu/wall=',timerel
totcpu = times(9)
totwal = times(10)
do i=1,4
  totcpu = totcpu+times(i)
  totwal = totwal+times(4+i)
end do
if (totwal <= Zero) then
  timerel = One
else
  timerel = totcpu/totwal
end if
write(u6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',totcpu,totwal,' (min);   cpu/wall=',timerel

if (nprocs > 1) then
  ! Parallel timing
  totcpu = times_parr(9)
  totwal = times_parr(10)
  times(9) = times_parr(9)
  do i=1,4
    totcpu = totcpu+times_parr(i)
    totwal = totwal+times_parr(4+i)
  end do
  if (totwal <= Zero) then
    timerel = One
  else
    timerel = totcpu/totwal
  end if
  write(u6,*)
  write(u6,*) 'Aggregate parallel timings:'
  write(u6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',totcpu,totwal,' (min);   cpu/wall=',timerel
end if

! Return
!mp call w_debug(.false.,.false.,'Triply done')
!? nprocs = nprocs0

call mma_deallocate(la)
call mma_deallocate(tmp)
call mma_deallocate(t1a)
call mma_deallocate(t1b)

return

9993 format(/1X,'T2-W-T3 contribution from current amplitudes ',ES18.10)
9991 format(1X,'T2-F-T3 contribution from current amplitudes ',ES18.10)
9990 format(1X,'T1-W-T3 contribution from current amplitudes ',ES18.10)
9994 format(/1X,'4th order MBPT tripleexcitation contribution ',ES18.10)
9995 format(/1X,'5th order noniterative E[5] ST correction ',ES18.10/)

end subroutine T3AMPL_BTI
