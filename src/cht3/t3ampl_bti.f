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

      subroutine barf(a)
      character*(*) a
      write(6,*) a
      call abend
      end


cmp      SUBROUTINE T3AMPL_BTI(W,OEH,OEP)
      SUBROUTINE T3AMPL_BTI(OEH,OEP)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    *** PROGRAM TRIPLY-UHF/RHF - REDUCED DIMENSION and timing - DRIVER
C
C    CALCULATION OF TRIPLE EXCITATION CONTRIBUTION TO
C    THE 4TH ORDER CORRELATION ENERGY OF MB RSPT AND CC-TRIPLES
C
C    Provides T3AMPL_BLOCKEDIMPROVED:  Accelerated algorithm.
c
c     Combines triply w/ trick by JN and limited I/Os by buffered
c     blocking scheme by PV. Blocks the virtual space to pieces
c     optimal either for parallelization or in using the available
c     memory
c     Involves also more robust and more reliable memory allocator.
c
c History
c     Buffered blocking scheme: PV, LAOG Grenoble, 16 april 2003.
c     First version w/ trick: JN, June 12, 2003
c     Parallel version: PV, 15 oct 2003.
c     Implemented integer offsets, PV, 14 may 2004.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer IUHF,LU(6)
      integer i,nuga,nugc,nga,ngb,ngc,vblock, it1
      integer NOAB,NNOAB,NUAB,NNUAB,iopt,iout,isp,krem
      real*8 OEH(*),OEP(*),ddot_,ccsdt,ccsdt4,energ(4),tccsd,
     $     ENSCF, RESULT,times(10),
     $     times_parr(10), totcpu, totwal, timerel
c     real*8 cpu0,cpu1,wall0,wall1
      character ICH*1, FN*6
      logical ifvo
      integer IHW
      common /hermit_addr/ IHW
      integer la,t1a,t1b
        logical lastcall,scored
      integer IT,ITLAST,NBF,NOMX,NU,MX2,NNO,NNU,NUO,NSO
cmp
c       integer maxspace
        integer nla,nlb
        real*8 enx1
cmp
      COMMON/PARAM/IT,ITLAST,NBF,NOMX,NU,MX2,NNO,NNU,NUO,NSO
      COMMON/IOIND/IOPT(96)
      COMMON/ENERGY/ENSCF, RESULT(102,5)
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      LOGICAL    RRT1,CSDAT,ORTHO,DETAIL
      COMMON/WHAT1/RRT1,CSDAT,ORTHO,DETAIL
cmp
#include "cht3_casy.fh"
#include "WrkSpc.fh"
cmp
#include "ndisk.fh"
#include "dupfiles.fh"
cmp!      include 'task_info_inc'
cmp!      include 'ws_conn_inc'
cmp
#include "para_info.fh"
#include "cht3_ccsd1.fh"
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
#endif
!???#ifdef _MOLCAS_MPP_
!???      integer mytype, krems(nhmx)
!???#endif
cmp!!      integer me, err, len_trim
!?      integer nprocs0
!?cmp      common /my_mpi_world_com/ me, nprocs
      integer LENPAR
      common/LENPAR_com/LENPAR
c
        integer itmp
        integer imy_tsk
        integer id,j
        logical rsv_tsk
        external rsv_tsk
c
        real*8 e_ccsd, e_scf
        real*8 cpu0_aaa,wall0_aaa,cpu0_aab,wall0_aab
c
!?      nprocs0=nprocs
c Uncomment the following to force sequential mode
!      nprocs=1

      IOUT=IOPT(14)

      if (nprocs.gt.1) then
         write(6,'(A,i4,A)') ' Parallel run on ',nprocs,' nodes'
      else
         write(6,'(A,i4,A)') ' Serial run'
      endif
      call xflush(6)
c
C t^ijk_abc . D^abc_ijk =
C
C alpha-alpha-alpha (K_ab^ir+K_ab^jr-K_ac^kr)*(L_rc^jk+L_rc^ki+L_rc^ij)
C or beta-beta-beta   + permutations bca, acb
C        - M_abc^iki -M_abc^iij -M_abc^jjk+M_abc^jij-M_abc^kki+M_abc^kjk
C         + permutations bca cab
C
C        where    M_abc^ijk=K_ab^ir*L_ra^jk
C
C  strategy:   blocks the number of virtuals to pieces by ng
C              then everything in the memory!
c
c free everything from W(1) including Hermit space
c (PV) done for slaves in calling routine slave_T3AMPL_BTI
c
C parallelization: distribute common UHF
c (PV) done for slaves in calling routine slave_T3AMPL_BTI
c
cmp!!      call gettim(cpu0,wall0)

      do I=1,6
         LU(I)=90+I
      enddo

      call zeroma(times,1,10)
      call zeroma(energ,1,4)
      IHW=0

c
C calculate the T2-T3 contraction?  ifvo is the answer ! open-shell stuff !
cmp         call get3dm('FOC-VO',w,noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
cmp         ifvo=(ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w,1,w,1).gt.1d-14)
        ifvo=.false.

C
C parallelization:  do this on each host >>>>>>> start
C
C allocates for two arrays of the size T1

        it1=NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)

        call GetMem('t3_ampl_t1a','Allo','Real',t1a,it1)
        call GetMem('t3_ampl_t1b','Allo','Real',t1b,it1)

        call zeroma(Work(t1a),1,it1)
        call zeroma(Work(t1b),1,it1)

         call GetMem('t3_ampl_t1','Allo','Real',itmp,
     & noab(1)*nuab(1))
         call GetMem('t3_ampl_la','Allo','Real',la,it1)

cmp!    read t1 amplitudes
        if (printkey.ge.10) then
           write (6,*) 'Reading ',it1,' t1_alpha, t1_beta amplitudes'
        end if
         call GetRest_t3 (Work(la),Work(itmp),e_ccsd)
c
cmp!         write (6,*) 'Ze t1 = ',
cmp!     & ddot_(it1,Work(la),1,Work(la),1)
c
C t1a and t1b always remain in the memory -- small fields
c
c parallelization:   on each host <<<<<<< end  - block done
c

C Remaining space can be used for the blocking in-core algorithms

        Call GetMem('t3_ampl_la','Max','Real',krem,krem)

        if (printkey.ge.10) then
           write(6,*)
           write(6,*) 'Available Memory before v_size_t3 = ', krem
           call xflush(6)
        end if
C
C determines the virtual block size
C vblock - the virtual orbitals block size
C look for optimal number with respect nprocs to v_size_t3
C
        write (6,*)
        write(6,'(2x,A)')
     &             'Starting triply with blocked virtuals algorithm'
        if (nprocs.gt.1) then
           write (6,*) ' Node Number ',MyRank
        end if
c
c Compute memory pattern.
c
        call v_size_t3(vblock,nprocs,krem,printkey)

        IUHF=1 ! open-shell stuff
        IF(IOPT(76).EQ.0)IUHF=0

        if (printkey.ge.10) then
           write (6,'(A,i1)') ' Closed-Shell calculation, IUHF = ',iuhf
        end if
c
C  create K(beta-alpha,beta-alpha) and L(alpha-beta,beta-alpha)
        ndup=0 !?????

cmp - print information on number of steps in loopa and loopb
c
        call check_loops(nuab(1),vblock,nla,nlb)
        write (6,*)
        write (6,'(A,i6)') 'Number of steps in loopa : ',nla
        write (6,'(A,i6)') 'Number of steps in loopb : ',nlb
c
c checking :
c
        if (t3_stopa.gt.nla) then
          write (6,*) 'Too many steps in t3_loopa requested in input'
          write (6,*) 'T3_LOOPA from input : ',
     & t3_starta,t3_stopa,t3_stopa-t3_startb+1
          write (6,*) 'total T3_LOOPA steps : ',nla
          call abend()
        end if
c
        if (t3_stopb.gt.nlb) then
          write (6,*) 'Too many steps in t3_loopb requested in input'
          write (6,*) 'T3_LOOPB from input : ',
     & t3_startb,t3_stopb,t3_stopb-t3_startb+1
          write (6,*) 'total T3_LOOPB steps : ',nlb
          call abend()
        end if
c
cmp
         if (gen_files) then
c
             write (6,*)
             write (6,*) 'Creating KMAT, LMAT Scratch Integral Files'//
     & '-----------------------------'
             write (6,*)
cmp
             call create_klvab_t3(vblock)
cmp
        Call CWTime(TCpu,TWall)
c
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if

        TCpu_l=TCpu
        TWall_l=TWall

cmp

             write (6,*) 'Create of Integrals done',
     & '-------------------------------------------------'
             write (6,*)
cmp
c
         else
             write (6,*)
             write (6,*) 'Skipping KMAT, LMAT scratch integral',
     & ' files generation as requested'
         end if
cmp
        cpu0_aaa=TCpu
        wall0_aaa=TWall
c
        if (.not.run_triples) then
        write (6,*)
        write (6,*) 'Exiting triples after scratch integral',
     & ' files generation as requested'
        return
        !? call abend()
        end if
cmp

C  creates K(alpha-alpha,alpha-alpha) (dummy for closed shell)
C  (master only)
cmp!
cmp!  open-shell stuff
cmp!
cmp!         if(iuhf.ne.0)then
cmp!            DO isp=1,1+iuhf
cmp!               call create_klvaa_t3(w(la),vblock,isp)
cmp!            enddo
cmp!         endif

C
C parallelization: files KMATxy LMATxy x=A,B  y=A,B (BA - closed shell)
C                  generated on each host

c Sorting and initialization timing

cmp!      call gettim(cpu1,wall1)  !  prerob na poriadne timingy
      Call CWTime(TCpu,TWall)
      times(9)=times(9)+TCpu-TCpu0
      times(10)=times(10)+TWall-TWall0

      DO isp=1,1+iuhf

         nuga=nuab(isp)/vblock
         if((nuga*vblock).lt.nuab(isp))nuga=nuga+1
         nugc=nuab(3-isp)/vblock
         if((nugc*vblock).lt.nuab(3-isp))nugc=nugc+1
C  creates K(alpha-alpha,alpha-alpha) (dummy for closed shell)
c!         if(iuhf.ne.0)then   ! open-shell stuff.
c!            call create_klvaa_t3(w(la),vblock,isp)
c!         endif
         FN(5:5)=ich(isp)
         FN(6:6)=ich(isp)
         FN(1:4)='KMAT'
         call multi_opendir(FN,LU(1))
         FN(1:4)='LMAT'
         call multi_opendir(FN,LU(2))
C
C parallelization: files KMATxx LMATxx are assumed to be available
C to all tasks on each host
C
C                  not for closed shell
C     alpha-alpha-alpha   or closed shell
c
c - skip t3loopa if needed
c
        if ((t3_starta.lt.0).and.(t3_startb.gt.0)) then
          write (6,*)
          write (6,*) 'Skipping t3loopa on user request '
          write (6,*)
          goto 494
        end if
c
         i=0
         do nga=1,nuga
         do ngb=1,nga
         do ngc=1,ngb
         i=i+1
         end do
         end do
         end do
c
        if (t3_starta.lt.0) then
           call GetMem ('Imy_tsk','Allo','Inte',imy_tsk,
     &                  3*i)
        else
           call GetMem ('Imy_tsk','Allo','Inte',imy_tsk,
     &                  (t3_stopa-t3_starta+1)*3)
        end if
c
        if (t3_starta.lt.0) then
         i=0
         do nga=1,nuga
          do ngb=1,nga
           do ngc=1,ngb
              i=i+1
              iWork(imy_tsk-3+3*i)=nga
              iWork(imy_tsk-3+3*i+1)=ngb
              iWork(imy_tsk-3+3*i+2)=ngc
           end do
          end do
         end do
        else
         i=0
         do nga=1,nuga
          do ngb=1,nga
           do ngc=1,ngb
              i=i+1
c
              if ((i.ge.t3_starta).and.(i.le.t3_stopa)) then
                 iWork(imy_tsk-3+3*(i-t3_starta+1))=nga
                 iWork(imy_tsk-3+3*(i-t3_starta+1)+1)=ngb
                 iWork(imy_tsk-3+3*(i-t3_starta+1)+2)=ngc
              end if
c
         end do
         end do
         end do
        end if
c
        if (t3_starta.gt.0) then ! for correct deallocation
            i=t3_stopa-t3_starta+1
        end if
c
        write (6,*)
        write (6,*) '# of tasks to be parallelized in t3loop a = ',i
        id=666
        lastcall=.false.
        scored=.false.
        call init_tsk(id,i)
98      if (.not. rsv_tsk(id,j)) goto 99
c
        nga=iWork(imy_tsk-3+3*j)
        ngb=iWork(imy_tsk-3+3*j+1)
        ngc=iWork(imy_tsk-3+3*j+2)
c
cmp
cmp        Call GetMem('(T)','Max','Real',maxspace,maxspace)
cmp        write (*,*) 'maxspace before ',maxspace
cmp
        call t3loopa(
     $               oeh(noab(1)*(isp-1)+1),
     $               oep(nuab(1)*(isp-1)+1),
     $               Work(t1a+noab(1)*nuab(1)*(isp-1)),
     $               Work(t1b+noab(1)*nuab(1)*(isp-1)),
     $               nga,ngb,ngc,vblock,energ,isp,LU,ifvo,
     &               lastcall,scored,j,enx1)
cmp
c update 5th order terms
c

cmp      call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,
cmp     $Work(t1a),1,noab(1)*nuab(1))

        call daxpy_((noab(1)*nuab(1)), 1.0d0,
     & Work(t1a+noab(1)*nuab(1)), 1,
     & Work(t1a), 1)
        ccsdt=2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
        call daxpy_((noab(1)*nuab(1)), -1.0d0,
     & Work(t1a+noab(1)*nuab(1)), 1,
     & Work(t1a), 1)
c
cmp
        Call CWTime(TCpu,TWall)
cmp
         if (t3_starta.lt.0) then
        write (6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))')
     & 'Tsk, nga, ngb, ngc, inc = ',
     & j,
     & nga,ngb,ngc,2.0d0*enx1,ccsdt,
     & TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
         else
        write (6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))')
     & 'Tsk, nga, ngb, ngc, inc = ',
     & j+t3_starta-1,
     & nga,ngb,ngc,2.0d0*enx1,ccsdt,
     & TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
         end if

         TCpu_l=TCpu
         TWall_l=TWall
cmp
cmp        Call GetMem('(T)','Max','Real',maxspace,maxspace)
cmp        write (*,*) 'maxspace after ',maxspace
cmp
c
        goto 98
99      continue
        write (6,*) 't3loopa finished'
        write (6,*)
        Call Free_tsk(id)
        call GetMem ('Imy_tsk','Free','Inte',imy_tsk,3*i)
c
cmp!         call gettim(cpu1,wall1)            ! dorob timingy !!!!!!
         Call CWTime(TCpu,TWall)
         times(isp)=times(isp)+TCpu-cpu0_aaa
         times(isp+4)=times(isp+4)+TWall-wall0_aaa

         write(6,'(1X,5A,D12.4,A,D12.4,A)')
     $        'Spin case ',ich(isp),ich(isp),ich(isp),' done:',
     $        (TCpu-cpu0_aaa),' CPU [s]',
     &        (TWall-wall0_aaa),' Wall [s]'
         call xflush(6)
cmp
        Call CWTime(TCpu,TWall)
c
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if

        TCpu_l=TCpu
        TWall_l=TWall
cmp

494        continue

cmp
        cpu0_aab=TCpu
        wall0_aab=TWall
cmp

        if ((t3_starta.gt.0).and.(t3_startb.lt.0)) then
          write (6,*)
          write (6,*) 'Skipping t3loopb on user request '
          write (6,*)
          goto 495
        end if
c

C     alpha-alpha-beta or beta-beta-alpha only UHF
         !!if(IUHF.ne.0)then
cmp!!            call gettim(cpu0,wall0)
            FN(5:5)=ich(3-isp)
            FN(6:6)=ich(isp)
            FN(1:4)='KMAT'
            call multi_opendir(FN,LU(3))
            FN(1:4)='LMAT'
            call multi_opendir(FN,LU(5))

            IF(IUHF.NE.0)THEN ! open-shell stuff
              FN(5:5)=ich(isp)
              FN(6:6)=ich(3-isp)
              FN(1:4)='KMAT'
              call multi_opendir(FN,LU(4))
              FN(1:4)='LMAT'
              call multi_opendir(FN,LU(6))
            ELSE
              LU(4)=LU(3)
              LU(6)=LU(5)
            ENDIF

cmp
         i=0
         do nga=1,nuga
         do ngb=1,nga
         do ngc=1,nugc
         i=i+1
         end do
         end do
         end do
c
        if (t3_startb.lt.0) then
        call GetMem ('Imy_tsk','Allo','Inte',imy_tsk,3*i)
        else
        call GetMem ('Imy_tsk','Allo','Inte',imy_tsk,
     & (t3_stopb-t3_startb+1)*3)
        end if
c
         i=0
         if (t3_startb.lt.0) then
          do nga=1,nuga
           do ngb=1,nga
            do ngc=1,nugc
               i=i+1
               iWork(imy_tsk-3+3*i)=nga
               iWork(imy_tsk-3+3*i+1)=ngb
               iWork(imy_tsk-3+3*i+2)=ngc
            end do
           end do
          end do
         else
          do nga=1,nuga
           do ngb=1,nga
            do ngc=1,nugc
               i=i+1
               if ((i.ge.t3_startb).and.(i.le.t3_stopb)) then
                  iWork(imy_tsk-3+3*(i-t3_startb+1))=nga
                  iWork(imy_tsk-3+3*(i-t3_startb+1)+1)=ngb
                  iWork(imy_tsk-3+3*(i-t3_startb+1)+2)=ngc
               end if
            end do
           end do
          end do
         end if
c
        if (t3_startb.gt.0) then
         i=t3_stopb-t3_startb+1
        end if
c
        write (6,*)
        write (6,*) '# of tasks to be parallelized in t3loopb = ',i
        id=667
        lastcall=.false.
        scored=.false.
        call init_tsk(id,i)
198      if (.not. rsv_tsk(id,j)) goto 199
c
        nga=iWork(imy_tsk-3+3*j)
        ngb=iWork(imy_tsk-3+3*j+1)
        ngc=iWork(imy_tsk-3+3*j+2)
cmp        write (6,'(A,4(i5,2x))') 'Tsk, nga, ngb, ngc = ',j,nga,ngb,ngc
c
cmp
cmp        Call GetMem('(T)','Max','Real',maxspace,maxspace)
cmp        write (*,*) 'maxspace before ',maxspace
cmp
        call t3loopb(oeh,oep,Work(t1a),Work(t1b),
     $                 nga,ngb,ngc,vblock,energ(3),isp,LU,ifvo,
     &                 lastcall,scored,j,enx1)
c
cmp???      call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,
cmp???     $Work(t1a),1,noab(1)*nuab(1))
        call daxpy_((noab(1)*nuab(1)), 1.0d0,
     & Work(t1a+noab(1)*nuab(1)), 1,
     & Work(t1a), 1)
        ccsdt=2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
        call daxpy_((noab(1)*nuab(1)), -1.0d0,
     & Work(t1a+noab(1)*nuab(1)), 1,
     & Work(t1a), 1)
c
cmp
        Call CWTime(TCpu,TWall)
cmp
         if (t3_startb.lt.0) then
        write (6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))')
     & 'Tsk, nga, ngb, ngc, inc = ',
     & j,
     & nga,ngb,ngc,2.0d0*enx1,ccsdt,
     & TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
         else
        write (6,'(A,i5,1x,3(i3,1x),2(f21.19,1x),2(f8.1,A,1x))')
     & 'Tsk, nga, ngb, ngc, inc = ',
     & j+t3_startb-1,
     & nga,ngb,ngc,2.0d0*enx1,ccsdt,
     & TCpu-TCpu_l,' CPU [s]',TWall-TWall_l,' Wall [s]'
         end if

         TCpu_l=TCpu
         TWall_l=TWall
c
cmp
cmp        Call GetMem('(T)','Max','Real',maxspace,maxspace)
cmp        write (*,*) 'maxspace before ',maxspace
cmp
c
        goto 198
199      continue
        write (6,*) 't3loopb finished'
        write (6,*)
        Call Free_tsk(id)
        call GetMem ('Imy_tsk','Free','Inte',imy_tsk,3*i)
cmp        write (6,*) ' energ = ',energ
cmp
c - deallocate arrays in t3loopb
c
cmp!            call gettim(cpu1,wall1)  ! urob poriadne timingy !!!!
            Call CWTime(TCpu,TWall)
            times(2+isp)=times(2+isp)+TCpu-cpu0_aab
            times(6+isp)=times(6+isp)+TWall-wall0_aab
            write(6,'(1X,5A,D12.4,A,D12.4,A)')
     $           'Spin case ',ich(isp),ich(isp),ich(3-isp),' done:',
     $           (TCpu-cpu0_aab),' CPU [s]',
     &           (TWall-wall0_aab),' Wall [s]'
            call xflush(6)
            do i=3,6
               close(LU(i))
            enddo
         !!endif   !! IUHF
         close (LU(1))
         close (LU(2))
      enddo                     ! isp
cmp!!      call gettim(cpu0,wall0) ! urob timingy!
cmp
        Call CWTime(TCpu,TWall)
c
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
cmp

495        continue

        do i=1,10
           times_parr(i)=times(i)
        end do

#ifdef _MOLCAS_MPP_

        it1=NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2)
        block
            real*8 :: real_buffer(1)
            real_buffer(1) = ccsdt
            call gadgop (real_buffer, size(real_buffer), '+')
            ccsdt = real_buffer(1)
        end block
        call gadgop (energ,4,'+')
        call gadgop (times_parr,10,'+')
c
#endif

      call xflush(6)

c Add MPI_Reduce contribution to sorting timing
cmp!cmp!      call gettim(cpu1,wall1)
cmp!      Call CWTime(TCpu,TWall)
cmp!      times(9)=times(9)+TCpu-TCpu0
cmp!      times(10)=times(10)+TWall-TWall0

C
C   T(CCSD)
C
      IF(IUHF.EQ.0)then
      energ(2)=energ(1)
      energ(4)=energ(3)
      endif
      tccsd=energ(1)+energ(2)+energ(3)+energ(4)

cmp!      IF(IUHF.EQ.0)then ! open-shell stuff
cmp!      call vadd(Work(t1a),1,Work(t1a+noab(1)*nuab(1)),1,
cmp!     $Work(t1a),1,noab(1)*nuab(1))
cmp!      ccsdt=2.d0*ddot_(noab(1)*nuab(1),Work(la),1,Work(t1a),1)
cmp!       else
cmp!      ccsdt=ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),Work(la),1,Work(t1a),1)
cmp!        write (*,*) 'ze co do ... ?'
cmp!        stop
cmp!        endif

cmp! for MOLCAS verify
        Call Get_dScalar('SCF energy',e_scf)
        Call Add_Info('CHT3ene', tccsd + ccsdt, 1, 6)
        Call Add_Info('E_CHT3', tccsd + ccsdt, 1, 6)
        Call Add_Info('E_HYPE', e_scf + e_ccsd + tccsd + ccsdt, 1, 6)
c    for NUMERICAL_GRADIENTS
        Call Put_cArray('Relax Method','CHT3    ',8)
        Call Store_Energies(1,tccsd+ccsdt+e_ccsd+e_scf,1)
cmp!

c Slaves have finished here...
      if (MyRank.ne.0) then
         return
      endif
        if (printkey.gt.1) then
        write (6,*)
        write (6,*) '--------------------------------------------------'
        write (6,*)
        write (6,*) 'GADGOPed energ & ccsdt values: '
        write (6,*)
        write (6,'(A,4(f15.12,1x))') 'energ (t2-w-t3 = 2*e1 + 2*e3) ',
     & energ
        write (6,*)
        write (6,'(A,f15.12)')      'Sum                           ',
     & 2.0d0*energ(1)+2.0d0*energ(3)
        write (6,'(A,f15.12)')      'ccsdt (e5th ord. ST)          ',
     & ccsdt
        write (6,*)
        write (6,*) '--------------------------------------------------'
        end if


c Contents of times array:
c
c                  CPU      WALL
c     sorting
c     aaa           1        4
c     bbb           2        6
c     aab           3        7
c     bba           4        8
c
c Master has only to sum up the results
      do i=1,10
         times(i)=times(i)/60.d0 ! now in minutes
      enddo
      if (nprocs.gt.1) then
         do i=1,10
            times_parr(i)=times_parr(i)/60.d0
         enddo
      endif

      if(ifvo)then
        write (6,*) 'ifvo correspond to open-shell system. Exiting'
        call abend()

cmp!         call get3dm('FOC-VO',w(la),noab(1)*nuab(1)+noab(2)*nuab(2),1,0)
cmp!         call transm(w(la),w(t1a),nuab(1),noab(1))
cmp!         call transm(w(noab(1)*nuab(1)+la),
cmp!     $        w(noab(1)*nuab(1)+t1a),nuab(2),noab(2))
c  E4 T2FT3
cmp!      IF(IUHF.EQ.0)then
cmp!      call vadd(w(t1b),1,w(t1b+noab(1)*nuab(1)),1,
cmp!     $w(t1b),1,noab(1)*nuab(1))
cmp!         ccsdt4=2.0*ddot_(noab(1)*nuab(1),w(t1a),1,w(t1b),1)
cmp!      else
cmp!         ccsdt4=ddot_(noab(1)*nuab(1)+noab(2)*nuab(2),w(t1a),1,w(t1b),1)
cmp!      endif
      else
         ccsdt4=0.d0
      endif

      RESULT(IT+3,5)=ccsdt4+ccsdt
      IF(IT.EQ.0)THEN
         RESULT(IT+1,5)=tccsd+ccsdt4+ccsdt
         if(ifvo)then
            WRITE(6,9993)TCCSD
            WRITE(6,9991)ccsdt4
            WRITE(6,9990)ccsdt
         endif
         WRITE(6,9994)TCCSD+ccsdt+ccsdt4
      ELSE
         RESULT(IT+2,5)=TCCSD
         WRITE(6,9993)TCCSD
         if(ifvo)then
            WRITE(6,9991)ccsdt4
            WRITE(6,9990)ccsdt
            WRITE(6,9995)CCSDT+ccsdt4
         else
            WRITE(6,9995)CCSDT
         endif
      ENDIF

        write (6,*) '*************************************************'
        write (6,*)
        write (6,*) 'Final Results : '
        write (6,*)
        write (6,'(A,f15.12)') '   (T)     corr. = ',tccsd + ccsdt
        write (6,'(A,f15.12)') '   CCSD(T) corr. = ',
     &                               e_ccsd + tccsd + ccsdt
        write (6,*)
        write (6,*) '*************************************************'
        write (6,*)

c Timing
      if (nprocs.gt.1) write(6,*) 'Master timings:'
      if (times(10).le.0.0d0) then
         timerel=1.0d0
      else
         timerel=times(9)/times(10)
      end if
      write(6,'(1x,a,2f13.3,a,f6.3)')'sorting cpu & wall',
     $     times(9),times(10),' (min);   cpu/wall=',timerel
      if (times(5).le.0.0d0) then
         timerel=1.0d0
      else
         timerel=times(1)/times(5)
      end if
      write(6,'(1x,a,2f13.3,a,f6.3)')'aaa     cpu & wall',
     $     times(1),times(5),' (min);   cpu/wall=',timerel
      if (times(6).le.0.0d0) then
         timerel=1.0d0
      else
         timerel=times(2)/times(6)
      end if
      IF(IUHF.NE.0)
     $write(6,'(1x,a,2f13.3,a,f6.3)')'bbb     cpu & wall',
     $     times(2),times(6),' (min);   cpu/wall=',timerel
      if (times(7).le.0.0d0) then
         timerel=1.0d0
      else
         timerel=times(3)/times(7)
      end if
      write(6,'(1x,a,2f13.3,a,f6.3)')'aab     cpu & wall',
     $     times(3),times(7),' (min);   cpu/wall=',timerel
      if (times(8).le.0.0d0) then
         timerel=1.0d0
      else
         timerel=times(4)/times(8)
      end if
      IF(IUHF.NE.0)
     $write(6,'(1x,a,2f13.3,a,f6.3)')'bba     cpu & wall',
     $     times(4),times(8),' (min);   cpu/wall=',timerel
      totcpu=times(9)
      totwal=times(10)
      do i=1,4
         totcpu=totcpu+times(i)
         totwal=totwal+times(4+i)
      enddo
      if (totwal.le.0.0d0) then
         timerel=1.0d0
      else
         timerel=totcpu/totwal
      end if
      write(6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',
     $     totcpu,totwal,' (min);   cpu/wall=',timerel
      if (nprocs.le.1) goto 2000

c Parallel timing
      totcpu=times_parr(9)
      totwal=times_parr(10)
      times(9)=times_parr(9)
      do i=1,4
         totcpu=totcpu+times_parr(i)
         totwal=totwal+times_parr(4+i)
      enddo
      if (totwal.le.0.0d0) then
         timerel=1.0d0
      else
         timerel=totcpu/totwal
      end if
      write(6,*)
      write(6,*) 'Aggregate parallel timings:'
      write(6,'(1x,a,2f13.3,a,f6.3)') 'total   cpu & wall',
     $     totcpu,totwal,' (min);   cpu/wall=',timerel

c Return
 2000         continue
cmp        call w_debug(.false.,.false.,'Triply done')
!?      nprocs=nprocs0

        call GetMem('t3_ampl_la','Free','Real',la,
     & NUAB(1)*NOAB(1)+NUAB(2)*NOAB(2))
        call GetMem('t3_ampl_t1','Free','Real',itmp,
     & noab(1)*nuab(1))
        call GetMem('t3_ampl_t1a','Free','Real',t1a,it1)
        call GetMem('t3_ampl_t1b','Free','Real',t1b,it1)
      return

 9993 FORMAT(/1X,'T2-W-T3 contribution from current amplitudes ',D18.10)
 9991 FORMAT(1X,'T2-F-T3 contribution from current amplitudes ',D18.10)
 9990 FORMAT(1X,'T1-W-T3 contribution from current amplitudes ',D18.10)
 9994 FORMAT (/1X,'4th order MBPT tripleexcitation contribution '
     $     ,D18.10)
 9995 FORMAT(/1X,'5th order noniterative E[5] ST correction '
     $     ,D18.10/)
c9997 FORMAT(/1X,'Total ST correction with ROHF reference '
c    $     ,D18.10/)
      end
c
