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
      subroutine v_size_t3(vblock,nprocs,krem,printkey)
      implicit none
      integer krem,N
      integer t3_size_a,t3_size
      integer isp,vblock,maxnu
      integer nuga,nugc,vblock_isp(2),nprocs,rest
      CHARACTER ich*1
cmp
        integer printkey, tmp
cmp
      INTEGER IOPT,NOAB,NNOAB,NUAB,NNUAB
      COMMON/UHF/NOAB(2),NNOAB(3),NUAB(2),NNUAB(3),ICH(3)
      COMMON/IOIND/IOPT(96)
C number of elementary subprocesses: nugc*nuga*(nuga+1)/2
C                                   + nuga*nugc*(nugc+1)/2
C                                   + nugc**3/6
C                                   + nuga**3/6
C  check this:
      maxnu=max(nuab(1),nuab(2))
      vblock_isp(1)=maxnu/nprocs

cmp
      tmp = 1
      if(maxnu.ge.100)tmp = int((2*nprocs)**(1.0d0/3.0d0))

      do while ((tmp*(tmp*(tmp+1))/2).lt.nprocs)
         tmp = tmp + 1
      enddo
      vblock_isp(1) = maxnu/tmp
cmp

cmp      if(vblock_isp(1).lt.40)then
cmp         if(maxnu.ge.160)vblock_isp(1)=40
cmp      endif

C adjusting to reasonably full last block
      vblock_isp(2)=vblock_isp(1)
      t3_size=krem+1
C
C brute force
C
      t3_size_a=0
      do isp=1,2
         vblock=vblock_isp(isp)+1
         N=noab(isp)+nuab(isp)
C this is a first entry - initialization (makes no harm if reapeated)
         do while (krem.lt.t3_size)
            vblock=vblock-1
!!      write(6,*)'whiblock',vblock,krem,t3_size
            t3_size=0
            nuga=nuab(isp)/vblock
            if((nuga*vblock).lt.nuab(isp))nuga=nuga+1
            nugc=nuab(3-isp)/vblock
            if((nugc*vblock).lt.nuab(3-isp))nugc=nugc+1
c  dummy allocations
            if(nuga.ne.1)then
C             call w_alloc(kab, noab(isp)*vblock*vblock*n ,'kaT3loopb')
               t3_size=t3_size+ noab(isp)*vblock*vblock*n +1
C             call w_alloc(kcb, noab(isp)*vblock*vblock*n ,'kbT3loopb')
               t3_size=t3_size+ noab(isp)*vblock*vblock*n +1

C             call w_alloc(kbc, vblock*vblock*n ,'kcT3loopb')
               t3_size=t3_size+ vblock*vblock*n +1
C             call w_alloc(kbc, vblock*vblock*n ,'kcT3loopb')
               t3_size=t3_size+ vblock*vblock*n +1

            else
C             call w_alloc(kab, noab(isp)*N*nnuab(isp) ,'kaT3loopb')
               t3_size=t3_size+ noab(isp)*N*nnuab(isp) +1

            endif
C          call w_alloc(kac, vblock*vblock*n ,'kcT3loopb')
            t3_size=t3_size+ vblock*vblock*n +1

C          call w_alloc(kca, noab(isp)*vblock*vblock*n ,'kbT3loopb')
            t3_size=t3_size+ noab(isp)*vblock*vblock*n +1

C           call w_alloc(kc,v block*vblock*n ,'kcT3loopb')
            t3_size=t3_size+v block*vblock*n +1

C           call w_alloc(la,n noab(isp)*vblock*n ,'laT3loopb')
            t3_size=t3_size+n noab(isp)*vblock*n +1

C          call w_alloc(lxa, nnoab(3)*vblock*n ,'lbaT3loopb')
            t3_size=t3_size+ nnoab(3)*vblock*n +1

            if(nuga.ne.1)then
C              call w_alloc(lb,n noab(isp)*vblock*n ,'lbT3loopb')
               t3_size=t3_size+n noab(isp)*vblock*n +1

C             call w_alloc(lxb, nnoab(3)*vblock*n ,'labT3loopb')
               t3_size=t3_size+ nnoab(3)*vblock*n +1
            endif
C          call w_alloc(lxc, nnoab(3)*vblock*n ,'lacT3loopb')
            t3_size=t3_size+ nnoab(3)*vblock*n +1

C          call w_alloc(t3a, vblock*vblock*vblock ,'t3aT3loopb')
            t3_size=t3_size+ vblock*vblock*vblock +1
C          call w_alloc(t3b, vblock*vblock*vblock ,'t3bT3loopb')
            t3_size=t3_size+ vblock*vblock*vblock +1

C          call w_alloc(vac, vblock*vblock*nnoab(3) ,'vbcT3loopb')
            t3_size=t3_size+ vblock*vblock*nnoab(3) +1
            if(nuga.ne.1) then
C             call w_alloc(vab, vblock*vblock*nnoab(isp) ,'vabT3loopb')
               t3_size=t3_size+ vblock*vblock*nnoab(isp) +1

C             call w_alloc(vbc, vblock*vblock*nnoab(3) ,'vacT3loopb')
               t3_size=t3_size+ vblock*vblock*nnoab(3) +1
            else
C             call w_alloc(vab, nnoab(isp)*nnuab(isp) ,'vabT3loopb')
               t3_size=t3_size+ nnoab(isp)*nnuab(isp) +1
            endif
C           call w_alloc(mi, noab(isp)*(vblock**3) ,'miT3loopb')
            t3_size=t3_size+ noab(isp)*(vblock**3) +1
C          call w_alloc(mij, N*vblock ,'mijT3loopb')
            t3_size=t3_size+ N*vblock +1
C
         enddo                  ! while
         vblock_isp(isp)=vblock
         if(isp.eq.1)t3_size_a=t3_size
      enddo
      vblock=min0(vblock_isp(1),vblock_isp(2))
      nuga=maxnu/vblock
      if(nuga*vblock.lt.maxnu)nuga=nuga+1
      if(mod(maxnu,vblock).ne.0)
     $     vblock=min0(vblock,maxnu/nuga+mod(maxnu,nuga))
C adjusting to reasonably full last block
      rest=mod(maxnu,vblock)
      do while (.not.((rest.eq.0).or.(rest.gt.(vblock-nuga))))
         vblock=vblock-1
         rest=mod(maxnu,vblock)
      enddo
      do isp=1,2
         t3_size=0
         nuga=nuab(isp)/vblock
         if((nuga*vblock).lt.nuab(isp))nuga=nuga+1
         nugc=nuab(3-isp)/vblock
         if((nugc*vblock).lt.nuab(3-isp))nugc=nugc+1
c  dummy allocations
         if(nuga.ne.1)then
C          call w_alloc(kab, noab(isp)*vblock*vblock*n ,'kaT3loopb')
            t3_size=t3_size+ noab(isp)*vblock*vblock*n +1
C          call w_alloc(kcb, noab(isp)*vblock*vblock*n ,'kbT3loopb')
            t3_size=t3_size+ noab(isp)*vblock*vblock*n +1

C          call w_alloc(kbc, vblock*vblock*n ,'kcT3loopb')
            t3_size=t3_size+ vblock*vblock*n +1
C          call w_alloc(kbc, vblock*vblock*n ,'kcT3loopb')
            t3_size=t3_size+ vblock*vblock*n +1

         else
C          call w_alloc(kab, noab(isp)*N*nnuab(isp) ,'kaT3loopb')
            t3_size=t3_size+ noab(isp)*N*nnuab(isp) +1

         endif
C       call w_alloc(kac, vblock*vblock*n ,'kcT3loopb')
         t3_size=t3_size+ vblock*vblock*n +1

C       call w_alloc(kca, noab(isp)*vblock*vblock*n ,'kbT3loopb')
         t3_size=t3_size+ noab(isp)*vblock*vblock*n +1

C        call w_alloc(kc,v block*vblock*n ,'kcT3loopb')
         t3_size=t3_size+v block*vblock*n +1

C        call w_alloc(la,n noab(isp)*vblock*n ,'laT3loopb')
         t3_size=t3_size+n noab(isp)*vblock*n +1

C       call w_alloc(lxa, nnoab(3)*vblock*n ,'lbaT3loopb')
         t3_size=t3_size+ nnoab(3)*vblock*n +1

         if(nuga.ne.1)then
C           call w_alloc(lb,n noab(isp)*vblock*n ,'lbT3loopb')
            t3_size=t3_size+n noab(isp)*vblock*n +1

C          call w_alloc(lxb, nnoab(3)*vblock*n ,'labT3loopb')
            t3_size=t3_size+ nnoab(3)*vblock*n +1
         endif
C       call w_alloc(lxc, nnoab(3)*vblock*n ,'lacT3loopb')
         t3_size=t3_size+ nnoab(3)*vblock*n +1

C      call w_alloc(t3a, vblock*vblock*vblock ,'t3aT3loopb')
        t3_size=t3_size+ vblock*vblock*vblock +1
C      call w_alloc(t3b, vblock*vblock*vblock ,'t3bT3loopb')
        t3_size=t3_size+ vblock*vblock*vblock +1

C      call w_alloc(vac, vblock*vblock*nnoab(3) ,'vbcT3loopb')
        t3_size=t3_size+ vblock*vblock*nnoab(3) +1
        if(nuga.ne.1) then
C         call w_alloc(vab, vblock*vblock*nnoab(isp) ,'vabT3loopb')
           t3_size=t3_size+ vblock*vblock*nnoab(isp) +1

C         call w_alloc(vbc, vblock*vblock*nnoab(3) ,'vacT3loopb')
           t3_size=t3_size+ vblock*vblock*nnoab(3) +1
        else
C         call w_alloc(vab, nnoab(isp)*nnuab(isp) ,'vabT3loopb')
           t3_size=t3_size+ nnoab(isp)*nnuab(isp) +1
        endif
C       call w_alloc(mi, noab(isp)*(vblock**3) ,'miT3loopb')
        t3_size=t3_size+ noab(isp)*(vblock**3) +1
C      call w_alloc(mij, N*vblock ,'mijT3loopb')
        t3_size=t3_size+ N*vblock +1
        if(isp.eq.1)t3_size_a=t3_size
      enddo
      write(6,*)
      write(6,'(2x,A,I5)')
     $     'Virtual orbitals will be treated in blocks of:',vblock
        if (printkey.ge.10) then
           write(6,'(2x,A,I11,A,I11,A)')'Memory requirement:',
     $           max(t3_size,t3_size_a),' Words;    remaining:',krem-
     $           max(t3_size,t3_size_a),' Words'
        end if
      call xflush(6)
      return
      end
