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
C
C----------------------------------------------------------------------|
C
      subroutine dkh_xpx(n,dkord,xord,vord,EL,ES,OL,OS,Ep,E0,dkcof,cc,
     &                   wr,rw,t1,t2,t3,t4,or,ro,e,rer,
     &                   or_,ro_,e_,rer_,s1,s2,wsav )
C
C Evaluate DKH transformation in moment space for property operator
C
      implicit none
C Input
      integer n,dkord,xord,vord
      Real*8 EL(n,n),ES(n,n),OL(n,n),OS(n,n),Ep(n),E0(n)
      Real*8 dkcof(*),cc(*)
      Real*8 wsav(n,n,*)
C Output : EL overwritten by the transformed property matrix
C Temp
      Real*8 wr(n,n),rw(n,n),t1(n,n),t2(n,n),t3(n,n),
     & t4(n,n),or(n,n,*),ro(n,n,*),e(n,n,*),rer(n,n,*),or_(n,n,*),
     & ro_(n,n,*),e_(n,n,*),rer_(n,n,*),s1(n,n,*),s2(n,n,*)
      integer i,j,k,ord,cou,ks,ioe
      logical ifodd
C
C Copy initial matrices
C
      do i=1,n
        do j=1,n
          e(j,i,1)   = EL(j,i)
          rer(j,i,1) = ES(j,i)
          or(j,i,1)  = OL(j,i)
          ro(j,i,1)  = OS(j,i)
        end do
      end do
C
      cou=0
      do 10 ord=1,xord
        do k=1,vord
          do i=1,n
            do j=1,n
              or_(j,i,k)  = 0.d0
              ro_(j,i,k)  = 0.d0
              e_(j,i,k)   = 0.d0
              rer_(j,i,k) = 0.d0
            end do
          end do
        end do
C
C Set up W(ord) matrix, copy from saved set generated from dkh_ham
C
        do i=1,n
          do j=1,n
            wr(j,i) = wsav(j,i,ord*2-1)
            rw(j,i) = wsav(j,i,ord*2  )
          end do
        end do
C
C Calculate [W(ord)->E(1-...)/O(1-...)]
C   note that no odd term will be eliminated
C
        do 30 ks=1,xord+1
C         ! W1 only apply to E1/O1
CDP          if(ord.eq.1.and.ks.ge.2) goto 30
          if(ord.gt.1.or.ks.lt.2)then
          do 40 ioe=1,2
            k=ks
            if(ioe.eq.1)then
              ifodd=.true.
            else
              ifodd=.false.
            end if
C           ! copy initial matrix
            if(ifodd)then
              do i=1,n
                do j=1,n
                  t1(j,i)    = or(j,i,k)
                  t2(j,i)    = ro(j,i,k)
                  or_(j,i,k) = or_(j,i,k) + t1(j,i)
                  ro_(j,i,k) = ro_(j,i,k) + t2(j,i)
                end do
              end do
            else
              do i=1,n
                do j=1,n
                  t1(j,i)     = e(j,i,k)
                  t2(j,i)     = rer(j,i,k)
                  e_(j,i,k)   = e_(j,i,k)   + t1(j,i)
                  rer_(j,i,k) = rer_(j,i,k) + t2(j,i)
                end do
              end do
            end if
            call dkh_wgene(n,ord,k,xord+1,ifodd,dkcof,wr,rw,t1,t2,
     &                     e_,rer_,or_,ro_,cou,s1,s2,t3,t4)
C         ! cycle for odd/even operator
40        end do
          end if
C       ! cycle for ks
30      end do
C       ! copy
        do k=1,xord+1
          do i=1,n
            do j=1,n
              or(j,i,k)    = or_(j,i,k)
              ro(j,i,k)    = ro_(j,i,k)
              e(j,i,k)     = e_(j,i,k)
              rer(j,i,k)   = rer_(j,i,k)
            end do
          end do
        end do
C     ! cycle for ord
10    end do
C
C Sum over all k<=xord+1 terms
C   +1 appears because X itself is not counted as an order
C
      do k=2,xord+1
        do i=1,n
          do j=1,n
            EL(j,i) = EL(j,i) + e(j,i,k)
          end do
        end do
      end do
C      write(6,*) "DKHX",xord," Total matmul",cou
C
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(dkord)
        call Unused_real_array(Ep)
        call Unused_real_array(E0)
        call Unused_real_array(cc)
      end if
      end
