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
      Subroutine UTMU( EXCH, N, Z, M1, M2 )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: EXCH, N
      Complex(kind=wp), intent(in)  :: M1(3,EXCH,EXCH)
      Complex(kind=wp), intent(in)  ::  Z(N,N)
      Complex(kind=wp), intent(out) :: M2(3,EXCH,EXCH)
c  local variables:
      Integer          :: L, I,J
      Logical          :: DBG
      Real(kind=wp)    :: dznrm2_,R1,R2
      External         :: dznrm2_
      Complex(kind=wp), allocatable :: TMP(:,:)

      Call qEnter('UTMU')

      DBG=.false.

      If((N<=0).OR.(EXCH<=0)) Then
         Write(6,'(A)') 'in UTMU:   EXCH or N<=0 !!!'
         Write(6,*) 'EXCH=', EXCH
         Write(6,*) 'N   =', N
         Call xFlush(6)
         Call xquit(128)
      End If

      If(N>EXCH)  Then
         Write(6,'(A)') 'in UTMU:   EXCH < N !!!'
         Write(6,*) 'EXCH=', EXCH
         Write(6,*) 'N   =', N
         Write(6,'(A)') 'Nothing is to be done >> Return'
         Call xFlush(6)
         Call xquit(128)
      End If

      R1=dznrm2_(3*EXCH*EXCH,M1,1)
      R2=dznrm2_(  N*N,Z,1)
      If((R1<1.0e-25_wp).or.(R2<1.0e-25_wp)) Then
         Write(6,'(A)') 'in UTMU:   M1 or Z are empty!!!'
         Write(6,*) 'norm(M1)=', R1
         Write(6,*) 'norm(Z )=', R2
        Return
      End If

      If (DBG) Then
         Write(6,'(A)') 'UTMU :: input moment'
         Do i=1,EXCH
            Do j=1,EXCH
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'|M1_l|',j,'>',(M1(l,i,j),l=1,3)
            End Do
         End Do
      End If




      Call mma_allocate(TMP,EXCH,EXCH,'TMP')
      Call zcopy_(3*EXCH*EXCH,(0.0_wp,0.0_wp),0,M2,1)

      If(N.eq.EXCH) Then

         Do L=1,3
            Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
            Call ZGEMM_('C', 'N', EXCH, EXCH, EXCH, (1.0_wp,0.0_wp),
     &                        Z, EXCH,
     &                M1(L,:,:), EXCH,             (0.0_wp,0.0_wp),
     &                      TMP, EXCH )
            Call ZGEMM_('N', 'N', EXCH, EXCH, EXCH, (1.0_wp,0.0_wp),
     &                      TMP, EXCH,
     &                        Z, EXCH,             (0.0_wp,0.0_wp),
     &                M2(L,:,:), EXCH )
         End Do !L

      Else

         Do L=1,3
           Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
           Call ZGEMM_('C', 'N',   N,   N,   N, (1.0_wp,0.0_wp),
     &              Z(1:N,1:N),   N,
     &           M1(L,1:N,1:N),   N,            (0.0_wp,0.0_wp),
     &                     TMP,   N )
           Call ZGEMM_('N', 'N',   N,    N,   N, (1.0_wp,0.0_wp),
     &                     TMP,   N,
     &              Z(1:N,1:N),   N,            (0.0_wp,0.0_wp),
     &           M2(L,1:N,1:N),   N )

           Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
           Call ZGEMM_('C', 'N',   N, EXCH,   N, (1.0_wp,0.0_wp),
     &              Z(1:N,1:N),   N,
     &        M1(L,1:N,1:EXCH),   N,            (0.0_wp,0.0_wp),
     &       TMP(  1:N,1:EXCH),   N )

           Do I=1,N
              Do J=N+1,EXCH
                 M2(L,I,J)=TMP(I,J)
                 M2(L,J,I)=CONJG(TMP(I,J))
              End Do
           End Do
           Do i=N+1,EXCH
              Do j=N+1,EXCH
                 M2(L,i,j)=M1(L,i,j)
              End Do
           End Do
         End Do !L

      End If !N.eq.exch

      If (DBG) Then
         Write(6,'(A)') 'UTMU :: input moment'
         Do i=1,EXCH
            Do j=1,EXCH
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'|M_l|',j,'>',(M1(l,i,j),l=1,3)
            End Do
         End Do
         Write(6,'(A)') 'UTMU :: unitary transformtion matrix'
         Do i=1,N
            Do j=1,N
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'| U |',j,'>',(Z(i,j),l=1,3)
            End Do
         End Do
         Write(6,'(A)') 'UTMU :: output moment'
         Do i=1,EXCH
            Do j=1,EXCH
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'|M_l|',j,'>',(M2(l,i,j),l=1,3)
            End Do
         End Do
      End If
      Call mma_deallocate(TMP)

      Call qExit('UTMU')
      Return
      End Subroutine utmu







      subroutine utmu2( exch, n, z, m )
      ! the same as utmu, except being that the input m is
      ! being transformed.
      implicit none
      integer, parameter           :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      integer, intent(in)            :: exch, n
      complex(kind=wp), intent(inout):: m(3,exch,exch)
      complex(kind=wp), intent(in)   ::  z(n,n)
c  local variables:
      integer          :: l, i,j, i1, j1
      logical          :: dbg
      real(kind=wp)    :: dznrm2_,r1,r2
      external         :: dznrm2_
      complex(kind=wp), allocatable :: tmp(:,:), mtmp(:,:,:)

      call qenter('utmu2')

      dbg=.false.

      if((n<=0).or.(exch<=0)) then
         write(6,'(a)') 'in utmu2:   exch or n<=0 !!!'
         write(6,*) 'exch=', exch
         write(6,*) 'n   =', n
         Call xflush(6)
         Call xquit(128)
      end if

      if(n>exch)  then
         write(6,'(a)') 'in utmu2:   exch < n !!!'
         write(6,*) 'exch=', exch
         write(6,*) 'n   =', n
         write(6,'(a)') 'nothing is to be done >> return'
         Call xflush(6)
         Call xquit(128)
      end if

      r1=dznrm2_(3*exch*exch,m,1)
      r2=dznrm2_(  n*n,z,1)
      if((r1<1.0e-25_wp).or.(r2<1.0e-25_wp)) then
         write(6,'(a)') 'in utmu2:   m or z are empty!!!'
         write(6,*) 'norm(m)=', r1
         write(6,*) 'norm(z)=', r2
         Return
      end if

      if (dbg) then
         write(6,'(a)') 'utmu2 :: input moment'
         do i=1,exch
            do j=1,exch
               write(6,'(a,i3,a,i3,a,3(2e16.8,2x))')
     &                         '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
            end do
         end do
      end if


      call mma_allocate(tmp,exch,exch,'tmp')
      if(n<exch) then
         call mma_allocate(mtmp,3,(exch-n),(exch-n),'mtmp')
         ! save the part which is not altered
         do i=n+1,exch
            do j=n+1,exch
               i1=i-n
               j1=j-n
               do l=1,3
                  mtmp(l,i1,j1)=m(l,i,j)
               end do
            end do
         end do
      end if


      if(n==exch) then

         do l=1,3
            call zcopy_( exch*exch,(0.0_wp,0.0_wp),0,tmp,1)
            call zgemm_('c', 'n', exch, exch, exch, (1.0_wp,0.0_wp),
     &                        z, exch,
     &                 m(l,:,:), exch,             (0.0_wp,0.0_wp),
     &                      tmp, exch )
            call zcopy_(exch*exch,(0.0_wp,0.0_wp),0,m(l,:,:),1)
            call zgemm_('n', 'n', exch, exch, exch, (1.0_wp,0.0_wp),
     &                      tmp, exch,
     &                        z, exch,             (0.0_wp,0.0_wp),
     &                m(l,:,:), exch )
         end do !l

      else

         do l=1,3
           call zcopy_( exch*exch,(0.0_wp,0.0_wp),0,tmp,1)
           call zgemm_('c', 'n',   n,   n,   n, (1.0_wp,0.0_wp),
     &              z(1:n,1:n),   n,
     &            m(l,1:n,1:n),   n,            (0.0_wp,0.0_wp),
     &                     tmp,   n )
            call zcopy_(n*n,(0.0_wp,0.0_wp),0,m(l,1:n,1:n),1)

           call zgemm_('n', 'n',   n,    n,   n, (1.0_wp,0.0_wp),
     &                     tmp,   n,
     &              z(1:n,1:n),   n,            (0.0_wp,0.0_wp),
     &            m(l,1:n,1:n),   n )

           call zcopy_( exch*exch,(0.0_wp,0.0_wp),0,tmp,1)
           call zgemm_('c', 'n',   n, exch,   n, (1.0_wp,0.0_wp),
     &              z(1:n,1:n),   n,
     &         m(l,1:n,1:exch),   n,            (0.0_wp,0.0_wp),
     &       tmp(  1:n,1:exch),   n )

           do i=1,n
              do j=n+1,exch
                  m(l,i,j)=tmp(i,j)
                  m(l,j,i)=conjg(tmp(i,j))
              end do
           end do
           do i=n+1,exch
              do j=n+1,exch
                 i1=i-n
                 j1=j-n
                 m(l,i,j)=mtmp(l,i1,j1)
              end do
           end do
         end do !l
      end if !n.eq.exch

      if (dbg) then
         write(6,'(a)') 'utmu2 :: unitary transformtion matrix'
         do i=1,n
            do j=1,n
               write(6,'(a,i3,a,i3,a,3(2e16.8,2x))')
     &                         '<',i,'| u |',j,'>',z(i,j)
            end do
         end do
         write(6,'(a)') 'utmu2 :: output moment'
         do i=1,exch
            do j=1,exch
               write(6,'(a,i3,a,i3,a,3(2e16.8,2x))')
     &                         '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
            end do
         end do
      end if

      if(n<exch) call mma_deallocate(mtmp)
      call mma_deallocate(tmp)

      if(dbg) write(6,*) 'at the end of utmu2'
      if(dbg) call prmom('utmu2, moment',m,n)

      call qexit('utmu2')
      return
      end subroutine utmu2






      Subroutine UTMUL( EXCH, N, Z, ML )
      ! the same as UTMU, except being that the input M is
      ! being transformed and only one projection is done at a time.
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)            :: EXCH, N
      Complex(kind=wp), intent(inout):: ML(EXCH,EXCH) ! one projection is done
      Complex(kind=wp), intent(in)   ::  Z(N,N)
c  local variables:
      Integer          :: I,J,i1,j1
      Logical          :: DBG
      Real(kind=wp)    :: dznrm2_,R1,R2
      External         :: dznrm2_
      Complex(kind=wp), allocatable :: TMP(:,:), MTMP(:,:)

      Call qEnter('UTMUL')

      DBG=.false.

      If((N<=0).OR.(EXCH<=0)) Then
         Write(6,'(A)') 'in UTMU2:   EXCH or N<=0 !!!'
         Write(6,*) 'EXCH=', EXCH
         Write(6,*) 'N   =', N
         Call xquit(128)
      End If


      If(N>EXCH)  Then
         Write(6,'(A)') 'in UTMU2:   EXCH < N !!!'
         Write(6,*) 'EXCH=', EXCH
         Write(6,*) 'N   =', N
         Write(6,'(A)') 'Nothing is to be done >> Return'
         Call xquit(128)
      End If

      R1=dznrm2_(EXCH*EXCH,ML,1)
      R2=dznrm2_(N*N,Z,1)
      If((R1<1.0e-25_wp).or.(R2<1.0e-25_wp)) Then
         Write(6,'(A)') 'in UTMU2:   M or Z are empty!!!'
         Write(6,*) 'norm(M)=', R1
         Write(6,*) 'norm(Z)=', R2
         Call xquit(128)
      End If

      If (DBG) Then
         Write(6,'(A)') 'UTMU :: input moment'
         Do i=1,EXCH
            Do j=1,EXCH
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'|ML|',j,'>',ML(i,j)
            End Do
         End Do
      End If


      Call mma_allocate(TMP,EXCH,EXCH,'TMP')
      If(N<EXCH) Then
         Call mma_allocate(MTMP,(EXCH-N),(EXCH-N),'MTMP')
         ! save the part which is not altered
         Do i=N+1,EXCH
            Do j=N+1,EXCH
               i1=i-N
               j1=j-N
               MTMP(i1,j1)=ML(i,j)
            End Do
         End Do
      End If


      If(N==EXCH) Then
         Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
         Call zgemm_('C', 'N', EXCH, EXCH, EXCH, (1.0_wp,0.0_wp),
     &                     Z, EXCH,
     &                    ML, EXCH,             (0.0_wp,0.0_wp),
     &                   TMP, EXCH )
         Call zcopy_(EXCH*EXCH,(0.0_wp,0.0_wp),0,ML(:,:),1)
         Call zgemm_('N', 'N', EXCH, EXCH, EXCH, (1.0_wp,0.0_wp),
     &                   TMP, EXCH,
     &                     Z, EXCH,             (0.0_wp,0.0_wp),
     &                    ML, EXCH )

      Else

         Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
         Call ZGEMM_('C', 'N',   N,   N,   N, (1.0_wp,0.0_wp),
     &            Z(1:N,1:N),   N,
     &           ML(1:N,1:N),   N,            (0.0_wp,0.0_wp),
     &                   TMP,   N )
         Call zcopy_(N*N,(0.0_wp,0.0_wp),0,ML(1:N,1:N),1)

         Call ZGEMM_('N', 'N',   N,    N,   N, (1.0_wp,0.0_wp),
     &                   TMP,   N,
     &            Z(1:N,1:N),   N,            (0.0_wp,0.0_wp),
     &          ML(1:N,1:N),   N )

         Call zcopy_( EXCH*EXCH,(0.0_wp,0.0_wp),0,TMP,1)
         Call ZGEMM_('C', 'N',   N, EXCH,   N, (1.0_wp,0.0_wp),
     &            Z(1:N,1:N),   N,
     &        ML(1:N,1:EXCH),   N,            (0.0_wp,0.0_wp),
     &     TMP(  1:N,1:EXCH),   N )

         Do I=1,N
            Do J=N+1,EXCH
                ML(I,J)=TMP(I,J)
                ML(J,I)=CONJG(TMP(I,J))
            End Do
         End Do
         Do i=N+1,EXCH
            Do j=N+1,EXCH
               i1=i-N
               j1=j-N
               ML(i,j)=MTMP(i1,j1)
            End Do
         End Do

      End If !N.eq.exch

      If (DBG) Then
         Write(6,'(A)') 'UTMU :: unitary transformtion matrix'
         Do i=1,N
            Do j=1,N
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'| U |',j,'>',Z(i,j)
            End Do
         End Do
         Write(6,'(A)') 'UTMU :: output moment'
         Do i=1,EXCH
            Do j=1,EXCH
               Write(6,'(A,i3,A,i3,A,3(2E16.8,2x))')
     &                         '<',i,'|ML|',j,'>',ML(i,j)
            End Do
         End Do
      End If

      If(N<EXCH) Call mma_deallocate(MTMP)
      Call mma_deallocate(TMP)

      Call qExit('UTMUL')
      Return
      End Subroutine utmul
