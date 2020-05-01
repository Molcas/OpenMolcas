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
      Subroutine mean_field( EXCH, N, H, X,Y,Z, zJ, T, W, thrs,
     &                       dM, sM, ST, dbg )

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in)          :: EXCH, N
      Real(kind=8), intent(in)    :: H, X,Y,Z, zJ, T
      Real(kind=8), intent(in)    :: thrs
      Real(kind=8), intent(in)    :: W(EXCH)
      Complex(kind=8), intent(in) :: DM(3,EXCH,EXCH)
      Complex(kind=8), intent(in) :: SM(3,EXCH,EXCH)
      Logical, intent(in)          :: dbg
      ! output
      Real(kind=8), intent(out)   :: ST(3)

      Integer                      :: iopt

      iopt=2

      If (iopt .eq. 1) Then

          If (dbg) Write(6,'(A)') 'mean_field: enter mean_field_exch'
          Call  mean_field_exch( N, H, X,Y,Z, zJ, T, thrs, W(1:N),
     &                           DM(1:3,1:N,1:N), SM(1:3,1:N,1:N), ST )
          If (dbg) Write(6,'(A)') 'mean_field: exit mean_field_exch'

      Else If (iopt .eq. 2) Then

          If (dbg) Write(6,'(A)') 'mean_field: enter mean_field_all'
          Call  mean_field_all( EXCH, N, H, X,Y,Z, zJ, T, thrs,
     &                          W, DM, SM, ST )
          If (dbg) Write(6,'(A)') 'mean_field: exit mean_field_all'

      Else

         Write(6,'(A)') 'MEAN_FIELD:  iopt is not defined:', iopt

      End If

      Return
      End subroutine mean_field







      Subroutine mean_field_exch( N, H, X,Y,Z, zJ, T, thrs, W,
     &                               dM, sM, ST )
!     this Subroutine computes the mean field of neighboring spins ST(3)
!     for zJ .ne. 0.0_wp
!     using ONLY Zeeman basis (N)

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in)       :: N
      Real(kind=8), intent(in) :: H, X,Y,Z, zJ, T, W(N)
      Complex(kind=8), intent(in) :: DM(3,N,N), SM(3,N,N)
      ! output
      Real(kind=8), intent(out) :: ST(3)

#include "stdalloc.fh"
      ! local variables:
      Integer :: i, l, mxIter, iter
      Logical :: DBG
      Real(kind=8)    :: WM(N), S(3), ZB, SCHK, THRS, SL(3)
      Complex(kind=8) :: ZM(N,N), SZ(3,N,N)
      Real(kind=8), allocatable :: RWORK(:)
      Complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)

      DBG=.false.

      MxIter=100
      THRS = 1.0D-12
      SCHK = 0.0_wp
      S = 0.0_wp
      ST = 0.0_wp

      ! temporary arrays used in ZEEM_SA:
      Call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
      Call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
      Call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
      Call mma_allocate(W_c,N,'ZEEM_W_c')

      ! zero everything:
      Call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
      Call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
      Call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
      Call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)
      ! determine first the average spin of neighboring
      ! molecules for each temperature point
      Do iter=1,mxIter
          WM(1:N) = 0.0_wp
          ZM(1:N,1:N) =(0.0_wp,0.0_wp)
          ! build and diagonalize the Zeeman Hamiltonian (size N x N)
          ! for the field direction (X,Y,Z) and strength (H)
          Call ZEEM_SA( N, H, X,Y,Z,  W(1:N), dM(1:3,1:N,1:N),
     &               sM(1:3,1:N,1:N), ST(1:3), zJ,
     &               WM(1:N), ZM(1:N,1:N),
     &               DBG, RWORK, HZEE, WORK, W_c )
          !Write(6,'(A,3ES16.8)') 'WM:',(WM(l),l=1,N)
          !Write(6,'(A,9(2ES16.8,2x))') 'ZM:',((ZM(l,i),l=1,N),i=1,N)
          ! transform the spin momenta to the Zeeman eigenstate basis
          SZ(1:3,1:N,1:N)=(0.0_wp,0.0_wp)
          Call UTMU( N, N, ZM(1:N,1:N), SM(1:3,1:N,1:N),
     &                                  SZ(1:3,1:N,1:N) )
!         compute the spin magnetization vector at this temperature (T):
          If(iter==mxIter) Then
            SL=0.0_wp
            SL(1)=S(1)
            SL(2)=S(2)
            SL(3)=S(3)
          End If
          Do l=1,3
            S(l)=0.0_wp
            ZB=0.0_wp
            Call calcmagn1( N, WM(1:N), SZ(l,1:N,1:N), T, S(l), ZB )
          End Do


          ! check if average spin is converged
          SCHK=0.0_wp
          Do L=1,3
            SCHK = SCHK + SQRT( (S(L)-ST(L))*(S(L)-ST(L)) )
          End Do

          If(DBG) Write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)')
     &                    'ST:   End of iteration',iter,':',
     &                    (ST(l),l=1,3),'DIFF:',(S(l)-ST(l) ,l=1,3)
          If(DBG) Write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)')
     &                    'SCHK: End of iteration',iter,':',SCHK

          ! decide to continue the iterative process or exit
          If( SCHK < THRS ) Then
             Go To 1001
          Else
             Do l=1,3
                ST(l)=0.0_wp
                ST(l)=S(l)
             End Do
          End If
        End Do ! iter


        Write(6,'(A, ES24.14)') 'This message shows that the '//
     &                          'average spin did NOT converge '//
     &                          'after 100 iterations. Temp.(in K)=',T
        Write(6,'(A,4ES24.14)') 'Field: (X, Y, Z), and (H):',X,Y,Z,H
        Write(6,'(A,4ES24.14)') 'Last values of the average spin: '//
     &                          '(Sx,Sy,Sz):',(ST(i),i=1,3)
        Write(6,'(A,4ES24.14)') 'Last values of the deviation:    '//
     &                          '          :',(SL(i)-ST(i),i=1,3)
        Write(6,'(A,4ES24.14)') 'Absolute value of the deviation: '//
     &                          '          :',SCHK
        Write(6,'(A,4ES24.14)') 'Convergence threshold:    THRS = '//
     &                          '          :',THRS
        Write(6,'(A         )') 'The program will continue, using '//
     &                          'the last value of the average spin'


        Return
1001  Continue

      ! deallocate temporary data:
      Call mma_deallocate(RWORK)
      Call mma_deallocate(HZEE)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W_c)

      Return
      End subroutine mean_field_exch











      Subroutine mean_field_all( EXCH, N, H, X,Y,Z, zJ, T, thrs, W,
     &                           dM, SM, ST)
!     this Subroutine computes the mean field of neighboring spins ST(3)
!     for zJ .ne. 0.0_wp
!     using ONLY Zeeman basis (N)

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in)          :: EXCH, N
      Real(kind=8), intent(in)    :: H, X,Y,Z, zJ, T, W(EXCH)
      Complex(kind=8), intent(in) :: DM(3,EXCH,EXCH), SM(3,EXCH,EXCH)
      ! output
      Real(kind=8), intent(out) :: ST(3)

#include "stdalloc.fh"
      ! local variables:
      Integer :: i, l, mxIter, iter
      Logical :: DBG
      Real(kind=8)    :: WM(EXCH), S(3), ZB, SCHK, THRS, SL(3)
      Complex(kind=8) :: ZM(N,N), SZ(3,EXCH,EXCH)
      Real(kind=8), allocatable :: RWORK(:)
      Complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)

      DBG=.false.

      MxIter=100
      THRS = 1.0D-12
      SCHK = 0.0_wp
      S = 0.0_wp
      ST = 0.0_wp

      ! temporary arrays used in ZEEM_SA:
      Call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
      Call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
      Call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
      Call mma_allocate(W_c,N,'ZEEM_W_c')

      ! zero everything:
      Call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
      Call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
      Call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
      Call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)
      ! determine first the average spin of neighboring
      ! molecules for each temperature point
      Do iter=1,mxIter
          WM(1:N) = 0.0_wp
          ZM(1:N,1:N) =(0.0_wp,0.0_wp)
          ! build and diagonalize the Zeeman Hamiltonian (size N x N)
          ! for the field direction (X,Y,Z) and strength (H)
          Call ZEEM_SA( N, H, X,Y,Z,  W(1:N), dM(1:3,1:N,1:N),
     &               sM(1:3,1:N,1:N), ST(1:3), zJ,
     &               WM(1:N), ZM(1:N,1:N),
     &               DBG, RWORK, HZEE, WORK, W_c  )
          If(N.ne.EXCH) Then
            Do i=N+1,EXCH
                WM(i)=W(i)
            End Do
          End If

!         transform the spin momenta to the Zeeman eigenstate basis
          Call zcopy_(3*EXCH*EXCH,[(0.0_wp,0.0_wp)],0,SZ,1)
          Call UTMU( EXCH, N, ZM(1:N,1:N), SM, SZ )
!         compute the spin magnetization vector at this temperature (T):
          If(iter==mxIter) Then
            SL=0.0_wp
            SL(1)=S(1)
            SL(2)=S(2)
            SL(3)=S(3)
          End If


          ZB=0.0_wp
          S(1:3)=0.0_wp
          If (N.eq.EXCH) Then
            Do L=1,3
                Call calcmagn1( EXCH, WM, SZ(l,1:EXCH,1:EXCH),
     &                          T, S(l), ZB)
            End Do
          Else
            Do L=1,3
               Call calcmagn2( EXCH, N, WM, T, H,
     &                         SZ, X, Y, Z,
     &                         L, S(l), ZB)
             End Do
          End If


          ! check if average spin is converged
          SCHK=0.0_wp
          Do L=1,3
            SCHK = SCHK + SQRT( (S(L)-ST(L))*(S(L)-ST(L)) )
          End Do

          If(DBG) Write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)')
     &                    'ST:   End of iteration',iter,':',
     &                    (ST(l),l=1,3),'DIFF:',(S(l)-ST(l) ,l=1,3)
          If(DBG) Write(6,'(A,i4,1x,A,3ES20.10,2x,A,3ES20.10)')
     &                    'SCHK: End of iteration',iter,':',SCHK

          ! decide to continue the iterative process or exit
          If( SCHK < THRS ) Then
             Go To 1001
          Else
             Do l=1,3
                ST(l)=0.0_wp
                ST(l)=S(l)
             End Do
          End If
        End Do ! iter


        Write(6,'(A, ES24.14)') 'This message shows that the '//
     &                          'average spin did NOT converge '//
     &                          'after 100 iterations. Temp.(in K)=',T
        Write(6,'(A,4ES24.14)') 'Field: (X, Y, Z), and (H):',X,Y,Z,H
        Write(6,'(A,4ES24.14)') 'Last values of the average spin: '//
     &                          '(Sx,Sy,Sz):',(ST(i),i=1,3)
        Write(6,'(A,4ES24.14)') 'Last values of the deviation:    '//
     &                          '          :',(SL(i)-ST(i),i=1,3)
        Write(6,'(A,4ES24.14)') 'Absolute value of the deviation: '//
     &                          '          :',SCHK
        Write(6,'(A,4ES24.14)') 'Convergence threshold:    THRS = '//
     &                          '          :',THRS
        Write(6,'(A         )') 'The program will continue, using '//
     &                          'the last value of the average spin'


        Return
1001  Continue

      ! deallocate temporary data:
      Call mma_deallocate(RWORK)
      Call mma_deallocate(HZEE)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W_c)


      Return
      End subroutine mean_field_all
