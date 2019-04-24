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
      Subroutine JKQPar(N1,N2,HEXCH,Jpar)
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: N1, N2
      Complex(kind=wp), intent(in)  :: HEXCH(N1,N1,N2,N2)
      Complex(kind=wp), intent(out) :: Jpar( (N1-1), -(N1-1):(N1-1),
     &                                       (N2-1), -(N2-1):(N2-1) )
! local variables
      Integer                       :: k1,k2,q1,q2,ipr,i1,i2,j1,j2,i
      Real(kind=wp)                 :: F, THRS, R, knm(12,0:12)
      Complex(kind=wp)              :: F1, F2, F12, CI
      Complex(kind=wp)              :: trace_exch, trace, fact
      Complex(kind=wp), allocatable :: O1(:,:), W1(:,:)
      Complex(kind=wp), allocatable :: O2(:,:), W2(:,:)
      Complex(kind=wp), allocatable :: O1_O2(:,:,:,:) !(N1,N1,N2,N2)
      Complex(kind=wp), allocatable :: O1_W2(:,:,:,:) !(N1,N1,N2,N2)
      Complex(kind=wp), allocatable :: W1_O2(:,:,:,:) !(N1,N1,N2,N2)
      Complex(kind=wp), allocatable :: W1_W2(:,:,:,:) !(N1,N1,N2,N2)

      Complex(kind=wp) :: Jcc( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )
      Complex(kind=wp) :: Jcs( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )
      Complex(kind=wp) :: Jsc( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )
      Complex(kind=wp) :: Jss( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )

      External         :: trace_exch, trace
      Logical          :: DBG
!      Real(kind=wp)    :: cm_to_MHz
!-----------------------------------------------------------------------
      knm=0.0_wp
      Call Set_knm(knm)
      Call mma_allocate(O1,N1,N1,'O1')
      Call mma_allocate(O2,N2,N2,'O2')
      Call mma_allocate(W1,N1,N1,'W1')
      Call mma_allocate(W2,N2,N2,'W2')
      Call mma_allocate(O1_O2,N1,N1,N2,N2,'O1_O2')
      Call mma_allocate(O1_W2,N1,N1,N2,N2,'O1_W2')
      Call mma_allocate(W1_O2,N1,N1,N2,N2,'W1_O2')
      Call mma_allocate(W1_W2,N1,N1,N2,N2,'W1_W2')

!-----------------------------------------------------------------------

c      cm_to_MHz=29979.2458_wp
      DBG=.false.
      ipr=1

c we need to project now the HEXCH: in products of ITOs
c  HEXCH = SUM(rank1,proj1,rank2,proj2)=
c         { B(rank1,proj1,rank2,proj2)* O1(rank1,proj1) * O2(rank2,proj2) }
c Naoya definition
c eq.40 in DoI:10.1103/PhysRevB.91.174438
      Jpar=(0.0_wp,0.0_wp)
      Do k1=1,N1-1
      Do q1=0,k1

        Do k2=1,N2-1
        Do q2=0,k2

        Call zcopy_(N1*N1,[(0.0_wp,0.0_wp)],0,O1,1)
        Call zcopy_(N2*N2,[(0.0_wp,0.0_wp)],0,O2,1)
        Call zcopy_(N1*N1,[(0.0_wp,0.0_wp)],0,W1,1)
        Call zcopy_(N2*N2,[(0.0_wp,0.0_wp)],0,W2,1)
        ! get the ITOs for each site:
        Call Stewens_matrixel(k1,q1, N1, O1,W1, ipr)
        Call Stewens_matrixel(k2,q2, N2, O2,W2, ipr)

        If((k1==1).and.(q1==1)) Then
          Call pa_prmat('JKQPar:  O1(k1=1,q1=1): ', O1, n1)
          Call pa_prmat('JKQPar:  W1(k1=1,q1=1): ', W1, n1)
        Else If((k1==1).and.(q1==0)) Then
          Call pa_prmat('JKQPar:  O1(k1=1,q1=0): ', O1, n1)
          Call pa_prmat('JKQPar:  W1(k1=1,q1=0): ', W1, n1)
        End If

c      If (DBG) Then
c         Write(6,*)
c         Write(6,'(A)') '--------------------------------'
c         Write(6,*)
c         Write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O1  k1 = ',k1,'q1 =',q1
c         Write(6,*)
c         Do i1=1,N1
c           Write(6,'(20(2ES14.7,1x))') ( (0.5_wp, 0.0_wp)*
c     &            (  CMPLX((-1)**q1)*W1(i1,j1) + O1(i1,j1)  ), j1=1,N1)
c         End Do
c         Write(6,*)
c         Write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W1  k1 = ',k1,'q1 =',q1
c         Write(6,*)
c         Do i1=1,N1
c           Write(6,'(20(2ES14.7,1x))') ( (0.0_wp,-0.5_wp)*
c     &            (  CMPLX((-1)**q1)*W1(i1,j1) - O1(i1,j1)  ), j1=1,N1)
c         End Do
c         Write(6,*)
c         Write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O2  k2 = ',k2,'q2 =',q2
c         Write(6,*)
c         Do i2=1,N2
c           Write(6,'(20(2ES14.7,1x))') ( (0.5_wp, 0.0_wp)*
c     &            (  CMPLX((-1)**q2)*W2(i2,j2) + O2(i2,j2)  ), j2=1,N2)
c         End Do
c         Write(6,*)
c         Write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W2  k2 = ',k2,'q2 =',q2
c         Write(6,*)
c         Do i2=1,N2
c           Write(6,'(20(2ES14.7,1x))') ( (0.0_wp,-0.5_wp)*
c     &            (  CMPLX((-1)**q2)*W2(i2,j2) - O2(i2,j2)  ), j2=1,N2)
c         End Do
c      End If

      ! Build 4 coupled tensor products:
      ! O1-O2
      ! O1-W2
      ! W1-O2
      ! W1-W2
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,O1_O2,1)
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,O1_W2,1)
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,W1_O2,1)
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,W1_W2,1)
      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2
              O1_O2(i1,j1,i2,j2) = O1(i1,j1) * O2(i2,j2)
              O1_W2(i1,j1,i2,j2) = O1(i1,j1) * W2(i2,j2)
              W1_O2(i1,j1,i2,j2) = W1(i1,j1) * O2(i2,j2)
              W1_W2(i1,j1,i2,j2) = W1(i1,j1) * W2(i2,j2)
            End Do
          End Do
        End Do
      End Do

      !    SP_HZFSO=trace(nDIMcf,HCF,DIP_O)
      !    SP_HZFSW=trace(nDIMcf,HCF,DIP_W)
      !    SP_MOW  =trace(nDIMcf,DIP_O,DIP_W)

      !    B(N,-M)=SP_HZFSO/SP_MOW
      !    B(N, M)=SP_HZFSW/SP_MOW


      FACT=(0.0_wp,0.0_wp)
      FACT=trace(N1,O1,W1)*trace(N2,O2,W2)
c      Write(6,'(A,4I3,2F20.13)') 'k1,q1,k2,q2,FACT=',k1,q1,k2,q2,FACT


      Jpar( k1,-q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, O1_O2) / FACT
      Jpar( k1, q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, W1_O2) / FACT
      Jpar( k1,-q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, O1_W2) / FACT
      Jpar( k1, q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, W1_W2) / FACT

c      Jpar( k1,-q1, k2,-q2 ) =  trace_exch(N1,N2,HEXCH, O1_O2)
c     &                        * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
c      Jpar( k1, q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, W1_O2)
c     &                        * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
c      Jpar( k1,-q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, O1_W2)
c     &                        * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
c      Jpar( k1, q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, W1_W2)
c     &                        * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
c



      If(DBG) Then
         If ((q1.eq.0).and.(q2.eq.0)) Then
            Write(6,'(A,2ES18.10)') 'FACT = ', FACT
            If(ABS(Jpar(k1,-q1,k2,-q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',q1,',',k2,',',q2,')=',
     &                  Jpar(  k1,     0,    k2,     0  )
            End If
         Else If((q1.eq.0).and.(q2.ne.0)) Then
            Write(6,'(A,2ES18.10)') 'FACT = ', FACT
            If(ABS(Jpar(k1,0,k2,-q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',q1,',',k2,',',-q2,')=',
     &                  Jpar(  k1,     0,    k2,    -q2)
            End If
            If(ABS(Jpar(k1,0,k2, q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',q1,',',k2,',',q2,')=',
     &                  Jpar(  k1,     0,    k2,    q2  )
            End If
         Else If((q1.ne.0).and.(q2.eq.0)) Then
            Write(6,'(A,2ES18.10)') 'FACT = ', FACT
            If(ABS(Jpar(k1,0,k2,-q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',-q1,',',k2,',',q2,')=',
     &                  Jpar(  k1,    -q1,    k2,     0  )
            End If
            If(ABS(Jpar(k1,0,k2, q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',', q1,',',k2,',',q2,')=',
     &                  Jpar(  k1,     q1,    k2,     0  )
            End If
         Else If((q1.ne.0).and.(q2.ne.0)) Then
            If(ABS(Jpar(k1,-q1,k2,-q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',-q1,',',k2,',',-q2,')=',
     &                  Jpar(  k1,    -q1,    k2,    -q2  )
            End If
            If(ABS(Jpar(k1, q1,k2,-q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',', q1,',',k2,',',-q2,')=',
     &                  Jpar(  k1,     q1,    k2,    -q2  )
            End If
            If(ABS(Jpar(k1,-q1,k2,q2)).gt. 0.5d-13) Then
               Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',',-q1,',',k2,',', q2,')=',
     &                  Jpar(  k1,    -q1,    k2,     q2  )
            End If
            If(ABS(Jpar(k1, q1,k2, q2)).gt. 0.5d-13) Then
            Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                 'Jpar(',k1,',', q1,',',k2,',', q2,')=',
     &                  Jpar(  k1,     q1,    k2,     q2  )
            End If
         End If
      End If ! DBG


        End Do
        End Do
      End Do
      End Do

      Jcc( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )=(0.0_wp, 0.0_wp)
      Jcs( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )=(0.0_wp, 0.0_wp)
      Jsc( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )=(0.0_wp, 0.0_wp)
      Jss( (N1-1),0:(N1-1), (N2-1),0:(N2-1) )=(0.0_wp, 0.0_wp)


      Do k1=1,N1-1
      Do k2=1,N2-1
      Do q1=0,k1
      Do q2=0,k2

        F1 =CMPLX((-1)**(   -q1),0,wp)
        F2 =CMPLX((-1)**(   -q2),0,wp)
        F12=CMPLX((-1)**(-q1-q2),0,wp)
        CI =(0.0_wp,-1.0_wp)

        Jcc(k1,q1,k2,q2) =          Jpar(k1, q1,k2, q2)
     &                     +   F2 * Jpar(k1, q1,k2,-q2)
     &                     +   F1 * Jpar(k1,-q1,k2, q2)
     &                     +  F12 * Jpar(k1,-q1,k2,-q2)

        Jss(k1,q1,k2,q2) =          Jpar(k1, q1,k2, q2)
     &                     -   F2 * Jpar(k1, q1,k2,-q2)
     &                     -   F1 * Jpar(k1,-q1,k2, q2)
     &                     +  F12 * Jpar(k1,-q1,k2,-q2)

        Jcs(k1,q1,k2,q2) =(         Jpar(k1, q1,k2, q2)
     &                     -   F2 * Jpar(k1, q1,k2,-q2)
     &                     +   F1 * Jpar(k1,-q1,k2, q2)
     &                     -  F12 * Jpar(k1,-q1,k2,-q2) ) * CI

        Jsc(k1,q1,k2,q2) =(         Jpar(k1, q1,k2, q2)
     &                     +   F2 * Jpar(k1, q1,k2,-q2)
     &                     -   F1 * Jpar(k1,-q1,k2, q2)
     &                     -  F12 * Jpar(k1,-q1,k2,-q2) ) * CI

        If(DBG) THEN
           If(ABS(Jcc(k1,q1,k2,q2)).gt. 0.5d-13) Then
           Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                'Jcc(',k1,',', q1,',',k2,',', q2,')=',
     &                 Jcc(  k1,     q1,    k2,     q2  )
           End If

           If(ABS(Jcs(k1,q1,k2,q2)).gt. 0.5d-13) Then
           Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                'Jcs(',k1,',', q1,',',k2,',', q2,')=',
     &                 Jcs(  k1,     q1,    k2,     q2  )
           End If

           If(ABS(Jsc(k1,q1,k2,q2)).gt. 0.5d-13) Then
           Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                'Jsc(',k1,',', q1,',',k2,',', q2,')=',
     &                 Jsc(  k1,     q1,    k2,     q2  )
           End If

           If(ABS(Jss(k1,q1,k2,q2)).gt. 0.5d-13) Then
           Write(6,'(A,2(i2,A,i3,A),2ES18.10)')
     &                'Jss(',k1,',', q1,',',k2,',', q2,')=',
     &                 Jss(  k1,     q1,    k2,     q2  )
           End If
        ENDIF

      End Do
      End Do
      End Do
      End Do



      IF(DBG) THEN
      THRS=1.0E-16_wp
      Write(6,'(A,ES18.7)') 'Real Exchange parameters with values '//
     &                      'larger than: ', THRS
      Write(6,'(120A)') ('-',i=1,119),'|'
      Write(6,'(A)')   ' k1 |  q1 || k2 |  q2 |'//
     &                 '-------- O1-O2 --------|'//
     &                 '-------- O1-W2 --------|'//
     &                 '-------- W1-O2 --------|'//
     &                 '-------- W1-W2 --------|'
      Do k1=1,N1-1
      Do q1=-k1,k1
        Write(6,'(A)') '----|-----||----|-----|'//
     &                 '-----------------------|'//
     &                 '-----------------------|'//
     &                 '-----------------------|'//
     &                 '-----------------------|'
      Do k2=1,N2-1
      Do q2=-k2,k2

      F = knm(k1,ABS(q1))*knm(k2,ABS(q2))
      R = ABS( Jcc(k1,ABS(q1),k2,ABS(q2)) )
     &   +ABS( Jcs(k1,ABS(q1),k2,ABS(q2)) )
     &   +ABS( Jsc(k1,ABS(q1),k2,ABS(q2)) )
     &   +ABS( Jss(k1,ABS(q1),k2,ABS(q2)) )

      If( R > THRS ) Then

      If ( (q1 <0).and.(q2 <0) ) Then

            Write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),4(ES21.14,1x,A))')
     &                k1,'|',q1,'||',k2,'|',q2,'| ',
     &                DBLE(Jcc(k1,ABS(q1),k2,ABS(q2)))*F ,'| ',
     &                DBLE(Jcs(k1,ABS(q1),k2,ABS(q2)))*F ,'| ',
     &                DBLE(Jsc(k1,ABS(q1),k2,ABS(q2)))*F ,'| ',
     &                DBLE(Jss(k1,ABS(q1),k2,ABS(q2)))*F ,'|'

      Else If ( (q1>=0).and.(q2 <0) ) Then

            Write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))')
     &                k1,'|',q1,'||',k2,'|',q2,'| ',
     &                DBLE(Jcs(k1,ABS(q1),k2,ABS(q2)))*F ,'|'

      Else If ( (q1 <0).and.(q2>=0) ) Then

            Write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))')
     &                k1,'|',q1,'||',k2,'|',q2,'| ',
     &                DBLE(Jsc(k1,ABS(q1),k2,ABS(q2)))*F ,'|'

      Else ! (q1>0).and.(q2>0)

            Write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))')
     &                k1,'|',q1,'||',k2,'|',q2,'| ',
     &                DBLE(Jcc(k1,ABS(q1),k2,ABS(q2)))*F ,'|'

      End if
      End if

      End Do
      End Do
      End Do
      End Do
      Write(6,'(120A)') ('-',i=1,119),'|'

      ENDIF



!-----------------------------------------------------------------------
      Call mma_deallocate(O1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W1)
      Call mma_deallocate(W2)
      Call mma_deallocate(O1_O2)
      Call mma_deallocate(O1_W2)
      Call mma_deallocate(W1_O2)
      Call mma_deallocate(W1_W2)

      Return
      End Subroutine JKQPar
