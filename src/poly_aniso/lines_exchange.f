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
      Subroutine Lines_Exchange( Jex, N1, N2, S1, S2, HAM )
!     this Subroutine calculates the Lines exchange interaction between
!     two sites, of the one interacting pair
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      ! input variables
      Integer, intent(in)           :: N1, N2
      Real(kind=8), intent(in)     :: Jex
      Complex(kind=8), intent(in)  :: S1(3,N1,N1)
      Complex(kind=8), intent(in)  :: S2(3,N2,N2)
      ! output variables
      Complex(kind=8), intent(out) ::  HAM( N1,N1, N2,N2 )
      ! local variables
      Integer          :: i1,i2,j1,j2,l
      Call qEnter('Lines_Exchange')

      If( (N1<=0).OR.(N2<=0) ) Return
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HAM,1)
      If (Jex==0.0_wp) Return

      ! kind=8, complex double precision
      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2

              Do l=1,3
                HAM(i1,j1, i2,j2) = HAM(i1,j1, i2,j2)
     &                         + cmplx(-Jex,0_wp,wp)
     &                         * S1(l,i1,j1)
     &                         * S2(l,i2,j2)
              End Do

            End Do
          End Do
        End Do
      End Do
      Call qExit('Lines_Exchange')
      Return
      End Subroutine Lines_Exchange


      Subroutine Aniso_Lines_Exchange3( Jex, N1, N2, S1, S2, HAM )
!     this Subroutine calculates the Lines exchange interaction between
!     two sites, of the one interacting pair
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      ! input variables
      Integer, intent(in)           :: N1, N2
      Real(kind=8), intent(in)     :: Jex(3)
      Complex(kind=8), intent(in)  :: S1(3,N1,N1)
      Complex(kind=8), intent(in)  :: S2(3,N2,N2)
      ! output variables
      Complex(kind=8), intent(out) ::  HAM( N1,N1, N2,N2 )
      ! local variables
      Integer          :: i1,i2,j1,j2,l
      Complex(kind=8) :: Jc(3)
      Real(kind=8)    :: dnrm2_
      External         :: dnrm2_

      Call qEnter('Aniso_Lines3')
      If( (N1<=0).OR.(N2<=0) ) Return
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HAM,1)
      If ( dnrm2_(3,Jex,1)==0.0_wp) Return

      Call zcopy_(3,[(0.0_wp,0.0_wp)],0,Jc,1)
      Do l=1,3
        Jc(l) = cmplx(-Jex(l),0.0_wp,wp)
      End Do

      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2

              Do l=1,3
                HAM(i1,j1, i2,j2) = HAM(i1,j1, i2,j2)
     &                            + Jc(l) * S1(l,i1,j1) * S2(l,i2,j2)
              End Do

            End Do
          End Do
        End Do
      End Do
      Call qExit('Aniso_Lines3')
      Return
      End Subroutine Aniso_Lines_Exchange3



      Subroutine Aniso_Lines_Exchange9( Jex, N1, N2, S1, S2, HAM )
!     this Subroutine calculates the Lines exchange interaction between
!     two sites, of the one interacting pair
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      ! input variables
      Integer, intent(in)           :: N1, N2
      Real(kind=8), intent(in)     :: Jex(3,3)
      Complex(kind=8), intent(in)  :: S1(3,N1,N1)
      Complex(kind=8), intent(in)  :: S2(3,N2,N2)
      ! output variables
      Complex(kind=8), intent(out) ::  HAM( N1,N1, N2,N2 )
      ! local variables
      Integer          :: i1,i2,j1,j2,l,m
      Complex(kind=8) :: Jc(3,3)
      Real(kind=8)    :: dnrm2_
      External         :: dnrm2_

      Call qEnter('Aniso_Lines9')
      If( (N1<=0).OR.(N2<=0) ) Return
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HAM,1)
      If ( dnrm2_(9,Jex,1)==0.0_wp) Return

      Call zcopy_(3*3,[(0.0_wp,0.0_wp)],0,Jc,1)
      Do l=1,3
        Do m=1,3
           Jc(l,m)=cmplx(-Jex(l,m),0.0_wp,wp)
        End Do
      End Do

      ! kind=8, complex double precision
      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2

              Do l=1,3
                Do m=1,3
                  HAM(i1,j1, i2,j2) = HAM(i1,j1, i2,j2)
     &                    + Jc(l,m) * S1( l, i1,j1 ) * S2( m, i2,j2 )
                End Do
              End Do

            End Do
          End Do
        End Do
      End Do
      Call qExit('Aniso_Lines9')
      Return
      End Subroutine Aniso_Lines_Exchange9



      Subroutine Dzyaloshinsky_Morya_Exchange(Jex, N1, N2, S1, S2, HAM)
!     this Subroutine calculates the Dzyaloshinsky-Morya exchange interaction between
!     two sites, of the one interacting pair
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      ! input variables
      Integer, intent(in)           :: N1, N2
      Real(kind=8), intent(in)     :: Jex(3)
      Complex(kind=8), intent(in)  :: S1(3,N1,N1)
      Complex(kind=8), intent(in)  :: S2(3,N2,N2)
      ! output variables
      Complex(kind=8), intent(out) ::  HAM( N1,N1, N2,N2 )
      ! local variables
      Integer          :: i1,i2,j1,j2,l
      Complex(kind=8) :: X,Y,Z, Jc(3)
      Real(kind=8)    :: dnrm2_
      External         :: dnrm2_

      Call qEnter('DM_exchange')
      If( (N1<=0).OR.(N2<=0) ) Return
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HAM,1)
      If ( dnrm2_(3,Jex,1)==0.0_wp) Return

      Jc=(0.0_wp,0.0_wp)
      Do l=1,3
        Jc(l)=cmplx( -Jex(l),0.0_wp,wp )
      End Do

      ! kind=8, complex double precision
      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2

              X=(0.0_wp,0.0_wp)
              Y=(0.0_wp,0.0_wp)
              Z=(0.0_wp,0.0_wp)

              X =  S1(2,i1,j1) * S2(3,i2,j2)
     &            -S1(3,i1,j1) * S2(2,i2,j2)

              Y =  S1(3,i1,j1) * S2(1,i2,j2)
     &            -S1(1,i1,j1) * S2(3,i2,j2)

              Z =  S1(1,i1,j1) * S2(2,i2,j2)
     &            -S1(2,i1,j1) * S2(1,i2,j2)

              HAM(i1,j1, i2,j2) = HAM(i1,j1, i2,j2)
     &                           + Jc(1)*X
     &                           + Jc(2)*Y
     &                           + Jc(3)*Z

            End Do
          End Do
        End Do
      End Do
      Call qExit('DM_exchange')
      Return
      End Subroutine Dzyaloshinsky_Morya_Exchange



      Subroutine JITO_Exchange_Int( MxR1, MxR2, imaxrank,
     &                              n1, n2, JR, JI, HAM )

!     this Subroutine calculates the anisotropic exchange interaction between
!     two sites, of the one interacting pair on the basis of input ITO parameters
      Implicit None
#include "stdalloc.fh"
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      ! input variables
      Integer, intent(in)           :: imaxrank(2)
      Integer, intent(in)           :: MxR1, MxR2
      Integer, intent(in)           :: n1, n2
      Real(kind=8), intent(in)     :: JR(MxR1,-MxR1:MxR1,
     &                                    MxR2,-MxR2:MxR2)
      Real(kind=8), intent(in)     :: JI(MxR1,-MxR1:MxR1,
     &                                    MxR2,-MxR2:MxR2)
      ! output variables
      Complex(kind=8), intent(out) ::  HAM( n1,n1, n2,n2 )

      ! lcoal variables:
      Integer                       :: ibuf, k1,q1,k2,q2,m1,m2,l1,l2
      Real(kind=8)                 :: jpar, RR, RI
!      Real(kind=8)                 :: rK1,rK2,rQ1,rQ2,rM1,rM2,rJ1,rJ2,
!     &                                 CGp1,CGp2,CGm1,CGm2,CG01,CG02
      Real(kind=8)                :: C01,C02
      Complex(kind=8)             :: J(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2)
      Complex(kind=8), allocatable :: O1(:,:), O2(:,:)
      Complex(kind=8), allocatable :: W1(:,:), W2(:,:)
      Complex(kind=8), allocatable :: OO(:,:,:,:), WW(:,:,:,:),
     &                                 OW(:,:,:,:), WO(:,:,:,:)
      Logical                       :: dbg
      Real(kind=8)                 :: dnrm2_
      External                      :: dnrm2_

      Call qEnter('JITO_Exchange_Int')
      dbg=.false.
! ----  initial checks
      If( (n1<=0).OR.(n2<=0) ) Return
      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,HAM,1)
      ibuf=0
      ibuf=MxR1*(2*MxR1+1)*MxR2*(2*MxR2+1)
      If(ibuf==0) Return
      jpar=0.0_wp
      jpar=dnrm2_(ibuf,JR(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)
     &    +dnrm2_(ibuf,JI(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)
      If (jpar==0.0_wp) Return
! ---- end initial checks
      Call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,
     &                 J(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)
      Do k1=1,MxR1,2
       Do k2=1,MxR2,2
        Do q1=-k1,k1
         Do q2=-k2,k2
          RR=JR(k1,q1,k2,q2)
          RI=JI(k1,q1,k2,q2)
          J(k1,q1,k2,q2)=cmplx(RR,RI,wp)
         End Do
        End Do
       End Do
      End Do
!----------------------------------------------------------------
      Call mma_allocate(O1,n1,n1,'operator O1')
      Call mma_allocate(W1,n1,n1,'operator W1')
      Call mma_allocate(O2,n2,n2,'operator O2')
      Call mma_allocate(W2,n2,n2,'operator W2')
      Call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
      Call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
      Call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
      Call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
!----------------------------------------------------------------
      Do k1=1,imaxrank(1),2
       Do q1=0,k1
        Do k2=1,imaxrank(2),2
         Do q2=0,k2
          ! generate the operator matrix K=ik, Q=iq, dimension=na
          Call ITO(n1,k1,q1,C01,O1,W1)
          Call ITO(n2,k2,q2,C02,O2,W2)
          !generate coupled operators:
          Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,OO,1)
          Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,OW,1)
          Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,WO,1)
          Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,WW,1)
          Do m1=1,n1
           Do m2=1,n1
            Do l1=1,n2
             Do l2=1,n2
               OO(m1,m2,l1,l2) = O1(m1,m2) * O2(l1,l2)
               OW(m1,m2,l1,l2) = O1(m1,m2) * W2(l1,l2)
               WO(m1,m2,l1,l2) = W1(m1,m2) * O2(l1,l2)
               WW(m1,m2,l1,l2) = W1(m1,m2) * W2(l1,l2)
             End Do
            End Do
           End Do
          End Do !m1
          ! compute the exchange Hamiltonian:
          If((q1==0).and.(q2==0)) Then
           Call zaxpy_(n1*n1*n2*n2,J(k1,  0,k2,  0),OO,1,HAM,1)
          Else If ((q1==0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,J(k1,  0,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,J(k1,  0,k2,-q2),OW,1,HAM,1)
          Else If ((q1.ne.0).and.(q2==0)) Then
           Call zaxpy_(n1*n1*n2*n2,J(k1, q1,k2,  0),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2,  0),WO,1,HAM,1)
          Else If ((q1.ne.0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,J(k1, q1,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,J(k1, q1,k2,-q2),OW,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2, q2),WO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2,-q2),WW,1,HAM,1)
          End If
         End Do
        End Do
       End Do
      End Do ! k1

      If(dbg) Then
       Write(6,'(A)') 'JITO_Exchange_Int: generated <m1,l1|HAM|m2,l2>'
       Do m1=1,n1
        Do m2=1,n1
         Do l1=1,n2
          Do l2=1,n2
           Write(6,'(A,i2,A,i2,A,i2,A,i2,A,2ES22.14)')
     &        '<',m1,',',l1,'| HAM |',m2,',',l2,'> = ',HAM(m1,m2,l1,l2)
          End Do
         End Do
        End Do
       End Do
      End If

      Call mma_deallocate(OO)
      Call mma_deallocate(OW)
      Call mma_deallocate(WO)
      Call mma_deallocate(WW)
      Call mma_deallocate(O1)
      Call mma_deallocate(W1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W2)

      Call qExit('JITO_Exchange_Int')
      Return
      End Subroutine JITO_Exchange_Int
