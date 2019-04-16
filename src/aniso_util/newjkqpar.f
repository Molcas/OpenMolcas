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
      Subroutine newjkqpar(n1,n2,H,J,B,S)

      Implicit none
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: n1, n2
      Complex(kind=wp),intent(in)   :: H(n1,n1,n2,n2)
      Complex(kind=wp), intent(out) ::
     &     J( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      Complex(kind=wp), intent(out) ::
     &     B( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      Complex(kind=wp), intent(out) ::
     &     S( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      ! local variables:
      Integer                       :: k1,k2,q1,q2,m1,m2,m12,i1,i2,j1,j2
      Real(kind=wp)                 :: cr1,cr2,C01,C02,r1,r2,F1,F2
      Complex(kind=wp)              :: cf1,cf2,c1,c2,c12,ci,cc1,cc2,
     &                                 trace_exch2
      Complex(kind=wp), allocatable :: O1(:,:), W1(:,:)
      Complex(kind=wp), allocatable :: O2(:,:), W2(:,:)
      Complex(kind=wp), allocatable :: OO(:,:,:,:), WO(:,:,:,:)
      Complex(kind=wp), allocatable :: OW(:,:,:,:), WW(:,:,:,:)
      Complex(kind=wp), allocatable :: HAM(:,:,:,:)
      Real(kind=wp)                 :: knm(12,0:12),dznrm2_
      External                      :: dznrm2_,trace_exch2
      Logical                       :: dbg

      Call qEnter('newJKQP')
!-------------------------------------------
      If( (n1<1).or.(n2<1)) Return
!-------------------------------------------
      dbg=.false.
      Call mma_allocate(O1,n1,n1,'operator O1')
      Call mma_allocate(W1,n1,n1,'operator W1')
      Call mma_allocate(O2,n2,n2,'operator O2')
      Call mma_allocate(W2,n2,n2,'operator W2')
      Call set_knm( knm )
!-------------------------------------------
      J(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1))=(0.0_wp,0.0_wp)
      B(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1))=(0.0_wp,0.0_wp)
      S(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1))=(0.0_wp,0.0_wp)
      Do k1=1,n1-1
       Do q1=0,k1
        Do k2=1,n2-1
         Do q2=0,k2
          cr1=0.0_wp
          cr2=0.0_wp
          m1=0
          m2=0
          r1=0.0_wp
          r2=0.0_wp
          cf1=(0.0_wp,0.0_wp)
          cf2=(0.0_wp,0.0_wp)

          m1 =(-1)**q1
          m2 =(-1)**q2
          m12=(-1)**(q1+q2)
          c1 =cmplx(m1,0.0_wp,wp)
          c2 =cmplx(m2,0.0_wp,wp)
          c12=cmplx(m12,0.0_wp,wp)
          ci =(0.0_wp,-1.0_wp)

          ! generate the operators for site 1 and 2:
          Call ITO(n1,k1,q1,C01,O1,W1)
          Call ITO(n2,k2,q2,C02,O2,W2)

          ! compute the scaling factors for various Operators
          Call coeff_redus_sub(n1,k1,cr1)
          Call coeff_redus_sub(n2,k2,cr2)
          cc1=cmplx(1.0_wp/(cr1*C01), 0.0_wp,wp)
          cc2=cmplx(1.0_wp/(cr2*C02), 0.0_wp,wp)
          r1=C01*C01*dble(2*k1+1)/dble(n1)
          r2=C02*C02*dble(2*k2+1)/dble(n2)
          cf1=c1*cmplx(r1,0.0_wp,wp)
          cf2=c2*cmplx(r2,0.0_wp,wp)

          If (dbg) Then
             ! use the old trace_exch function
             Call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
             Call mma_allocate(WO,n1,n1,n2,n2,'operator WO')
             Call mma_allocate(OW,n1,n1,n2,n2,'operator OW')
             Call mma_allocate(WW,n1,n1,n2,n2,'operator WW')

             Do i1=1,n1
              Do j1=1,n1
               Do i2=1,n2
                Do j2=1,n2
                 OO(i1,j1,i2,j2)=O1(i1,j1)*O2(i2,j2)
                 OW(i1,j1,i2,j2)=O1(i1,j1)*W2(i2,j2)
                 WO(i1,j1,i2,j2)=W1(i1,j1)*O2(i2,j2)
                 WW(i1,j1,i2,j2)=W1(i1,j1)*W2(i2,j2)
                End Do
               End Do
              End Do
             End Do
             ! find the parameters:
             ! in the Naoya's operators
             J(k1,-q1, k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,O2)
             J(k1,-q1, k2, q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,W2)
             J(k1, q1, k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,O2)
             J(k1, q1, k2, q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,W2)

             Call mma_deallocate(OO)
             Call mma_deallocate(WO)
             Call mma_deallocate(OW)
             Call mma_deallocate(WW)
          Else
             ! find the parameters:
             ! in the Naoya's operators
             J(k1,-q1, k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,O2)
             J(k1,-q1, k2, q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,W2)
             J(k1, q1, k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,O2)
             J(k1, q1, k2, q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,W2)
          End If

          ! in the Liviu operators
          B(k1,-q1, k2,-q2) = J(k1,-q1, k2,-q2)*cc1*cc2
          B(k1,-q1, k2, q2) = J(k1,-q1, k2, q2)*cc1*cc2
          B(k1, q1, k2,-q2) = J(k1, q1, k2,-q2)*cc1*cc2
          B(k1, q1, k2, q2) = J(k1, q1, k2, q2)*cc1*cc2

!         generate real parameters for Extended Stevens Operators formalism:
          If((q1>0).and.(q2>0)) Then
! BB = (Jmm + (-1)^q2 Jmp + (-1)^q1 Jpm + (-1)^(q1 + q2) Jpp);
               S(k1, q1, k2, q2) =        B( k1,-q1, k2,-q2)
     &                             + c2 * B( k1,-q1, k2, q2)
     &                             + c1 * B( k1, q1, k2,-q2)
     &                             +c12 * B( k1, q1, k2, q2)
! BC = (Jmm - (-1)^q2 Jmp + (-1)^q1 Jpm - (-1)^(q1 + q2) Jpp) (-I);
               S(k1, q1, k2,-q2) = (      B( k1,-q1, k2,-q2)
     &                             - c2 * B( k1,-q1, k2, q2)
     &                             + c1 * B( k1, q1, k2,-q2)
     &                             -c12 * B( k1, q1, k2, q2) )*ci
! CB = (Jpp + (-1)^q2 Jpm - (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm) (-I);
               S(k1,-q1, k2, q2) = (      B( k1,-q1, k2,-q2)
     &                             + c2 * B( k1,-q1, k2, q2)
     &                             - c1 * B( k1, q1, k2,-q2)
     &                             -c12 * B( k1, q1, k2, q2) )*ci
! CC = (-Jpp + (-1)^q2 Jpm + (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm);
               S(k1,-q1, k2,-q2) =       -B( k1,-q1, k2,-q2)
     &                             + c2 * B( k1,-q1, k2, q2)
     &                             + c1 * B( k1, q1, k2,-q2)
     &                             -c12 * B( k1, q1, k2, q2)
          Else If ((q1==0).and.(q2>0))  Then
               S(k1, q1, k2, q2) = (      B( k1, q1, k2,-q2)
     &                             + c2 * B( k1, q1, k2, q2) )
               S(k1, q1, k2,-q2) = (      B( k1, q1, k2,-q2)
     &                             - c2 * B( k1, q1, k2, q2) )*ci
          Else If ((q1>0).and.(q2==0))  Then
               S(k1, q1, k2, q2) = (      B( k1,-q1, k2, q2)
     &                             + c1 * B( k1, q1, k2, q2) )
               S(k1,-q1, k2, q2) = (      B( k1,-q1, k2, q2)
     &                             - c1 * B( k1, q1, k2, q2) )*ci
          Else If ((q1==0).and.(q2==0))  Then
               S(k1, q1, k2, q2) =        B( k1, q1, k2, q2)
          End If

         End Do
        End Do
       End Do
      End Do

      ! scale the Stevens parameters to match the ESO operators:
      Do k1=1,n1-1
       Do q1=-k1,k1
        Do k2=1,n2-1
         Do q2=-k2,k2
           F1=knm(k1,abs(q1))
           F2=knm(k2,abs(q2))
           S(k1,q1,k2,q2) = S(k1,q1,k2,q2)*cmplx(F1*F2,0.0_wp,wp)
         End Do
        End Do
       End Do
      End Do

      ! in case verification is needed:
      If(dbg) Then
       Call mma_allocate(HAM,n1,n1,n2,n2,'recovered HAM_S')
       Write(6,'(A)') 'Extracted exchange parameters: J(k1,q1,k2,q2):'
       Write(6,'(A)') 'using Naoya ITO:'
       Do k1=1,n1-1
        Do q1=-k1,k1
         Do k2=1,n2-1
          Do q2=-k2,k2
!           If( ABS(S(k1,q1,k2,q2)).gt.1.0e-20_wp) Then
           Write(6,'(4(A,i2),A,2ES22.13)')
     &       'J(',k1,',',q1,',',k2,',',q2,') = ',J(k1,q1,k2,q2)
!           End If
          End Do
         End Do
        End Do
       End Do
       Call recover_exch_HAM_from_Naoya_ITO(n1,n2,J,HAM)
       Write(6,'(A,ES20.10)') 'recover from Naoya Jkqkq parameters'
       Write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',
     &                         dznrm2_(n1*n1*n2*n2,HAM-H,1)

       Write(6,'(A)') 'Extracted exchange parameters: B(k1,q1,k2,q2):'
       Write(6,'(A)') 'using Liviu ITO:'
       Do k1=1,n1-1
        Do q1=-k1,k1
         Do k2=1,n2-1
          Do q2=-k2,k2
!           If( ABS(S(k1,q1,k2,q2)).gt.1.0e-20_wp) Then
           Write(6,'(4(A,i2),A,2ES22.13)')
     &       'B(',k1,',',q1,',',k2,',',q2,') = ',B(k1,q1,k2,q2)
!           End If
          End Do
         End Do
        End Do
       End Do
       Call recover_exch_HAM_from_Liviu_ITO(n1,n2,B,HAM)
       Write(6,'(A,ES20.10)') 'recover from Liviu Bkqkq parameters'
       Write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',
     &                         dznrm2_(n1*n1*n2*n2,HAM-H,1)

       Write(6,'(A)') 'Extracted exchange parameters: S(k1,q1,k2,q2):'
       Write(6,'(A)') 'using Stevens ESO:'
       Do k1=1,n1-1
        Do q1=-k1,k1
         Do k2=1,n2-1
          Do q2=-k2,k2
!           If( ABS(S(k1,q1,k2,q2)).gt.1.0e-20_wp) Then
           Write(6,'(4(A,i2),A,2ES22.13)')
     &       'S(',k1,',',q1,',',k2,',',q2,') = ',S(k1,q1,k2,q2)
!           End If
          End Do
         End Do
        End Do
       End Do
       Call recover_exch_HAM_from_Stevens_ESO(n1,n2,S,HAM)
       Write(6,'(A,ES20.10)') 'recover from Stevens Skqkq parameters'
       Write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',
     &                         dznrm2_(n1*n1*n2*n2,HAM-H,1)
       Call mma_deallocate(HAM)
      End If

!==================================================================
      Call mma_deallocate(O1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W1)
      Call mma_deallocate(W2)

      Call qExit('newJKQP')

      Return
      End subroutine newjkqpar
