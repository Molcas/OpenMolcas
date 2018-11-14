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
      Subroutine recover_exch_HAM_from_Naoya_ITO(n1,n2,J,HAM)
      Implicit none
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: n1, n2
      Complex(kind=wp), intent(in)  ::
     &     J( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      Complex(kind=wp),intent(out)  :: HAM(n1,n1,n2,n2)
      ! local variables:
      Integer                       :: k1,k2,q1,q2,m1,m2,l1,l2
      Real(kind=wp)                 :: C01,C02
      Complex(kind=wp), allocatable :: O1(:,:), O2(:,:),
     &                                 W1(:,:), W2(:,:)
      Complex(kind=wp), allocatable :: OO(:,:,:,:), WW(:,:,:,:),
     &                                 OW(:,:,:,:), WO(:,:,:,:)
!---------------------------------------------------------------------
!  recover the original HAMILTONIAN using the J parameters
!==================================================================
      Call mma_allocate(O1,n1,n1,'operator O1')
      Call mma_allocate(O2,n2,n2,'operator O2')
      Call mma_allocate(W1,n1,n1,'operator W1')
      Call mma_allocate(W2,n2,n2,'operator W2')
      Call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
      Call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
      Call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
      Call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
      Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,HAM,1)
      Do k1=1,n1-1
       Do q1=0,k1
        Do k2=1,n2-1
         Do q2=0,k2
          ! generate the operator matrix K=ik, Q=iq, dimension=na
          Call ITO(n1,k1,q1,C01,O1,W1)
          Call ITO(n2,k2,q2,C02,O2,W2)
          !generate coupled operators:
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OW,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WW,1)
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
         End Do !q
        End Do !k
       End Do !q
      End Do !k

      Call mma_deallocate(O1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W1)
      Call mma_deallocate(W2)
      Call mma_deallocate(OO)
      Call mma_deallocate(OW)
      Call mma_deallocate(WO)
      Call mma_deallocate(WW)

      Return
      End Subroutine recover_exch_HAM_from_Naoya_ITO





      Subroutine recover_exch_HAM_from_Liviu_ITO(n1,n2,B,HAM)
      Implicit none
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: n1, n2
      Complex(kind=wp), intent(in)  ::
     &     B( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      Complex(kind=wp),intent(out)  :: HAM(n1,n1,n2,n2)
      ! local variables:
      Integer                       :: k1,k2,q1,q2,m1,m2,l1,l2
      Complex(kind=wp)              :: redME1,redME2
      Complex(kind=wp), allocatable :: O1(:,:), O2(:,:),
     &                                 W1(:,:), W2(:,:)
      Complex(kind=wp), allocatable :: OO(:,:,:,:), WW(:,:,:,:),
     &                                 OW(:,:,:,:), WO(:,:,:,:)

!---------------------------------------------------------------------
!  recover the original HAMILTONIAN using the J parameters
!==================================================================
      Call mma_allocate(O1,n1,n1,'operator O1')
      Call mma_allocate(O2,n2,n2,'operator O2')
      Call mma_allocate(W1,n1,n1,'operator W1')
      Call mma_allocate(W2,n2,n2,'operator W2')
      Call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
      Call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
      Call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
      Call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
      Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,HAM,1)
      Do k1=1,n1-1
       Do q1=0,k1
        Do k2=1,n2-1
         Do q2=0,k2
          ! generate the operator matrix K=ik, Q=iq, dimension=na
          Call Liviu_ITO(n1,k1,q1,O1,W1,redME1)
          Call Liviu_ITO(n2,k2,q2,O2,W2,redME2)
          !generate coupled operators:
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OW,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WW,1)
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
           Call zaxpy_(n1*n1*n2*n2,B(k1,  0,k2,  0),OO,1,HAM,1)
          Else If ((q1==0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,B(k1,  0,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,B(k1,  0,k2,-q2),OW,1,HAM,1)
          Else If ((q1.ne.0).and.(q2==0)) Then
           Call zaxpy_(n1*n1*n2*n2,B(k1, q1,k2,  0),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,B(k1,-q1,k2,  0),WO,1,HAM,1)
          Else If ((q1.ne.0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,B(k1, q1,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,B(k1, q1,k2,-q2),OW,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,B(k1,-q1,k2, q2),WO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,B(k1,-q1,k2,-q2),WW,1,HAM,1)
          End If
         End Do !q
        End Do !k
       End Do !q
      End Do !k

      Call mma_deallocate(O1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W1)
      Call mma_deallocate(W2)
      Call mma_deallocate(OO)
      Call mma_deallocate(OW)
      Call mma_deallocate(WO)
      Call mma_deallocate(WW)

      Return
      End Subroutine recover_exch_HAM_from_Liviu_ITO






      Subroutine recover_exch_HAM_from_Stevens_ESO(n1,n2,S,HAM)
      Implicit none
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: n1, n2
      Complex(kind=wp), intent(in)  ::
     &     S( (n1-1), -(n1-1):(n1-1), (n2-1), -(n2-1):(n2-1) )
      Complex(kind=wp),intent(out)  :: HAM(n1,n1,n2,n2)
      ! local variables:
      Integer                       :: k1,k2,q1,q2,m1,m2,l1,l2
      Complex(kind=wp)              :: redME1,redME2
      Complex(kind=wp), allocatable :: O1(:,:), O2(:,:),
     &                                 W1(:,:), W2(:,:)
      Complex(kind=wp), allocatable :: OO(:,:,:,:), WW(:,:,:,:),
     &                                 OW(:,:,:,:), WO(:,:,:,:)

!---------------------------------------------------------------------
!  recover the original HAMILTONIAN using the S parameters
!==================================================================
      Call mma_allocate(O1,n1,n1,'operator O1')
      Call mma_allocate(O2,n2,n2,'operator O2')
      Call mma_allocate(W1,n1,n1,'operator W1')
      Call mma_allocate(W2,n2,n2,'operator W2')
      Call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
      Call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
      Call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
      Call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
      Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,HAM,1)
      Do k1=1,n1-1
       Do q1=0,k1
        Do k2=1,n2-1
         Do q2=0,k2
          ! generate the operator matrix K=ik, Q=iq, dimension=na
          Call ESO(n1,k1,q1,O1,W1,redME1)
          Call ESO(n2,k2,q2,O2,W2,redME2)
          !generate coupled operators:
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,OW,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WO,1)
          Call zcopy_(n1*n1*n2*n2,(0.0_wp,0.0_wp),0,WW,1)
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
           Call zaxpy_(n1*n1*n2*n2,S(k1,  0,k2,  0),OO,1,HAM,1)
          Else If ((q1==0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,S(k1,  0,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,S(k1,  0,k2,-q2),OW,1,HAM,1)
          Else If ((q1.ne.0).and.(q2==0)) Then
           Call zaxpy_(n1*n1*n2*n2,S(k1, q1,k2,  0),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,S(k1,-q1,k2,  0),WO,1,HAM,1)
          Else If ((q1.ne.0).and.(q2.ne.0)) Then
           Call zaxpy_(n1*n1*n2*n2,S(k1, q1,k2, q2),OO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,S(k1, q1,k2,-q2),OW,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,S(k1,-q1,k2, q2),WO,1,HAM,1)
           Call zaxpy_(n1*n1*n2*n2,S(k1,-q1,k2,-q2),WW,1,HAM,1)
          End If
         End Do !q
        End Do !k
       End Do !q
      End Do !k

      Call mma_deallocate(O1)
      Call mma_deallocate(O2)
      Call mma_deallocate(W1)
      Call mma_deallocate(W2)
      Call mma_deallocate(OO)
      Call mma_deallocate(OW)
      Call mma_deallocate(WO)
      Call mma_deallocate(WW)

      Return
      End Subroutine recover_exch_HAM_from_Stevens_ESO

