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
      Subroutine transHam( n1, n2, rot1, rot2, MM1, MM2, typ1, typ2,
     &                     H, HT, iopt)
c purpose:  transform the exchange Hamiltonian
c           matrices to the basis of their local pseudospins
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "stdalloc.fh"
!     exchange basis of both sites
      Integer, intent(in)           :: n1, n2, iopt
      Real(kind=8), intent(in)     :: rot1(3,3)
      Real(kind=8), intent(in)     :: rot2(3,3)
      Complex(kind=8), intent(in)  :: MM1(3,n1,n1)
      Complex(kind=8), intent(in)  :: MM2(3,n2,n2)
      Complex(kind=8), intent(in)  :: H( n1,n1,n2,n2)
      Complex(kind=8), intent(out) :: HT(n1,n1,n2,n2)
      Character(Len=1), intent(in)  :: typ1, typ2
      ! local variables
      Integer                       :: i1,i2,j1,j2,i
      Real(kind=8), allocatable    :: gt1(:), gt2(:)
      Real(kind=8), allocatable    :: ax1(:,:), ax2(:,:)
      Complex(kind=8), allocatable :: M1(:,:,:),   M2(:,:,:)
      Complex(kind=8), allocatable :: MR1(:,:,:), MR2(:,:,:)
      Complex(kind=8), allocatable ::   Z1(:,:), Z2(:,:)
      Complex(kind=8), allocatable :: TMP1(:,:), TMP2(:,:)
      Complex(kind=8), allocatable :: HI(:,:,:,:) !HI(n1,n1,n2,n2)
      Logical ::  DBG

      DBG= .false.

!-----------------------------------------------------------------------
      Call mma_allocate(gt1,3,'gt1')
      Call mma_allocate(gt2,3,'gt2')
      Call mma_allocate(ax1,3,3,'ax1')
      Call mma_allocate(ax2,3,3,'ax2')
      If(n1>0) Then
         Call mma_allocate(M1,3,n1,n1,'M1')
         Call mma_allocate(MR1,3,n1,n1,'MR1')
         Call mma_allocate(Z1,n1,n1,'Z1')
         Call mma_allocate(TMP1,n1,n1,'TMP1')
      End If
      If(n2>0) Then
         Call mma_allocate(M2,3,n2,n2,'M2')
         Call mma_allocate(MR2,3,n2,n2,'MR2')
         Call mma_allocate(Z2,n2,n2,'Z2')
         Call mma_allocate(TMP2,n2,n2,'TMP2')
      End If
      If((n1>0).AND.(n2>0)) Then
         Call mma_allocate(HI,n1,n1,n2,n2,'HI')
      End If
      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,HT,1)
!-----------------------------------------------------------------------


      ! rotate the magnetic moments to their general coordinate system
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,M1,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,M2,1)
      Call rotmom( MM1, n1, rot1, M1 )
      Call rotmom( MM2, n2, rot2, M2 )

      If (iopt.eq.1) Then ! local coordinate system
         Call dcopy_(3,[0.0_wp],0,gt1,1)
         Call dcopy_(3,[0.0_wp],0,gt2,1)
         Call dcopy_(3*3,[0.0_wp],0,ax1,1)
         Call dcopy_(3*3,[0.0_wp],0,ax2,1)
         If(DBG) Then
            Write(6,'(A)') 'TRANSHAM:: local g tensors and axes:'
            Call atens( M1, n1, gt1, ax1, 2)
            Call atens( M2, n2, gt2, ax2, 2)
         Else
            Call atens( M1, n1, gt1, ax1, 1)
            Call atens( M2, n2, gt2, ax2, 1)
         End If

      Else If (iopt==2) Then ! general coordinate system
         Call dcopy_(3*3,[0.0_wp],0,ax1,1)
         Call dcopy_(3*3,[0.0_wp],0,ax2,1)
         Do i1=1,3
           ax1(i1,i1)=1.0_wp
           ax2(i1,i1)=1.0_wp
         End Do
      End If



c----------------------------------------------------------------------
      ! rotate magnetic moments to their local magnetic axes:
      Call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,MR1,1)
      Call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,MR2,1)
      If (  ((typ1.eq.'A').AND.(typ2.eq.'A')) .OR.
     &      ((typ1.eq.'B').AND.(typ2.eq.'B')) .OR.
     &      ((typ1.eq.'B').AND.(typ2.eq.'C')) .OR.
     &      ((typ1.eq.'C').AND.(typ2.eq.'B')) .OR.
     &      ((typ1.eq.'C').AND.(typ2.eq.'C')) )  Then
         ! both sites have been computed ab initio
         Call rotmom2( M1, n1, ax1, MR1 )
         Call rotmom2( M2, n2, ax2, MR2 )
      Else If( ((typ1.eq.'A').AND.(typ2.eq.'B')).OR.
     &         ((typ1.eq.'A').AND.(typ2.eq.'C')) ) Then
         !  site 2 is generated isotropic
         !  use axes of site 1
         Call rotmom2( M1, n1, ax1, MR1 )
         Call rotmom2( M2, n2, ax1, MR2 )
      Else If( ((typ1.eq.'B').AND.(typ2.eq.'A')).OR.
     &         ((typ1.eq.'C').AND.(typ2.eq.'A')) ) Then
         !  site 1 is generated isotropic
         !  use axes of site 2
         Call rotmom2( M1, n1, ax2, MR1 )
         Call rotmom2( M2, n2, ax2, MR2 )
      End If

      If(DBG) Then
        If(iopt.eq.1) Then
          Write(6,'(A,i3)') 'site 1'
          Call prMom('transHam:: magnetic moment, coordinate '//
     &               'system of LOCAL main magnetic axes, site 1',
     &                MR1, n1)
          Write(6,'(A,i3)') 'site 2'
          Call prMom('transHam:: magnetic moment, coordinate '//
     &               'system of LOCAL main magnetic axes, site 2',
     &                MR2, n2)
        Else
          Write(6,'(A,i3)') 'site 1'
          Call prMom('transHam:: magnetic moment, coordinate '//
     &               'system of GENERAL main magnetic axes, site 1',
     &                MR1, n1)
          Write(6,'(A,i3)') 'site 2'
          Call prMom('transHam:: magnetic moment, coordinate '//
     &               'system of GENERAL main magnetic axes, site 2',
     &                MR2, n2)
        End If
      End If



c----------------------------------------------------------------------
      ! find local pseudospin on each site:
      Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,Z1,1)
      Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,Z2,1)
      Call pseudospin( MR1, n1, Z1, 3,1, 1 )
      Call pseudospin( MR2, n2, Z2, 3,1, 1 )


      If (DBG) Then
         Call pa_prmat('Matrix Z1',Z1,n1)
         Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,TMP1,1)
         Call zgemm_('C','N', n1, n1, n1, (1.0_wp,0.0_wp),
     &               Z1, n1,
     &               Z1, n1,(0.0_wp,0.0_wp),
     &               TMP1, n1 )
         Write(6,'(A)') 'transHam:  Verify unitarity of Z1'
         Do i=1,n1
            Write(6,'(A,i2,A,i2,A,2ES20.10)')
     &               'conjg(Z1)*Z1:  (',i,',',i,')=', TMP1(i,i)
         End Do
         Call pa_prmat('Matrix Z2',Z2,n2)
         Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,TMP2,1)
         Call zgemm_('C','N', n2, n2, n2, (1.0_wp,0.0_wp),
     &               Z2, n2,
     &               Z2, n2,(0.0_wp,0.0_wp),
     &               TMP2, n2 )
         Write(6,'(A)') 'transHam:  Verify unitarity of Z2'
         Do i=1,n2
            Write(6,'(A,i2,A,i2,A,2ES20.10)')
     &               'conjg(Z2)*Z2:  (',i,',',i,')=', TMP2(i,i)
         End Do
         Call UTMU2( n1, n1, Z1, MR1 )
         Call UTMU2( n2, n2, Z2, MR2 )
         Call prMom('transHam:: magnetic moment, coordinate '//
     &              'MR1, site 1',
     &               MR1, n1)
         Call prMom('transHam:: magnetic moment, coordinate '//
     &              'MR2, site 2',
     &               MR2, n2)
      End If


      ! save a local copy:
      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,HI,1)
      Call zcopy_(n1*n1*n2*n2,H,1,HI,1)

      Do i1=1,n1
        Do j1=1,n1
          Do i2=1,n2
            Do j2=1,n2
!              If (DBG.AND.(ABS(H(i1,j1,i2,j2)).gt.0.5d-13)) Then
                 Write(6,'(A,4(i2,A),2ES22.14)')
     &              'H (',i1,',',j1,',',i2,',',j2,')=',H(i1,j1,i2,j2)
!              End If
            End Do
          End Do
        End Do
      End Do

      ! transform the Hamiltonian to local pseudospins:
      Do i2=1,n2
        Do j2=1,n2
          Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,TMP1,1)
          Call zgemm_('C','N', n1, n1, n1, (1.0_wp,0.0_wp),
     &                Z1, n1,
     &                HI(1:n1,1:n1, i2,j2), n1,(0.0_wp,0.0_wp),
     &                TMP1, n1 )


          Call zcopy_(n1*n1,[(0.0_wp,0.0_wp)],0,HI(1:n1,1:n1, i2,j2),1)
          Call zgemm_('N','N', n1, n1, n1, (1.0_wp,0.0_wp),
     &                TMP1, n1,
     &                Z1, n1, (0.0_wp,0.0_wp),
     &                HI(1:n1,1:n1, i2,j2), n1 )
        End Do
      End Do

      Do i1=1,n1
        Do j1=1,n1
          Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,TMP2,1)
          Call zgemm_( 'C','N', n2, n2, n2, (1.0_wp,0.0_wp),
     &                Z2, n2,
     &                HI(i1,j1, 1:n2,1:n2), n2, (0.0_wp,0.0_wp),
     &                TMP2, n2 )

          Call zcopy_(n2*n2,[(0.0_wp,0.0_wp)],0,HI(i1,j1,1:n2,1:n2),1)
          Call zgemm_('N','N', n2, n2, n2, (1.0_wp,0.0_wp),
     &                TMP2, n2,
     &                Z2, n2, (0.0_wp,0.0_wp),
     &                HI(i1,j1,1:n2,1:n2), n2 )
        End Do
      End Do

      Call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,HT,1)
      Call zcopy_(n1*n1*n2*n2,HI,1,HT,1)


      Do i1=1,n1
        Do j1=1,n1
          Do i2=1,n2
            Do j2=1,n2
!              If (DBG.AND.(ABS(HT(i1,j1,i2,j2)).gt.0.5d-13)) Then
                 Write(6,'(A,4(i2,A),2ES22.14)')
     &              'HT(',i1,',',j1,',',i2,',',j2,')=',HT(i1,j1,i2,j2)
!              End If
            End Do
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
      Call mma_deallocate(gt1)
      Call mma_deallocate(gt2)
      Call mma_deallocate(ax1)
      Call mma_deallocate(ax2)
      If(n1>0) Then
         Call mma_deallocate(M1)
         Call mma_deallocate(MR1)
         Call mma_deallocate(Z1)
         Call mma_deallocate(TMP1)
      End If
      If(n2>0) Then
         Call mma_deallocate(M2)
         Call mma_deallocate(MR2)
         Call mma_deallocate(Z2)
         Call mma_deallocate(TMP2)
      End If
      If((n1>0).AND.(n2>0)) Then
         Call mma_deallocate(HI)
      End If

      Return
      End Subroutine transHam
