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
      Subroutine SPIN_PHASE(MM,dim,Zinp,Zout)
C
C     The RASSI program gives a random phase to the spin-orbit functions.
C
C     This routine performs a simple check with the obtained spin functions,
C     in order to determine the phase of the spin functions.
C     If the phase is not the same, Then the spin functions will be multiplied
C     with the correspondind coefficient that sets the same phase to all spin
C     eigenfunctions
C
      Implicit None
#include "stdalloc.fh"
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: dim
      Complex(kind=wp), intent(in)  :: mm(3,dim,dim)
      Complex(kind=wp), intent(in)  :: Zinp(dim,dim)
      Complex(kind=wp), intent(out) :: Zout(dim,dim)
! ------------------------------------------------------------
      Integer                       :: i, j, i1, i2, l
      Real(kind=wp), allocatable    :: rxr(:) !dim)
      Real(kind=wp), allocatable    :: rxi(:) !dim)
      Complex(kind=wp), allocatable :: r(:) !(dim)
      Complex(kind=wp), allocatable :: phs(:,:,:)  !3,dim,dim)
      Complex(kind=wp), allocatable :: tmp(:,:) !dim,dim
      Logical :: dbg

      Call qEnter('s_phase')

      dbg=.false.

      Call mma_allocate(rxr,dim,'rxr')
      Call mma_allocate(rxi,dim,'rxi')
      Call mma_allocate(r,dim,'r')
      Call mma_allocate(phs,3,dim,dim,'phs')
      Call mma_allocate(tmp,dim,dim,'tmp')
! ------------------------------------------------------------
      Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
      Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
      Call zcopy_(dim,(0.0_wp,0.0_wp),0,r,1)
      Call dcopy_(dim,0.0_wp,0,rxr,1)
      Call dcopy_(dim,0.0_wp,0,rxi,1)
      rxr(1)=1.0_wp
      rxi(1)=0.0_wp

      Do i=1,dim-1
        j = i+1
        r(j)=(0.0_wp,0.0_wp)
        Do i1=1,dim
          Zout(i1,1)=Zinp(i1,1)
        End Do

        Do i1=1,dim
          Do i2=1,dim
            phs(1,i,j)=phs(1,i,j)+ MM(1,i1,i2)*CONJG(Zout(i1,i))*
     &                                               Zinp(i2,j)
          End Do
        End Do

        If(ABS(phs(1,i,j)).gt.1.0e-14_wp ) Then
          rxr(j)= DBLE(phs(1,i,j))/ABS(phs(1,i,j))
          rxi(j)=AIMAG(phs(1,i,j))/ABS(phs(1,i,j))
        Else
          rxr(j)=1.0_wp
          rxi(j)=0.0_wp
        End If
        r(1)=(1.0_wp,0.0_wp)

        ! kind=8, complex double precision
        r(j)=CMPLX( rxr(j), rxi(j), kind=wp )

        Do i1=1,dim
          Zout(i1,j)=CONJG(r(j))*Zinp(i1,j)
        End Do

        If(dbg) Then
          Write(6,'(A,i2,A,2ES24.14)') 'SPIN-PHASE:'//
     &                        ' R(',j,') = ', CONJG(r(j))
        End If
      End Do ! i


      Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
      Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
      Call zgemm_('C','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &            Zout(  1:dim,1:dim), dim,
     &              mm(1,1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim )
      Call zgemm_('N','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim,
     &            Zout(  1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             phs(1,1:dim,1:dim), dim )
cc convention:
cc    mX(i,i+1) => Real, negative
cc    mY(i,i+1) => imag, positive
cc    mZ(i,i)   => diagonal
      Do i=1,dim-1,2
        j=i+1
        If(DBLE( phs(1,i,j) ).gt.0.0_wp ) Then
          Do i1=1,dim
            Zout(i1,j)=-Zout(i1,j)
          End Do
        End If
      End Do

      If(dbg) Then

        Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
        Do l=1,3
          Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
          Call zgemm_('C','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &                Zout(  1:dim,1:dim), dim,
     &                  mm(l,1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &                 TMP(  1:dim,1:dim), dim )
          Call zgemm_('N','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &                 TMP(  1:dim,1:dim), dim,
     &                Zout(  1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &                 phs(l,1:dim,1:dim), dim )
        End Do

       Do i=1,dim
         Do j=1,dim
           Write(6,'(a,i2,a,i2,a,2ES24.14)') 'SPIN-PHASE:'//
     &      '  Zout(',i,',',j,') = ',Zout(i,j)
         End Do
       End Do

       Write(6,'(//)')
       Do i=1,dim
        Do j=1,dim
         If((j==(i-1)).or.(j==(i+1))) Then
          Write(6,'(A,i2,A,i2,A, 3(2ES24.14,3x))') 'SPIN-PHASE:'//
     &             ' PHS(',i,',',j,') = (x,y,z) =',(phs(l,i,j),l=1,3)
         Else
           cycle
         End If
        End Do
       End Do
      End If



! ------------------------------------------------------------
      Call mma_deallocate(rxr)
      Call mma_deallocate(rxi)
      Call mma_deallocate(r)
      Call mma_deallocate(phs)
      Call mma_deallocate(tmp)
      Call qExit('s_phase')

      Return
      End
c





      Subroutine SPIN_PHASE2(MM,dim,Zinp,Zout)
C
C     The RASSI program gives a ranDom phase to the spin-orbit functions.
C
C     This routine performs a simple check with the obtained spin functions,
C     in order to determine the phase of the spin functions.
C     If the phase is not the same, Then the spin functions will be multiplied
C     with the correspondind coefficient that sets the same phase to all spin
C     eigenfunctions
C
      Implicit None
#include "stdalloc.fh"
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: dim
      Complex(kind=wp), intent(in)  :: mm(3,dim,dim)
      Complex(kind=wp), intent(in)  :: Zinp(dim,dim)
      Complex(kind=wp), intent(out) :: Zout(dim,dim)
! ------------------------------------------------------------
      Integer                       :: i, j, i1, l
      Complex(kind=wp)              :: t
      Real(kind=wp), allocatable    :: rxr(:) !dim)
      Real(kind=wp), allocatable    :: rxi(:) !dim)
      Complex(kind=wp), allocatable :: r(:) !(dim)
      Complex(kind=wp), allocatable :: phs(:,:,:)  !3,dim,dim)
      Complex(kind=wp), allocatable :: tmp(:,:) !dim,dim
      Logical :: dbg

      Call qEnter('s_phase')
      dbg=.true.

      Call mma_allocate(rxr,dim,'rxr')
      Call mma_allocate(rxi,dim,'rxi')
      Call mma_allocate(r,dim,'r')
      Call mma_allocate(phs,3,dim,dim,'phs')
      Call mma_allocate(tmp,dim,dim,'tmp')
! ------------------------------------------------------------
      Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
      Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
      Call zcopy_(dim,(0.0_wp,0.0_wp),0,r,1)
      Call dcopy_(dim,0.0_wp,0,rxr,1)
      Call dcopy_(dim,0.0_wp,0,rxi,1)
      rxr(1)=1.0_wp
      rxi(1)=0.0_wp

      ! compute magnetic moment, X
      Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
      Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
      Call zgemm_('C','N', dim,  dim,  dim, (1.0_wp,0.0_wp),
     &            Zinp(  1:dim,1:dim), dim,
     &              mm(1,1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim )
      Call zgemm_('N','N', dim,  dim,  dim, (1.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim,
     &            Zinp(  1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             phs(1,1:dim,1:dim), dim )

      r(1)=(1.0_wp,0.0_wp)
      Do i=1,dim-1
        j=i+1

        If(ABS(phs(1,i,j)).gt.1.0e-14_wp ) Then
          rxr(j)= DBLE(phs(1,i,j))/ABS(phs(1,i,j))
          rxi(j)=AIMAG(phs(1,i,j))/ABS(phs(1,i,j))
        Else
          rxr(j)=1.0_wp
          rxi(j)=0.0_wp
        End If

        ! kind=8, complex double precision
        r(j)=CMPLX( rxr(j),-rxi(j), kind=wp )

        If(dbg) Then
          Write(6,'(A,i2,A,2ES24.14)') 'SPIN-PHASE:'//
     &                        ' R(',j,') = ', r(j)*r(i)
        End If
      End Do

      t=(1.0_wp,0.0_wp)
      Do j=1,dim
        t=t*r(j)
        Do i1=1,dim
          Zout(i1,j)=t*Zinp(i1,j)
        End Do
      End Do

      ! compute the momentum using the ZOUT functions:
      Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
      Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
      Call zgemm_('C','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &            Zout(  1:dim,1:dim), dim,
     &              mm(1,1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim )
      Call zgemm_('N','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &             TMP(  1:dim,1:dim), dim,
     &            Zout(  1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &             phs(1,1:dim,1:dim), dim )
cc convention:
cc    mX(i,i+1) => Real, negative
cc    mY(i,i+1) => imag, positive
cc    mZ(i,i)   => diagonal
      Do i=1,dim-1,2
        j=i+1
        If(DBLE( phs(1,i,j) ).gt.0.0_wp ) Then
          Do i1=1,dim
            Zout(i1,j)=-Zout(i1,j)
          End Do
        End If
      End Do


      If(dbg) Then

        Call zcopy_(3*dim*dim,(0.0_wp,0.0_wp),0,phs,1)
        Do l=1,3
           Call zcopy_(  dim*dim,(0.0_wp,0.0_wp),0,tmp,1)
           CALL ZGEMM_('C','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &                 Zout(  1:dim,1:dim), dim,
     &                   mm(l,1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &                  TMP(  1:dim,1:dim), dim )
           CALL ZGEMM_('N','N',  dim,  dim,  dim, (1.0_wp,0.0_wp),
     &                  TMP(  1:dim,1:dim), dim,
     &                 Zout(  1:dim,1:dim), dim, (0.0_wp,0.0_wp),
     &                  phs(l,1:dim,1:dim), dim )
        End Do

       Do i=1,dim
         Do j=1,dim
           Write(6,'(a,i2,a,i2,a,2ES24.14)') 'SPIN-PHASE:'//
     &      '  Zout(',i,',',j,') = ',Zout(i,j)
         End Do
       End Do

       Write(6,'(//)')
       Do i=1,dim
        Do j=1,diM
         If((j.eq.(i-1)).or.(j.eq.(i+1))) Then
          Write(6,'(A,i2,A,i2,A, 3(2ES24.14,3x))') 'SPIN-PHASE:'//
     &             ' PHS(',i,',',j,') = (x,y,z) =',(phs(l,i,j),l=1,3)
         Else
           cycle
         End If
        End Do
       End Do
      End If


! ------------------------------------------------------------
      Call mma_deallocate(rxr)
      Call mma_deallocate(rxi)
      Call mma_deallocate(r)
      Call mma_deallocate(phs)
      Call mma_deallocate(tmp)
      Call qExit('s_phase')

      Return
      End
c
