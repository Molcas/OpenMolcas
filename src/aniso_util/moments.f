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
      Subroutine moments(N,MS,MM,iprint)

      Implicit None
#include "stdalloc.fh"
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: N,iprint
      Complex(kind=wp), intent(in) :: MM(3,N,N), MS(3,N,N)

      Integer                       :: I,J,L,i1,i2,iDir
      Real(kind=wp)                 :: g_e
      Complex(kind=wp), allocatable :: Z(:,:) ! N,N
      Complex(kind=wp), allocatable :: AMS(:,:,:)
      Complex(kind=wp), allocatable :: AML(:,:,:)
      Complex(kind=wp), allocatable :: AMM(:,:,:) !(3,N,N),
      Complex(kind=wp), allocatable :: Mf(:,:), Sf(:,:), Lf(:,:) !(3,3)
!-----------------------------------------------------------------------
      Call qEnter('moments')

      g_e=2.0023193043718_wp

      If(N<1) Return
      Call mma_allocate(Z,N,N,'Z')
      Call mma_allocate(AMS,3,N,N,'AMS')
      Call mma_allocate(AML,3,N,N,'AML')
      Call mma_allocate(AMM,3,N,N,'AMM')
      Call mma_allocate(Mf,3,3,'Mf')
      Call mma_allocate(Sf,3,3,'Sf')
      Call mma_allocate(Lf,3,3,'Lf')
!-----------------------------------------------------------------------
      Call zcopy_(3*3,[(0.0_wp,0.0_wp)],0,Lf,1)
      Call zcopy_(3*3,[(0.0_wp,0.0_wp)],0,Mf,1)
      Call zcopy_(3*3,[(0.0_wp,0.0_wp)],0,Sf,1)
      Do iDir=1,3
        Call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,Z,1)
        Call zcopy_(3*N*N,[(0.0_wp,0.0_wp)],0,AMM,1)
        Call zcopy_(3*N*N,[(0.0_wp,0.0_wp)],0,AML,1)
        Call zcopy_(3*N*N,[(0.0_wp,0.0_wp)],0,AMS,1)

        Call pseudospin(MM,N,Z,iDir,1,1)

        Do l=1,3
          Do i=1,N
            Do j=1,N
              Do i1=1,N
                Do i2=1,N
        AMM(l,i,j) = AMM(l,i,j) + MM(l,i1,i2) * conjg(Z(i1,i))*Z(i2,j)
        AMS(l,i,j) = AMS(l,i,j) + MS(l,i1,i2) * conjg(Z(i1,i))*Z(i2,j)
                End Do
              End Do
        AML(l,i,j)= -AMM(l,i,j) - g_e * AMS(l,i,j)
            End Do
          End Do
          Mf(iDir,l)=AMM(l,1,1)
          Sf(iDir,l)=AMS(l,1,1)
          Lf(iDir,l)=AML(l,1,1)
        End Do

        If(iprint.gt.3) Then
          Write(6,*)
          Write(6,'(2x,a)') 'MOMENTS:  AMM(l,:,:)'
          Do l=1,3
            Write(6,*)
            Write(6,'(a,i3)') 'PROJECTION2 =' , l
            Do i=1,N
              Write(6,'(20(2F12.8,2x))') (AMM(l,i,j), j=1,N)
            End Do
          End Do
          Write(6,*)
          Write(6,'(2x,a)') 'MOMENTS:  AMS(l,:,:)'
          Do l=1,3
            Write(6,*)
            Write(6,'(a,i3)') 'PROJECTION2 =' , l
            Do i=1,N
              Write(6,'(20(2F12.8,2x))') (AMS(l,i,j), j=1,N)
            End Do
          End Do
          Write(6,*)
          Write(6,'(2x,a)') 'MOMENTS:  AML(l,:,:)'
          Do l=1,3
            Write(6,*)
            Write(6,'(a,i3)') 'PROJECTION2 =' , l
            Do i=1,N
              Write(6,'(20(2F12.8,2x))') (AML(l,i,j), j=1,N)
            End Do
          End Do
        End If
      End Do !iDir

      Write(6,'(A)') '--------------|--- H || Xm --|'//
     &                              '--- H || Ym --|'//
     &                              '--- H || Zm --|'
      Write(6,'(65A)') ('-',i=1,14),'|', ( ('-',i=1,14),'|',j=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') ' <1| mu_X |1> |',
     & ( DBLE(Mf(iDir,1)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') ' <1| mu_Y |1> |',
     & ( DBLE(Mf(iDir,2)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') ' <1| mu_Z |1> |',
     & ( DBLE(Mf(iDir,3)),'|',iDir=1,3 )
      Write(6,'(65a)') ('-',i=1,14),'|', ( ('-',i=1,14),'|',j=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| L_X |1> |',
     & ( DBLE(Lf(iDir,1)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| L_Y |1> |',
     & ( DBLE(Lf(iDir,2)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| L_Z |1> |',
     & ( DBLE(Lf(iDir,3)),'|',iDir=1,3 )
      Write(6,'(65a)') ('-',i=1,14),'|', ( ('-',i=1,14),'|',j=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| S_X |1> |',
     & ( DBLE(Sf(iDir,1)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| S_Y |1> |',
     & ( DBLE(Sf(iDir,2)),'|',iDir=1,3 )
      Write(6,'(A,3(F13.9,1x,A))') '  <1| S_Z |1> |',
     & ( DBLE(Sf(iDir,3)),'|',iDir=1,3 )
      Write(6,'(65a)') ('-',i=1,59), '|'

!-----------------------------------------------------------------------
      Call mma_deallocate(Z)
      Call mma_deallocate(AMS)
      Call mma_deallocate(AML)
      Call mma_deallocate(AMM)
      Call mma_deallocate(Mf)
      Call mma_deallocate(Sf)
      Call mma_deallocate(Lf)

      Call qExit('moments')
      Return
      End
