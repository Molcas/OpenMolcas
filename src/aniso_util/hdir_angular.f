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
      Subroutine hdir2(nP,L,dX,dY,dZ,Ang,iprint)
!  angstep == steps in the angular distribution. It defines the numer of points
!             in which M will be computed;
!       nP == number of points, ( nP = 360/angstep )
!        L == cartesian component of the magnetisation torque
!             If L=1 (i.e.X), rotation of the M occurs in the YZ plane
!             If L=2 (i.e.Y), rotation of the M occurs in the XZ plane
!             If L=3 (i.e.Z), rotation of the M occurs in the XY plane

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer       :: nP, L, iprint
      Real(kind=8) :: dX(nP), dY(nP), dZ(nP), Ang(nP)
      !local variables
      Integer       :: i
      Real(kind=8) :: AngStep,AngRad
      Real(kind=8) :: pi

      Call qEnter('hdir_angular')
      pi=3.1415926535897932384626433832795028841971_wp
      dX(:)=0.0_wp
      dY(:)=0.0_wp
      dZ(:)=0.0_wp
      Ang(:)=0.0_wp
      AngStep=0.0_wp
      AngRad =0.0_wp
      AngStep=360.0_wp/dble(nP-1)

      If (L.eq.1) Then
        dY(1)=1.0_wp
        dZ(1)=0.0_wp
        Do i=1,nP
          AngRad=dble(i-1)*AngStep*Pi/180.0_wp
          Ang(i)=dble(i-1)*AngStep
           dY(i)=cos(AngRad)
           dZ(i)=sin(AngRad)
        End Do
      Else If(L.eq.2) Then
        dX(1)=1.0_wp
        dZ(1)=0.0_wp
        Do i=1,nP
          AngRad=dble(i-1)*AngStep*Pi/180.0_wp+122.625_wp*Pi/180.0_wp
          Ang(i)=dble(i-1)*AngStep
           dX(i)=cos(AngRad)
           dZ(i)=sin(AngRad)
        End Do
      Else If(L.eq.3) Then
        dX(1)=1.0_wp
        dY(1)=0.0_wp
        Do i=1,nP
          AngRad=dble(i-1)*AngStep*Pi/180.0_wp
          Ang(i)=dble(i-1)*AngStep
           dX(i)=cos(AngRad)
           dY(i)=sin(AngRad)
        End Do
      Else
        Write(6,'(A   )') 'Error. Parametr L can take only Integer '//
     &                    'values 1, 2 or 3.'
        Write(6,'(A,I5)') 'Current value: L = ', L
      End If

      If(iprint.gt.2) Then
        Write(6,'(A,I5)') 'Angular grid for Magnetization Torque, '//
     &                    'Cartesian Component =',L
        Write(6,'(2x,A,4x,A,5x,3(10X,A,10x))') 'Nr.','Angle','X','Y','Z'
        Do i=1,nP
          Write(6,'(I4,F10.3,3x,3F21.14)') i,Ang(i),dX(i),dY(i),dZ(i)
        End Do
      End If

      Call qExit('hdir_angular')
      Return
      End
