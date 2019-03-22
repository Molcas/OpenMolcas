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
      Subroutine rotmom(MOM, N, R, MOMR )
      ! inverse rotation
      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: N
      Real(kind=wp), intent(in) :: R(3,3) !rotation matrix
!     initial momentum matrix
      Complex(kind=wp), intent(in) :: MOM(3,N,N)
!     rotated momentum matrix
      Complex(kind=wp), intent(out) :: MOMR(3,N,N)
c  local variables
      Integer :: i,j,l,k
      Complex(kind=wp) :: RC(3,3)

      Call qEnter('rotmom')
c rotate the matrix

      Call zcopy_(3*N*N,(0.0_wp,0.0_wp),0,MOMR,1)

      Do l=1,3
         Do k=1,3
            RC(l,k) = (0.0_wp,0.0_wp)
            RC(l,k) = cmplx(R(l,k),0.0_wp,kind=wp)
         End Do
      End Do

      Do i=1,N
         Do j=1,N
            Do l=1,3
               Do k=1,3
                  MOMR(l,i,j) = MOMR(l,i,j) + RC(l,k) * MOM(k,i,j)
               End Do
            End Do
         End Do
      End Do

      Call qExit('rotmom')
      Return
      End
c------------------------------------------------------------------------
      Subroutine rotmom2(MOM, N, R, MOMR)
      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: N
      Real(kind=wp),intent(in) :: R(3,3) !rotation matrix
      Complex(kind=wp),intent(in) :: MOM(3,N,N) !initial momentum matrix
!     rotated momentum matrix
      Complex(kind=wp),intent(out) :: MOMR(3,N,N)
c  local variables
      Integer i,j,l,k
      Complex(kind=wp) :: RC(3,3)

      Call qEnter('rotmom2')
c rotate the matrix
      Call zcopy_(3*N*N,(0.0_wp,0.0_wp),0,MOMR,1)

      Do l=1,3
         Do k=1,3
            RC(l,k) = (0.0_wp,0.0_wp)
            RC(l,k) = cmplx(R(l,k),0.0_wp,kind=wp)
         End Do
      End Do

      Do i=1,N
         Do j=1,N
            Do l=1,3
               Do k=1,3
                  MOMR(l,i,j) = MOMR(l,i,j) + RC(k,l) * MOM(k,i,j)
               End Do
            End Do
         End Do
      End Do

      Call qExit('rotmom2')
      Return
      End


