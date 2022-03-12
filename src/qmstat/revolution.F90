!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!-- Construct rotation matrix.
!
      Subroutine Revolution(v,Rinv,Rotte)
      Implicit Real*8 (a-h,o-z)

      Dimension v(3),Rotte(3,3)
      Dimension u(3),w(3),t(3)

!
!-- Obtain base-vectors for the plane to which v in the normal vector.
!
      Call PlaneVectors(u,w,v,Rinv)

!
!-- Normalize v.
!
      t(1)=Rinv*v(1)
      t(2)=Rinv*v(2)
      t(3)=Rinv*v(3)

!
!-- Assemble rotation matrix
!
      Rotte(1,1)=u(1)
      Rotte(1,2)=u(2)
      Rotte(1,3)=u(3)
      Rotte(2,1)=w(1)
      Rotte(2,2)=w(2)
      Rotte(2,3)=w(3)
      Rotte(3,1)=t(1)
      Rotte(3,2)=t(2)
      Rotte(3,3)=t(3)

      Return
      End
