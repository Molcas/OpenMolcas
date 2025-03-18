!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine G2qtoG2r(G2r,G2q,nG2,nG2r)
      use Constants, only: One, Two
      use input_mclr, only: ntAsh
      Implicit None
      INTEGER nG2,nG2r
      Real*8,DIMENSION(nG2 )::G2q
      Real*8,DIMENSION(nG2r)::G2r
      INTEGER iB,jB,kB,lB,iDij,iRij,iDkl,iRkl,iijkl,iRijkl
      Real*8 Fact
      Integer i,j,itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      Do iB=1,ntash
       Do jB=1,ntash
        iDij=iTri(ib,jB)
        iRij=jb+(ib-1)*ntash
        Do kB=1,ntash
         Do lB=1,ntash
          iDkl=iTri(kB,lB)
          iRkl=lb+(kb-1)*ntash
          fact=One
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
          iijkl=itri(iDij,iDkl)
          iRijkl=itri(iRij,iRkl)
          G2r(iRijkl)=Fact*G2q(iijkl)
         End Do
        End Do
       End Do
      End Do
      End Subroutine G2qtoG2r
