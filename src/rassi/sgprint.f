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
      Subroutine SGPrint(SGS)
      use definitions, only: iwp, u6
      use gugx, only: SGStruct
      implicit None
      Type (SGStruct), intent(in):: SGS

      integer(kind=iwp) nLev, nVert, MidLev, MVSta, MVEnd, i, iv, ic

C Unpack structure SGS:
      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
      MVEnd  =SGS%MVEnd

      Write(u6,*)' Split-Graph UGA. Graph description:'
      Write(u6,*)' Nr of levels:',nLev
      Write(u6,*)' Orbital symmetry labels:'
      Write(u6,'(1x,30i2)')(SGS%ISm(i),i=1,nLev)
      Write(u6,*)' Nr of vertices:',nVert
      Write(u6,*)
      Write(u6,*)' Vertex    L  N    A  B  C      '//
     &               'Downchain table        Upchain table'
      Write(u6,*)
      Do iv=1,nVert
        Write(u6,'(1x,i4,5x,2i3,2x,3i3,5x,4i4,5x,4i4)')iv,
     &  (SGS%DRT(iv,i-1),i=1,5),
     &  (SGS%Down(iv,ic),ic=0,3),
     &  (SGS%Up  (iv,ic),ic=0,3)
      End Do
      Write(u6,*)
      Write(u6,*)' Mid Level:',MidLev
      Write(u6,*)' Mid Vertices:',MVSta,'...',MVEnd
      Write(u6,*)
      Write(u6,*)' Modified Arc Weight table:'
      Write(u6,*)'           Coupling case number'
      Write(u6,*)' Vertex      0    1    2    3'
      Write(u6,*)
      Do iv=1,nVert
        Write(u6,'(1x,i4,5x,4i5)')
     &             iv,(SGS%MAW(iv,ic),ic=0,3)
      End Do
      End Subroutine SGPrint
