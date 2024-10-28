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
      use gugx, only: SGStruct
      implicit real*8 (a-h,o-z)
      Type (SGStruct) SGS

C Unpack structure SGS:
      nLev   =SGS%nLev
      nVert  =SGS%nVert
      MidLev =SGS%MidLev
      MVSta  =SGS%MVSta
      MVEnd  =SGS%MVEnd

      Write(6,*)' Split-Graph UGA. Graph description:'
      Write(6,*)' Nr of levels:',nLev
      Write(6,*)' Orbital symmetry labels:'
      Write(6,'(1x,30i2)')(SGS%ISm(i),i=1,nLev)
      Write(6,*)' Nr of vertices:',nVert
      Write(6,*)
      Write(6,*)' Vertex    L  N    A  B  C      '//
     &               'Downchain table        Upchain table'
      Write(6,*)
      Do iv=1,nVert
        Write(6,'(1x,i4,5x,2i3,2x,3i3,5x,4i4,5x,4i4)')iv,
     &  (SGS%DRT(iv,i-1),i=1,5),
     &  (SGS%Down(iv,ic),ic=0,3),
     &  (SGS%Up  (iv,ic),ic=0,3)
      End Do
      Write(6,*)
      Write(6,*)' Mid Level:',MidLev
      Write(6,*)' Mid Vertices:',MVSta,'...',MVEnd
      Write(6,*)
      Write(6,*)' Modified Arc Weight table:'
      Write(6,*)'           Coupling case number'
      Write(6,*)' Vertex      0    1    2    3'
      Write(6,*)
      Do iv=1,nVert
        Write(6,'(1x,i4,5x,4i5)')
     &             iv,(SGS%MAW(iv,ic),ic=0,3)
      End Do
      end Subroutine SGPrint
