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
      Subroutine SGPrint(iSGStruct)
      implicit real*8 (a-h,o-z)
#include "Struct.fh"
      Dimension iSGStruct(nSGSize)
#include "WrkSpc.fh"

C Unpack structure iSGStruct:
      nLev   =iSGStruct(2)
      lISm   =iSGStruct(3)
      nVert  =iSGStruct(4)
      lDRT   =iSGStruct(5)
      lDown  =iSGStruct(6)
      lUp    =iSGStruct(7)
      MidLev =iSGStruct(8)
      MVSta  =iSGStruct(9)
      MVEnd  =iSGStruct(10)
      lMAW   =iSGStruct(11)
      Write(6,*)' Split-Graph UGA. Graph description:'
      Write(6,*)' Nr of levels:',nLev
      Write(6,*)' Orbital symmetry labels:'
      Write(6,'(1x,30i2)')(iWork(lISm+i),i=0,nLev-1)
      Write(6,*)' Nr of vertices:',nVert
      Write(6,*)
      Write(6,*)' Vertex    L  N    A  B  C      '//
     &               'Downchain table        Upchain table'
      Write(6,*)
      Do iv=1,nVert
        Write(6,'(1x,i4,5x,2i3,2x,3i3,5x,4i4,5x,4i4)')iv,
     &  (IWork(lDRT-1+iv+nVert*(i-1)),i=1,5),
     &  (IWork(lDown-1+iv+nVert*ic),ic=0,3),
     &  (IWork(  lUp-1+iv+nVert*ic),ic=0,3)
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
     &             iv,(IWork(lMAW-1+iv+nVert*ic),ic=0,3)
      End Do
      return
      end
