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
      Subroutine MkCom
*
************************************************************************
*                                                                      *
*     Objective: Initialize the Command tabels                         *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*

#include "input.fh"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Initialize command tables                                        *
*                                                                      *
*----------------------------------------------------------------------*
*
      Do iCom=1,nCom
         Do iSub1=0,MxSub1
            Do iSub2=0,MxSub2
               ComCtl(iCom,iSub1,iSub2)=0
               Do iParm=0,MxParm
                  ComTab(iCom,iSub1,iSub2,iParm)="????"
                  ComStk(iCom,iSub1,iSub2,iParm)=.false.
                  ComVal(iCom,iSub1,iSub2,iParm)=0.0
               End Do
            End Do
         End Do
      End Do
*
*----------------------------------------------------------------------*
*     Define vocabulary                                                *
*----------------------------------------------------------------------*
*
      ComTab(1,0,0,0)='TITL'
      ComTab(2,0,0,0)='FFPT'
         ComTab(2,1,0,0)='DIPO'
            ComTab(2,1,1,0)='COMP'
               ComTab(2,1,1,1)=' X= '
               ComTab(2,1,1,2)=' Y= '
               ComTab(2,1,1,3)=' Z= '
         ComTab(2,2,0,0)='QUAD'
            ComTab(2,2,1,0)='COMP'
               ComTab(2,2,1,1)='XX= '
               ComTab(2,2,1,2)='XY= '
               ComTab(2,2,1,3)='XZ= '
               ComTab(2,2,1,4)='YY= '
               ComTab(2,2,1,5)='YZ= '
               ComTab(2,2,1,6)='ZZ= '
               ComTab(2,2,1,7)='RR= '
            ComTab(2,2,2,0)='ORIG'
               ComTab(2,2,2,1)=' X= '
               ComTab(2,2,2,2)=' Y= '
               ComTab(2,2,2,3)=' Z= '
               ComTab(2,2,2,4)=' N= '
         ComTab(2,6,0,0)='OCTU'
            ComTab(2,6,1,0)='COMP'
               ComTab(2,6,1,1)='XXX='
               ComTab(2,6,1,2)='XXY='
               ComTab(2,6,1,3)='XXZ='
               ComTab(2,6,1,4)='XYY='
               ComTab(2,6,1,5)='XYZ='
               ComTab(2,6,1,6)='XZZ='
               ComTab(2,6,1,7)='YYY='
               ComTab(2,6,1,8)='YYZ='
               ComTab(2,6,1,9)='YZZ='
               ComTab(2,6,1,10)='ZZZ='
            ComTab(2,6,2,0)='ORIG'
               ComTab(2,6,2,1)=' X= '
               ComTab(2,6,2,2)=' Y= '
               ComTab(2,6,2,3)=' Z= '
               ComTab(2,6,2,4)=' N= '
         ComTab(2,3,0,0)='EFLD'
            ComTab(2,3,1,0)='COMP'
               ComTab(2,3,1,1)=' X= '
               ComTab(2,3,1,2)=' Y= '
               ComTab(2,3,1,3)=' Z= '
            ComTab(2,3,2,0)='ORIG'
               ComTab(2,3,2,1)=' X= '
               ComTab(2,3,2,2)=' Y= '
               ComTab(2,3,2,3)=' Z= '
               ComTab(2,3,2,4)=' N= '
         ComTab(2,4,0,0)='EFGR'
            ComTab(2,4,1,0)='COMP'
               ComTab(2,4,1,1)='XX= '
               ComTab(2,4,1,2)='XY= '
               ComTab(2,4,1,3)='XZ= '
               ComTab(2,4,1,4)='YY= '
               ComTab(2,4,1,5)='YZ= '
               ComTab(2,4,1,6)='ZZ= '
            ComTab(2,4,2,0)='ORIG'
               ComTab(2,4,2,1)=' X= '
               ComTab(2,4,2,2)=' Y= '
               ComTab(2,4,2,3)=' Z= '
               ComTab(2,4,2,4)=' N= '
         ComTab(2,5,0,0)='RELA'
               ComTab(2,5,0,1)=' W= '
      ComTab(3,0,0,0)='GLBL'
      ComTab(4,0,0,0)='EXTR'
      ComTab(5,0,0,0)='END '
*
*----------------------------------------------------------------------*
*     Set control tables                                               *
*----------------------------------------------------------------------*
*
      ComCtl(1,0,0)=0
      ComCtl(2,0,0)=5
         ComCtl(2,1,0)=1
            ComCtl(2,1,1)=3
         ComCtl(2,2,0)=2
            ComCtl(2,2,1)=7
            ComCtl(2,2,2)=4
         ComCtl(2,6,0)=2
            ComCtl(2,6,1)=10
            ComCtl(2,6,2)=4
         ComCtl(2,3,0)=2
            ComCtl(2,3,1)=3
            ComCtl(2,3,2)=4
         ComCtl(2,4,0)=2
            ComCtl(2,4,1)=6
            ComCtl(2,4,2)=4
         ComCtl(2,5,0)=0
            ComCtl(2,5,1)=1
      ComCtl(3,0,0)=0
      ComCtl(4,0,0)=0
      ComCtl(5,0,0)=0
*
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*
      Return
      End
