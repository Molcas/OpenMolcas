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
************************************************************************
*                                                                      *
      Subroutine Do_NInt1X(AOInt,mGrid,
     &                    TabAO1,TabAO2,nBfn,nGrid_Tot,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nBfn,nBfn,nD),
     &       TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)

      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
            Do iGrid = 1, mGrid
               ToAdd1 = ToAdd1
     &                + TabAO1(1,iGrid,iCB,iD)*TabAO2(1,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
*
      Return
*
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt2X(AOInt,mGrid,
     &                    TabAO1,TabAO2,nBfn,nGrid_Tot,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nBfn,nBfn,nD),
     &       TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*
*                                                                      *
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt3X(AOInt,mGrid,
     &                    TabAO1,TabAO2,nBfn,nGrid_Tot,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nBfn,nBfn,nD),
     &       TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*                                                                      *
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
     &               + TabAO1(5,iGrid,iCB,iD) *(TabAO2(5,iGrid,jCB)
     &                                         +TabAO2(8,iGrid,jCB)
     &                                        +TabAO2(10,iGrid,jCB))
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
*
      Return
      End
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Subroutine Do_NInt4X(AOInt,mGrid,
     &                    TabAO1,TabAO2,nBfn,nGrid_Tot,nD,mAO,nFn)
*                                                                      *
************************************************************************
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(nBfn,nBfn,nD),
     &       TabAO1(nFn,mGrid,nBfn,nD),
     &       TabAO2(mAO,mGrid,nBfn)
*
      nGrid_Tot=nGrid_Tot+mGrid*nBfn**2
*                                                                      *
************************************************************************
*                                                                      *
      Do iD = 1, nD
      Do iCB = 1, nBfn
*                                                                      *
************************************************************************
*                                                                      *
*
         Do jCB = 1, iCB
*                                                                      *
************************************************************************
*                                                                      *
            ToAdd1=Zero
*
            Do iGrid = 1, mGrid
               ToAdd1= ToAdd1
     &               + TabAO1(1,iGrid,iCB,iD) * TabAO2(1,iGrid,jCB)
     &               + TabAO1(2,iGrid,iCB,iD) * TabAO2(2,iGrid,jCB)
     &               + TabAO1(3,iGrid,iCB,iD) * TabAO2(3,iGrid,jCB)
     &               + TabAO1(4,iGrid,iCB,iD) * TabAO2(4,iGrid,jCB)
            End Do
            AOInt(iCB,jCB,iD) = ToAdd1
            AOInt(jCB,iCB,iD) = ToAdd1
*                                                                      *
************************************************************************
*                                                                      *
         End Do
      End Do
      End Do
*
      Return
      End
