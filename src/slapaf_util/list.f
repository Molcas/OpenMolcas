************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993, Roland Lindh                                     *
*               Giovanni Ghigo                                         *
************************************************************************
      SubRoutine List(Line,Lbl,gq,mInt,nIter)
************************************************************************
*                                                                      *
* Object: to print gradient or internal coordinate lists               *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             1993                                                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 gq(mInt,nIter)
      Character Lbl(mInt)*8, Format*72, Line*(*)
#include "print.fh"
#include "real.fh"
*
      iRout = 118
      iPrint = nPrint(iRout)
      Lu=6
*
      Write (Lu,*)
      Write (Lu,*)
      Write (Lu,*) Line
*
      MxWdth=132
      nLbl=8+1
      nRow=10
      inc = Min((MxWdth-nLbl)/nRow,nIter)
*
      Do 10 ii = 1, nIter, inc
         Write (Lu,*)
         Write(Format,'(A,I2,A)') '(A,1X,',inc,'(I5,5X))'
         Write (Lu,Format) 'Iter.   ',(i,i=ii,Min(ii+inc-1,nIter))
         Write (Lu,*)
         Write(Format,'(A,I2,A)') '(A,1X,',inc,'(F9.5,1X))'
         Do 20 igq = 1, mInt
            Write (Lu,Format) Lbl(igq),
     &            (gq(igq,i),i=ii,Min(ii+inc-1,nIter))
 20      Continue
         Write (Lu,*)
         Write (Lu,*)
 10   Continue
      Write (Lu,*)
*
      Return
      End

      SubRoutine ListU(Lu,Lbl,gq,mInt,nIter)
************************************************************************
*                                                                      *
* Object: to print last gradient or internal coordinate list           *
*                                                                      *
* Called from: RlxCtl                                                  *
*                                                                      *
*   Author: Giovanni Ghigo                                             *
*   Adapted from Lis by Roland Lindh, Dept. of Theoretical Chemistry,  *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 gq(mInt,nIter)
      Character Lbl(mInt)*8
#include "print.fh"
#include "real.fh"

      Write (Lu,*)
      Write (Lu,*) '****************************'
      Write (Lu,*) '* Value of internal forces *'
      Write (Lu,*) '----------------------------'
      Do igq = 1, mInt
         Write (Lu,'(1X,A8,1X,F9.5)') Lbl(igq),gq(igq,nIter)
      EndDo
      Write (Lu,*)

      Return
      End
