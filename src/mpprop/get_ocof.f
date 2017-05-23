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
       Subroutine Get_OCOF(nPrim,nBas,Work,nVec_p,OCOF)
       Implicit REAL*8 (a-h,o-z)

#include "Ocof.fh"
       DIMENSION  Work(nVec_p)

       iStdOut = 6 ! Added by EB
       iVec_p = 0
       Do i=1,nBas
          Do j=1,nPrim
             OCOF(i,j) = Work((i-1)*nPrim+j)
             iVec_p = iVec_p +1
             If (iVec_p.gt.nVec_p) Then
                Write(iStdOut,*) 'iVec_p.gt.nVec_p'
                Write(iStdOut,*) iVec_p,'.gt.',nVec_p
                Write(iStdOut,*) 'nPrim=',nPrim
                Call ABEND()
             End If
          End Do
       End Do
       Return
       End
