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
       Integer Function nMemAm(nShBF,nIrrep,nS,jS,iOffA,Out_of_Core)
       Integer nShBf(0:nIrrep-1,nS), iOffA(4,0:nIrrep-1)
       Logical Out_of_Core
*
       If (Out_Of_Core) Then
*
*         Only a subblock of the upper trianular storage.
*
          nMemAm=0
          Do iIrrep = 0, nIrrep-1
             nj = nShBf(iIrrep,jS)
             nl = 0
             Do lS = 1, jS
                nl = nl+ nShBf(iIrrep,lS)
             End Do
*            Offset to where this block starts
             iOffA(1,iIrrep)=nMemAm
*            # of basis functions for the jSs shell
             iOffA(2,iIrrep)=nj
*            max # number of basis functions for this block
             iOffA(4,iIrrep)=nl
*            Update nMemAm with the sise of this subblock.
             nn = nl-nj
             nMemAm = nMemAm + nl*(nl+1)/2 - nn*(nn+1)/2
          End Do
       Else
*
*         The whole A matrix is in core! Upper triangular storage.
*
          nMemAm=0
          Do iIrrep = 0, nIrrep-1
             nj = nShBf(iIrrep,jS)
             ni = 0
             Do lS = 1, jS-1
                ni = ni + nShBf(iIrrep,lS)
             End Do
             nl = ni + nj
*            Offset to where this block starts
             iOffA(1,iIrrep)=nMemAm+ni*(ni+1)/2
*            # of basis functions for the jSs shell
             iOffA(2,iIrrep)=nj
             iOffA(4,iIrrep)=nl
*            Update nMemAm with the size of the whole block for this irrep.
             Do lS = jS+1,nS
                nl = nl + nShBf(iIrrep,lS)
             End Do
             nMemAm = nMemAm + nl*(nl+1)/2
          End Do
       End If
*
       Return
       End
