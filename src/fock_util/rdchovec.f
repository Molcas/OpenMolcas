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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine RdChoVec(Vec,nDim,nVec,iVec1,Lunit)
************************************************************
*
*   <DOC>
*     <Name>RdChoVec</Name>
*     <Syntax>Call RdChoVec(Vec,nDim,nVec,iVec1,Lunit)</Syntax>
*     <Arguments>
*       \Argument{Vec}{Array containing the vectors}{Real*8}
*                {out}
*       \Argument{nDim}{Vector dimension}{Integer}{in}
*       \Argument{nVec}{Number of vectors}{Integer}{in}
*       \Argument{iVec1}{Index of first vector}{Integer}{in}
*       \Argument{Lunit}{Logical unit number}{Integer}{in}
*     </Arguments>
*     <Purpose>
*        Read dense Cholesky vectors in SO basis
*     </Purpose>
*     <Dependencies>
*        Cholesky vectors must have been sorted into the
*        dense representation (by the Cholesky utility)
*        and the storage file must be open
*     </Dependencies>
*     <Author>
*        T. B. Pedersen
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Read dense SO Cholesky vectors L$^J_{pq}$ where
*        J=iVec1,iVec1+1,...,iVec1+nVec-1 and indices
*        p >= q.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Real*8 Vec(*)
      Integer nDim, nVec, iVec1, Lunit

      Integer iOpt, iAdr, nTot

      iOpt = 2
      iAdr = nDim*(iVec1 - 1) + 1
      nTot = nDim*nVec
      Call dDaFile(Lunit,iOpt,Vec,nTot,iAdr)

      End
