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
      Subroutine Cho_X_VecRd(Scr,lScr,jVec1,jVec2,iSym,jNum,iRedC,mUsed)
************************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_VecRd</Name>
*     <Syntax>Call Cho\_X\_VecRd(Scr,lScr,jVec1,jVec2,iSym,jNum,iRedC,
*                                mUsed)
*     </Syntax>
*     <Arguments>
*       \Argument{Scr}{contains the vectors on exit}{Real*8}{out}
*       \Argument{lScr}{dimension of Scr}{Integer}{in}
*       \Argument{jVec1}{first vector to read}{Integer}{in}
*       \Argument{jVec2}{last vector allowed to read}{Integer}{in}
*       \Argument{iSym}{vector symmetry}{Integer}{in}
*       \Argument{jNum}{number of vectors actually read}{Integer}{out}
*       \Argument{iRedC}{reduced set stored at location 3 on entry as
*                        well as exit}{Integer}{inout}
*       \Argument{mUsed}{amount of memory actually used
*                       (in real*8 words)}{Integer}{out}
*     </Arguments>
*     <Purpose>Read as many Cholesky vectors as possible in the range
*              jVec1,jVec1+1,jVec1+2,...,jVec2.
*     </Purpose>
*     <Dependencies>cholesky.fh and choptr.fh must be initialized
*     </Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The vectors are returned in their native storage (reduced set).
*        Starting with vector jVec1, this routine will read as many
*        vectors as possible, although at most (jVec2-jVec1+1) vectors
*        are read.
*        On entry as well as exit, iRedC identifies the
*        reduced set stored at "location 3" (-1 if none or unknown).
*        On exit, jNum is the number
*        of vectors actually read and mUsed is the memory
*        (in real*8 words) actually used.
*     </Description>
*    </DOC>
*
************************************************************************
      Implicit None
      Integer lScr, jVec1, jVec2, iSym, jNum, iRedC, mUsed
      Real*8  Scr(lScr)
#include "cholesky.fh"

      Integer l_jVec2

      If (iSym.lt.1 .or. iSym.gt.8) Then
         jNum = 0
         mUsed = 0
      Else
         l_jVec2 = min(NumCho(iSym),jVec2)
         Call Cho_VecRd(Scr,lScr,jVec1,l_jVec2,iSym,jNum,iRedC,mUsed)
      End If

      End
