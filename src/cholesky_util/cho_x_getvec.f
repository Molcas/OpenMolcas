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
      Subroutine Cho_X_GetVec(ChoVec,LenVec,NumVec,iVec1,iSym,Scr,lScr)
************************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_GetVec</Name>
*     <Syntax>Call Cho\_X\_GetVec(ChoVec,LenVec,NumVec,iVec1,iSym,
*                                 Scr,lScr)</Syntax>
*     <Arguments>
*       \Argument{ChoVec}{Cholesky vectors}{Real*8}{out}
*       \Argument{LenVec}{vector dimension}{Integer}{in}
*       \Argument{NumVec}{number of vectors to read}{Integer}{in}
*       \Argument{iVec1}{first vector to read}{Integer}{in}
*       \Argument{iSym}{vector symmetry}{Integer}{in}
*       \Argument{Scr}{scratch array for read}{Real*8}{in}
*       \Argument{lScr}{dimension of Scr array}{Integer}{in}
*     </Arguments>
*     <Purpose>Read Cholesky vectors</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        read Cholesky vectors iVec=iVec1,....,iVec1+NumVec-1
*        of symmetry iSym from file. The vectors are returned
*        in the reduced set storage defined by reduced set index array
*        location "2".
*     </Description>
*    </DOC>
*
************************************************************************
      Implicit None
      Integer LenVec, NumVec, iVec1, iSym, lScr
      Real*8  ChoVec(LenVec,NumVec), Scr(lScr)

      Call Cho_GetVec(ChoVec,LenVec,NumVec,iVec1,iSym,Scr,lScr)

      End
