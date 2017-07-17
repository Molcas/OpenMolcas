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
*  Cho_X_GetVec
*
*> @brief
*>   Read Cholesky vectors
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Read Cholesky vectors \c iVec = \p iVec1, ...., \p iVec1+NumVec-1
*> of symmetry \p  iSym from file. The vectors are returned
*> in the reduced set storage defined by reduced set index array
*> location ``2``.
*>
*> @param[out] ChoVec Cholesky vectors
*> @param[in]  LenVec vector dimension
*> @param[in]  NumVec number of vectors to read
*> @param[in]  iVec1  first vector to read
*> @param[in]  iSym   vector symmetry
*> @param[in]  Scr    scratch array for read
*> @param[in]  lScr   dimension of \p Scr array
************************************************************************
      Subroutine Cho_X_GetVec(ChoVec,LenVec,NumVec,iVec1,iSym,Scr,lScr)
      Implicit None
      Integer LenVec, NumVec, iVec1, iSym, lScr
      Real*8  ChoVec(LenVec,NumVec), Scr(lScr)

      Call Cho_GetVec(ChoVec,LenVec,NumVec,iVec1,iSym,Scr,lScr)

      End
