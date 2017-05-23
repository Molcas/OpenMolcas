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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine DmpOne
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Print the TOC of the one electron integral file                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "OneRc.fh"
#include "OneFlags.fh"

#include "OneDat.fh"
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Write (6,*)
      Write (6,*) ' Auxiliary info on the ONEINT file'
      Write (6,*) ' ---------------------------------'
      Write (6,*)
      Write (6,*) 'pLu     =',pLu  ,' AuxOne(pLu)   =',AuxOne(pLu)
      Write (6,*) 'pOpen   =',pOpen,' AuxOne(pOpen) =',AuxOne(pOpen)
      Write (6,*)
      Write (6,*) ' TOC of the ONEINT file'
      Write (6,*) ' ----------------------'
      Write (6,*)
      Write (6,*) ' pFID  =',pFID  ,' TocOne(pFID)  =',TocOne(pFID)
      Write (6,*) ' pVersN=',pVersN,' TocOne(pVersN)=',TocOne(pVersN)
      Write (6,*) ' pTitle=',pTitle,' TocOne(pTitle)=',TocOne(pTitle)
      Write (6,*) ' pOp   =',pOp   ,' TocOne(pOp)   =',TocOne(pOp)
      Write (6,*) ' pSym  =',pSym  ,' TocOne(pSym)  =',TocOne(pSym)
      Write (6,*) ' pSymOp=',pSymOp,' TocOne(pSymOp)=',TocOne(pSymOp)
      Write (6,*) ' pBas  =',pBas  ,' TocOne(pBas)  =',TocOne(pBas)
      Write (6,*) ' pAtom =',pAtom ,' TocOne(pAtom) =',TocOne(pAtom)
      Write (6,*) ' pCoord=',pCoord,' TocOne(pCoord)=',TocOne(pCoord)
      Write (6,*) ' pPot  =',pPot  ,' TocOne(pPot)  =',TocOne(pPot)
      Write (6,*) ' pCoM  =',pCoM  ,' TocOne(pCoM)  =',TocOne(pCoM)
      Write (6,*) ' pCoC  =',pCoC  ,' TocOne(pCoC)  =',TocOne(pCoC)
      Write (6,*) ' pALbl =',pALbl ,' TocOne(pALbl) =',TocOne(pALbl)
      Write (6,*) ' pType =',pType ,' TocOne(pType) =',TocOne(pType)
      Write (6,*) ' pChrge=',pChrge,' TocOne(pChrge)=',TocOne(pChrge)
      Write (6,*) ' pIndex=',pIndex,' TocOne(pIndex)=',TocOne(pIndex)
      Write (6,*) ' pNext =',pNext ,' TocOne(pNext) =',TocOne(pNext)
      Write (6,*) ' pEnd  =',pEnd  ,' TocOne(pEnd)  =',TocOne(pEnd)
      Write (6,*)
*----------------------------------------------------------------------*
*     Terminate                                                         *
*----------------------------------------------------------------------*
      Return
      End
