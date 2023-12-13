!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine DmpOne()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Print the TOC of the one electron integral file                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: AuxOne, pALbl, pAtom, pBas, pChrge, pCoC, pCoM, pCoord, pEnd, pFID, pIndex, pNext, pOp, pOption, pPot, pSym, &
                  pSymOp, pTitle, pType, pVersN, TocOne
use Definitions, only: u6

implicit none

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,*) ' Auxiliary info on the ONEINT file'
write(u6,*) ' ---------------------------------'
write(u6,*)
write(u6,*) ' AuxOne%Lu =',AuxOne%Lu
write(u6,*) ' AuxOne%Opn=',AuxOne%Opn
write(u6,*)
write(u6,*) ' TOC of the ONEINT file'
write(u6,*) ' ----------------------'
write(u6,*)
write(u6,*) ' pFID   =',pFID,' TocOne(pFID)   =',TocOne(pFID)
write(u6,*) ' pVersN =',pVersN,' TocOne(pVersN) =',TocOne(pVersN)
write(u6,*) ' pTitle =',pTitle,' TocOne(pTitle) =',TocOne(pTitle)
write(u6,*) ' pOp    =',pOp,' TocOne(pOp)    =',TocOne(pOp)
write(u6,*) ' pSym   =',pSym,' TocOne(pSym)   =',TocOne(pSym)
write(u6,*) ' pSymOp =',pSymOp,' TocOne(pSymOp) =',TocOne(pSymOp)
write(u6,*) ' pBas   =',pBas,' TocOne(pBas)   =',TocOne(pBas)
write(u6,*) ' pAtom  =',pAtom,' TocOne(pAtom)  =',TocOne(pAtom)
write(u6,*) ' pCoord =',pCoord,' TocOne(pCoord) =',TocOne(pCoord)
write(u6,*) ' pPot   =',pPot,' TocOne(pPot)   =',TocOne(pPot)
write(u6,*) ' pCoM   =',pCoM,' TocOne(pCoM)   =',TocOne(pCoM)
write(u6,*) ' pCoC   =',pCoC,' TocOne(pCoC)   =',TocOne(pCoC)
write(u6,*) ' pALbl  =',pALbl,' TocOne(pALbl)  =',TocOne(pALbl)
write(u6,*) ' pType  =',pType,' TocOne(pType)  =',TocOne(pType)
write(u6,*) ' pChrge =',pChrge,' TocOne(pChrge) =',TocOne(pChrge)
write(u6,*) ' pIndex =',pIndex,' TocOne(pIndex) =',TocOne(pIndex)
write(u6,*) ' pNext  =',pNext,' TocOne(pNext)  =',TocOne(pNext)
write(u6,*) ' pOption=',pOption,' TocOne(pOption)=',TocOne(pOption)
write(u6,*) ' pEnd   =',pEnd,' TocOne(pEnd)   =',TocOne(pEnd)
write(u6,*)

!----------------------------------------------------------------------*
! Terminate                                                            *
!----------------------------------------------------------------------*
return

end subroutine DmpOne
