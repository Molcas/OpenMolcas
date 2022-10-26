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

implicit integer(A-Z)
#include "OneDat.fh"

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
write(6,*)
write(6,*) ' Auxiliary info on the ONEINT file'
write(6,*) ' ---------------------------------'
write(6,*)
write(6,*) 'pLu     =',pLu,' AuxOne(pLu)   =',AuxOne(pLu)
write(6,*) 'pOpen   =',pOpen,' AuxOne(pOpen) =',AuxOne(pOpen)
write(6,*)
write(6,*) ' TOC of the ONEINT file'
write(6,*) ' ----------------------'
write(6,*)
write(6,*) ' pFID  =',pFID,' TocOne(pFID)  =',TocOne(pFID)
write(6,*) ' pVersN=',pVersN,' TocOne(pVersN)=',TocOne(pVersN)
write(6,*) ' pTitle=',pTitle,' TocOne(pTitle)=',TocOne(pTitle)
write(6,*) ' pOp   =',pOp,' TocOne(pOp)   =',TocOne(pOp)
write(6,*) ' pSym  =',pSym,' TocOne(pSym)  =',TocOne(pSym)
write(6,*) ' pSymOp=',pSymOp,' TocOne(pSymOp)=',TocOne(pSymOp)
write(6,*) ' pBas  =',pBas,' TocOne(pBas)  =',TocOne(pBas)
write(6,*) ' pAtom =',pAtom,' TocOne(pAtom) =',TocOne(pAtom)
write(6,*) ' pCoord=',pCoord,' TocOne(pCoord)=',TocOne(pCoord)
write(6,*) ' pPot  =',pPot,' TocOne(pPot)  =',TocOne(pPot)
write(6,*) ' pCoM  =',pCoM,' TocOne(pCoM)  =',TocOne(pCoM)
write(6,*) ' pCoC  =',pCoC,' TocOne(pCoC)  =',TocOne(pCoC)
write(6,*) ' pALbl =',pALbl,' TocOne(pALbl) =',TocOne(pALbl)
write(6,*) ' pType =',pType,' TocOne(pType) =',TocOne(pType)
write(6,*) ' pChrge=',pChrge,' TocOne(pChrge)=',TocOne(pChrge)
write(6,*) ' pIndex=',pIndex,' TocOne(pIndex)=',TocOne(pIndex)
write(6,*) ' pNext =',pNext,' TocOne(pNext) =',TocOne(pNext)
write(6,*) ' pEnd  =',pEnd,' TocOne(pEnd)  =',TocOne(pEnd)
write(6,*)

!----------------------------------------------------------------------*
! Terminate                                                            *
!----------------------------------------------------------------------*
return

end subroutine DmpOne
