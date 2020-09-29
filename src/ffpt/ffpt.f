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
* Copyright (C) 1992, Markus P. Fuelscher                              *
*               1999, Roland Lindh                                     *
*               2005, Anders Ohrn                                      *
************************************************************************
      Subroutine FFPT(ireturn)
************************************************************************
*                                                                      *
*                    ######  ######  #####    #####                    *
*                    #       #       #    #     #                      *
*                    #####   #####   #    #     #                      *
*                    #       #       #####      #                      *
*                    #       #       #          #                      *
*                    #       #       #          #                      *
*                                                                      *
*     A utility to perform finite field perturbation calculations      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*     Modified for dynamic memory allocation by                        *
*     R. Lindh March 1999                                              *
*                                                                      *
*     Added local perturbation by                                      *
*     A. Ohrn October 2005                                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (A-H,O-Z)
#include "input.fh"
#include "WrkSpc.fh"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
C     Call Hello
      Call MkCom
      Call Rd1Int_FFPT
*----------------------------------------------------------------------*
      nSize=0
      nTemp=0
      Do i = 1, nSym
         nSize=nSize+nBas(i)*(nBas(i)+1)/2
         nTemp=nTemp+nBas(i)
      End Do
      nTemp=nTemp**2+4
      nSize=nSize+4
      Call GetMem('H0','Allo','Real',ipH0,nSize)
      Call GetMem('Ovlp','Allo','Real',ipOvlp,nSize)
      Call GetMem('RR','Allo','Real',ipRR,nSize)
      Call GetMem('Temp','Allo','Real',ipTemp,nTemp)
*----------------------------------------------------------------------*
      Call RdInp_FFPT
      Call PrInp_FFPT
      Call PtAdd(Work(ipH0),Work(ipOvlp),Work(ipRR),nSize,
     &           Work(ipTemp),nTemp)
*----------------------------------------------------------------------*
      Call GetMem('Temp','Free','Real',ipTemp,nTemp)
      Call GetMem('RR','Free','Real',ipRR,nSize)
      Call GetMem('Ovlp','Free','Real',ipOvlp,nSize)
      Call GetMem('H0','Free','Real',ipH0,nSize)
*----------------------------------------------------------------------*
      Call FastIO('STATUS')
*
      ireturn=0
      return
      End
