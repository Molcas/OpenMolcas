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
* Copyright (C) 2000, Markus P. Fuelscher                              *
************************************************************************
      Subroutine ClnMO(CMO)

************************************************************************
*                                                                      *
*     In order to preserve symmetry of orbitals which belong           *
*     to a point group of higher order than those allowed in           *
*     MOLCAS, it is sometimes necessary to enforce the MO coefficients *
*     to remain zero. This subroutine does that job in every           *
*     iteration.                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 2000                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     This subroutine replaces the one written earlier by L. Serrano   *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 2000                                 *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"

      Dimension CMO(*)

* Prelude

      Call qEnter ('Clen')
* Body

      iOff = ipCleanMask-1
      iiOff=0
      Do iSym = 1,nSym
        mBas = nBas(iSym)
        Do i = 1,mBas
          ii = (i-1)*mBas+iOff
          iii = (i-1)*mBas+iiOff

          Do j = 1,mBas
            ij = j+ii
            iij= j+iii
            If ( iWork(ij).eq.1 ) CMO(iij) = 0.0D0
          End Do
        End Do
        iOff = iOff+mBas*mBas
        iiOff = iiOff+mBas*mBas
 
      End Do

* Epilogue

      Call qExit('Clen')

      Return
      End
