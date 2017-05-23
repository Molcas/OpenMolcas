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
* Copyright (C) 1998, Markus P. Fuelscher                              *
************************************************************************
      SubRoutine Guess_m(CMO)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Diagonalize core Hamiltonian to get starting orbitals.           *
*                                                                      *
*     calling arguments:                                               *
*     CMO     : real*8, output                                         *
*               starting vectors                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1998                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)

*     global definitions

#include "rasdim.fh"
#include "warnings.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='GUESS   ')
#include "rasscf.fh"
#include "WrkSpc.fh"

*     calling arguments

      Real*8 CMO(*)

*     local definitions

      Character*8 Label
      Parameter ( zero = 0.0d0 , one = 1.0d0 )

*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*

*     allocate work space
      Call GetMem('scr1','Allo','Real',iTmp1,nTot1)

*     load bare nuclei Hamiltonian

      iRc    = -1
      iOpt   =  6
      iComp  =  1
      iSyLbl =  1
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp1),iSyLbl)
      If ( iRc.ne.0 ) Then
        Write(LF,*)' RASSCF tried to construct start orbitals from'
        Write(LF,*)' diagonalization of core Hamiltonian, but ran into'
        Write(LF,*)' a severe error: Failed to read the Hamiltonian'
        Write(LF,*)' from the ONEINT file. Something may be wrong with'
        Write(LF,*)' the file.'
        Call Quit(_RC_IO_ERROR_READ_)
      End If

*     diagonalize bare nuclei Hamiltonian

      i1 = iTmp1
      i2 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        Call dCopy_(iBas*iBas,zero,0,CMO(i2),1)
        Call dCopy_(iBas,one,0,CMO(i2),iBas+1)
        Call Jacob(Work(i1),CMO(i2),iBas,iBas)
        Call JacOrd(Work(i1),CMO(i2),iBas,iBas)
        i1 = i1+(iBas*iBas+iBas)/2
        i2 = i2+iBas*iBas
      End Do

*     deallocate work space
      Call GetMem('scr1','Free','Real',iTmp1,nTot1)

      Return
      End
