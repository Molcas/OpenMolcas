************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Freq_Molden(Freq,nFreq,Vectors,nVectors,nSym,
     &                       Intens,mDisp,RedMas)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Freq(nFreq), Vectors(nVectors),Intens(nFreq),RedMas(nFreq)
      Integer mDisp(nSym)
#include "Molcas.fh"
      Character*2 Element(MxAtom*8)
*
      Call QEnter('Freq_Molden')
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('Freq',' ',Freq,1,nFreq)
      Call RecPrt('Intens',' ',Intens,1,nFreq)
      Call RecPrt('Vectors',' ',Vectors,1,nVectors)
      Write (6,*) 'mDisp=',mDisp
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Open input file for MOLDEN
*
      Lu_9=9
      Lu_9=isFreeUnit(Lu_9)
      call molcas_open(Lu_9,'MD_FREQ')
*                                                                      *
************************************************************************
*                                                                      *
      Write (Lu_9,*) '[MOLDEN FORMAT]'
*                                                                      *
************************************************************************
*                                                                      *
*---- Write frequecnies to Molden input file
*
      Write(Lu_9,*) '[N_FREQ]'
      Write(Lu_9,*) nFreq
      Write (Lu_9,*) '[FREQ]'
      Do iFreq = 1, nFreq
         Write (Lu_9,*) Freq(iFreq)
      End Do
      Write (Lu_9,*) '[INT]'
      Do iFreq = 1, nFreq
         Write (Lu_9,*) Intens(iFreq)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Write coordinates of all centers
*
      Call Get_nAtoms_All(nCoord)
      Call Allocate_Work(ipCoord,3*nCoord)
      Call Get_Coord_All(Work(ipCoord),nCoord)
#ifdef _DEBUG_
      Call RecPrt('Coord(all)',' ',Work(ipCoord),3,nCoord)
#endif
      Call Get_Name_All(Element)
      ipTemp=ipCoord
      Write (Lu_9,*) '[NATOM]'
      Write (Lu_9,*) nCoord

      Write (Lu_9,*) '[FR-COORD]'
      Do iCoord = 1, nCoord
         Write (Lu_9,*) Element(iCoord),(Work(ipTemp+i),i=0,2)
         ipTemp=ipTemp+3
      End Do
      Call Free_Work(ipCoord)
*                                                                      *
************************************************************************
*                                                                      *
*---- Write normal modes, observe that the order here is the same as
*     for the coordinates.
*
      Call Get_iScalar('Unique atoms',nUnique_Atoms)
      Call Get_nAtoms_All(nAll_Atoms)
      Call GetMem('NMode','Allo','Real',ipNMode,3*nAll_Atoms*nFreq)
      Call FZero(Work(ipNMode),3*nAll_Atoms*nFreq)
      Call Get_NMode_All(Vectors,nVectors,nFreq,nUnique_Atoms,
     &                   Work(ipNMode),nAll_Atoms,mDisp)
      Write (Lu_9,*) '[FR-NORM-COORD]'
#ifdef _DEBUG_
      Call RecPrt('Normal Modes',' ',Work(ipNMode),3*nAll_Atoms,nFreq)
#endif
      ipTemp=ipNMode
      Do iFreq = 1, nFreq
         Write (Lu_9,*) 'vibration ', iFreq
         Do iCoor = 1, nAll_Atoms
            Write (Lu_9,*) (Work(ipTemp+i),i=0,2)
            ipTemp=ipTemp+3
         End Do
      End Do
      Call GetMem('NMode','Free','Real',ipNMode,3*nAll_Atoms*nFreq)
*                                                                      *
************************************************************************
*                                                                      *
*     Alessio Valentini 2018 - add reduced masses to freq.molden file
*     in order to use this file for computation of initial conditions
*     in semiclassical molecular dynamics
*
      Write (Lu_9,*) '[RMASS]'
      Do iFreq = 1, nFreq
         Write (Lu_9,*) RedMas(iFreq)
      End do


      Close(Lu_9)
      Call QExit('Freq_Molden')
      Return
      End
