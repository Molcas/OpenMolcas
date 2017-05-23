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
      Subroutine TdmTrans(nBas)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension nBas(MxSym)

      Character TDMchar*6

      Logical Exist

*
*--- Sag hej till publiken.
*
      Write(6,*)
      Write(6,*)'     Transforming the transition density matrices.'

*
*--- Inquire if the ToFile is in WorkDir.
*
      Call f_Inquire(RassiM,Exist)
      If(.not.Exist) then
        Write(6,*)
        Write(6,*)'No Transition density matrix file found.'
        Write(6,*)'Did you use the TOFIle keyword in RASSI?'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif
      Call f_Inquire(EigV,Exist)
      If(.not.Exist) then
        Write(6,*)
        Write(6,*)'No Rassi eigenvectors found.'
        Write(6,*)'Did you use the TOFIle keyword in RASSI?'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

*
*--- Compute number of 'primitive' states.
*
      nStatePrim=0
      Do 111, i=1,NrFiles
        nStatePrim=nStatePrim+NrStates(i)
111   Continue

*
*--- Open EigV file and read information.
*
      Lu=92
      Call DaName(Lu,EigV)
      iDisk=0

*
*--- Read RASSCF overlap and H-matrix.
*
      nSize=nStatePrim*(nStatePrim+1)/2
      Call GetMem('NonOrtH','Allo','Real',iNonH,nSize)
      Call GetMem('NonOrtS','Allo','Real',iNonS,nSize)
      kaunt=0
      Do 201, i=1,nStatePrim
        Do 203, j=1,i
          Call dDaFile(Lu,2,Work(iNonH+kaunt),1,iDisk)
          kaunt=kaunt+1
203     Continue
201   Continue
      kaunt=0
      Do 202, i=1,nStatePrim
        Do 204, j=1,i
          Call dDaFile(Lu,2,Work(iNonS+kaunt),1,iDisk)
          kaunt=kaunt+1
204     Continue
202   Continue
      If(iPrint.ge.10) then
        Call TriPrt('RASSCF Hamiltonian',' ',Work(iNonH),nStatePrim)
        Call TriPrt('RASSCF Overlaps',' ',Work(iNonS),nStatePrim)
      Endif
      Call DaClos(Lu)

*
*--- Construct CASSI state basis.
*
      Call ContRASBas(nBas,nStatePrim,iNonH,iNonS,iEig2)
      Call GetMem('NonOrtH','Free','Real',iNonH,nSize)
      Call GetMem('NonOrtS','Free','Real',iNonS,nSize)

*
*--- Now transform from 'primitive' RASSCF to 'contracted' RASSI states.
*
      Call RasRasTrans(nBas(1),nStatePrim,iEig2,iPrint)

*
*--- If requested, obtain reduced MO-basis, otherwise just go as
*    usual.
*
      If(MoAveRed) then
        Call MoReduce(nBas,nRedMO,ipAvRed)
        Write(TDMchar,'(A)')'TDMSCR'
        Call FetchTDM(nRedMO,nState,iBigT,TDMchar)
      Else
        Write(6,*)'     ----- Use AO-representation of the transition'
     &//' density matrix.'
        nRedMO=0 !Only a dummy.
      Endif

*
*--- Finished!
*
      Write(6,*)'     ...Done!'

      Return
      End
