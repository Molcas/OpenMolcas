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
      Subroutine Expectus(QMMethod,HmatOld,Vmat,VpolMat,Smat,MxDim
     &                   ,iVEC,nDim,lEig,iEig,ip_ExpVal)
      Implicit Real*8 (a-h,o-z)

#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension HmatOld(MxDim),Vmat(MxDim),VpolMat(MxDim),Smat(MxDim)
      Character QMMethod*5
      Logical lEig

*
*-- Take different path for different QM-method.
*
      If(QMMethod(1:5).eq.'RASSI') then
*
*---- For how many roots are the eigenvalues to be computed.
*
        If(lEig) then
          nRoots=iEig
        Else
          nRoots=nDim
        Endif
*
*---- Loop over roots and compute expectation values according to
*     well-known formulas.
*
        nDTri=nDim*(nDim+1)/2
        Call GetMem('DenTemp','Allo','Real',iDTmp,nDTri)
        Call GetMem('ExpVals','Allo','Real',ip_ExpVal,4*nRoots)
        Do 801, iRoot=1,nRoots
*
*------ Generate density matrix for relevant root.
*
          Call DensiSt(Work(iDTmp),Work(iVEC),iRoot,nDim,nDim)
*
*------ Expectation values.
*
          Work(ip_ExpVal+4*(iRoot-1)+0)=Ddot_(nDTri,Work(iDTmp),iOne
     &                                      ,HmatOld,iOne)
          Work(ip_ExpVal+4*(iRoot-1)+1)=Ddot_(nDTri,Work(iDTmp),iOne
     &                                      ,Vmat,iOne)
          Work(ip_ExpVal+4*(iRoot-1)+2)=Ddot_(nDTri,Work(iDTmp),iOne
     &                                      ,Vpolmat,iOne)
          Work(ip_ExpVal+4*(iRoot-1)+3)=Ddot_(nDTri,Work(iDTmp),iOne
     &                                      ,Smat,iOne)
801     Continue
        Call GetMem('DenTemp','Free','Real',iDTmp,nDTri)

*
*-- If its SCF we are running.
*
      Elseif(QMMethod(1:5).eq.'SCF  ') then
        nDTri=nDim*(nDim+1)/2
        Call GetMem('DenTemp','Allo','Real',iDTmp,nDTri)
        Call GetMem('ExpVals','Allo','Real',ip_ExpVal,4)
        Call Densi_MO(Work(iDTmp),Work(iVEC),1,iEig,nDim,nDim)
*
*------ Expectation values.
*
        Work(ip_ExpVal+0)=Ddot_(nDTri,Work(iDTmp),iOne,HmatOld,iOne)
        Work(ip_ExpVal+1)=Ddot_(nDTri,Work(iDTmp),iOne,Vmat,iOne)
        Work(ip_ExpVal+2)=Ddot_(nDTri,Work(iDTmp),iOne,Vpolmat,iOne)
        Work(ip_ExpVal+3)=Ddot_(nDTri,Work(iDTmp),iOne,Smat,iOne)
        Call GetMem('DenTemp','Free','Real',iDTmp,nDTri)

*
*-- Shit happens.
*
      Else
        Write(6,*)
        Write(6,*)' Now how did this happen, says Expectus!'
        Call Quit(_RC_INTERNAL_ERROR_ )
      Endif

*
*-- What's you major malfunction, numb nuts!
*
      Return
      End
