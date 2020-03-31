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
************************************************************************
*                                                                      *
*  Subroutine Shell_MxSchwz:  gets max integral estimates for each     *
*                             shell pair...                            *
*                                                                      *
************************************************************************
      SubRoutine Shell_MxSchwz(nSkal,Schwz_Shl)
c----------------------------------------------------------------------
      use k2_setup
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      Integer nSkal
      Real*8 Schwz_Shl(nSkal,nSkal)
*
#include "ndarray.fh"
#include "real.fh"
#include "itmax.fh"
#include "nsd.fh"
#include "setup.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "k2.fh"
*
      nElem(i)=(i+1)*(i+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
c     Call GetMem('_scf','List','Real',iDum,iDum)
*     loop over shell pair...
      call dcopy_(nSkal*nSkal,[Zero],0,Schwz_Shl,1)
      nSDp1=nSD+1
      Do iS = 1, nSkal
        iShll= iSD( 0,iS)
        If (AuxShell(iShll) .and. iS.ne.nSkal) Go To 100
        iShell=iSD(11,iS)
        iPrimi=iSD( 5,iS)
        iCmp=iSD(2,iS)
        iAng=iSD(1,iS)
        iCnttp=iSD(13,iS)
        Do jS = 1, iS
          jShll= iSD( 0,jS)
          If (AuxShell(iShll).and..Not.AuxShell(jShll)) Go To 200
          If (AuxShell(jShll) .and. jS.eq.nSkal) Go To 200
C         Write (*,*) 'Shell_..:iS,jS=',iS,jS
          jShell=iSD(11,jS)
          jPrimj=iSD( 5,jS)
          jCmp=iSD(2,jS)
          jAng=iSD(1,jS)
          jCnttp=iSD(13,jS)
          nZeta=iPrimi*jPrimj
          If (iShell.ge.jShell) Then
            ijS = iShell*(iShell-1)/2 + jShell
          Else
            ijS = jShell*(jShell-1)/2 + iShell
          End If
          k2ij  = Indk2(1,ijS)
          nDCRR = Indk2(2,ijS)
C         Write (*,*) 'nDCRR=',nDCRR
          ijCmp=nElem(iAng)*nElem(jAng)
          If (.Not.DoGrad_) ijCmp=0
          nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
          nHm=nHm*nIrrep
          If (DoHess_) nHm=0
*         i9    = ip_abMax(nZeta)-1
          i9    = ip_EstI(nZeta)-1
          i10   = nZeta*(nDArray+2*ijCmp)+nDScalar+nHm
*         now loop over  R operator...
          If (fmass(iCnttp).eq.fmass(jCnttp)) Then
             Schwz_tmp=Data_k2(k2ij+i9)
             Do lDCRR = 1, nDCRR-1
                Schwz_tmp=Max(Schwz_tmp,Data_k2(k2ij+i10*lDCRR+i9))
             End Do
          Else
             Schwz_tmp=0.0D0
          End If
          Schwz_Shl(jS,iS)=Schwz_tmp
          Schwz_Shl(iS,jS)=Schwz_tmp
 200      Continue
        End Do
 100    Continue
      End Do
*     Call RecPrt('Schwz_shl',' ',Schwz_Shl,nSkal,nSkal)
      return
      end
