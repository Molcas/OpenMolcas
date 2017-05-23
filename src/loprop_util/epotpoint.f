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
      Subroutine EPotPoint(iPotte,nPick,ipPick,ipDPick,nEPP
     &                    ,ipT,ipTi,NucNr
     &                    ,nB,iAtom,jAtom,ip_Center)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.fh"

      Character*10 Label

*
*-- Loop through all points, pick out the relevant ones, obtain the
*   partial expectation value from the appropriate basis and return.
*
      nB2=nB*(nB+1)/2
      nB22=nB**2
      Call GetMem('DSq','Allo','Real',iDSq,nB22)
      Call Get_D1ao(ip_D,nDens)
      Call Dsq(Work(ip_D),Work(iDSq),1,nB,nB)
      Call GetMem('TEMP','Allo','Real',iTEMP,nB22)
      Call GetMem('DTrans','Allo','Real',iDTrans,nB22)
*
*-- Contravariant transformation of density matrix.
*
      Call DGEMM_('N','N',nB,nB,nB,1.0d0,Work(ipTi),nB,Work(iDSq),nB
     &          ,0.0d0,Work(iTEMP),nB)
      Call DGEMM_('N','T',nB,nB,nB,1.0d0,Work(iTEMP),nB,Work(ipTi),nB
     &          ,0.0d0,Work(iDTrans),nB)
      Call GetMem('Points','Allo','Real',iPP,nB2+4)
      Call GetMem('PointsSq','Allo','Real',iPSq,nB22)
      Call GetMem('PointsTr','Allo','Real',iPTr,nB22)
      Do iPoint=1,nPick
        iPo=iWork(ipPick+iPoint-1)
        Write(Label,'(A3,I5)')'EF0',iPo
        irc=-1
        iOpt=0
        iSmLbl=0
        iComp=1
        Call RdOne(irc,iOpt,Label,iComp,Work(iPP),iSmLbl)
        Call Square(Work(iPP),Work(iPSq),1,nB,nB)
*
*-- Covariant transformation of the matrix for the potential in this
*   particular point.
*
        Call DGEMM_('T','N',nB,nB,nB,1.0d0,Work(ipT),nB,Work(iPSq),nB
     &            ,0.0d0,Work(iTEMP),nB)
        Call DGEMM_('N','N',nB,nB,nB,1.0d0,Work(iTEMP),nB,Work(ipT),nB
     &            ,0.0d0,Work(iPTr),nB)
        dEx=0.0d0
        kaunter=0
*
*-- The usual stuff to get the localized value.
*
        Do iB1=1,nB
          Do iB2=1,nB
            iC1=iWork(ip_Center+iB1-1)
            iC2=iWork(ip_Center+iB2-1)
            If((iC1.eq.iAtom.and.iC2.eq.jAtom).or.
     &         (iC1.eq.jAtom.and.iC2.eq.iAtom)) then
              dEx=dEx+Work(iDTrans+kaunter)*Work(iPTr+kaunter)
            Endif
            kaunter=kaunter+1
          Enddo
        Enddo
*
*-- Accumulate in the return vector.
*
        If(iAtom.eq.jAtom) then
          Work(iPotte+iPoint-1)=-1.0d0*dEx
     &                         +dble(NucNr)/Work(ipDPick+iPoint-1)
        Else
          Work(iPotte+iPoint-1)=-1.0d0*dEx
        Endif
      Enddo

*
*-- Deallocate.
*
      Call GetMem('DSq','Free','Real',iDSq,nB22)
      Call GetMem('TEMP','Free','Real',iTEMP,nB22)
      Call GetMem('DTrans','Free','Real',iDTrans,nB22)
      Call GetMem('Points','Free','Real',iPP,nB2+4)
      Call GetMem('PointsSq','Free','Real',iPSq,nB22)
      Call GetMem('PointsTr','Free','Real',iPTr,nB22)
      Call GetMem('Dens','Free','Real',ip_D,nDens)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nEPP)
      End
