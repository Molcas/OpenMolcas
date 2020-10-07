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
*
* Not properly worked through. Do not use!
*
      Subroutine Mbpt2Corr(nBas,Cmo)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "qm1.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "warnings.fh"
      Real*8, Allocatable:: Diff(:)
      Dimension Cmo(MxBas**2)

      Write(6,*)
      Write(6,*)'MP2 density correction is requested.'
      Write(6,*)' -- perturbative correlation correction to the solute '
     &//'density.'
*
*-- No-no zone!
*
      Write(6,*)
      Write(6,*)'THIS OPTION IS NOT PROPERLY WORKED THROUGH! SHOULD NOT'
     &//' BE USED!'
      Call Quit(_RC_GENERAL_ERROR_)
*--- Check that the density difference is sound.
      iT=nBas*(nBas+1)/2
      Call mma_allocate(Diff,iT,Label='Diff')
      Call Get_D1ao(Diff,iT)
      If(iPrint.ge.10) then
        Call TriPrt('Non-reduced difference density matrix',' '
     &             ,Diff,nBas)
      Endif
*--- Transform density difference to orbital basis.
      Call GetMem('SqDenA','Allo','Real',ipSqD,nBas**2)
      Call GetMem('SqDenM','Allo','Real',ipSqE,nBas**2)
      Call GetMem('TEMP','Allo','Real',ipTEMP,nBas**2)
      Call GetMem('Inv','Allo','Real',iI,nBas**2)
      Call GetMem('RedSq','Allo','Real',iRedSq,nBas**2)
      call dcopy_(nBas**2,[ZERO],iZERO,Work(ipSqD),iONE)
      call dcopy_(iOrb(1)**2,[ZERO],iZERO,Work(ipSqE),iONE)
      call dcopy_(nBas*iOrb(1),[ZERO],iZERO,Work(ipTEMP),iONE)
*--- Do not forget the density matrix convention in Molcas.
      Call Dsq(Diff,Work(ipSqD),iONE,nBas,nBas)
*--- Inverse of orbital file and transformation.
      Call Minv(Cmo,Work(iI),Ising,Det,nBas)
      Call Dgemm_('N','N',nBas,nBas,nBas,ONE,Work(iI),nBas,Work(ipSqD)
     &          ,nBas,ZERO,Work(ipTEMP),nBas)
      Call Dgemm_('N','T',nBas,nBas,nBas,ONE,Work(ipTEMP),nBas,Work(iI)
     &          ,nBas,ZERO,Work(ipSqE),nBas)
*--- Remove all except the suck-out orbitals.
      kaunt1=0
      Do i=1,nBas
        Do j=1,nBas
          If(i.le.iOrb(1).and.j.le.iOrb(1)) then
            Work(iRedSq+kaunt1)=Work(ipSqE+kaunt1)
          Else
            Work(iRedSq+kaunt1)=0.0d0
          Endif
          kaunt1=kaunt1+1
        Enddo
      Enddo
*--- Make a check of the trace. Should be small.
      kaunter=0
      Trace_MP2=0
      Do 108, iB1=1,nBas
        do jjj=1,nBas
          If(iB1.eq.jjj)Trace_MP2=Trace_MP2+Work(iRedSq+kaunter)
          kaunter=kaunter+1
        enddo
108   Continue
      If(iPrint.ge.10) then
        Write(6,*)'Trace: ',Trace_MP2
      Endif
*--- Make things a bit more tidy.
      kaunt1=0
      kaunt2=0
      Do i=1,iOrb(1)
        Do j=1,nBas
          If(j.le.iOrb(1)) then
            Work(ipSqE+kaunt1)=Work(iRedSq+kaunt2)
            kaunt1=kaunt1+1
          Endif
          kaunt2=kaunt2+1
        Enddo
      Enddo
      Call SqToTri_q(Work(ipSqE),DenCorrD,iOrb(1))

*--- Transform back if we want to keep things in AO-basis. Not
*    used in QMSTAT at the present. If you wish, comment away the
*    code below 'make things a bit more tidy' and you are in
*    ready to rumble.
*      Call Dgemm_('N','N',nBas,nBas,nBas,ONE,Cmo,nBas
*     &          ,Work(iRedSq),nBas,ZERO,Work(ipTEMP),nBas)
*      Call Dgemm_('N','T',nBas,nBas,nBas,ONE,Work(ipTEMP),nBas
*     &          ,Cmo,nBas,ZERO,Work(ipSqE),nBas)
*      k=0
*      Do i=1,nBas
*        Do j=1,nBas
*          If(i.ne.j)Work(ipSqE+k)=Work(ipSqE+k)*2
*            k=k+1
*        Enddo
*      Enddo
*      Call SqToTri_q(Work(ipSqE),Work(ipTrDiffD),nBas)
*      If(iPrint.ge.10) then
*        Call TriPrt('Reduced difference density matrix',' '
*     &             ,Work(ipTrDiffD),nBas)
*      Endif

      Call mma_deallocate(Diff)
      Call GetMem('SqDenA','Free','Real',ipSqD,nBas**2)
      Call GetMem('SqDenM','Free','Real',ipSqE,nBas**2)
      Call GetMem('TEMP','Free','Real',ipTEMP,nBas**2)
      Call GetMem('Inv','Free','Real',iI,nBas**2)
      Call GetMem('RedSq','Free','Real',iRedSq,nBas**2)

      Return
      End
