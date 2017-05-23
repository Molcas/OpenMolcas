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
      Subroutine Interf(i_root,Ene,isuseene,iscasvb)
************************************************************************
*                                                                      *
*     Object: Driver toward MOLDEN interface                           *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "casvb.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='INTERF  ')
      Character*10 Filename
      Character*80 Note
      Dimension Ene(*)
*
* Local print level:
      IPRLEV=IPRLOC(7)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute memory requirements and allocate memory
*
      nB =0
      nB2=0
      Do iS=1,nSym
         nB =nB +nBas(iS)
         nB2=nB2+nBas(iS)**2
      End Do
*
      Call GetMem('OCCA','Allo','Real',ipOccA,nB)
      Call GetMem('OCCB','Allo','Real',ipOccB,nB)
      Call GetMem('ENERGY','Allo','Real',ipEA,2*nB)
      ipEB=ipEA+nB
      Call GetMem('CMOA','Allo','Real',ipCA,nB**2)
      Call GetMem('CMOB','Allo','Real',ipCB,nB**2)
      Call GetMem('AdCMOA','Allo','Real',mAdCMOA,nB2)
      Call GetMem('AdCMOB','Allo','Real',mAdCMOB,nB2)
*                                                                      *
************************************************************************
*                                                                      *
C -For the moment: Orbital energies just zero
      If (isuseene.ne.0) then
         Do i=1,nB
            Work(ipEA+i-1)=Ene(i)
            Work(ipEB+i-1)=Ene(i)
         End Do
      Else
         Call FZero(Work(ipEA),2*nB)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Get the coeff. of sym adapted basis functions (ipCA, ipCB) and
*     the spin orbital occupations (ipOccA, ipOccB)
*
      Call Dens_IF(i_root,Work(ipCA),Work(ipCB),Work(ipOccA),
     &                                          Work(ipOccB))
      Call Dens_IF_SCF(Work(ipCA),Work(mAdCMOA),'B')
      Call Dens_IF_SCF(Work(ipCB),Work(mAdCMOB),'B')
*                                                                      *
************************************************************************
*                                                                      *
*     Write out info on a temporary vector file.
*
      Note='Temporary orbital file for the MOLDEN interface.'
      LuTmp=50
      LuTmp=IsFreeUnit(LuTmp)
      iUHF=IFVB
      If (i_root.ne.0) iUHF=1
      Call WrVec_('TMPORB',LuTmp,'COE',iUHF,nSym,nBas,nBas,
     &            Work(mAdCMOA),Work(mAdCMOB),
     &            Work(ipOccA),Work(ipOccB),
     &            Work(ipEA),Work(ipEB),
     &            iDum,Note,0)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('AdCMOB','Free','Real',mAdCMOB,nB2)
      Call GetMem('AdCMOA','Free','Real',mAdCMOA,nB2)
      Call GetMem('CMOA','Free','Real',ipCA,nB**2)
      Call GetMem('CMOB','Free','Real',ipCB,nB**2)
      Call GetMem('ENERGY','Free','Real',ipEA,nB)
      Call GetMem('OCCA','Free','Real',ipOccA,nB)
      Call GetMem('OCCB','Free','Real',ipOccB,nB)
*                                                                      *
************************************************************************
*                                                                      *
      If (i_root.ne.0) Then
         If (i_root.le.9) Then
            Write(filename,'(A7,I1)') 'MD_CAS.',i_root
         Else If (i_root.le.99) Then
            Write(filename,'(A7,I2)') 'MD_CAS.',i_root
         Else If (i_root.le.999) Then
            Write(filename,'(A7,I3)') 'MD_CAS.',i_root
         Else
            filename='MD_CAS.x'
         End If
      Else
         filename='MD_CAS'
      End If
      if(iscasvb.eq.1) filename='MD_VB'
*                                                                      *
************************************************************************
*                                                                      *
*     Call the generic MOLDEN interface
*
      Call Molden_Interface(iUHF,'TMPORB',filename,.False.)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
