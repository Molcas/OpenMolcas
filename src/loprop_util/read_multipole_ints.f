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
      Subroutine Read_Multipole_Int(lMax,ip_sq_mu,nBas,ip_mu,Ttot,Temp,
     &                              Origin,rMPq,nElem,nBas1,nBas2,
     &                           nBasMax,nTemp,nSym,ipP,Restart,Utility)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 Temp(nTemp), Ttot(nBas2), Origin(3,0:lMax),
     &       rMPq(0:nElem-1)
      Integer ip_mu(0:nElem-1), ip_sq_mu(0:nElem-1), nBas(0:nSym-1)
      Character*8 Label
      Character*16 RunFile_dLabel, RunFile_iLabel, RunFile_iLabel2
      Logical Restart, Found, Utility
      Dimension idum(1)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Note that we will store two sets of the integrals.
*
*     Set 1: used to compute the localized properties.
*     Set 2: used in the finite differenc calculation of the
*            polarizability.
*
*     Set 1 will be desymmetrized in case of symmetry.
*
*     Pointers to Set 1 are stored in ip_sq_mu.
*     Pointers to Set 2 are stored in ip_mu.
*
      RunFile_dLabel  = 'LoProp Integrals'
      RunFile_iLabel  = 'LoProp nInts'
      RunFile_iLabel2 = 'LoProp iSyLbl'
      nInts_Tot = 0
      Call Allocate_iWork(ip_nComp,nElem)
      Call Allocate_iWork(ip_iSyLbl,nElem)
      If (Restart) Then
         Call Qpg_dArray(RunFile_dLabel,Found,nInts_tot)
         If (.Not. Found) Then
            Write(6,*) 'LoProp Integrals not available on the RunFile.'
            Call Abend()
         End If
         Call Allocate_Work(ip_all_ints,nInts_Tot)
         Call Get_dArray(RunFile_dLabel,Work(ip_all_ints),nInts_tot)
         Call Get_iArray(RunFile_iLabel,iWork(ip_nComp),nElem)
         Call Get_iArray(RunFile_iLabel2,iWork(ip_iSyLbl),nElem)
      End If
      nInts=0
      iOpt0=0
      iOpt1=1
      Label='Mltpl  X'
      mu = -1
      iOff = 0
      Do l = 0, lMax
         nComp = (l+1)*(l+2)/2
         Write (Label(8:8),'(I1)') l
         Do iComp = 1, nComp
#ifdef _DEBUG_
            Write (6,*) 'l,iComp=',l,iComp
#endif
            mu = mu + 1
            If (Restart) Then
               Call Allocate_Work(ip_mu(mu),iWork(ip_nComp+mu))
               nInts = iWork(ip_nComp+mu) - 4
               Call dCopy_(iWork(ip_nComp+mu),Work(ip_all_ints+iOff),1,
     &                                       Work(ip_mu(mu)),1)
               iSyLbl = iWork(ip_iSyLbl+mu)
               iOff = iOff + iWork(ip_nComp+mu)
            Else
               iRc=-1
               iSyLbl=0
               Call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
               nInts=idum(1)
               Call Allocate_Work(ip_mu(mu),nInts+4)
               If ( iRc.ne.0 ) Then
                  Write (6,*) 'Polar: error reading length of mu!'
                  Write (6,*) 'Mu=',mu
                  Call QTrace
                  Call Abend()
               End If
               Call RdOne(iRc,iOpt0,Label,iComp,Work(ip_mu(mu)),iSyLbl)
               If ( iRc.ne.0 ) Then
                  Write (6,*) 'Polar: error reading mu!'
                  Write (6,*) 'Mu=',mu
                  Call QTrace
                  Call Abend()
               End If
               iWork(ip_iSyLbl+mu) = iSyLbl
               iWork(ip_nComp+mu)  = nInts+4
               nInts_Tot = nInts_Tot + iWork(ip_nComp+mu)
            End If
*
*           Transform multipole moment integrals to new basis
*
            Call Allocate_Work(ip_Tmp,nBas1**2)
            Call Fzero(Work(ip_Tmp),nBas1**2)
            If (nSym.eq.1) Then
               ip_sq_mu(mu) = ip_Tmp
            Else
               Call GetMem('SMatrix','Allo','Real',ip_sq_mu(mu),
     &                     nBas1**2)
            End If
*
*           Now I square the dipole moment integral matrix because it
*           has been stored as a low triangle
*
            iOfft = ip_mu(mu)
            iOffs = ip_Tmp
            Do iSym = 0, nSym-1
               Do jSym = 0, iSym
                  ijSym=iEor(iSym,jSym)
                  If (iAnd(iSyLbl,2**ijSym).eq.0) Then
                     Go To 20
                  End If
                  If (nBas(iSym)*nBas(jSym).eq.0) Then
                     Go To 20
                  End If
                  If (iSym.eq.jSym) Then
*
*                    Now I square the overlap matrix because it has been
*                    stored as a lower triangle
*
                     Call Square(Work(iOfft),Work(iOffs),1,nBas(iSym),
     &                           nBas(iSym))
*
                     iOfft = iOfft + nBas(iSym)*(nBas(iSym)+1)/2
                  Else
                     call dcopy_(nBas(iSym)*nBas(jSym),Work(iOfft),1,
     &                                                Work(iOffs),1)
                     iOfft = iOfft + nBas(iSym)*nBas(jSym)
                  End If
                  iOffs = iOffs + nBas(iSym)*nBas(jSym)
 20               Continue
               End Do
            End Do
*
            If (nSym.ne.1) Then
*
*------------- Desymmetrize
*
               nScr=nBasMax*nBas1
               Call Allocate_Work(ipScr,nScr)
               Call FZero(Work(ip_sq_mu(mu)),nBas1**2)
               Call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,
     &                           Work(ip_sq_mu(mu)),nBas,nBas1,
     &                           Work(ipP),nSym,iSyLbl)
               Call Free_Work(ipScr)
               Call Free_Work(ip_Tmp)
*
            End If
*
*           Now I transform with my Transformation Matrix Ttot the
*           multipole moment integral matrices
*
#ifdef _DEBUG_
            Call RecPrt('Multipole Integrals in AO Basis',' ',
     &                  Work(ip_sq_mu(mu)),nBas1,nBas1)
#endif
            Call TransMu(Work(ip_sq_mu(mu)),nBas1,Ttot,Temp)
#ifdef _DEBUG_
            Call RecPrt('Multipole Integrals in LoProp Basis',' ',
     &                  Work(ip_sq_mu(mu)),nBas1,nBas1)
#endif
*
*           Pick up the nuclear value for this operator
*
            rMPq(mu)=Work(ip_mu(mu)+nInts+3)
         End Do
*
*        Pick up the origin of this operator
*
         call dcopy_(3,Work(ip_mu(mu)+nInts),1,Origin(1,l),1)
*
      End Do
      If ((.Not. Restart) .AND. (.Not. Utility)) Then
         Call Allocate_Work(ip_all_ints,nInts_Tot)
         mu = -1
         iOff = 0
         Do l = 0, lMax
            nComp = (l+1)*(l+2)/2
            Do iComp = 1, nComp
               mu = mu + 1
               Call dCopy_(iWork(ip_nComp+mu),Work(ip_mu(mu)),1,
     &                                       Work(ip_all_ints+iOff),1)
               iOff = iOff + iWork(ip_nComp+mu)
            End Do
         End Do
         Call Put_dArray(RunFile_dLabel,Work(ip_all_ints),nInts_Tot)
         Call Put_iArray(RunFile_iLabel,iWork(ip_nComp),nElem)
         Call Put_iArray(RunFile_iLabel2,iWork(ip_iSyLbl),nElem)
      End If
      If (.Not. Utility) Then
         Call Free_Work(ip_all_ints)
      End If
      Call Free_iWork(ip_nComp)
      Call Free_iWork(ip_iSyLbl)
#ifdef _DEBUG_
      Call RecPrt('Origin',' ',Origin,3,lMax+1)
      Call RecPrt('rMPq',' ',rMPq,1,nElem)
      Call xSpot('Exit  Read_Multipole_Int')
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
