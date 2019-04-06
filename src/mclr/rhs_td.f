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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine RHS_td(Temp1,Temp2,Temp3,Temp4,
     &               Temp5,Temp6,temp7,
     &               rKappa,ipst,iDisp,lOper,
     &               CMO,jdisp,jspin,CI)
********************************************************************
*                                                                  *
*                                                                  *
*    Purpose:                                                      *
*            Read the perturbed fock operator and one electron     *
*            hamiltonian from disk and add the connection part     *
*            to the hessian.                                       *
*                                                                  *
*     In :                                                         *
*                loper : Symmetry operator for perurbation         *
*                idisp : Perturbation component                    *
*     Out                                                          *
*                rKappa: Preconditioned RHS for the perturbation   *
*                                                                  *
*     Temporary                                                    *
*                Temp1,Temp2,Temp3                                 *
*                                                                  *
*      Calling:                                                    *
*              qenter                                              *
*              rdmck                                               *
*              oneindtra                                           *
*              dgemm                                               *
*              dcopy                                               *
*              dgetmo                                              *
*              qexit                                               *
*                                                                  *
* Author: Anders Bernhardsson, 1995                                *
*         Theoretical Chemistry, University of Lund                *
********************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"
* for the integrals needed in sigma gen
#include "glbbas_mclr.fh"
#include "lbbas1.fh"

      Character*8 Label
      Logical CI
      Real*8 Temp1(nDens),rKappa(nDens),Temp4(nDens),
     &      Temp2(nDens),Temp3(nDens),CMO(nCMO),Temp5(nDens),
     &      Temp6(nDens),temp7(ndens)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
      one=1.0d0
      debug=.true.
      iRC=-1
      idsym=loper+1
      iOpt=0
      iOp=2**loper
*
*-------------------------------------------------------------------*
*
*     Read in connection matrix
*     and transform it to MO basis
*
*
*
*      Call Getmem('rhs1','CHECK','REAL',idum,idum)
      If (iAnd(ntpert(idisp),2**3).eq.8) Then
      iRC=-1
      iOpt=0
      iOp=2**loper
      Label='OVRGRD  '
      Call dRdMCK(iRC,iOpt,Label,DspVec(iDisp),Temp7,iop)
      If (iRc.ne.0) Goto 992
      ip=1
      Do iS=1,nSym
       Do jS=1,is
        If (iEOr(iS-1,jS-1).eq.loper) Then
        If (nBas(is)*nBas(js).ne.0) Then
           If (is.eq.js) Then
             Call Square(Temp7(ipMatLT(is,js)),
     &                   Temp6,
     &                   1,nBas(is),nBas(is))
             ip=ip+nBas(is)*(nBas(iS)+1)/2
           Else
             call dcopy_(nBas(iS)*nBas(jS),
     &                  Temp7(ipMatLt(is,js)),1,
     &                  Temp6,1)
           End If
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(jS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Temp6,nBas(iS),
     &                 0.0d0,Temp5,nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(jS),nBAs(jS),
     &                 1.0d0,Temp5,nBas(iS),
     &                 Work(ipCMO+ipCM(jS)-1),nBas(jS),
     &                 0.0d0,Temp1(ipMat(iS,jS)),nBas(is))
           If (is.ne.js) Then
           Call DGEMM_('T','T',
     &                 nBas(jS),nBas(iS),nBAs(jS),
     &                 1.0d0,Work(ipCMO+ipCM(jS)-1),nBas(js),
     &                 Temp6,nBas(iS),
     &                 0.0d0,Temp5,nBas(jS))
           Call DGEMM_('N','N',
     &                 nBas(js),nBas(iS),nBas(iS),
     &                 1.0d0,Temp5,nBas(jS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Temp1(ipMat(jS,iS)),nBas(jS))
          End If

        End If
        End If
       End Do
      End Do
      End If
*      Call Getmem('rhs2','CHECK','REAL',idum,idum)
*
*-------------------------------------------------------------------*
*
*     Read in derivative of hamiltonian
*
      rone=0.0d0
      If (iMethod.eq.2) Then
        Call GetMem('MOOIX','ALLO','REAL',ipMOX,n2dens)
        call dcopy_(n2dens,[0.0d0],0,Work(ipMOX),1)
      Else
        ipMOX = ip_Dummy
      End If
      Call GetMem('FIX','ALLO','REAL',ipFiX,nDens2)
      Call IntX(Work(ipFIX),temp7,temp6,temp5,temp4,rkappa,
     &          work(ipMOX),loper,idisp,rone)
*      Call Getmem('rhs3','CHECK','REAL',idum,idum)
*
*-------------------------------------------------------------------*
*
#ifdef _GRAD_
       If (idSym.eq.1) Then
       r=0.0d0
       Do iS=1,nSym
       Do k=1,nBas(iS)
        Do j=1,nAsh(is)+nish(is)
         Do i=1,nAsh(is)+nIsh(is)
          If (i.eq.j.and.i.le.nish(is).and.j.le.nish(is))
     &    Then
           rde=2.0d0
          Else If  (i.gt.nish(is).and.j.gt.nish(is)) Then
           rde=Work(ipG1-1+itri(i-nish(is)+nA(is),
     &               j-nIsh(is)+nA(is)))
          Else
           rde=0.0d0
          end if

          r=r+Temp1(ipMat(is,is)-1+nBas(is)*(j-1)+k)*
     &        Work(kint1-1+ipmat(is,is)-1+nBas(is)*(i-1)+k)*rDe
         End Do
        End Do
       End Do
       End Do
       Write(6,*) 'Connect one',r
       End If
#endif
*
*                       C O N N E C T
*
*
*                  {kappa MO}           F({D,k}){F,k}
*                                                 ~
*     Area for one index transformed integrals (pj|kl)
*
*      Call Getmem('rhs4','CHECK','REAL',idum,idum)
*
      If (iAnd(ntpert(idisp),2**3).eq.8)  Then
       If (iMethod.eq.2) Then
        Call GetMem('MOOIT','ALLO','REAL',ipMOT ,nmba)
        Call GetMem('MOOIT2','ALLO','REAL',ipMOT2,nmba)
        call dcopy_(nmba,[0.0d0],0,work(ipMOT),1)
        call dcopy_(nmba,[0.0d0],0,work(ipMOT2),1)
       End If
*
*      kappa rmo Fi Fa
*
       Call r2ElInt(Temp1,Work(ipMOT),  Work(ipMOT2),
     &             Temp4,Temp5,ndens2,iDSym,1.0d0,-0.5d0,0)

       If (imethod.eq.2) Then
        Call DaXpY_(nmba,1.0d0,Work(ipmot2),1,Work(ipmot),1 )
        Call GetMem('MOOIT2','FREE','REAL',ipMOT2,nmba)
       End If
*----- ix  ix  ~i
       call dcopy_(ndens2,[0.0d0],0,temp7,1)
*----- F  =F  + F
       Call DaXpY_(nDens2,One,Temp4,1,Work(ipFIX),1)
*
*       Call Getmem('rhs5','CHECK','REAL',idum,idum)
*

       If (iMethod.eq.2)
     &  Call CreQ_td(Temp6,Work(ipMOT),Work(ipG2sq),loper+1)
*
       Do iS=1,nSym
        jS=iEOr(iS-1,loper)+1
*------ F~=2*Fi~
        Call DaXpY_(nIsh(is)*nBas(js),2.0d0,
     &            Temp4(ipMat(js,is)),1,Temp7(ipMat(js,is)),1)
        If (iMethod.eq.2) Then
*------- F~=F~+2*FA~
         Call DaXpY_(nIsh(is)*nBas(js),2.0d0,
     &            Temp5(ipMat(js,is)),1,Temp7(ipMat(js,is)),1)
         Do iAsh=1,nAsh(iS)
          Do jAsh=1,nAsh(is)
           Dij=Work(ipg1t+itri(iash+nA(is),jAsh+nA(is))-1)
*
*           F~=F~+DFi~
*
           Call DaXpY_(nBas(jS),Dij,
     &               Temp4(ipMat(js,is)+nBas(js)*(nish(is)+iAsh-1)),1,
     &               Temp7(ipMat(js,is)+nBas(js)*(nish(is)+jAsh-1)),1)
          End Do
         End Do
*------- F~=F~+Q~
         Call DaXpY_(nAsh(is)*nBas(js),1.0d0,
     &            Temp6(ipMatba(js,is)),1,
     &            Temp7(ipMat(js,is)+nBas(js)*nIsh(is)),1)
        End If
       End Do ! is
*
*       Call Getmem('rhs6','CHECK','REAL',idum,idum)
*
#ifdef _GRAD_
       If (idsym.eq.1) then
       r2=0.0d0
       Do i=1,nSym
        Do j=1,nIsh(i)+nAsh(i)
         r2=Temp7(ipMat(i,i)+Nbas(i)*(j-1)+j-1)+r2
        End Do
       End Do
       Write(6,*) 'Connect 2',-r2
       Write(6,*) 'Connect' , 0.5d0*(r-r2)
       renc=-0.5d0*(r-r2)
       end if
#endif
*
*
*     Calculate connection contribution to hessian
*
      End If ! ntpert
      Call Hess(Temp7,rkappa,Temp1,temp4,Temp5,temp6,
     &            Temp3,loper+1,jdisp,idisp)
*
*      Call Getmem('rhs7','CHECK','REAL',idum,idum)
*
*---- F=F~+Fx
      If (iAnd(ntpert(idisp),2**3).eq.8)
     & call daxpy_(nDens,One,Temp7,1,rKappa,1)
#ifdef _GRAD_
       if (idsym.eq.1) Then
       r22=0.0d0
       Do i=1,nSym
        Do j=1,nIsh(i)+nAsh(i)
         r22=rKappa(ipMat(i,i)+Nbas(i)*(j-1)+j-1)+r22
        End Do
       End Do
       Write(6,*) 'Total 2',r22
       Call Getmem('JJJJ','ALLO','REAL',ipGGG,ndisp)
       Call Getmem('JJJJ','ALLO','REAL',ipGG,ndisp)
       iopt=0
       irc=-1
       Call RdMCK(iRC,iOpt,'NUCGRAD',1,Work(ipGGG),1)
       rone2=rone-r
       Work(ipGG-1+idisp)=0.5d0*(rone2+r22)+Work(ipGGG+idisp-1)
*      write(*,*)  Work(ipGG-1+idisp)
       Call Getmem('JJJJ','Free','REAL',ipGGG,ndisp)
       Call Getmem('JJJJ','Free','REAL',ipGG,ndisp)
        end if
#endif
*
*      Call Getmem('rhs7b','CHECK','REAL',idum,idum)
*
*     Add connection to 2el MO integrals
*
*---- Adds (pb|cd) to triangular (ab|cd)
      If (iMethod.eq.2.and.iAnd(ntpert(idisp),2**2).eq.4)
     &Call ABXpY(Work(ipMOT),Work(ipMOX),idsym)
*
      If (CI) Then
       ipMX=0
       If (iAnd(ntPert(idisp),2**3).ne.0) ipMX=ipMOX
*
*      Call Getmem('rhs7c','CHECK','REAL',idum,idum)
*
       Call CiSigma_td(0,State_Sym,iEor(State_sym-1,idsym-1)+1,
     &             ipFix,ipMx,idum,ipCI,ipst,'N')
C
*       Call RECPRT('IpST',' ',Work(ipin(ipST)),nConf1*2,1)
*       Call RECPRT('ipFix',' ',Work(ipFix),ndens2,1)
*       Call RECPRT('ipmox',' ',Work(ipmox),ndens2,1)
C
*      Call Getmem('rhs7d','CHECK','REAL',idum,idum)
*
       If (idsym.eq.1) Then
        EnA=E2_td(Work(ipFix),Work(ipmox),idsym-1,idisp)
        Call DaXpY_(nConf1,-Ena,Work(ipin(ipCI))
     &   ,1,Work(ipin(ipst)),1)
       End If
       call dscal_(nconf1,2.0d0,Work(ipin(ipst)),1)
      End If
*
      Call DYAX(ndens2,2.0d0,rkappa,1,Temp1,1)
C
*      Call Getmem('rhs8','CHECK','REAL',idum,idum)
*      Call RECPRT('IpST',' ',Work(ipin(ipST)),nConf1*2,1)
*      Stop 10
C
      Do iS=1,nSym
        js=iEOR(is-1,loper)+1
        If (nbas(is)*nBas(js).ne.0)
     &  Call DGESUB(Temp1(ipMat(is,js)),nBas(is),'N',
     &              Temp1(ipMat(js,is)),nBas(js),'T',
     %              rKappa(ipMat(is,js)),nBas(is),
     &              nBas(is),nBas(js))
      End Do
*
      Call GetMem('FIX','FREE','REAL',ipFix,ndens2)
      If (iMethod.eq.2)
     & Call GetMem('MOOIX','FREE','REAL',ipMOX,n2dens)
      If (iMethod.eq.2.and.iAnd(ntpert(idisp),2**3).eq.8)
     & CaLL GetMem('MOOIT','FREE','REAL',ipMOT,nmba)
      return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Temp2)
        Call Unused_real_array(CMO)
        Call Unused_integer(jspin)
      End If
*
 992  Write (6,*)
      Write (6,*) ' *** Error in subroutine RHS_TD ***'
      Write (6,*) ' Error when reading OVRGRD from MCKINT '
      Write (6,*)
*
      End
