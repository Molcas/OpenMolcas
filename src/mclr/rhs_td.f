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
* Author: Anders Bernhardsson, 1995                                *
*         Theoretical Chemistry, University of Lund                *
********************************************************************
      use ipPage, only: W
      use Arrays, only: G2sq, G1t
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Pointers.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
* for the integrals needed in sigma gen
#include "glbbas_mclr.fh"
#include "lbbas1.fh"
      Character(LEN=8) Label
      Logical CI
      Real*8 Temp1(nDens),rKappa(nDens),Temp4(nDens),
     &      Temp2(nDens),Temp3(nDens),CMO(nCMO),Temp5(nDens),
     &      Temp6(nDens),temp7(ndens)
      Real*8 rDum(1)
      Real*8, Allocatable:: FiX(:),MOX(:),MOT(:),MOT2(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
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
     &                 1.0d0,CMO(ipCM(iS)),nBas(iS),
     &                 Temp6,nBas(iS),
     &                 0.0d0,Temp5,nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(jS),nBAs(jS),
     &                 1.0d0,Temp5,nBas(iS),
     &                 CMO(ipCM(jS)),nBas(jS),
     &                 0.0d0,Temp1(ipMat(iS,jS)),nBas(is))
           If (is.ne.js) Then
           Call DGEMM_('T','T',
     &                 nBas(jS),nBas(iS),nBAs(jS),
     &                 1.0d0,CMO(ipCM(jS)),nBas(js),
     &                 Temp6,nBas(iS),
     &                 0.0d0,Temp5,nBas(jS))
           Call DGEMM_('N','N',
     &                 nBas(js),nBas(iS),nBas(iS),
     &                 1.0d0,Temp5,nBas(jS),
     &                 CMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,Temp1(ipMat(jS,iS)),nBas(jS))
          End If

        End If
        End If
       End Do
      End Do
      End If
*
*-------------------------------------------------------------------*
*
*     Read in derivative of hamiltonian
*
      rone=0.0d0
      If (iMethod.eq.2) Then
        Call mma_allocate(MOX,n2dens,Label='MOX')
      Else
        Call mma_allocate(MOX,1,Label='MOX')  ! Dummy allocation
      End If
      MOX(:)=Zero
      Call mma_allocate(FiX,nDens2,Label='FiX')
      Call IntX(FIX,temp7,temp6,temp5,temp4,rkappa,
     &          MOX,loper,idisp,rone)
*
*-------------------------------------------------------------------*
*
*                       C O N N E C T
*
*
*                  {kappa MO}           F({D,k}){F,k}
*                                                 ~
*     Area for one index transformed integrals (pj|kl)
*
*
      If (iAnd(ntpert(idisp),2**3).eq.8)  Then
       If (iMethod.eq.2) Then
        Call mma_allocate(MOT, nmba,Label='MOT')
        Call mma_allocate(MOT2,nmba,Label='MOT2')
        MOT (:)=Zero
        MOT2(:)=Zero
*
*       kappa rmo Fi Fa
*
*       IFG: this was outside "if (imethod.eq.2)",
*            probably a bug? ipmot & ipmot2 would be uninitialized
        Call r2ElInt(Temp1,MOT,MOT2,
     &               Temp4,Temp5,ndens2,iDSym,1.0d0,-0.5d0,0)
        Call DaXpY_(nmba,One,MOT2,1,MOT,1 )
        Call mma_deallocate(MOT2)
       End If
*----- ix  ix  ~i
       call dcopy_(ndens2,[Zero],0,temp7,1)
*----- F  =F  + F
       Call DaXpY_(nDens2,One,Temp4,1,FIX,1)
*

       If (iMethod.eq.2) Call CreQ_td(Temp6,MOT,G2sq,loper+1)
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
           Dij=G1t(itri(iash+nA(is),jAsh+nA(is)))
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
*     Calculate connection contribution to hessian
*
      End If ! ntpert
      Call Hess(Temp7,rkappa,Temp1,temp4,Temp5,temp6,
     &            Temp3,loper+1,jdisp,idisp)
*
*
*---- F=F~+Fx
      If (iAnd(ntpert(idisp),2**3).eq.8)
     & call daxpy_(nDens,One,Temp7,1,rKappa,1)
*
*     Add connection to 2el MO integrals
*
*---- Adds (pb|cd) to triangular (ab|cd)
      If (iMethod.eq.2.and.iAnd(ntpert(idisp),2**2).eq.4)
     &Call ABXpY(MOT,MOX,idsym)
*
      If (CI) Then
       If (iAnd(ntPert(idisp),2**3).ne.0) Then
          Call CiSigma_td(0,State_Sym,iEor(State_sym-1,idsym-1)+1,
     &                    Fix,nDens2,MOX,SIZE(MOX),rdum,1,ipCI,ipst,'N',
     &                    .True.)
       Else
          Call CiSigma_td(0,State_Sym,iEor(State_sym-1,idsym-1)+1,
     &                    Fix,nDens2,rdum,1,rdum,1,ipCI,ipst,'N',
     &                    .False.)
       End If
C
       irc=ipin(ipST)
       If (idsym.eq.1) Then
        EnA=E2_td(Fix,MOX,idsym-1,idisp)
        irc=ipin(ipCI)
        Call DaXpY_(nConf1,-Ena,W(ipCI)%Vec,1,W(ipST)%Vec,1)
       End If
       call dscal_(nconf1,2.0d0,W(ipST)%Vec,1)
      End If
*
      Call DYAX(ndens2,2.0d0,rkappa,1,Temp1,1)
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
      Call mma_deallocate(FIX)
      Call mma_deallocate(MOX)
      If (Allocated(MOT)) Call mma_deallocate(MOT)

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
