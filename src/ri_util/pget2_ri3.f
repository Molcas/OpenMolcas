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
* Copyright (C) 1992,2007, Roland Lindh                                *
************************************************************************
      SubRoutine PGet2_RI3(iCmp,iBas,jBas,kBas,lBas,
     &                  Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                  DSO,DSSO,nDSO,ExFac,CoulFac,PMax,V_k,mV_k,
     &                  ZpK,nSA,nAct)
************************************************************************
*  Object: to assemble the 2nd order density matrix of a SCF wave      *
*          function from the 1st order density matrix.                 *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             January '92.                                             *
*                                                                      *
*             Modified for 3-center RI gradients, March 2007           *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use pso_stuff, only: lPSO, nnp, Thpkl, ipAorb
      use Basis_Info, only: nBas, nBas_Aux
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "real.fh"
#include "print.fh"
#include "exterm.fh"
      Real*8 PSO(nijkl,nPSO), DSO(nDSO,nSA), DSSO(nDSO), V_k(mV_k,nSA),
     &       Zpk(*)
      Integer iCmp(4), iAO(4), iAOst(4)
      Logical Shijij
*     Local Array
      Integer jSym(0:7), kSym(0:7), lSym(0:7), nAct(0:7)
      Integer iTwoj(0:7),nCumnnP(0:7),nCumnnP2(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
#include "ymnij.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      kYmnij(l)=iWork(ipYmnij(1)-1+l)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 39
      iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
      iPrint=99
      If (iPrint.ge.99) Then
         iComp = 1
         Call PrMtrx(' In PGET_RI3:DSO ',[iD0Lbl],iComp,1,D0)
         Call RecPrt('V_K',' ',V_K,1,mV_K)
         Write (6,*)
         Write (6,*) 'Distribution of Ymnij'
         Do iSym = 1, nIrrep
           If (nYmnij(iSym,1).gt.0) Then
           Write (6,*) 'iSym=',iSym
           Do i= iOff_Ymnij(iSym,1)+1,iOff_Ymnij(iSym,1)+nYmnij(iSym,1)
              Write (6,*) 'kYmnij=',kYmnij(i)
           End Do
           End If
         End Do
         Write (6,*) 'jbas,kbas,lbas=',jBas,kBas,lBas
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(Cpu1,Wall1)
*
      Fac = One / Four
      lOper=1
      PMax=Zero
      iSO = 1
      Call FZero(PSO,nijkl*nPSO)
*
      If (lPSO) Then
         nCumnnP(0)=0
         nBas_Aux(0)=nBas_Aux(0)-1
         Do i=1,nIrrep-1
           nCumnnP(i)=nCumnnP(i-1)+nnP(i-1)*nBas_Aux(i-1)
         End Do
         nBas_Aux(0)=nBas_Aux(0)+1
      EndIf
*
*     i2, j2, jBas: auxiliary basis
*     i3, j3, kBas: valence basis
*     i4, j4, lBas: valence basis
*
*     Note when j2 is symmetric then we can have both Coulomb and
*     exchange contributions, while for j2 asymmetric we will
*     only have exchange contributions
*
      MemSO2 = 0
      Do i2 = 1, iCmp(2)
         njSym = 0
         Do j = 0, nIrrep-1
           If (iAOtSO(iAO(2)+i2,j)>0) Then
               jSym(njSym) = j
               njSym = njSym + 1
            End If
         End Do
         Do i3 = 1, iCmp(3)
            nkSym = 0
            Do 301 j = 0, nIrrep-1
               If (iAOtSO(iAO(3)+i3,j)>0) Then
                  kSym(nkSym) = j
                  nkSym = nkSym + 1
               End If
301         Continue
            Do i4 = 1, iCmp(4)
               nlSym = 0
               Do 401 j = 0, nIrrep-1
                  If (iAOtSO(iAO(4)+i4,j)>0) Then
                     lSym(nlSym) = j
                     nlSym = nlSym + 1
                  End If
 401           Continue
*              Write (*,*)
*              Write (*,*) 'i2,i3,i4=',i2,i3,i4
*              Write (*,*) '=================='
*              Write (*,*)
*
*------Loop over irreps which are spanned by the basis function.
*
          Do js = 0, njSym-1
             j2 = jSym(js)
*            nJ = nChOrb(j2,iSO)
             nJ = jBas
*
             If (lPSO) Then
               ntmp=0
               Do j4=0,nIrrep-1
                 j3=iEOR(j4,j2)
                 If (j3.le.j4) nCumnnP2(j3)=ntmp
                 If (j3.eq.j4) ntmp=ntmp+nAct(j3)*(nAct(j3)+1)/2
                 If (j3.lt.j4) ntmp=ntmp+nAct(j3)*nAct(j4)
               End Do
               Do j4=0,nIrrep-1
                 j3=iEOR(j4,j2)
                 If (j3.gt.j4) nCumnnP2(j3)=nCumnnP2(j4)
               End Do
             EndIf
*
*
             Do 310 ks = 0, nkSym-1
                j3 = kSym(ks)
                j23 = iEor(j2,j3)
                nk = nYmnij(j3+1,1)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
*
*               Pointers to the full list of the X_mu,i elements
*               Note this list runs over all basis functions mu
*               (kBas*iCmp(3)). Here we only want to pick up the
*               subblock for a fixed iCmp(3) value.
*
                If (nk.lt.nChOrb(j3,iSO).and.ExFac.ne.Zero.and.
     &              nk.gt.0) Then
*
*                  Offset to where the block starts (jbas,kbas,i3)
*
                   jp_Xki=ip_CMOi(1) + iOff_CMOi(j3+1,1)
     &                   + nChOrb(j3,iSO)*(kSO-1)-1
*
*                  Loop over the auxiliary basis functions which has
*                  significant contributions to the k shell.
*
                   imo=0
                   Do k = 1, nk
                      kmo=kYmnij(k+iOff_Ymnij(j3+1,1))
*
                      jCMOk=jp_Xki+kmo
                      jr   =jr_Xki(1)+imo
                      call dcopy_(kBas,Work(jCMOk),nChOrb(j3,iSO),
     &                                 Work(jr),  nk)
*
                      imo = imo +1
                   End Do
*                  Reset pointers
                   jp_Xki=jr_Xki(1)
*                  Write (*,*) 'i3,j3=',i3,j3
*                  Call RecPrt('X(i,mu)C',' ',Work(jp_Xki),nk,kBas)
                Else If (ExFac.ne.Zero.and.nk.gt.0) Then
                   jp_Xki=ip_CMOi(1) + iOff_CMOi(j3+1,1)
     &                   + nChOrb(j3,iSO)*(kSO-1)
*                  Write (*,*) 'i3,j3=',i3,j3
*                  Call RecPrt('X(i,mu)R',' ',Work(jp_Xki),nk,kBas)
                Else
                   jp_Xki=ip_Dummy
                End If
*
                Do 410 ls = 0, nlSym-1
                   j4 = lSym(ls)
*                  Write (*,*)
*                  Write (*,*) 'j2,j3,j4=',j2,j3,j4
*                  Write (*,*)
                   If (j23.ne.j4) Go To 410
*                  Write (*,*) 'Pass symmetry check'
*                  Write (*,*) 'j2,j3,j4=',j2,j3,j4
                   nl = nYmnij(j4+1,1)
*                  Write (*,*) 'j4,nl=',j4,nl
                   lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
*                  Pointers to the full list of the X_mu,i elements
*
                   If (nl.lt.nChOrb(j4,iSO).and.ExFac.ne.Zero.and.
     &                 nl.gt.0) Then
                      jp_Xli=ip_CMOi(1) + iOff_CMOi(j4+1,1)
     &                      + nChOrb(j4,iSO)*(lSO-1)-1
                      imo=0
                      Do l = 1, nl
                         lmo=kYmnij(l+iOff_Ymnij(j4+1,1))
*
                         jCMOl=jp_Xli+lmo
                         jr   =jr_Xli(1)+imo
                         call dcopy_(lBas,Work(jCMOl),nChOrb(j4,iSO),
     &                                    Work(jr),   nl)
*
                         imo = imo +1
                      End Do
*                     Reset pointers
                      jp_Xli=jr_Xli(1)
*                     Write (*,*) 'i4,j4=',i4,j4
*                     Call RecPrt('X(j,nu)C',' ',Work(jp_Xli),nl,lBas)
                   Else If (ExFac.ne.Zero.and.nl.gt.0) Then
                      jp_Xli=ip_CMOi(1) + iOff_CMOi(j4+1,1)
     &                      + nChOrb(j4,iSO)*(lSO-1)
*                     Write (*,*) 'i4,j4=',i4,j4
*                     Call RecPrt('X(j,nu)R',' ',Work(jp_Xli),nl,lBas)
                   Else
                      jp_Xli=ip_Dummy
                   End If
*
                   MemSO2 = MemSO2 + 1
*                  Write (6,*) 'MemSO2=',MemSO2
*
                   jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                   jSO_off = jSO - nBas(j2)
*
                   ExFac_ = ExFac
                   If (nJ*nk*nl.eq.0) ExFac=Zero
*                  Write (6,*) 'i2,i3,i4=',i2,i3,i4
*                  Write (6,*) 'nJ,nk,nl=',nJ,nk,nl
*                  Write (6,*) 'ExFac=',ExFac
*                                                                      *
************************************************************************
*                                                                      *
*                  Read a block of C_kl^J and transform it to
*                  AO basis.
*
                   If (ExFac.ne.Zero) Then
*
*                  Read C(i,j,J) for a fix i2 value
*
                   lCVec = nIJR(j3+1,j4+1,iSO)*jBas
                   iAdr =  iAdrCVec(j2+1,j3+1,1)
     &                  +  nIJR(j3+1,j4+1,iSO)*(jSO_Off-1)
                   Call dDaFile(LuCVector(j2+1,1),2,Work(ip_Cijk),
     &                          lCVec,iAdr)
*                  Call RecPrt('C(ij,K)',' ',Work(ip_CijK),
*    &                         nIJR(j3+1,j4+1,iSO),jBas)
*
*                  Extract only those C_kl^Js for which we deem k and l
*                  to belong to the shell-pair and to be of
*                  significance. Use temporary memory location at
*                  ip_CilK.
*
                   If (nk*nl.lt.nChOrb(j3,iSO)*nChOrb(j4,iSO)) Then
                      ij=0
                      Do j=1,nl
                         jmo=kYmnij(j+iOff_Ymnij(j4+1,1))
                         Do i=1,nk
                            imo=kYmnij(i+iOff_Ymnij(j3+1,1))
*
                            jC=ip_Cijk-1+nChOrb(j3,iSO)*(jmo-1)+imo
                            jr=ip_CilK+ij
*
*                           For this particular ij combination pick
*                           the whole row.
*
                            call dcopy_(jBas,Work(jC),nChOrb(j3,iSO)*
     &                                               nChOrb(j4,iSO),
     &                                      Work(jr),nk*nl)
                            ij=ij+1
                         End Do
                      End Do
*
*                     Copy back to original memory position.
*
                      call dcopy_(nk*nl*jBas,Work(ip_CilK),1,
     &                                      Work(ip_CijK),1)
                   End If
*
*                  Transform according to Eq. 16 (step 4) and
*                  generate B_kl^J. This is a transformation from
*                  the MO basis, ij, to the AO basis mn.
*
*                  E(jK,m) = Sum_i C(i,jK)' * X(i,m)
*
*                  Call RecPrt('C(ij,K)',' ',Work(ip_CijK),
*    &                         nk*nl,jBas)
*                  Call RecPrt('X(i,mu)',' ',Work(jp_Xki),
*    &                          nk,kBas)
                   Call dGEMM_('T','N',nl*jBas,kBas,nk,
     &                        1.0D0,Work(ip_CijK),nk,
     &                              Work(jp_Xki),nk,
     &                        0.0D0,Work(ip_CilK),nl*jBas)
*
*                  B(Km,n) = Sum_j E(j, Km)' * X(j,n)
*
                   Call dGEMM_('T','N',jBas*kBas,lBas,nl,
     &                         1.0D0,Work(ip_CilK),nl,
     &                               Work(jp_Xli),nl,
     &                         0.0D0,Work(ip_BklK),jBas*kBas)
*                  Write (*,*) 'i2,i3,i4=',i2,i3,i4
*                  Write (*,*) 'MemSO2=',MemSO2
*                  Call RecPrt('B_kl^K',' ',Work(ip_BklK),
*    &                         jBas,kBas*lBas)

                   End If
*
*                  Write (6,*)
*                  Write (6,*) 'j2,j3,j4=',j2,j3,j4
*                  Write (6,*)
*                                                                      *
************************************************************************
*                                                                      *
*                  Active term (CASSCF and SA-CASSCF)
*
                   If (lPSO) Then
                     Call dzero(Thpkl,jBas*kBas*lBas)
                     If (nAct(j3)*nAct(j4).eq.0) Go to 21
                     Do iVec=1,nAVec
                       iMO1=1
                       iMO2=1
                       iVec_=iVec
                       fact=1.0d0
                       If (iVec.eq.2) iMO2=2
                       If (iVec.eq.3) fact=2.0d0
                       If (iVec.eq.4) Then
                         iMO1=2
                         iVec_=2
                       EndIf
*
                       Do jAOj = 0, jBas-1
                         jSOj = jSO + jAOj - nBas(j2)
                         jp=nCumnnP(j2)+(jSOj-1)*nnP(j2)+nCumnnP2(j3)
                         Do lAOl = 0, lBas-1
                           lSOl = lSO + lAOl
                           lp=(lSOl-1)*nAct(j4)+ipAOrb(j4,iMO1)

                           If (j3.eq.j4) Then
                             Do kAct=1,nAct(j3)
*Zpk(*,iVec)
                               tmp=ddot_(kAct,Zpk(jp+kAct*(kAct-1)/2+1),
     &                             1,Work(lp),1)
                               Do lAct=kAct+1,nAct(j4)
                                 tmp=tmp+Zpk(jp+lAct*(lAct-1)/2+kAct)*
     &                               Work(lp+lAct-1)
                               End Do
                               Work(ip_Cilk+kAct-1)=tmp
                             End Do
                           Else
                             If (j3.lt.j4) Then
                               Call dGeMV_('N',nAct(j3),nAct(j4),1.0d0,
     &                                    Zpk(jp+1),nAct(j3),Work(lp),1,
     &                                    0.0d0,Work(ip_CilK),1)
                             Else
                               Call dGeMV_('T',nAct(j4),nAct(j3),1.0d0,
     &                                    Zpk(jp+1),nAct(j4),Work(lp),1,
     &                                    0.0d0,Work(ip_CilK),1)
                             EndIf
                           EndIf
*
                           iThpkl= jAOj+ lAOl*kBas*jBas+1
                           kp=ipAOrb(j3,iMO2)+(kSO-1)*nAct(j3)
                           Call dGeMV_('T',nAct(j3),kBas,1.0d0,
     &                         Work(kp),
     &                         nAct(j3),Work(ip_Cilk),1,1.0d0,
     &                         Thpkl(iThpkl),jBas)

                         End Do
                       End Do
                     End Do
 21                  Continue
                   EndIf
*
*                                                                      *
************************************************************************
*                                                                      *
                   If (ExFac .ne. Zero) Then
*                                                                      *
************************************************************************
*                                                                      *
*
#define _EXCHANGE_
                     If (j3.ne.j4) Then
                       If (lPSO) Then
*                     Exchange and active contributions
#define _ACTIVE_
#include "pget2_ri3.fh"
#undef _ACTIVE_
                       Else
*                     Exchange contribution
#include "pget2_ri3.fh"
                       EndIf
                     Else
#define _COULOMB_
                       If (lPSO) Then
*                     Coulomb, Exchange and active contributions
#define _ACTIVE_
#include "pget2_ri3.fh"
#undef _ACTIVE_
                       Else
*                     Coulomb and Exchange contributions
#include "pget2_ri3.fh"
                       EndIf
#undef _COULOMB_
                     End If
#undef _EXCHANGE_
*                                                                      *
************************************************************************
*                                                                      *
                   Else If (ExFac.eq.Zero .and. j3.eq.j4) Then
*                                                                      *
************************************************************************
*                                                                      *
*
#define _COULOMB_
                       If (lPSO) Then
*                     Coulomb and active contributions
#define _ACTIVE_
#include "pget2_ri3.fh"
#undef _ACTIVE_
                       Else
*                     Coulomb only contribution
#include "pget2_ri3.fh"
                       EndIf
#undef _COULOMB_
*                                                                      *
************************************************************************
*                                                                      *
                   End If
*                                                                      *
************************************************************************
*                                                                      *
*
                   ExFac = ExFac_
*
 410            Continue
 310         Continue
          End Do
*
            End Do
         End Do
      End Do
      If (nPSO.ne.MemSO2) Then
        Write (6,*) ' PGET_RI3: nPSO.ne.MemSO2'
        Write (6,*) nPSO, MemSO2
        Call Abend
      End If
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Call RecPrt(' In PGET_RI3:PSO ',' ',PSO,nijkl,nPSO)
      End If
      Call GetMem(' Exit PGET_RI3','CHECK','REAL',iDum,iDum)
#endif

      Call CWTime(Cpu2,Wall2)
      Cpu = Cpu2 - Cpu1
      Wall = Wall2 - Wall1
      tbvec(1) = tbvec(1) + Cpu
      tbvec(2) = tbvec(2) + Wall
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(iBas)
         Call Unused_logical(Shijij)
         Call Unused_real_array(DSSO)
      End If
      End
