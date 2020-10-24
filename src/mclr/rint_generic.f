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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine RInt_Generic(rkappa,rmos,rmoa,Fock,Q,Focki,Focka,
     &               idsym,reco,jspin)
      use Arrays, only: CMO_Inv
*
*                              ~
*     Constructs  F  = <0|[E  ,H]|0> ( + <0|[[E  , Kappa],H]|0> )
*                  pq       pq                 pq
*
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"
#include "standard_iounits.fh"
      Real*8 Fock(nDens2),focka(nDens2),rkappa(nDens2),
     &       Focki(ndens2),Q(ndens2),rMOs(*),rmoa(*)
      Integer ipAsh(2)
      Logical Fake_CMO2,DoAct
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Focki=',DDot_(nDens2,Focki,1,Focki,1)
      Write (LuWr,*) 'Focka=',DDot_(nDens2,Focka,1,Focka,1)
#endif
*
      Fact=-1.0d0
      One=1.0d0
      call dcopy_(ndens2,[0.0d0],0,Fock,1)
*
      If (.not.newCho) Then
        Call GetMem('MOTemp2','ALLO','REAL',ipMT1,nmba)
        Call GetMem('MOTemp1','ALLO','REAL',ipMT2,nmba)
        call dcopy_(nmba,[0.0d0],0,Work(ipMT1),1)
        call dcopy_(nmba,[0.0d0],0,Work(ipMT2),1)

        Call R2ElInt(rkappa,work(ipMT1),work(ipMT2),
     &               focki,focka,
     &               nDens2,idSym,ReCo,Fact,jspin)
#ifdef _DEBUGPRINT_
        Write (LuWr,*) 'MT1=',DDot_(nmba,Work(ipMT1),1,Work(ipMT1),1)
        Write (LuWr,*) 'MT2=',DDot_(nmba,Work(ipMT2),1,Work(ipMT2),1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
        If (iMethod.eq.2) Then
           Call GetMem('QTemp','ALLO','REAL',ipq,ndens2)
           Call CreQ(Q,Work(ipMT1),Work(ipG2tpp),idsym)
           Call CreQ(Work(ipQ),Work(ipMT2),Work(ipG2tmm),idsym)
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'Q=',DDot_(nDens2,Q,1,Q,1)
           Write (LuWr,*) 'ipQ=',DDot_(nDens2,Work(ipQ),1,Work(ipQ),1)
#endif
           call daxpy_(ndens2,1.0d0,Work(ipQ),1,Q,1)
           Call GetMem('QTemp','Free','REAL',ipq,ndens2)
        End If
*************************************************************************
*                                                                       *
*        Cholesky code                                                  *
*                                                                       *
*************************************************************************
      Else
        Fake_CMO2=.false.
        DoAct=.true.
*
        ipkappa  =ip_of_work(rkappa(1))
        ipFockI  =ip_of_work(FockI(1))
        ipFockA  =ip_of_work(FockA(1))
        ipMO1    =ip_of_work(rmos(1))
        ipQ      =ip_of_work(Q(1))
*
**      Form inactive density
*
        Call GetMem('DI ','ALLO','REAL',ipDI,nDens2)
        call dcopy_(nDens2,[0.0d0],0,Work(ipDI),1)

        Do iS=1,nsym
         Do iB=1,nIsh(iS)
          ip=ipCM(iS)+(ib-1)*nOrb(is)+ib-1
          Work(ipDI+ip-1)=2.0d0
         End Do
        End Do
*
**      Form AO 1-index transform inactive density
*
        Call GetMem('TmpDns2','ALLO','REAL',ipDens2,nDens2)
        Call GetMem('TDensI','ALLO','REAL',ipDLT,nDens2)
        Do iS=1,nSym
          If (nOrb(iS).ne.0) Then
            Do jS=1,nSym
              If (iEOr(iS-1,jS-1)+1.eq.idsym.and.nOrb(jS).ne.0) Then
                   Call DGEMM_('N','N',
     &                         nOrb(iS),nOrb(jS),nOrb(jS),
     &                         1.0d0,rkappa(ipMat(is,js)),nOrb(iS),
     &                         Work(ipDI+ipCM(js)-1),nOrb(jS),0.0d0,
     &                         Work(ipDens2+ipMat(iS,jS)-1),nOrb(iS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         1.0d0,Work(ipDens2+ipMat(iS,jS)-1),
     &                         nOrb(iS),Work(ipCMO+ipCM(is)-1),
     &                         nOrb(iS),0.0d0,
     &                         Work(ipDLT+ipMat(jS,iS)-1),nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),1.0d0,
     &                         Work(ipDLT+ipMat(jS,iS)-1),nOrb(iS),
     &                         Work(ipCMO+ipCM(js)-1),nOrb(jS),0.0d0,
     &                         Work(ipDens2+ipMat(iS,jS)-1),nOrb(jS))

                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         1.0d0,Work(ipDI+ipCM(js)-1),
     &                         nOrb(iS),Work(ipCMO+ipCM(is)-1),
     &                         nOrb(iS),0.0d0,
     &                         Work(ipDLT+ipMat(jS,iS)-1),nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),1.0d0,
     &                         Work(ipDLT+ipMat(jS,iS)-1),nOrb(iS),
     &                         Work(ipCMO+ipCM(js)-1),nOrb(jS),0.0d0,
     &                         Work(ipDI+ipCM(js)-1),nOrb(jS))
              EndIf
            End Do
          EndIf
        End Do
        call Fold_Mat(nSym,nOrb,Work(ipDens2),Work(ipDLT))
        Call DScal_(ndens2,ReCo,Work(ipDLT),1)
*
**      Form active density and MO coefficients
*
        nVB=0
        If (iMethod.eq.2) Then
          na2=0
          nAct=0
          nG2=0
          Do iSym=1,nSym
            nVB = nVB + nAsh(iSym)*nOrb(iSym)
            na2 = na2+nAsh(iSym)**2
            nact=nact+nAsh(iSym)
            nAG2=0
            Do jSym=1,nSym
              kSym=iEOr(jsym-1,isym-1)+1
              nAG2=nAg2+nAsh(jSym)*nAsh(kSym)
            End Do
            nG2=nG2+nAG2**2
          End Do
          nAtri=nact*(nact+1)/2
          nAtri=nAtri*(nAtri+1)/2
          Call GetMem('Cva','Allo','Real',ipAsh(1),2*nVB)
          ipAsh(2)=ipAsh(1)+nVB
          Call GetMem('DA','Allo','Real',ipDA,na2)
*
          ioff=0
          ioff1=0
          ioff4=1
          ioffD=0
          Do iS=1,nSym
            ioff2 = ioff + nOrb(iS)*nIsh(iS)
            ioff5 = ioff4+ nOrb(iS)*nIsh(iS)
            Do iB=1,nAsh(iS)
              ioff3=ioff2+nOrb(iS)*(iB-1)
              call dcopy_(nOrb(iS),Work(ipCMO+ioff3),1,
     &                  Work(ipAsh(1)+ioff1+iB-1),nAsh(iS))
             Do jB=1,nAsh(iS)
              ip=ioffD+ib+(jB-1)*nAsh(is)
              iA=nA(is)+ib
              jA=nA(is)+jb
              ip2=itri(iA,jA)
              Work(ipDA+ip-1)=Work(ipG1t+ip2-1)
             End Do
            End Do
*MGD to check
            Call DGEMM_('T','T',nAsh(iS),nOrb(iS),nOrb(iS),1.0d0,
     &                  rkappa(ioff5),nOrb(iS),
     &                  Work(ipCMO+ioff),nOrb(iS),
     &                  0.0d0,Work(ipAsh(2)+ioff1),nAsh(iS))
            ioff=ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
            ioff1=ioff1+nAsh(iS)*nOrb(iS)
            ioffD=ioffD+nAsh(iS)**2
            ioff4=ioff4+nOrb(iS)**2
          End Do
*
**      Expand 2-body density matrix
*
          Call GetMem('G2x','ALLO','REAL',ipG2x,nG2)
          ipGx=ipG2x
          Do ijS=1,nSym
            Do iS=1,nSym
              jS=iEOR(is-1,ijS-1)+1
              Do kS=1,nSym
                lS=iEOR(kS-1,ijS-1)+1
                Do kAsh=1,nAsh(ks)
                  Do lAsh=1,nAsh(ls)
                    ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
                    Do iAsh=1,nAsh(is)
                      Do jAsh=1,nAsh(js)
                        iij =itri(iAsh+nA(is),jAsh+nA(jS))
                        ipG=ipG2+itri(iij,ikl)-1
                        Work(ipGx)=Work(ipG)
                        ipGx=ipGx+1
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End If
*
**      Allocate temp arrays and zero Fock matrices
*
        Call GetMem('Coul+Exc','ALLO','REAL',ipJI,nDens2*4)
        ipKI=ipJI+nDens2
        ipJA=ipKI+nDens2
        ipKA=ipJA+nDens2
        call dcopy_(ndens2*4,[0.0d0],0,Work(ipJI),1)
        call dcopy_(ndens2,[0.0d0],0,FockA,1)
        call dcopy_(ndens2,[0.0d0],0,FockI,1)
        call dcopy_(nATri,[0.0d0],0,Work(ipMO1),1)
*
**      Compute the whole thing
*
        iread=2 ! Asks to read the half-transformed Cho vectors
        ip_CMO_Inv = ip_of_work(CMO_Inv(1))
        Call CHO_LK_MCLR(ipDLT,ipDI,ipDA,ipG2x,ipkappa,
     &                   ipJI,ipKI,ipJA,ipKA,ipFockI,ipFockA,
     &                   ipMO1,ipQ,ipAsh,ipCMO,ip_CMO_inv,
     &                   nIsh, nAsh,nIsh,DoAct,Fake_CMO2,
     &                   LuAChoVec,LuIChoVec,iread)
*
**      Calculate contribution from uncontracted indexes
*
        Do iS=1,nSym
          jS=iEOr(iS-1,iDSym-1)+1
          If (nOrb(iS)*nOrb(jS).ne.0) Then
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Reco*Fact,
     &                  Work(ipFIMO+ipCM(iS)-1),nOrb(is),
     &                  rkappa(ipMat(iS,jS)),nOrb(iS),
     &                  One,FockI(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,
     &                  rkappa(ipMat(iS,jS)),nOrb(is),
     &                  Work(ipFIMO+ipCM(jS)-1),nOrb(jS),
     &                  One,FockI(ipMat(iS,jS)),nOrb(is))
            If (iMethod.eq.2) Then
               Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),
     &                     Reco*Fact,
     &                     Work(ipFAMO+ipCM(iS)-1),nOrb(is),
     &                     rkappa(ipMat(iS,jS)),nOrb(iS),
     &                     One,FockA(ipMat(iS,jS)),nOrb(iS))
               Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,
     &                     rkappa(ipMat(iS,jS)),nOrb(is),
     &                     Work(ipFAMO+ipCM(jS)-1),nOrb(jS),
     &                     One,FockA(ipMat(iS,jS)),nOrb(is))
            End If
          End If
        End Do
*
**      Deallocate
*
        Call GetMem('TmpDns2','FREE','REAL',ipDens2,nDens2)
        Call GetMem('Coul+Exc','FREE','REAL',ipJI,nDens2*4)
        If (iMethod.eq.2) Then
          Call GetMem('G2x','Free','REAL',ipG2x,na2*na2)
          Call GetMem('Cva','Free','Real',ipAsh(1),2*nVB)
          Call GetMem('DA','Free','Real',ipDA,na2)
        EndIf
        Call GetMem('TDensI','FREE','REAL',ipDLT,nDens2)
        Call GetMem('DI','FREE','REAL',ipDI,nDens2)
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*
      Do iS=1,nSym
         jS=iEOr(iS-1,idsym-1)+1
*
*                I    A
*        F  =2( F  + F  )
*         pi     pi   pi
*
         If (nIsh(iS)*nOrb(jS).gt.0) Then
            Call DaXpY_(nIsh(iS)*nOrb(jS),2.0d0,
     &                 Focki(ipMat(js,is)),1,
     &                 Fock(ipMat(js,is)),1)
            If (iMethod.eq.2)  Call DaXpY_(nIsh(iS)*nOrb(jS),2.0d0,
     &                                    Focka(ipMat(js,is)),1,
     &                                    Fock(ipMat(js,is)),1)
         End If
*
         If (nOrb(iS).gt.0) Then
            Do iAsh=1,nAsh(jS)
               Do jAsh=1,nAsh(js)
                  ipF=ipMat(js,is)+nIsh(js)+jAsh-1
                  ipFI=ipMat(is,js)+(nIsh(js)+iAsh-1)*nOrb(is)
                  Dij=Work(ipg1t+itri(iash+nA(js),jAsh+nA(js))-1)

*                I
*        F  = F - F  D
*         ap   ap  ap ba
*
                  ipF= ipMat(is,js)+(Nish(js)+iAsh-1)*nOrb(is)
                  ipFI=ipMat(is,js)+(Nish(js)+jAsh-1)*nOrb(is)
*
*                I
*        F  = F + F  D
*         pa   pa  pb ab
*
                  Call DaXpY_(nOrb(is),Dij,
     &                       Focki(ipFI),1,
     &                       Fock(ipF),1)
               End Do
            End Do
         End If
*
*
*      F  = F  + Q
*       pa   pa   pa
*
         If (nAsh(iS)*nOrb(jS).gt.0)
     &      Call DaXpY_(nAsh(is)*nOrb(js),1.0d0,
     &                 Q(ipMatba(js,is)),1,
     &                 Fock(ipMat(js,is)+nOrb(js)*nIsh(is)),1)
*
*
*      F  = F  - Q
*       ap   ap   ap
*
      End Do
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif
*
      Call DYAX(ndens2,2.0d0,Fock,1,Focka,1)
      Do iS=1,nSym
        js=iEOR(is-1,idsym-1)+1
        If (nOrb(is)*nOrb(js).ne.0)
     &  Call DGESUB(Focka(ipMat(is,js)),nOrb(is),'N',
     &              Focka(ipMat(js,is)),nOrb(js),'T',
     %              Fock(ipMat(is,js)),nOrb(is),
     &              nOrb(is),nOrb(js))
      End Do
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif
*
      Call AddGrad(rKappa,Fock,idsym,2.0d0*fact)
      If (.not.newCho) Then
        Call GetMem('MOTemp3','ALLO','REAL',ipMT3,nmba)
        Call DZAXPY(nmba,1.0d0,Work(ipMT1),1,
     &             Work(ipMT2),1,Work(ipMT3),1)
        Call PickMO_MCLR(Work(ipMT3),rmos,idsym)
*
        If (ispop.ne.0) Then
           Call DZAXPY(nmba,-1.0d0,Work(ipMT1),1,
     &                 Work(ipMT2),1,Work(ipMT3),1)
           Call PickMO_MCLR(Work(ipMT3),rmoa,idsym)
        End If
        Call GetMem('MOTemp3','FREE','REAL',ipMT3,nmba)
        Call GetMem('MOTemp2','FREE','REAL',ipMT2,nmba)
        Call GetMem('MOTemp1','FREE','REAL',ipMT1,nmba)
      EndIf
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Exit RInt_Generic'
#endif
      Return
      End
