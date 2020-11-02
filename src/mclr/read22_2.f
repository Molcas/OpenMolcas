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
      SubRoutine Read22_2(MO1,Fock,Q,FockI,FockA,Temp2,Scr,Temp3)
********************************************************************
*                                                                  *
*   Constructs         everything                                  *
*                                                                  *
*                                                                  *
*                                                                  *
*   Output:MO     :MO integrals                                    *
*          Fock   :Fock matrix (one index transformed integrals)   *
*          MOtilde:MO (one index transformed integrals)            *
*                                                                  *
********************************************************************
      use Arrays, only: CMO, Int1, G1t
      Implicit Real*8(a-h,o-z)
#include "real.fh"
#include "Pointers.fh"
#include "standard_iounits.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
      Real*8 Fock(nDens2),FockI(nDens2),FockA(nDens2),
     &       Temp2(nDens2),Temp3(ndens2),Q(nDens2),
     &       MO1(*), Scr(*)
      Logical Fake_CMO2,DoAct
      Real*8, Allocatable:: DLT(:), JA(:), KA(:), Ash(:), DA(:), G2x(:)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(ndens2,[0.0d0],0,focki,1)
      call dcopy_(ndens2,[0.0d0],0,focka,1)
      If(TwoStep.and.(StepType.eq.'RUN2')) Then
        iaddressQDAT=0
        Call dcopy_(ndens2,[0.0d0],0,fock,1)
        Call dcopy_(ndens2,[0.0d0],0,Q,1)
        Call ddafile(LuQDAT,2,FockA,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,FockI,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,Fock ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,Q    ,nDens2,iaddressQDAT)
        goto 101
      End If

*
      nas=0
      Do is=1,nSym
       nAS=nAS+nAsh(is)
      end do
*
      If (newCho) Go to 50
      Do iS=1,nSym
         Do jS=1,iS
            ijS=iEOr(iS-1,jS-1)+1
            Do kS=1,nSym
               Do lS=1,ks
                  If (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS).eq.0)
     &                Go To 100
                  If (iEOr(kS-1,lS-1).ne.ijS-1) Goto 100
*                                                                      *
************************************************************************
*                                                                      *
                  ijB1=0
                  Do iB=1,nB(iS)
                     nnB=nB(jS)
                     If (iS.eq.jS) nnB=iB
                     Do jB=1,nnB
*
*                       Read a symmetry block
*
                        Call COUL(kS,lS,iS,jS,iB,jB,Temp2,Scr)
                        ijB1=ijB1+1
*                                                                      *
************************************************************************
*                                                                      *
*                       Add to unrotated inactive fock matrix
*
*                       Fkl sum(i)     2(ii|kl) -(ki|li)
*
*
*                       Coulomb term: F  =2(ii|kl)
*                                      kl
*
                        If (iS.eq.jS.and.iB.eq.jB .and.
     &                     (iB.le.nIsh(iS)))
     &                     Call DaXpY_(nOrb(kS)*nOrb(lS),Two,
     &                                Temp2,1,Focki(ipCM(kS)),1)

*                                                                      *
************************************************************************
*                                                                      *
*                       Add to unrotated active fock matrix
*
*                       Coulomb term: F  =2(ij|kl)d   i<j
*                                      kl          ij
*
                        If (iMethod.eq.iCASSCF) Then
*
                           If (iS.eq.jS) Then
                              If (((iB.gt.nIsh(is)).and.(nAsh(iS).ne.0))
     &                                .and.
     &                            ((jB.gt.nIsh(js)).and.(nAsh(jS).ne.0))
     &                           ) Then
                                 ipD=iTri(jB-nIsh(jS)+nA(jS),
     &                               iB-nIsh(is)+nA(iS))
                                 Fact=Two
                                 If (iB.eq.jB) Fact=one
                                 Call DaXpY_(nOrb(kS)*nOrb(lS),
     &                                      Fact*G1t(ipD),
     &                                      Temp2,1,FockA(ipCM(kS)),1)
                              End If
                           End If
*
                        End If
*                                                                      *
************************************************************************
*                                                                      *
                     End Do  ! jB
                  End Do     ! iB
*                                                                      *
************************************************************************
*                                                                      *
 100              Continue
               End Do          ! lS
            End Do             ! kS
         End Do                ! jS
      End Do                   ! iS
*                                                                      *
************************************************************************
*                                                                      *
*    Construct Q matrix: Q = sum(jkl)(pj|kl)d
*                         pi                 ijkl
      If (iMethod.eq.iCASSCF) Then
         Call CreQ2(Q,Work(ipG2),1,Temp2,Scr,nDens2)
*
*        Sort out MO (ij|kl)
*
         Do iS=1,nSym
            Do jS=1,iS
               Do kS=1,iS
                  lS=iEOR(iEOr(iS-1,jS-1),kS-1)+1
                  If (lS.gt.kS.or.(iS.eq.kS.and.lS.gt.jS)) Goto 123
*
                  Do iB=1,nAsh(iS)
                     iib=ib+nA(iS)
                     nnb=nAsh(jS)
                     If (iS.eq.jS) nnb=ib
                     Do jB=1,nnB
                        jjb=jb+nA(jS)
*
                        Call Coul(kS,lS,iS,jS,iB+nIsh(iS),jB+nIsh(jS),
     &                            Temp2,Scr)
*
                        nnK=nAsh(kS)
                        If (iS.eq.kS)  nnK=iB
                        Do kB=1,nnk
                           kkb=kb+nA(kS)
                           nnL=nAsh(lS)
                           If (kS.eq.lS)  nnL=kB
                           If (iib.eq.kkb) nnL=jB
                           Do lB=1,nnL
                              llb=lb+nA(lS)
*
                              ip2 = (lB+nIsh(lS)-1)*nBas(kS)
     &                            + kB+nIsh(kS)
*
                              ip1=iTri(iTri(iib,jjb),iTri(kkb,llb))
                              MO1(ip1) = Temp2(ip2)
                           End Do
                        End Do
                     End Do
                  End Do
*
 123              Continue
               End Do
            End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         kS=iS
         Do js=1,nSym
            lS=jS
            If (iEor(iEor(is-1,js-1),iEor(ks-1,ls-1)).ne.0) Go To 200
            If (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS).eq.0)   Go To 200
            JLB=1
            JLB1=0
            JLBas=0
            jlBB=0
            Do LB=1,nB(LS)
               Do JB=1,nB(JS)
*                                                                      *
************************************************************************
*                                                                      *
                  Call EXCH(is,js,ks,ls,jb,lb,Temp2,Scr)
*                                                                      *
************************************************************************
*                                                                      *
*                 Add to unrotated inactive fock matrix
*
*                 Exchange term: F  =-(ij|kj)
*                                 ik
*
                  If (jS.eq.lS.and.jB.eq.lB.and.
     &                (jB.le.nIsh(jS)))
     &               Call DaXpY_(nOrb(iS)*nOrb(kS),-One,
     &                          Temp2,1,Focki(ipCM(iS)),1)

*                                                                      *
************************************************************************
*                                                                      *
*                 Add to unrotated active fock matrix
*
*                 Exchange term: F  =-1/2(ij|kl)d
*                                 ik             jl
*
                  If (iMethod.eq.iCASSCF) Then
                     If (jS.eq.lS) Then
                        If (((jB.gt.nIsh(js)).and.(nAsh(jS).ne.0)).and.
     &                      ((lB.gt.nIsh(ls)).and.(nAsh(lS).ne.0))
     &                     ) Then
                          ipD=iTri(lB-nIsh(lS)+nA(lS),
     &                             jB-nIsh(js)+nA(jS))
                          Call DaXpY_(nOrb(iS)*nOrb(kS),-half*G1t(ipD),
     &                               Temp2,1,FockA(ipCM(iS)),1)
                        End If
                     End If
                  End If
*                                                                      *
************************************************************************
*                                                                      *
               End Do
            End Do
 200        Continue
         End Do
      End Do
*
************************************************************************
*                                                                      *
*     Cholesky code                                                    *
*                                                                      *
************************************************************************
 50   Continue
      If (NewCho) Then
        Fake_CMO2=.true.
        DoAct=.true.
*
**      Construct inactive density matrix
*
        call dcopy_(nDens2,[0.0d0],0,temp2,1)
        Do is=1,nSym
          Do iB=1,nIsh(is)
            ip=ipCM(iS)+(ib-1)*nOrb(is)+ib-1
            Temp2(ip)=2.0d0
          End Do
        End Do
*
**      Transform to AO basis
*
        Do iS=1,nSym
           If (nIsh(iS).ne.0) Then
              jS=iS
              Call DGEMM_('T','T',nIsh(jS),nOrb(iS),nIsh(iS),
     &                    1.0d0,Temp2(ipCM(iS)),nOrb(iS),
     &                    CMO(ipCM(is)),nOrb(iS),
     &                    0.0d0,Temp3(ipMat(jS,iS)),nOrb(jS))
              Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nIsh(iS),
     &                    1.0d0,Temp3(ipMat(jS,iS)),nOrb(iS),
     &                    CMO(ipCM(js)),nOrb(jS),
     &                    0.0d0,Temp2(ipCM(iS)),nOrb(jS))
           EndIf
        End Do
*
        Call mma_allocate(DLT,nDens2,Label='DLT')
        call Fold_Mat(nSym,nOrb,Temp2,DLT)
*
**      Form active CMO and density
*
        nAct=0
        If (iMethod.eq.iCASSCF) Then
          nVB=0
          na2=0
          nG2=0
          Do iSym=1,nSym
            nVB = nVB + nAsh(iSym)*nOrb(iSym)
            na2=na2+nAsh(iSym)**2
            nAct=nAct+nAsh(iSym)
            nAG2=0
            Do jSym=1,nSym
              kSym=iEOr(jsym-1,isym-1)+1
              nAG2=nAg2+nAsh(jSym)*nAsh(kSym)
            End Do
            nG2=nG2+nAG2**2
          End Do
          Call mma_allocate(Ash,nVB,Label='Ash')
          Call mma_allocate(DA,na2,Label='DA')
*
          ioff=0
          ioff1=0
          ioffA=0
          Do iSym=1,nSym
            ioff2 = ioff + nOrb(iSym)*nIsh(iSym)
            do ikk=1,nAsh(iSym)
               ioff3=ioff2+nOrb(iSym)*(ikk-1)
               call dcopy_(nOrb(iSym),CMO(1+ioff3),1,
     &                   Ash(ioff1+ikk),nAsh(iSym))
               ik=ikk+nA(iSym)
               Do ill=1,ikk-1
                 il=ill+nA(iSym)
                 ikl=ik*(ik-1)/2+il
                 DA(ioffA+(ikk-1)*nAsh(iSym)+ill)=G1t(ikl)
                 DA(ioffA+(ill-1)*nAsh(iSym)+ikk)=G1t(ikl)
               End Do
               ikl=ik*(ik-1)/2+ik
               DA(ioffA+(ikk-1)*nAsh(iSym)+ikk)=G1t(ikl)
            End Do
            ioff=ioff+nOrb(iSym)**2
            ioff1=ioff1+nAsh(iSym)*nOrb(iSym)
            ioffA=ioffA+nAsh(iSym)*nAsh(iSym)
*            iofftA=ioffA+nAsh(iSym)*(nAsh(iSym)+1)/2
          End Do
          Call DScal_(na2,half,DA,1)
*
**      Expand 2-body density matrix
*
          Call mma_allocate(G2x,nG2,Label='G2x')
          ipGx=0
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
                        ipGx=ipGx+1
                        G2x(ipGx)=Work(ipG)
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        EndIf
*
**      Let's go
*
        ipDI    = ip_of_work(Temp2(1))
        ipkappa = ip_Dummy
        ipJI    = ip_of_work(Temp3(1))
        ipKI    = ip_of_work(Scr(1))
        ipFockI = ip_of_work(FockI(1))
        ipFockA = ip_of_work(FockA(1))
        ipMO1   = ip_of_work(MO1(1))
        ipQ     = ip_of_work(Q(1))
        Call mma_allocate(JA,nDens2,Label='JA')
        Call mma_allocate(KA,nDens2,Label='KA')
*
        call dcopy_(nDens2,[0.0d0],0,Temp3,1)
        call dcopy_(nDens2,[0.0d0],0,Scr,1)
        call dcopy_(nDens2,[0.0d0],0,FockI,1)
        call dcopy_(nDens2,[0.0d0],0,FockA,1)
        JA(:)=0.0D0
        KA(:)=0.0D0
        call dcopy_(nDens2,[0.0d0],0,Q,1)
*
        istore=1 ! Ask to store the half-transformed vectors
! BIGOT FIXME
        Call WarningMessage(2,
     &     'There is probably a bug here, ipAsh should have two '//
     &     'elements.')
        Call Abend()
!       ip_CMO_inv = ip_of_work(CMO_Inv(1))
!       ipCMO      = ip_of_work(CMO(1))
        ipDLT      = ip_of_Work(DLT(1))
        ipJA       = ip_of_Work(JA(1))
        ipKA       = ip_of_Work(KA(1))
        ipAsh      = ip_of_Work(Ash(1))
        ipDA       = ip_of_Work(DA(1))
        ipG2X      = ip_of_Work(G2X(1))
!       Call CHO_LK_MCLR(ipDLT,ipDI,ipDA,ipG2x,ipkappa,
!    &                   ipJI,ipKI,ipJA,ipKA,ipFockI,ipFockA,
!    &                   ipMO1,ipQ,ipAsh,ipCMO,ip_CMO_inv,
!    &                   nIsh, nAsh,nIsh,DoAct,Fake_CMO2,
!    &                   LuAChoVec,LuIChoVec,istore)
        nAtri=nAct*(nAct+1)/2
        nAtri=nAtri*(nAtri+1)/2
        Call DScal_(nAtri,0.25D0,MO1,1)
        Call DScal_(nDens2,-0.5d0,FockI,1)
*
        Call mma_deallocate(JA)
        Call mma_deallocate(KA)
        Call mma_deallocate(DLT)
        If (iMethod.eq.iCASSCF) Then
           Call mma_deallocate(G2x)
           Call mma_deallocate(Ash)
           Call mma_deallocate(DA)
        EndIf
      EndIf
************************************************************************
*                                                                      *
*     End of new Cho                                                   *
*                                                                      *
************************************************************************
*
      Call DaXpY_(ndens2,One,Int1,1,FockI,1)
      call dcopy_(ndens2,[0.0d0],0,Fock,1)
*
      Do iS=1,nSym
         If (nOrb(iS).eq.0) Go To 300
*
         ip_1 = ipCM(iS)
         If (nIsh(iS).gt.0)
     &      Call DYaX(nOrb(iS)*nIsh(is),2.0d0,
     &                FockI(ipCM(iS)),1,
     &                Fock (ipCM(iS)),1)
         If (iMethod.eq.iCASSCF) Then
            If (nIsh(iS).gt.0)
     &         Call DaXpY_(nOrb(iS)*nIsh(is),2.0d0,
     &                    FockA(ipCM(iS)),1,
     &                    Fock (ipCM(iS)),1)
            If (nAsh(iS).gt.0)
     &         Call DYaX(nOrb(iS)*nAsh(is),1.0d0,
     &                   Q(ipMatba(iS,is)),1,
     &                   Fock(ipCM(iS)+nIsh(is)*nOrb(is)),1)
            Do iAsh=1,nAsh(is)
               ipi=ipCM(iS)+nOrb(is)*(nIsh(is)+iAsh-1)
               Do jAsh=1,nAsh(is)
                  ipj=ipCM(iS)+nOrb(is)*(nIsh(is)+jAsh-1)
                  ni=nA(is)+iAsh
                  nj=nA(is)+jAsh
                  ipD=iTri(ni,nj)
                  call daxpy_(nOrb(is),G1t(ipD),
     &                       FockI(ipi),1,
     &                       Fock (ipj),1)
               End Do
            End Do
         End If
*
 300     Continue
      End Do
 101  Continue
      renergy=0.0d0
      rcora=0.0d0
      Do iS=1,nSym
      Do iB=1,nAsh(is)+nIsh(is)
      rEnergy=rEnergy+Fock(ipCM(is)+nOrb(iS)*(ib-1)+ib-1)
      End Do
      End Do
      rcorei=0.0d0
      rcorea=0.0d0
      rcor=0.0d0
      Do iS=1,nSym
       iptmp=ipCM(iS)
       Do iB=1,nIsh(is)
       rcorei=rcorei+2.0d0*Int1(iptmp)
       rcor=rcor+2.0d0*Focki(iptmp)
       iptmp=iptmp+nOrb(iS)+1
       End Do

       Do iB=1,nAsh(iS)
        Do jB=1,nAsh(iS)
         iiB=nA(iS)+ib
         ijB=nA(iS)+jb
         iij=iTri(iib,ijb)
         iiB=nIsh(iS)+ib
         ijB=nIsh(iS)+jb
         rcorea=rcorea+G1t(iij)*Int1(ipCM(is)-1+nOrb(is)*(iib-1)+ijB)

         rcora=rcora+G1t(iij)*Focki(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
        End Do
       End Do
      End Do
      rin_ene=0.5d0*(rcor+rcorei)
      rcore=rCorei+rcoreA
      If (debug) Then
         Write(6,*) 'Checking energy',0.5d0*renergy+potnuc+half*rcore
         Write(6,*) 'Checking energy',0.5d0*renergy,potnuc,half*rcore
         write(6,*)
      End if

      If(TwoStep.and.(StepType.eq.'RUN1')) Then
        iaddressQDAT=0
        Call ddafile(LuQDAT,1,FockA,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,FockI,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,Fock ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,Q    ,nDens2,iaddressQDAT)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
