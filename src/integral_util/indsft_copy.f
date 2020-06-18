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
* Copyright (C) 1990,1998, Roland Lindh                                *
*               1998, Martin Schuetz                                   *
************************************************************************
* integral sifting/copying complex module
* Martin G. Schuetz/ University of Stuttgart/ March' 98
*
      SubRoutine Integral_Copy(iCmp,iShell,MapOrg,
     &                         iBas,jBas,kBas,lBas,kOp,
     &                         Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                         AOInt,SOInt,nSOint,
     &                         iSOSym,nSkal,nSOs,
     &                         TInt,nTInt,FacInt,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         Dens,Fock,LDens,ExFac,NDens,
     &                         ind,nind,FckNoClmb,FckNoExch)
*     calls the proper routines IndSft_Copy/IndSft_Copy_Spc/PLF
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        iTOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4),
     &        nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        nShOffi(0:7), nShOffj(0:7), nShOffk(0:7), nShOffl(0:7)
      Logical Shijij,IJeqKL,FckNoClmb,FckNoExch
*
      External SO_BAddr_Inc_ijkl,PL_BAddr_Inc_ijkl
      Logical EqShls
*
      If (Petite) Then
        Call PLF_Copy(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &               iShell,MapOrg,iAO,iAOst,Shijij.and.IJeqKL,
     &               iBas,jBas,kBas,lBas,kOp,TInt,nTInt,FacInt,
     &               nShi(0),nShj(0),nShk(0),nShl(0),
     &               nShOffi(0),nShOffj(0),nShOffk(0),nShOffl(0),
     &               PL_BAddr_Inc_ijkl)
      Else
        EqShls=iShell(1).eq.iShell(2) .or. iShell(3).eq.iShell(4)
     &         .or. iShell(1).eq.iShell(3).and.iShell(2).eq.iShell(4)
     &         .or. iShell(1).eq.iShell(4).and.iShell(2).eq.iShell(3)
        If (EqShls) Then
          Call IndSft_Copy_Spc(iCmp,iShell,MapOrg,
     &                         iBas,jBas,kBas,lBas,Shijij,
     &                         iAO,iAOst,ijkl,
     &                         SOInt,nSOint,
     &                         iSOSym,nSkal,nSOs,
     &                         TInt,nTInt,FacInt,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         SO_BAddr_Inc_ijkl)
        Else
          Call IndSft_Copy(iCmp,iShell,MapOrg,
     &                     iBas,jBas,kBas,lBas,Shijij,
     &                     iAO,iAOst,ijkl,
     &                     SOInt,nSOint,
     &                     iSOSym,nSkal,nSOs,
     &                     TInt,nTInt,FacInt,iTOffs,nSym,
     &                     nShi,nShj,nShk,nShl,
     &                     nShOffi,nShOffj,nShOffk,nShOffl,
     &                     SO_BAddr_Inc_ijkl)
        End If
      End If
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(Dens)
         Call Unused_real(Fock)
         Call Unused_integer(LDens)
         Call Unused_real(ExFac)
         Call Unused_integer(NDens)
         Call Unused_integer(ind)
         Call Unused_integer(nind)
         Call Unused_logical(FckNoClmb)
         Call Unused_logical(FckNoExch)
      End If
      End


      SubRoutine Integral_Copy_jikl(iCmp,iShell,MapOrg,
     &                         iBas,jBas,kBas,lBas,kOp,
     &                         Shijij,IJeqKL,iAO,iAOst,ijkl,
     &                         AOInt,SOInt,nSOint,
     &                         iSOSym,nSkal,nSOs,
     &                         TInt,nTInt,FacInt,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl)
*     calls the proper routines IndSft_Copy/IndSft_Copy_Spc/PLF
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
*
      Real*8 AOInt(*), SOInt(*), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4),
     &        iAOst(4), kOp(4), iSOSym(2,nSOs),
     &        iTOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4),
     &        nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        nShOffi(0:7), nShOffj(0:7), nShOffk(0:7), nShOffl(0:7)
      Logical Shijij,IJeqKL
*
      External SO_BAddr_Inc_jikl,PL_BAddr_Inc_jikl
      Logical EqShls
*
      If (Petite) Then
        Call PLF_Copy(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &               iShell,MapOrg,iAO,iAOst,Shijij.and.IJeqKL,
     &               iBas,jBas,kBas,lBas,kOp,TInt,nTInt,FacInt,
     &               nShi(0),nShj(0),nShk(0),nShl(0),
     &               nShOffi(0),nShOffj(0),nShOffk(0),nShOffl(0),
     &               PL_BAddr_Inc_jikl)
      Else
        EqShls=iShell(1).eq.iShell(2) .or. iShell(3).eq.iShell(4)
     &         .or. iShell(1).eq.iShell(3).and.iShell(2).eq.iShell(4)
     &         .or. iShell(1).eq.iShell(4).and.iShell(2).eq.iShell(3)
        If (EqShls) Then
          Call IndSft_Copy_Spc(iCmp,iShell,MapOrg,
     &                         iBas,jBas,kBas,lBas,Shijij,
     &                         iAO,iAOst,ijkl,
     &                         SOInt,nSOint,
     &                         iSOSym,nSkal,nSOs,
     &                         TInt,nTInt,FacInt,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         SO_BAddr_Inc_jikl)
        Else
          Call IndSft_Copy(iCmp,iShell,MapOrg,
     &                     iBas,jBas,kBas,lBas,Shijij,
     &                     iAO,iAOst,ijkl,
     &                     SOInt,nSOint,
     &                     iSOSym,nSkal,nSOs,
     &                     TInt,nTInt,FacInt,iTOffs,nSym,
     &                     nShi,nShj,nShk,nShl,
     &                     nShOffi,nShOffj,nShOffk,nShOffl,
     &                     SO_BAddr_Inc_jikl)
        End If
      End If
      Return
      End


      SubRoutine IndSft_Copy(iCmp,iShell,MapOrg,
     &                       iBas,jBas,kBas,lBas,
     &                       Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                       iSOSym,nSkal,nSOs,
     &                       TInt,nTInt,FacInt,iTOffs,nSym,
     &                       nShi,nShj,nShk,nShl,
     &                       nShOffi,nShOffj,nShOffk,nShOffl,
     &                       SO_BAddr_Inc)
************************************************************************
*  object: to sift and index the SO integrals.                         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          april '90                                                   *
*          Hacked by M&R during end of January in Stuttgart 1998.      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
*
      External SO_BAddr_Inc
      Real*8 SOint(ijkl,nSOint), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4), jj(4),
     &        iAOst(4), iSOSym(2,nSOs),
     &        iTOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4),
     &        nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        nShOffi(0:7), nShOffj(0:7), nShOffk(0:7), nShOffl(0:7)
      Logical Shijij
*     local array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwo(0:7)
      Data iTwo/1,2,4,8,16,32,64,128/
*
*     Call qEnter('IndSft_Copy')
      irout = 39
      iprint = nprint(irout)
      If (iPrint.ge.49) Then
         r1=DDot_(ijkl*nSOInt,SOInt,1,[One],0)
         r2=DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
         Write (6,*) ' Sum=',r1
         Write (6,*) ' Dot=',r2
      End If
      If (iprint.ge.99)
     &   Call RecPrt(' In Indsft_Copy:SOint ',' ',SOint,ijkl,nSOint)
      memSO2 = 0
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      Do 100 i1 = 1, iCmp(1)
         niSym=0
         Do 101 j = 0, nIrrep-1
            If (iand(IrrCmp(inds(iShell(1))+i1),2**j).ne.0) Then
               iSym(niSym) = j
               niSym=niSym+1
            End If
101      Continue
         Do 200 i2 = 1, iCmp(2)
            njSym=0
            Do 201 j = 0, nIrrep-1
               If (iand(IrrCmp(inds(iShell(2))+i2),2**j).ne.0) Then
                  jSym(njSym) = j
                  njSym=njSym+1
               End If
201         Continue
            Do 300 i3 = 1, iCmp(3)
               nkSym=0
               Do 301 j = 0, nIrrep-1
                  If (iand(IrrCmp(inds(iShell(3))+i3),2**j).ne.0) Then
                     kSym(nkSym) = j
                     nkSym=nkSym+1
                  End If
301            Continue
               Do 400 i4 = 1, iCmp(4)
                  Do 401 j = 0, nIrrep-1
                     lSym(j)=iand(IrrCmp(inds(iShell(4))+i4),iTwo(j))
401               Continue
*
*      loop over Irreps which are spanned by the basis function.
*      again, the loop structure is restricted to ensure unique
*      integrals.
*
       Do 110 is = 0, niSym-1
          j1 = iSym(is)
          jj(1) = j1
          j2max = nIrrep-1
          Do 210 js = 0, njSym-1
             j2 = jSym(js)
             jj(2) = j2
             j12 = ieor(j1,j2)
             Do 310 ks = 0, nkSym-1
                j3 = kSym(ks)
                j4 = ieor(j12,j3)
                If (lSym(j4).eq.0) go to 310
                jj(3) = j3
                jj(4) = j4
*
                memSO2 = memSO2 + 1
*
*               Compute offset to Symmetry Block
                Call SO_BAddr_Inc(jj,MapOrg,iTOffs,nSym,
     &                        nShi,nShj,nShk,nShl,
     &                        nShOffi,nShOffj,nShOffk,nShOffl,
     &                        Inci,Incj,Inck,Incl,IOff4)
*
*               Compute absolute starting SO index
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1                  !   in jj(4)
                  iOff4l=iOff4+lSOl*Incl
                  Do kSOk = kSO, kSO+kBas-1               !   in jj(3)
                    iOff4kl=iOff4l+kSOk*Inck
                    Do jSOj = jSO, jSO+jBas-1            !   in jj(2)
                      iOff4jkl=iOff4kl+jSOj*Incj
                      Do iSOi = iSO, iSO+iBas-1         !   in jj(1)
                        nijkl = nijkl + 1
                        iOff4ijkl=iOff4jkl+iSOi*Inci
                        TInt(iOff4ijkl)=SOint(nijkl,memSO2)*FacInt
                      End Do
                    End Do
                  End Do
                End Do
*
310          Continue
210       Continue
110    Continue
*
400            Continue
300         Continue
200      Continue
100   Continue
*
*     Call qExit('IndSft_Copy')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Shijij)
         Call Unused_integer_array(iSOSym)
         Call Unused_integer(nSkal)
      End If
      End


      SubRoutine IndSft_Copy_SpC(iCmp,iShell,MapOrg,
     &                           iBas,jBas,kBas,lBas,
     &                           Shijij, iAO, iAOst, ijkl,SOint,nSOint,
     &                           iSOSym,nSkal,nSOs,
     &                           TInt,nTInt,FacInt,iTOffs,nSym,
     &                           nShi,nShj,nShk,nShl,
     &                           nShOffi,nShOffj,nShOffk,nShOffl,
     &                           SO_BAddr_Inc)
************************************************************************
*  object: to sift and index the SO integrals (special case for shell  *
*          degeneracy...)                                              *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*  Author: Martin G. Schuetz, Theochem, Uni Stuttgart, March 1998      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
*
      External SO_BAddr_Inc
      Real*8 SOint(ijkl,nSOint), TInt(nTInt)
      Integer iCmp(4), iShell(4), iAO(4), jj(4), jjj(4),
     &        iAOst(4), iSOSym(2,nSOs),
     &        iTOffs(0:nSym-1,0:nSym-1,0:nSym-1), MapOrg(4),
     &        nShi(0:7), nShj(0:7), nShk(0:7), nShl(0:7),
     &        nShOffi(0:7), nShOffj(0:7), nShOffk(0:7), nShOffl(0:7)
      Logical Shijij, Shij, Shkl, qijij, qij, qkl
*     local array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer iTwo(0:7)
      Data iTwo/1,2,4,8,16,32,64,128/
      Logical usShij, usShkl, usShik, usShjl, usShil, usShjk
*
*     Call qEnter('IndSft_Copy_SpC')
      irout = 39
      iprint = nprint(irout)
*
      Inci_2=0
      Incj_2=0
      Inck_2=0
      Incl_2=0
      IOff4_2=0
      Inci_3=0
      Incj_3=0
      Inck_3=0
      Incl_3=0
      IOff4_3=0
      Inci_2a=0
      Incj_2a=0
      Inck_2a=0
      Incl_2a=0
      IOff4_2a=0
      Inci_3a=0
      Incj_3a=0
      Inck_3a=0
      Incl_3a=0
      IOff4_3a=0
*
      k12=0
      k34=0
      If (iPrint.ge.49) Then
         r1=DDot_(ijkl*nSOInt,SOInt,1,[One],0)
         r2=DDot_(ijkl*nSOInt,SOInt,1,SOInt,1)
         Write (6,*) ' Sum=',r1
         Write (6,*) ' Dot=',r2
      End If
      If (iprint.ge.99)
     &   Call RecPrt(' In Indsft_Copy:SOint ',' ',SOint,ijkl,nSOint)
      memSO2 = 0
*
*     quadruple loop over elements of the basis functions angular
*     description. loops are reduced to just produce unique SO integrals
*     observe that we will walk through the memory in AOint in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
*     map shell indices to basis functions mapping...
*     and set up corresponding logicals...
      usShij = iShell(1).eq.iShell(2)
      usShkl = iShell(3).eq.iShell(4)
      usShik = iShell(1).eq.iShell(3)
      usShjl = iShell(2).eq.iShell(4)
      usShil = iShell(1).eq.iShell(4)
      usShjk = iShell(2).eq.iShell(3)
*
      Do 100 i1 = 1, iCmp(1)
         Do 101 j = 0, nIrrep-1
            iSym(j) = iand(IrrCmp(inds(iShell(1))+i1),iTwo(j))
101      Continue
         jCmpMx = iCmp(2)
         If (Shij) jCmpMx = i1
         Do 200 i2 = 1, jCmpMx
            Do 201 j = 0, nIrrep-1
               jSym(j) = iand(IrrCmp(inds(iShell(2))+i2),iTwo(j))
201         Continue
            qij = i1.eq.i2
            If (iShell(2).gt.iShell(1)) then
               i12 = iCmp(2)*(i1-1) + i2
            else
               i12 = iCmp(1)*(i2-1) + i1
            End If
            Do 300 i3 = 1, iCmp(3)
               Do 301 j = 0, nIrrep-1
                  kSym(j)=iand(IrrCmp(inds(iShell(3))+i3),iTwo(j))
301            Continue
               lCmpMx = iCmp(4)
               If (Shkl) lCmpMx = i3
               Do 400 i4 = 1, lCmpMx
                  Do 401 j = 0, nIrrep-1
                     lSym(j)=iand(IrrCmp(inds(iShell(4))+i4),2**j)
401               Continue
                  qkl = i3.eq.i4
                  If (iShell(4).gt.iShell(3)) then
                     i34 = iCmp(4)*(i3-1) + i4
                  else
                     i34 = iCmp(3)*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) go to 400
                  qijij = Shijij .and. i12.eq.i34
*
*      loop over Irreps which are spanned by the basis function.
*      again, the loop structure is restricted to ensure unique
*      integrals.
*
       Do 110 j1 = 0, nIrrep-1
          If (iSym(j1).eq.0) Go To 110
          jj(1) = j1
          j2max = nIrrep-1
          If (Shij .and. qij) j2max = j1
          Do 210 j2 = 0, j2max
             If (jSym(j2).eq.0) Go To 210
             jj(2) = j2
             j12 = ieor(j1,j2)
             If (qijij) then
                If (Shij .and. qij) then
                    k12 = j1*(j1+1)/2 + j2+1
                else If (Shij) then
                    k12 = nIrrep*j1 + j2+1
                else If (iShell(1).gt.iShell(2)) then
                    k12 = nIrrep*j1 + j2+1
                else
                    k12 = nIrrep*j2 + j1+1
                End If
             End If
*
             Do 310 j3 = 0, nIrrep-1
                If (kSym(j3).eq.0) Go To 310
                j4 = ieor(j12,j3)
                If (lSym(j4).eq.0) go to 310
                If (Shkl .and. qkl .and. j4.gt.j3) go to 310
                If (qijij) then
                   If (Shkl .and. qkl) then
                      k34 = j3*(j3+1)/2 + j4+1
                   else If (Shkl) then
                      k34 = nIrrep*j3 + j4+1
                   else If (iShell(3).gt.iShell(4)) then
                      k34 = nIrrep*j3 + j4+1
                   else
                      k34 = nIrrep*j4 + j3+1
                   End If
                   If (k34.gt.k12) go to 310
                End If
                jj(3) = j3
                jj(4) = j4
*
*               Compute offset to Symmetry Block (ij|kl)
                Call SO_BAddr_Inc(jj,MapOrg,iTOffs,nSym,
     &                            nShi,nShj,nShk,nShl,
     &                            nShOffi,nShOffj,nShOffk,nShOffl,
     &                            Inci,Incj,Inck,Incl,IOff4)
                If (Shijij) Then
*                 Compute offset to Symmetry Block (kl|ij)
                  jjj(1)=jj(3)
                  jjj(2)=jj(4)
                  jjj(3)=jj(1)
                  jjj(4)=jj(2)
                  Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         Inci_2,Incj_2,Inck_2,Incl_2,IOff4_2)
*                 Compute offset to Symmetry Block (lk|ji)
                  jjj(1)=jj(4)
                  jjj(2)=jj(3)
                  jjj(3)=jj(2)
                  jjj(4)=jj(1)
                  Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         Inci_3,Incj_3,Inck_3,Incl_3,IOff4_3)
*                 Compute offset to Symmetry Block (kl|ji)
                  jjj(1)=jj(3)
                  jjj(2)=jj(4)
                  jjj(3)=jj(2)
                  jjj(4)=jj(1)
                  Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         Inci_2a,Incj_2a,Inck_2a,Incl_2a,IOff4_2a)
*                 Compute offset to Symmetry Block (lk|ij)
                  jjj(1)=jj(4)
                  jjj(2)=jj(3)
                  jjj(3)=jj(1)
                  jjj(4)=jj(2)
                  Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                         nShi,nShj,nShk,nShl,
     &                         nShOffi,nShOffj,nShOffk,nShOffl,
     &                         Inci_3a,Incj_3a,Inck_3a,Incl_3a,IOff4_3a)
                End If
*               Compute offset to Symmetry Block (ji|kl)
                jjj(1)=jj(2)
                jjj(2)=jj(1)
                jjj(3)=jj(3)
                jjj(4)=jj(4)
                Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                       nShi,nShj,nShk,nShl,
     &                       nShOffi,nShOffj,nShOffk,nShOffl,
     &                       Inci_4,Incj_4,Inck_4,Incl_4,IOff4_4)
*               Compute offset to Symmetry Block (ji|lk)
                jjj(1)=jj(2)
                jjj(2)=jj(1)
                jjj(3)=jj(4)
                jjj(4)=jj(3)
                Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                       nShi,nShj,nShk,nShl,
     &                       nShOffi,nShOffj,nShOffk,nShOffl,
     &                       Inci_4a,Incj_4a,Inck_4a,Incl_4a,IOff4_4a)
*               Compute offset to Symmetry Block (ij|lk)
                jjj(1)=jj(1)
                jjj(2)=jj(2)
                jjj(3)=jj(4)
                jjj(4)=jj(3)
                Call SO_BAddr_Inc(jjj,MapOrg,iTOffs,nSym,
     &                       nShi,nShj,nShk,nShl,
     &                       nShOffi,nShOffj,nShOffk,nShOffl,
     &                       Inci_5,Incj_5,Inck_5,Incl_5,IOff4_5)
*
*               Compute absolute starting SO index
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
*
                memSO2 = memSO2 + 1
                nijkl = 0
                Do lSOl = lSO, lSO+lBas-1                  !   in jj(4)
                  iOff4l=iOff4+lSOl*Incl
                  iOff4l_2=IOff4_2+lSOl*Incj_2
                  iOff4l_3=IOff4_3+lSOl*Inci_3
                  iOff4l_2a=IOff4_2a+lSOl*Incj_2a
                  iOff4l_3a=IOff4_3a+lSOl*Inci_3a
                  iOff4l_4=iOff4_4+lSOl*Incl_4
                  iOff4l_4a=iOff4_4a+lSOl*Inck_4a
                  iOff4l_5=iOff4_5+lSOl*Inck_5
                  Do kSOk = kSO, kSO+kBas-1               !   in jj(3)
                    iOff4kl=iOff4l+kSOk*Inck
                    iOff4kl_2=IOff4l_2+kSOk*Inci_2
                    iOff4kl_3=IOff4l_3+kSOk*Incj_3
                    iOff4kl_2a=IOff4l_2a+kSOk*Inci_2a
                    iOff4kl_3a=IOff4l_3a+kSOk*Incj_3a
                    iOff4kl_4=iOff4l_4+kSOk*Inck_4
                    iOff4kl_4a=iOff4l_4a+kSOk*Incl_4a
                    iOff4kl_5=iOff4l_5+kSOk*Incl_5
                    Do jSOj = jSO, jSO+jBas-1            !   in jj(2)
                      iOff4jkl=iOff4kl+jSOj*Incj
                      iOff4jkl_2=iOff4kl_2+jSOj*Incl_2
                      iOff4jkl_3=iOff4kl_3+jSOj*Inck_3
                      iOff4jkl_2a=iOff4kl_2a+jSOj*Inck_2a
                      iOff4jkl_3a=iOff4kl_3a+jSOj*Incl_3a
                      iOff4jkl_4=iOff4kl_4+jSOj*Inci_4
                      iOff4jkl_4a=iOff4kl_4a+jSOj*Inci_4a
                      iOff4jkl_5=iOff4kl_5+jSOj*Incj_5
                      Do iSOi = iSO, iSO+iBas-1         !   in jj(1)
                        nijkl = nijkl + 1
                        iOff4ijkl=iOff4jkl+iSOi*Inci
                        iOff4ijkl_2=iOff4jkl_2+iSOi*Inck_2
                        iOff4ijkl_3=iOff4jkl_3+iSOi*Incl_3
                        iOff4ijkl_2a=iOff4jkl_2a+iSOi*Incl_2a
                        iOff4ijkl_3a=iOff4jkl_3a+iSOi*Inck_3a
                        iOff4ijkl_4=iOff4jkl_4+iSOi*Incj_4
                        iOff4ijkl_4a=iOff4jkl_4a+iSOi*Incj_4a
                        iOff4ijkl_5=iOff4jkl_5+iSOi*Inci_5
                        Val=SOint(nijkl,memSO2)*FacInt
                        TInt(iOff4ijkl)=Val
*                       duplicate, since there is shell degeneracy...
                        If (Shijij) Then
                          If (usShik .and. usShjl) Then
*                           (ij|kl) -> (kl|ij)
                            TInt(iOff4ijkl_2)=Val
                            If (usShij)
*                             (ij|kl) -> (kl|ji)
     &                        TInt(iOff4ijkl_2a)=Val
                          End If
                          If (usShil .and. usShjk) Then
*                           (ij|kl) -> (lk|ji)
                            TInt(iOff4ijkl_3)=Val
                            If (usShij)
*                             (ij|kl) -> (lk|ij)
     &                        TInt(iOff4ijkl_3a)=Val
                          End If
                        End If
                        If (usShij) Then
*                         (ij|kl) -> (ji|kl)
                          TInt(iOff4ijkl_4)=Val
                          If (usShkl)
*                           (ij|kl) -> (ji|lk)
     &                      TInt(iOff4ijkl_4a)=Val
                        End If
                        If (usShkl)
*                         (ij|kl) -> (ij|lk)
     &                    TInt(iOff4ijkl_5)=Val
                      End Do
                    End Do
                  End Do
                End Do
*
310          Continue
210       Continue
110    Continue
*
400            Continue
300         Continue
200      Continue
100   Continue
*
*     Call qExit('IndSft_Copy_SpC')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iSOSym)
         Call Unused_integer(nSkal)
      End If
      End


      SubRoutine SO_BAddr_Inc_ijkl(iSym_scr,MapOrg,iTOffs,nSym,
     &                             nShi,nShj,nShk,nShl,
     &                             nShOffi,nShOffj,nShOffk,nShOffl,
     &                             Inci,Incj,Inck,Incl,iBAddr)
*     compute addresses for integrals in order ijkl...
      Implicit Real*8 (A-H,O-Z)
*     declaration of subr parameters...
      Integer iSym_scr(4),MapOrg(4),iTOffs(0:nSym-1,0:nSym-1,0:nSym-1),
     &        nShi(0:7),nShj(0:7),nShk(0:7),nShl(0:7),
     &        nShOffi(0:7),nShOffj(0:7),nShOffk(0:7),nShOffl(0:7),
     &        Inci,Incj,Inck,Incl,iBAddr
*     declaration of local variables...
      Integer Inc(4),IncP(4),IrrOrg(4)
*
      Do i = 1, 4
        IrrOrg(MapOrg(i))=iSym_scr(i)
      End Do
      iSy=IrrOrg(1)
      jSy=IrrOrg(2)
      kSy=IrrOrg(3)
      lSy=IrrOrg(4)
      iBAddr=iTOffs(ksy,jSy,iSy)+1
      Inc(4) = 1
      Inc(3) = nShl(lSy)
      Inc(2) = nShl(lSy)*nShk(kSy)
      Inc(1) = nShl(lSy)*nShk(kSy)*nShj(jSy)
      iBAddr = iBAddr - Inc(4)*(nShOffl(lSy)+1)
     &                - Inc(3)*(nShOffk(kSy)+1)
     &                - Inc(2)*(nShOffj(jSy)+1)
     &                - Inc(1)*(nShOffi(iSy)+1)
      Do i = 1, 4
        IncP(i) = Inc(MapOrg(i))
      End Do
      Inci = IncP(1)
      Incj = IncP(2)
      Inck = IncP(3)
      Incl = IncP(4)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nShi)
      End


      SubRoutine SO_BAddr_Inc_jikl(iSym_scr,MapOrg,iTOffs,nSym,
     &                             nShi,nShj,nShk,nShl,
     &                             nShOffi,nShOffj,nShOffk,nShOffl,
     &                             Inci,Incj,Inck,Incl,iBAddr)
*     compute addresses for integrals in order jikl...
      Implicit Real*8 (A-H,O-Z)
*     declaration of subr parameters...
      Integer iSym_scr(4),MapOrg(4),iTOffs(0:nSym-1,0:nSym-1,0:nSym-1),
     &        nShi(0:7),nShj(0:7),nShk(0:7),nShl(0:7),
     &        nShOffi(0:7),nShOffj(0:7),nShOffk(0:7),nShOffl(0:7),
     &        Inci,Incj,Inck,Incl,iBAddr
*     declaration of local variables...
      Integer Inc(4),IncP(4),IrrOrg(4)
*
      Do i = 1, 4
        IrrOrg(MapOrg(i))=iSym_scr(i)
      End Do
      iSy=IrrOrg(1)
      jSy=IrrOrg(2)
      kSy=IrrOrg(3)
      lSy=IrrOrg(4)
      iBAddr=iTOffs(ksy,jSy,iSy)+1
      Inc(4) = 1
      Inc(3) = nShl(lSy)
      Inc(2) = nShl(lSy)*nShk(kSy)*nShi(iSy)
      Inc(1) = nShl(lSy)*nShk(kSy)
      iBAddr = iBAddr - Inc(4)*(nShOffl(lSy)+1)
     &                - Inc(3)*(nShOffk(kSy)+1)
     &                - Inc(2)*(nShOffj(jSy)+1)
     &                - Inc(1)*(nShOffi(iSy)+1)
      Do i = 1, 4
        IncP(i) = Inc(MapOrg(i))
      End Do
      Inci = IncP(1)
      Incj = IncP(2)
      Inck = IncP(3)
      Incl = IncP(4)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(nShj)
      End


      SubRoutine CpTpNdShlB(iS,jS,kS,lS,nShi,nShj,nShk,nShl,
     &                      nShOffi,nShOffj,nShOffk,nShOffl,TInt,
     &                      iTOffs,nSym,IntOrd_jikl)
************************************************************************
*  object: copy/transpose non diagonal shell block in "diagonal"       *
*          merged shell block                                          *
*          if IntOrd_jikl==.TRUE. integral order within symblk: jikl   *
*                           else  integral order within symblk: ijkl   *
*                                                                      *
*          Martin Schuetz, Theoretische Chemie Stuttgart               *
*          version march'98                                            *
************************************************************************
      use index_arrays, only: iShOff
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "shinf.fh"
*
*     declaration of subroutine parameters ...
      Integer iS,jS,kS,lS, nShi(0:7),nShj(0:7),nShk(0:7),nShl(0:7),
     &        nShOffi(0:7),nShOffj(0:7),nShOffk(0:7),nShOffl(0:7),
     &        iTOffs(0:nSym-1,0:nSym-1,0:nSym-1)
      Real*8 TInt(*)
      Logical IntOrd_jikl
*
      Do 100 Kirp = 0, nIrrep-1
        nShBFk=iWork(ipShBF+(kS-1)*nIrrep+Kirp)
        If (nShBFk.eq.0) Go To 100
        kSOb=iShOff(Kirp,kS)
        Do 110 Jirp = 0, nIrrep-1
          nShBFj=iWork(ipShBF+(jS-1)*nIrrep+Jirp)
          If (nShBFj.eq.0) Go To 110
          jSOb=iShOff(Jirp,jS)
          Do 120 Iirp = 0, nIrrep-1
            nShBFi=iWork(ipShBF+(iS-1)*nIrrep+Iirp)
            Lirp=iEor(iEor(Jirp,Iirp),Kirp)
            nShBFl=iWork(ipShBF+(lS-1)*nIrrep+Lirp)
            If (nShBFi*nShBFl.eq.0) Go To 120
            iSOb=iShOff(Iirp,iS)
            lSOb=iShOff(Lirp,lS)
            ipSO3G_n=iTOffs(Kirp,Jirp,Iirp)+1
            ipSO3G_t=iTOffs(Lirp,Jirp,Iirp)+1
            Incl_n = 1
            Inck_n = nShl(Lirp)
            Inck_t = 1
            Incl_t = nShl(Kirp)
            If (IntOrd_jikl) Then
              Incj   = nShl(Lirp)*nShk(Kirp)*nShi(Iirp)
              Inci   = nShl(Lirp)*nShk(Kirp)
            Else
              Incj   = nShl(Lirp)*nShk(Kirp)
              Inci   = nShl(Lirp)*nShk(Kirp)*nShj(Jirp)
            End If
            IOff4_n = ipSO3G_n - Incl_n*(nShOffl(Lirp)+1)
     &                         - Inck_n*(nShOffk(Kirp)+1)
     &                         - Incj*(nShOffj(Jirp)+1)
     &                         - Inci*(nShOffi(Iirp)+1)
            IOff4_t = ipSO3G_t - Inck_t*(nShOffk(Kirp)+1)
     &                         - Incl_t*(nShOffl(Lirp)+1)
     &                         - Incj*(nShOffj(Jirp)+1)
     &                         - Inci*(nShOffi(Iirp)+1)
            nklf_Grp=nShk(Kirp)*nShl(Lirp)
            Do iSO = iSOb, iSOb+nShBFi-1
              IOff4_n_i=IOff4_n+iSO*Inci
              IOff4_t_i=IOff4_t+iSO*Inci
              Do jSO = jSOb, jSOb+nShBFj-1
                IOff4_n_ij=IOff4_n_i+jSO*Incj
                IOff4_t_ij=IOff4_t_i+jSO*Incj
                Do kSO = kSOb, kSOb+nShBFk-1
                  IOff4_n_ijk=IOff4_n_ij+kSO*Inck_n
                  IOff4_t_ijk=IOff4_t_ij+kSO*Inck_t
                  Do lSO = lSOb, lSOb+nShBFl-1
                    IOff4_n_ijkl=IOff4_n_ijk+lSO*Incl_n
                    IOff4_t_ijkl=IOff4_t_ijk+lSO*Incl_t
                    TInt(IOff4_t_ijkl)=TInt(IOff4_n_ijkl)
                  End Do
                End Do
              End Do
            End Do
  120     Continue
  110   Continue
  100 Continue
      Return
      End
