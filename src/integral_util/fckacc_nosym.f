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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      Subroutine FckAcc_NoSymq(iAng, iCmp, jCmp, kCmp, lCmp, Shijij,
     &                  iShll, iShell, nijkl,
     &                  AOInt,FMat,DMat,nDens,
     &                  iAO,iAOst,iBas,jBas,kBas,lBas,
     &                  DoCoul,DoExch,Dij,Dkl,Dik,Dil,Djk,Djl,ExFac)
************************************************************************
*                                                                      *
*  Object: to accumulate contibutions from the AO integrals directly   *
*          to the symmetry adapted Fock matrix.                        *
*                                                                      *
*          This version uses square density and fock matrices          *
*          Modifications by HJW, 28.12.98                              *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          In addition to this complication we have that the order of  *
*          indices in the integrals are not ordered canonically but    *
*          rather in an order such that the contraction step will be   *
*          optimal. Hence, special care has to be taken when tracing   *
*          the density with the integrals so that both entities have   *
*          the same order.                                             *
*                                                                      *
*          The Fock matrix is computed in lower triangular form.       *
*                                                                      *
*          The density matrix is not folded if the shell indices and   *
*          the angular indices are identical.                          *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. February '93                            *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Real_Info, only: ThrInt
      Implicit Real*8 (A-H,O-Z)
      Intrinsic Max, Min
#include "real.fh"
#include "print.fh"
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), FMat(nDens),
     &       DMat(nDens)
      Logical Shij, Shkl, Shijij, DoCoul, DoExch
      Integer iAng(4), iShell(4), iShll(4),
     &        iAO(4), iAOst(4), nCmpx(4), nBasx(4)
*
      Parameter (nCBMax=200)
      Integer Indx(3,nCBMax,4)
*
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.DoExch.and..Not.DoCoul) Return
      iRout = 38
      iPrint = nPrint(iRout)
*
#ifdef _DEBUG_
*     Write (*,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),
*    &            DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One  ,0)
      If (iPrint.ge.49) Then
         Write (6,*) iCmp,jCmp,kCmp,lCmp
         Write (6,*) 'iAO=',iAO
         Write (6,*) 'iAOst=',iAOst
         Write (6,*) 'iShell=',iShell
         Write (6,*) 'iShll=',iShll
         Write (6,*) 'iAng=',iAng
         Write (6,*) DoCoul,DoExch,Shijij
         Write (6,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,[One],0),
     &                            DDot_(nDens,DMat,1,[One],0)
         Write (6,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
     &               AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
     &               AOInt,1,[One],0)
      End If
#endif
*
      DMax=Max(Dij,Dkl,Dik,Dil,Djk,Djl)
      If (dmax.gt.Zero) Then
        thr=ThrInt/DMax
      Else
        Return
      End If
c
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      ntg=nbas(0)
      Fac=One
      If (Shij) Fac=Fac*Half
      If (Shkl) Fac=Fac*Half
      If (Shijij) Fac=Fac*Half
      Fac_C=Four*Fac
      Fac_E=-Fac*Exfac
C     If (nijkl.ne.ibas*jbas*kbas*lbas) Call SysHalt( 'fckacc_nosym' )
*
      nCmpx(1)=iCmp
      nCmpx(2)=jCmp
      nCmpx(3)=kCmp
      nCmpx(4)=lCmp
      nBasx(1)=iBas
      nBasx(2)=jBas
      nBasx(3)=kBas
      nBasx(4)=lBas
      nCB_Max=Max(iCmp*iBas,jCmp*jBas,kCmp*kBas,lCmp*lBas)
      If (nCB_Max.gt.nCBMax) Then
         Call WarningMessage(2,'FckAcc_NoSym: nCB_Max.gt.nCBMax')
         Write (6,*) 'nCB_Max=',nCB_Max
         Call Abend()
      End If
      Do ii = 1, 4
         iCB = 0
         Do iC = 1, nCmpx(ii)
            iSO=iAOtSO(iAO(ii)+iC,0)+iAOSt(ii)
            Do i = 1, nBasx(ii)
               iCB = iCB + 1
               Indx(1,iCB,ii)=iC
               Indx(2,iCB,ii)=i
               Indx(3,iCB,ii)=iSO
               iSO = iSO + 1
            End Do
         End Do
      End Do
*
      nkl=kBas*lBas
      nij=iBas*jBas
*
      If (.Not.DoExch) Then
         Do lCB = 1, lCmp*lBas
            i4 =Indx(1,lCB,4)
            l  =Indx(2,lCB,4)
            lSO=Indx(3,lCB,4)
            ll=(lSO-1)*ntg
*
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
               kl=ll+kSO
*
               kl_=(l-1)*kBas+k
*
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jj=(jSO-1)*ntg
*
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ij=jj+iSO
*
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = Fac_C*AOInt(ijkl,i1,i2,i3,i4)
*
                     If (Abs(AOijkl).gt.Thr) Then
*
                        FMat(ij) = FMat(ij)+AOijkl*DMat(kl)
                        FMat(kl) = FMat(kl)+AOijkl*DMat(ij)
*
                     End If
                  End Do
               End Do
            End Do
         End Do
      Else If (.Not.DoCoul) Then
         Do lCB = 1, lCmp*lBas
            i4 =Indx(1,lCB,4)
            l  =Indx(2,lCB,4)
            lSO=Indx(3,lCB,4)
            ll=(lSO-1)*ntg
*
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
*
               kl_=(l-1)*kBas+k
*
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jk=kk+jSO
                  jl=ll+jSO
*
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ik=kk+iSO
                     il=ll+iSO
*
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = Fac_E*AOInt(ijkl,i1,i2,i3,i4)
*
                     If (Abs(AOijkl).gt.Thr) Then
*
                        FMat(ik) = FMat(ik)+AOijkl*DMat(jl)
                        FMat(il) = FMat(il)+AOijkl*DMat(jk)
                        FMat(jl) = FMat(jl)+AOijkl*DMat(ik)
                        FMat(jk) = FMat(jk)+AOijkl*DMat(il)
*
                     End If
                  End Do
               End Do
            End Do
         End Do
      Else
         Do lCB = 1, lCmp*lBas
            i4 =Indx(1,lCB,4)
            l  =Indx(2,lCB,4)
            lSO=Indx(3,lCB,4)
            ll=(lSO-1)*ntg
*
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
               kl=ll+kSO
*
               kl_=(l-1)*kBas+k
*
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jj=(jSO-1)*ntg
                  jk=kk+jSO
                  jl=ll+jSO
*
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ik=kk+iSO
                     il=ll+iSO
                     ij=jj+iSO
*
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = AOInt(ijkl,i1,i2,i3,i4)
                     If (Abs(AOijkl).gt.Thr) Then
*
                        FMat(ij) = FMat(ij)+Fac_C*AOijkl*DMat(kl)
                        FMat(kl) = FMat(kl)+Fac_C*AOijkl*DMat(ij)
                        FMat(ik) = FMat(ik)+Fac_E*AOijkl*DMat(jl)
                        FMat(il) = FMat(il)+Fac_E*AOijkl*DMat(jk)
                        FMat(jl) = FMat(jl)+Fac_E*AOijkl*DMat(ik)
                        FMat(jk) = FMat(jk)+Fac_E*AOijkl*DMat(il)
*
                     End If
                  End Do
               End Do
            End Do
         End Do
      End If
*
C     If (iPrint.ge.99) Then
C        Write (*,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,One,0),
C    &                            DDot_(nDens,DMat,1,One,0)
c     End If
*
*     Call GetMem(' Exit FckAcc','CHECK','REAL',iDum,iDum)
      Return
#ifndef _DEBUG_
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iAng)
         Call Unused_integer_array(iShll)
      End If
#endif
      End

      Subroutine FckAcc_NoSym(iAng, iCmp, jCmp, kCmp, lCmp, Shijij,
     &                  iShll, iShell, nijkl,
     &                  AOInt,FMat,DMat,nDens,
     &                  iAO,iAOst,iBas,jBas,kBas,lBas,ExFac)
************************************************************************
*                                                                      *
*  Object: to accumulate contibutions from the AO integrals directly   *
*          to the symmetry adapted Fock matrix.                        *
*                                                                      *
*          The indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*          In addition to this complication we have that the order of  *
*          indices in the integrals are not ordered canonically but    *
*          rather in an order such that the contraction step will be   *
*          optimal. Hence, special care has to be taken when tracing   *
*          the density with the integrals so that both entities have   *
*          the same order.                                             *
*                                                                      *
*          The Fock matrix is computed in lower triangular form.       *
*                                                                      *
*          The density matrix is not folded if the shell indices and   *
*          the angular indices are identical.                          *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden. February '93                            *
************************************************************************
      use SOAO_Info, only: iAOtSO
      use Real_Info, only: CutInt
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), FMat(nDens),
     &       DMat(nDens)
      Logical Shij, Shkl, Shijij
      Integer iAng(4), iShell(4), iShll(4),
     &        in(4), iAO(4), iAOst(4)
*     Local Arrays
*
*     Statement Function
*
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      iRout = 38
      iPrint = nPrint(iRout)
*
*     Write (*,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),
*    &            DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One  ,0)
*     If (iPrint.ge.49) Then
*        Write (*,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
*    &               AOInt,1,One,0)
*     End If
*
*     Quadruple loop over elements of the basis functions angular
*     description. Loops are reduced to just produce unique SO integrals
*     Observe that we will walk through the memory in AOInt in a
*     sequential way.
*
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
*
      Fac=One
      If (Shij) Fac=Fac*Half
      If (Shkl) Fac=Fac*Half
      If (Shijij) Fac=Fac*Half
      Fac_C=Four*Fac
      Fac_E=-Fac*ExFac
*
      Do 100 i1 = 1, iCmp
         in(1)=iAOtSO(iAO(1)+i1,0)+iAOst(1)
         Do 200 i2 = 1, jCmp
            in(2)=iAOtSO(iAO(2)+i2,0)+iAOst(2)
            Do 300 i3 = 1, kCmp
               in(3)=iAOtSO(iAO(3)+i3,0)+iAOst(3)
               Do 400 i4 = 1, lCmp
                  in(4)=iAOtSO(iAO(4)+i4,0)+iAOst(4)
*
*................ Get the starting SO indices.
*
                  iSO = in(1)
                  jSO = in(2)
                  kSO = in(3)
                  lSO = in(4)
*
                  nijkl=0
                  Do l = lSO, lSO+lBas-1
                     Do k = kSO, kSO+kBas-1
                        kl=iTri(k,l)
                        D_kl=DMat(kl)*Fac_C
                        F_kl=Zero
*
                        Do j = jSO, jSO+jBas-1
                           jk=iTri(j,k)
                           jl=iTri(j,l)
                           D_jl=DMat(jl)*Fac_E
                           D_jk=DMat(jk)*Fac_E
                           F_jl=Zero
                           F_jk=Zero
                           Do i = iSO, iSO+iBas-1
                              nijkl=nijkl+1
                              AOijkl = AOInt(nijkl,i1,i2,i3,i4)
                              If (Abs(AOijkl).lt.CutInt) Go To 10
                              ij=iTri(i,j)
                              ik=iTri(i,k)
                              il=iTri(i,l)
*
                              FMat(ij) = FMat(ij)+AOijkl*D_kl
                              F_kl     = F_kl    +AOijkl*DMat(ij)
*
                              FMat(ik) = FMat(ik)+AOijkl*D_jl
                              F_jl     = F_jl    +AOijkl*DMat(ik)
                              FMat(il) = FMat(il)+AOijkl*D_jk
                              F_jk     = F_jk    +AOijkl*DMat(il)
*
 10                           Continue
                           End Do
                           FMat(jl)=FMat(jl)+F_jl*Fac_E
                           FMat(jk)=FMat(jk)+F_jk*Fac_E
                        End Do
                        FMat(kl)=FMat(kl)+F_kl*Fac_C
                     End Do
                  End Do
*
 400           Continue
 300        Continue
 200     Continue
 100  Continue
*
*     Call GetMem(' Exit FckAcc','CHECK','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iAng)
         Call Unused_integer_array(iShll)
      End If
      End
