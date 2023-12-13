!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine FckAcc_NoSymq(iCmp, jCmp, kCmp, lCmp, Shijij,
     &                         iShell, nijkl,AOInt,FMat,DMat,nDens,
     &                         iAO,iAOst,iBas,jBas,kBas,lBas,
     &                  DoCoul,DoExch,Dij,Dkl,Dik,Dil,Djk,Djl,ExFac)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmetry adapted Fock matrix.                        *
!                                                                      *
!          This version uses square density and fock matrices          *
!          Modifications by HJW, 28.12.98                              *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!          In addition to this complication we have that the order of  *
!          indices in the integrals are not ordered canonically but    *
!          rather in an order such that the contraction step will be   *
!          optimal. Hence, special care has to be taken when tracing   *
!          the density with the integrals so that both entities have   *
!          the same order.                                             *
!                                                                      *
!          The Fock matrix is computed in lower triangular form.       *
!                                                                      *
!          The density matrix is not folded if the shell indices and   *
!          the angular indices are identical.                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. February '93                            *
!***********************************************************************
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Gateway_Info, only: ThrInt
      use Constants, only: Zero, One, Four, Half
      Implicit None
      Integer nijkl, iCmp, jCmp, kCmp, lCmp, nDens
      Integer iBas, jBas, kBas, lBas
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), FMat(nDens),
     &       DMat(nDens)
      Logical Shij, Shkl, Shijij, DoCoul, DoExch
      Integer iShell(4), iAO(4), iAOst(4)
      Real*8 Dij, Dkl, Dik, Dil, Djk, Djl, ExFac
!
      Intrinsic Max
      Integer, Parameter:: nCBMax=200
      Integer Indx(3,nCBMax,4)
      Integer nCmpx(4), nBasx(4)
      Integer ntg, ii, jj, kk, ll, i1, i2, i3, i4, i, j, k, l,
     &        ij, kl, ik, il, jk, jl, ij_, kl_, ijkl, ncb_Max, ic,
     &        iSO, jSO, kSO, lSO, icb, jcb, kcb, lcb, nij
      Real*8 DMax, Thr, Fac, Fac_C, Fac_E, AOijkl
#ifdef _DEBUGPRINT_
      Real*8, External :: DDot_
#endif
!
!                                                                      *
!***********************************************************************
!                                                                      *
      If (.Not.DoExch.and..Not.DoCoul) Return
#ifdef _DEBUGPRINT_
!     Write (*,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),
!    &            DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One  ,0)
      Write (6,*) iCmp,jCmp,kCmp,lCmp
      Write (6,*) 'iAO=',iAO
      Write (6,*) 'iAOst=',iAOst
      Write (6,*) 'iShell=',iShell
      Write (6,*) DoCoul,DoExch,Shijij
      Write (6,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,[One],0),
     &                         DDot_(nDens,DMat,1,[One],0)
      Write (6,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
     &            AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
     &            AOInt,1,[One],0)
#endif
!
      DMax=Max(Dij,Dkl,Dik,Dil,Djk,Djl)
      If (dmax.gt.Zero) Then
        thr=ThrInt/DMax
      Else
        Return
      End If
!
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
      ntg=nbas(0)
      Fac=One
      If (Shij) Fac=Fac*Half
      If (Shkl) Fac=Fac*Half
      If (Shijij) Fac=Fac*Half
      Fac_C=Four*Fac
      Fac_E=-Fac*Exfac
!     If (nijkl.ne.ibas*jbas*kbas*lbas) Call SysHalt( 'fckacc_nosym' )
!
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
!
      nij=iBas*jBas
!
      If (.Not.DoExch) Then
         Do lCB = 1, lCmp*lBas
            i4 =Indx(1,lCB,4)
            l  =Indx(2,lCB,4)
            lSO=Indx(3,lCB,4)
            ll=(lSO-1)*ntg
!
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
               kl=ll+kSO
!
               kl_=(l-1)*kBas+k
!
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jj=(jSO-1)*ntg
!
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ij=jj+iSO
!
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = Fac_C*AOInt(ijkl,i1,i2,i3,i4)
!
                     If (Abs(AOijkl).gt.Thr) Then
!
                        FMat(ij) = FMat(ij)+AOijkl*DMat(kl)
                        FMat(kl) = FMat(kl)+AOijkl*DMat(ij)
!
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
!
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
!
               kl_=(l-1)*kBas+k
!
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jk=kk+jSO
                  jl=ll+jSO
!
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ik=kk+iSO
                     il=ll+iSO
!
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = Fac_E*AOInt(ijkl,i1,i2,i3,i4)
!
                     If (Abs(AOijkl).gt.Thr) Then
!
                        FMat(ik) = FMat(ik)+AOijkl*DMat(jl)
                        FMat(il) = FMat(il)+AOijkl*DMat(jk)
                        FMat(jl) = FMat(jl)+AOijkl*DMat(ik)
                        FMat(jk) = FMat(jk)+AOijkl*DMat(il)
!
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
!
            Do kCB = 1, kCmp*kBas
               i3 =Indx(1,kCB,3)
               k  =Indx(2,kCB,3)
               kSO=Indx(3,kCB,3)
               kk=(kSO-1)*ntg
               kl=ll+kSO
!
               kl_=(l-1)*kBas+k
!
               Do jCB = 1, jCmp*jBas
                  i2 =Indx(1,jCB,2)
                  j  =Indx(2,jCB,2)
                  jSO=Indx(3,jCB,2)
                  jj=(jSO-1)*ntg
                  jk=kk+jSO
                  jl=ll+jSO
!
                  Do iCB = 1, iCmp*iBas
                     i1 =Indx(1,iCB,1)
                     i  =Indx(2,iCB,1)
                     iSO=Indx(3,iCB,1)
                     ik=kk+iSO
                     il=ll+iSO
                     ij=jj+iSO
!
                     ij_=(j-1)*iBas+i
                     ijkl=(kl_-1)*nij+ij_
                     AOijkl = AOInt(ijkl,i1,i2,i3,i4)
                     If (Abs(AOijkl).gt.Thr) Then
!
                        FMat(ij) = FMat(ij)+Fac_C*AOijkl*DMat(kl)
                        FMat(kl) = FMat(kl)+Fac_C*AOijkl*DMat(ij)
                        FMat(ik) = FMat(ik)+Fac_E*AOijkl*DMat(jl)
                        FMat(il) = FMat(il)+Fac_E*AOijkl*DMat(jk)
                        FMat(jl) = FMat(jl)+Fac_E*AOijkl*DMat(ik)
                        FMat(jk) = FMat(jk)+Fac_E*AOijkl*DMat(il)
!
                     End If
                  End Do
               End Do
            End Do
         End Do
      End If
!
!     Write (6,*) 'FMAT,DMAT=',DDot_(nDens,FMat,1,One,0),
!    &                         DDot_(nDens,DMat,1,One,0)
!
      Return
      End Subroutine FckAcc_NoSymq
