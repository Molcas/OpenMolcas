**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993,1998, Roland Lindh                                *
!***********************************************************************
      Subroutine FckAcc(iAng,iCmp_, Shijij,
     &                  iShll, iShell, kOp, nijkl,
     &                  AOInt,TwoHam,nDens,Scrt,nScrt,
     &                  iAO,iAOst,iBas,jBas,kBas,lBas,
     &                  Dij,ij1,ij2,ij3,ij4,
     &                  Dkl,kl1,kl2,kl3,kl4,
     &                  Dik,ik1,ik2,ik3,ik4,
     &                  Dil,il1,il2,il3,il4,
     &                  Djk,jk1,jk2,jk3,jk4,
     &                  Djl,jl1,jl2,jl3,jl4,
     &                  FT,nFT,DoCoul,DoExch,ExFac)
!***********************************************************************
!                                                                      *
!  Object: to accumulate contributions from the AO integrals directly  *
!          to the symmetry adapted Fock matrix.                        *
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
!                                                                      *
!     Modified July '98 in Tokyo by R. Lindh                           *
!***********************************************************************
      use Basis_Info, only: Shells
      use SOAO_Info, only: iAOtSO
      use Real_Spherical, only: iSphCr
      use Symmetry_Info, only: iChBas, iOper, nIrrep, Prmt
      use Gateway_Info, only: ThrInt, CutInt
      use Constants, only: Zero, One, Two, Half, Quart
      Implicit None
!
!     Since I sometimes use Scrt as an anchor to reach into the
!     density matrix a check if I'm out of bounds does not work!
!
      Integer nijkl, nDens, nScrt, nFT,
     &       ij1,ij2,ij3,ij4,
     &       kl1,kl2,kl3,kl4,
     &       ik1,ik2,ik3,ik4,
     &       il1,il2,il3,il4,
     &       jk1,jk2,jk3,jk4,
     &       jl1,jl2,jl3,jl4
      Integer iBas, jBas, kBas, lBas
      Real*8, target:: Scrt(nScrt)
      Integer iCmp_(4)
      Real*8 AOInt(nijkl,iCmp_(1),iCmp_(2),iCmp_(3),iCmp_(4)),
     &       TwoHam(nDens)
      Real*8, target::
     &       Dij(ij1*ij2+1,ij3,ij4),
     &       Dkl(kl1*kl2+1,kl3,kl4),
     &       Dik(ik1*ik2+1,ik3,ik4),
     &       Dil(il1*il2+1,il3,il4),
     &       Djk(jk1*jk2+1,jk3,jk4),
     &       Djl(jl1*jl2+1,jl3,jl4)
      Real*8, Target:: FT(nFT)
      Logical Shijij, Qijij, DoCoul, DoExch,
     &        iShij, iShkl, iQij, iQkl,
     &        iQik, iShik, iQil, iShil, iQjk, iShjk, iQjl, iShjl,
     &        lFij, lFkl, lFik, lFjl, lFil, lFjk
      Integer iAng(4), iShell(4), iShll(4), kOp(4), kOp2(4),
     &        iAO(4), iAOst(4), iCmpa(4)
      Real*8 ExFac

!     Local Arrays
      Integer mijkl, jCmpMx, lCmpMx, j, ix, ijkl, ip, iOpt, iIrrep
      Integer iChBs, jChBs, kChBs, lChBs
      Integer ixyz, iOff
      Integer iSym(4)
      Integer iCmp, jCmp, kCmp, lCmp
      Integer ii, jj, kk, ll
      Integer i1, i2, i3, i4, i12, i34
      Real*8  pEa, pRb, pTc, pTSd
      Real*8 Fac_ij, Fac_kl, Fac_ik, Fac_jl, Fac_il, Fac_jk
      Real*8 D_ij, D_kl, D_ik, D_jl, D_il, D_jk
      Real*8 Vij, Vkl, Vik, Vjl, Vil, Vjk, Vijkl
      Real*8 Fact, pFctr
      Real*8, pointer:: Fij(:,:,:), Fkl(:,:,:), Fik(:,:,:),
     &                  Fil(:,:,:), Fjk(:,:,:), Fjl(:,:,:)
      Real*8, pointer:: pDij(:), pDkl(:), pDik(:),
     &                  pDil(:), pDjk(:), pDjl(:)
      Integer nF, ipF
#ifdef _DEBUGPRINT_
      Real*8, External:: XDot
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement Function
!
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
!                                                                      *
!***********************************************************************
!                                                                      *
      If (.Not.DoCoul.and..Not.DoExch) Return

      iCmp=iCmp_(1)
      jCmp=iCmp_(2)
      kCmp=iCmp_(3)
      lCmp=iCmp_(4)
#ifdef _DEBUGPRINT_
      Call RecPrt('FckAcc:AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
      Write (6,*)'Dij=',XDot(Dij,ij1,ij2,ij3,ij4)
      Write (6,*)'Dkl=',XDot(Dkl,kl1,kl2,kl3,kl4)
      Write (6,*)'Dik=',XDot(Dik,ik1,ik2,ik3,ik4)
      Write (6,*)'Dil=',XDot(Dil,il1,il2,il3,il4)
      Write (6,*)'Djk=',XDot(Djk,jk1,jk2,jk3,jk4)
      Write (6,*)'Djl=',XDot(Djl,jl1,jl2,jl3,jl4)
      Call RecPrt('AOInt',' ',AOInt,nijkl,iCmp*jCmp*kCmp*lCmp)
#endif
!
      If (iBas*jBas*kBas*lBas.gt.nScrt) Then
         Call WarningMessage(2,'FckAcc: nScrt too small!')
         Call Abend()
      End If
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))

      kOp2(1) = iOper(kOp(1))
      kOp2(2) = iOper(kOp(2))
      kOp2(3) = iOper(kOp(3))
      kOp2(4) = iOper(kOp(4))
      iCmpa(1) = iCmp
      iCmpa(2) = jCmp
      iCmpa(3) = kCmp
      iCmpa(4) = lCmp
      lFij = .False.
      lFkl = .False.
      lFik = .False.
      lFjl = .False.
      lFil = .False.
      lFjk = .False.
!
      ipF   = 1
      nF    = iBas*jBas*iCmpa(1)*iCmpa(2)
      Fij(1:iBas*jBas,1:iCmpa(1),1:iCmpa(2)) => FT(ipF:ipF+nF-1)

!
      ipF   = ipF   + nF
      nF    = kBas*lBas*iCmpa(3)*iCmpa(4)
      Fkl(1:kBas*lBas,1:iCmpa(3),1:iCmpa(4)) => FT(ipF:ipF+nF-1)
!
      ipF   = ipF   + nF
      nF    = iBas*kBas*iCmpa(1)*iCmpa(3)
      Fik(1:iBas*kBas,1:iCmpa(1),1:iCmpa(3)) => FT(ipF:ipF+nF-1)
!
      ipF   = ipF   + nF
      nF    = jBas*lBas*iCmpa(2)*iCmpa(4)
      Fjl(1:jBas*lBas,1:iCmpa(2),1:iCmpa(4)) => FT(ipF:ipF+nF-1)
!
      ipF   = ipF   + nF
      nF    = iBas*lBas*iCmpa(1)*iCmpa(4)
      Fil(1:iBas*lBas,1:iCmpa(1),1:iCmpa(4)) => FT(ipF:ipF+nF-1)
!
      ipF   = ipF   + nF
      nF    = jBas*kBas*iCmpa(2)*iCmpa(3)
      Fjk(1:jBas*kBas,1:iCmpa(2),1:iCmpa(3)) => FT(ipF:ipF+nF-1)
!
      ipF   = ipF   + nF
      FT(1:ipF-1)=Zero
!
!     Quadruple loop over elements of the basis functions angular
!     description. Loops are reduced to just produce unique SO integrals
!     Observe that we will walk through the memory in AOInt in a
!     sequential way.
!
      iShij = iShell(1).eq.iShell(2)
      iShkl = iShell(3).eq.iShell(4)
      iShik = iShell(1).eq.iShell(3)
      iShil = iShell(1).eq.iShell(4)
      iShjk = iShell(2).eq.iShell(3)
      iShjl = iShell(2).eq.iShell(4)
      mijkl=iBas*jBas*kBas*lBas
      Do i1 = 1, iCmp
         ix = 0
         Do j = 0, nIrrep-1
            If (iAOtSO(iAO(1)+i1,j)>0) ix = iEor(ix,2**j)
         End Do
         iSym(1)=ix
         jCmpMx = jCmp
         If (iShij) jCmpMx = i1
         iChBs = iChBas(ii+i1)
         If (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
         pEa = Prmt(iOper(kOp(1)),iChBs)
         Do i2 = 1, jCmpMx
            ix = 0
            Do j = 0, nIrrep-1
               If (iAOtSO(iAO(2)+i2,j)>0) ix = iEor(ix,2**j)
            End Do
            iSym(2)=ix
            jChBs = iChBas(jj+i2)
            If (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
            pRb = Prmt(iOper(kOp(2)),jChBs)
            If (iShell(2).gt.iShell(1)) Then
               i12 = jCmp*(i1-1) + i2
            Else
               i12 = iCmp*(i2-1) + i1
            End If
            Do i3 = 1, kCmp
               ix = 0
               Do j = 0, nIrrep-1
                  If (iAOtSO(iAO(3)+i3,j)>0) ix = iEor(ix,2**j)
               End Do
               iSym(3)=ix
               lCmpMx = lCmp
               If (iShkl) lCmpMx = i3
               kChBs = iChBas(kk+i3)
               If (Shells(iShll(3))%Transf)
     &            kChBs = iChBas(iSphCr(kk+i3))
               pTc = Prmt(iOper(kOp(3)),kChBs)
               Do i4 = 1, lCmpMx
                  ix = 0
                  Do j = 0, nIrrep-1
                     If (iAOtSO(iAO(4)+i4,j)>0) ix = iEor(ix,2**j)
                  End Do
                  iSym(4)=ix
                  lChBs = iChBas(ll+i4)
                  If (Shells(iShll(4))%Transf)
     &               lChBs = iChBas(iSphCr(ll+i4))
                  pTSd= Prmt(iOper(kOp(4)),lChBs)
                  If (iShell(4).gt.iShell(3)) Then
                     i34 = lCmp*(i3-1) + i4
                  Else
                     i34 = kCmp*(i4-1) + i3
                  End If
                  If (Shijij .and. i34.gt.i12) Cycle
                  vijkl=Zero
                  do ijkl=1,mijkl
                    vijkl=max(vijkl,abs(AOInt(ijkl,i1,i2,i3,i4)))
                  end do
!                 vijkl = DNrm2_(iBas*jBas*kBas*lBas,
!    &                          AOInt(1,i1,i2,i3,i4),1)
                  If (vijkl.lt.CutInt) Cycle
!
                  Qijij = Shijij .and. i12.eq.i34
                  iQij = iShij.and.i1.eq.i2
                  iQkl = iShkl.and.i3.eq.i4
                  iQik = iShik.and.i1.eq.i3
                  iQil = iShil.and.i1.eq.i4
                  iQjk = iShjk.and.i2.eq.i3
                  iQjl = iShjl.and.i2.eq.i4
                  pFctr=pEa*pRb*pTc*pTSd
!***********************************************************************
!
                  Fac_ij = pFctr
                  Fac_kl = pFctr
                  Fac_ik =-Quart * pFctr
                  Fac_jl =-Quart * pFctr
                  Fac_il =-Quart * pFctr
                  Fac_jk =-Quart * pFctr
                  If (iQij) Fac_ij = Fac_ij * Half
                  If (iQkl) Fac_kl = Fac_kl * Half
                  If (iQjl.and. .Not.iQik) Fac_ik = Fac_ik * Two
                  If (iQik.and. .Not.iQjl) Fac_jl = Fac_jl * Two
                  If (iQjk.and. .Not.iQil) Fac_il = Fac_il * Two
                  If (iQil.and. .Not.iQjk) Fac_jk = Fac_jk * Two
!
                  D_ij=Two
                  If (iQij) D_ij=One
                  D_kl=Two
                  If (iQkl) D_kl=One
                  D_ik=Two
                  If (iQik) D_ik=One
                  D_jl=Two
                  If (iQjl) D_jl=One
                  D_il=Two
                  If (iQil) D_il=One
                  D_jk=Two
                  If (iQjk) D_jk=One
!
                  Fac_ij=Fac_ij*D_kl
                  Fac_kl=Fac_kl*D_ij
                  Fac_jl=Fac_jl*D_ik
                  Fac_ik=Fac_ik*D_jl
                  Fac_jk=Fac_jk*D_il
                  Fac_il=Fac_il*D_jk
!
                  If (Qijij) Then
                     Fac_kl = Zero
                     Fac_jk = Zero
                  End If
                  If (iQij.and.iQkl) Then
                     Fac_jl = Zero
                     Fac_il = Zero
                     Fac_jk = Zero
                  End If
                  If (iQij.or.iQkl) Then
                     Fac_il = Zero
                     Fac_jk = Zero
                  End If
!***********************************************************************
!
                  iOpt=0
                  ip = 1
                  If (DoCoul .and.
     &                iAnd(iSym(1),iSym(2)).ne.0 .and.
     &                iAnd(iSym(3),iSym(4)).ne.0) Then
                     iOpt = iOpt + 1
!
                     If (iShell(3).lt.iShell(4)) Then
                        vkl = Dkl(kl1*kl2+1,i4,i3)
                        pDkl(1:kl1*kl2) => Scrt(ip:ip+kl1*kl2-1)
                        Call DGetMO(Dkl(1,i4,i3),kl1,kl1,kl2,pDkl,kl2)
                        ip = ip + kl1*kl2
                     Else
                        vkl = Dkl(kl1*kl2+1,i3,i4)
                        pDkl(1:kl1*kl2) => Dkl(1:kl1*kl2,i3,i4)
                     End If
                     If (iShell(1).lt.iShell(2)) Then
                        vij = Dij(ij1*ij2+1,i2,i1)
                        pDij(1:ij1*ij2) => Scrt(ip:ip+ij1*ij2-1)
                        Call DGeTMO(Dij(1,i2,i1),ij1,ij1,ij2,pDij,ij2)
                        ip = ip + ij1*ij2
                     Else
                        vij = Dij(ij1*ij2+1,i1,i2)
                        pDij(1:ij1*ij2) => Dij(1:ij1*ij2,i1,i2)
                     End  If
                     If (vkl*vijkl*Abs(Fac_ij).lt.ThrInt .and.
     &                   vij*vijkl*Abs(Fac_kl).lt.ThrInt) Then
                        iOpt = iOpt -1
                     Else
                        lFij=.True.
                        lFkl=.True.
                     End If
                  End If
                  If (DoExch .and.
     &                iAnd(iSym(1),iSym(3)).ne.0 .and.
     &                iAnd(iSym(2),iSym(4)).ne.0) Then
                     iOpt = iOpt + 2
!
                     If (iShell(2).lt.iShell(4)) Then
                        vjl=Djl(jl1*jl2+1,i4,i2)
                        pDjl(1:jl1*jl2) => Scrt(ip:ip+jl1*jl2-1)
                        Call DGeTMO(Djl(1,i4,i2),jl1,jl1,jl2,pDjl,jl2)
                        ip = ip + jl1*jl2
                     Else
                        vjl=Djl(jl1*jl2+1,i2,i4)
                        pDjl(1:jl1*jl2) => Djl(1:jl1*jl2,i2,i4)
                     End If
                     If (iShell(1).lt.iShell(3)) Then
                        vik=Dik(ik1*ik2+1,i3,i1)
                        pDik(1:ik1*ik2) => Scrt(ip:ip+ik1*ik2-1)
                        Call DGeTMO(Dik(1,i3,i1),ik1,ik1,ik2,pDik,ik2)
                        ip = ip + ik1*ik2
                     Else
                        vik=Dik(ik1*ik2+1,i1,i3)
                        pDik(1:ik1*ik2) => Dik(1:ik1*ik2,i1,i3)
                     End If
                     If (vik*vijkl*Abs(Fac_jl) .lt. ThrInt .and.
     &                   vjl*vijkl*Abs(Fac_ik) .lt. ThrInt) Then
                        iOpt = iOpt - 2
                     Else
                        lFik = .True.
                        lFjl = .True.
                     End If
                  End If
                  If (DoExch .and.
     &                iAnd(iSym(1),iSym(4)).ne.0 .and.
     &                iAnd(iSym(2),iSym(3)).ne.0) Then
                     iOpt = iOpt + 4
!
                     If (iShell(2).lt.iShell(3)) Then
                        vjk = Djk(jk1*jk2+1,i3,i2)
                        pDjk(1:jk1*jk2) => Scrt(ip:ip+jk1*jk2-1)
                        Call DGeTMO(Djk(1,i3,i2),jk1,jk1,jk2,pDjk,jk2)
                        ip = ip + jk1*jk2
                     Else
                        vjk = Djk(jk1*jk2+1,i2,i3)
                        pDjk(1:jk1*jk2) => Djk(1:jk1*jk2,i2,i3)
                     End If
                     If (iShell(1).lt.iShell(4)) Then
                        vil = Dil(il1*il2+1,i4,i1)
                        pDil(1:il1*il2) => Scrt(ip:ip+il1*il2-1)
                        Call DGeTMO(Dil(1,i4,i1),il1,il1,il2,pDil,il2)
                        ip = ip + il1*il2
                     Else
                        vil = Dil(il1*il2+1,i1,i4)
                        pDil(1:il1*il2) => Dil(1:il1*il2,i1,i4)
                     End If
                     If (vil*vijkl*Abs(Fac_jk) .lt. ThrInt .and.
     &                   vjk*vijkl*Abs(Fac_il) .lt. ThrInt) Then
                        iOpt = iOpt -4
                     Else
                        lFil = .True.
                        lFjk = .True.
                     End If
                  End If
                  If (ip-1.gt.nScrt) Then
                     Call WarningMessage(2,'FckAcc: nScrt too small!')
                     Call Abend()
                  End If
                  Select  Case (iOpt)
!
!***********************************************************************
                  Case (1)
                  Call Fck1(AOInt(:,i1,i2,i3,i4),
     &                      pDij,Fij(:,i1,i2),Fac_ij,
     &                      pDkl,Fkl(:,i3,i4),Fac_kl)
                  Case (2)
                  Call Fck2(AOInt(:,i1,i2,i3,i4),
     &                      pDik,Fik(:,i1,i3),Fac_ik,
     &                      pDjl,Fjl(:,i2,i4),Fac_jl)
                  Case (3)
                  Call Fck3(AOInt(:,i1,i2,i3,i4),
     &                      pDij,Fij(:,i1,i2),Fac_ij,
     &                      pDkl,Fkl(:,i3,i4),Fac_kl,
     &                      pDik,Fik(:,i1,i3),Fac_ik,
     &                      pDjl,Fjl(:,i2,i4),Fac_jl)
                  Case (4)
                  Call Fck4(AOInt(:,i1,i2,i3,i4),
     &                      pDil,Fil(:,i1,i4),Fac_il,
     &                      pDjk,Fjk(:,i2,i3),Fac_jk)
                  Case (5)
                  Call Fck5(AOInt(:,i1,i2,i3,i4),
     &                      pDij,Fij(:,i1,i2),Fac_ij,
     &                      pDkl,Fkl(:,i3,i4),Fac_kl,
     &                      pDil,Fil(:,i1,i4),Fac_il,
     &                      pDjk,Fjk(:,i2,i3),Fac_jk)
                  Case (6)
                  Call Fck6(AOInt(:,i1,i2,i3,i4),
     &                      pDik,Fik(:,i1,i3),Fac_ik,
     &                      pDjl,Fjl(:,i2,i4),Fac_jl,
     &                      pDil,Fil(:,i1,i4),Fac_il,
     &                      pDjk,Fjk(:,i2,i3),Fac_jk)
                  Case (7)
                  Call Fck7(AOInt(:,i1,i2,i3,i4),
     &                      pDij,Fij(:,i1,i2),Fac_ij,
     &                      pDkl,Fkl(:,i3,i4),Fac_kl,
     &                      pDik,Fik(:,i1,i3),Fac_ik,
     &                      pDjl,Fjl(:,i2,i4),Fac_jl,
     &                      pDil,Fil(:,i1,i4),Fac_il,
     &                      pDjk,Fjk(:,i2,i3),Fac_jk)
                  Case Default
                  End Select
                  Nullify(pDij,pDkl,pDik,pDil,pDjk,pDjl)
!***********************************************************************
!
               End Do
            End Do
         End Do
      End Do
!
      iIrrep=0
      Fact=One
      If (lFij)
     &Call FckDst(TwoHam,nDens,Fij,iBas,jBas,iCmpa(1),iCmpa(2),
     &            kOp2(1),kOp2(2),iIrrep,
     &            iShij,
     &            iAO(1),iAO(2),iAOst(1),iAOst(2),
     &            Fact)
!
      If (lFkl)
     &Call FckDst(TwoHam,nDens,Fkl,kBas,lBas,iCmpa(3),iCmpa(4),
     &            kOp2(3),kOp2(4),iIrrep,
     &            iShkl,
     &            iAO(3),iAO(4),iAOst(3),iAOst(4),
     &            Fact)
!
      If (lFik)
     &Call FckDst(TwoHam,nDens,Fik,iBas,kBas,iCmpa(1),iCmpa(3),
     &            kOp2(1),kOp2(3),iIrrep,
     &            iShik,
     &            iAO(1),iAO(3),iAOst(1),iAOst(3),
     &            Fact)
!
      If (lFjl)
     &Call FckDst(TwoHam,nDens,Fjl,jBas,lBas,iCmpa(2),iCmpa(4),
     &            kOp2(2),kOp2(4),iIrrep,
     &            iShjl,
     &            iAO(2),iAO(4),iAOst(2),iAOst(4),
     &            Fact)
      If (lFil)
     &Call FckDst(TwoHam,nDens,Fil,iBas,lBas,iCmpa(1),iCmpa(4),
     &            kOp2(1),kOp2(4),iIrrep,
     &            iShil,
     &            iAO(1),iAO(4),iAOst(1),iAOst(4),
     &            Fact)
      If (lFjk)
     &Call FckDst(TwoHam,nDens,Fjk,jBas,kBas,iCmpa(2),iCmpa(3),
     &            kOp2(2),kOp2(3),iIrrep,
     &            iShjk,
     &            iAO(2),iAO(3),iAOst(2),iAOst(3),
     &            Fact)

      Nullify(Fij,Fkl,Fik,Fil,Fjk,Fjl)
      Contains

      Subroutine Fck1(AOInt,Dij,Fij,Cij,Dkl,Fkl,Ckl)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dij(iBas,jBas), Fij(iBas,jBas),
     &       Dkl(kBas,lBas), Fkl(kBas,lBas)
      Real*8 Cij, Ckl
      Integer i, j, k, l
      Real*8 F_kl, D_kl, Vijkl
!
      Do l = 1, lBas
         Do k = 1, kBas
            F_kl = Zero
            D_kl = Dkl(k,l)*Cij
            Do j = 1, jBas
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)
!
                  Fij(i,j) = Fij(i,j) +D_kl*Vijkl
                  F_kl = F_kl + Dij(i,j)*Vijkl
!
               End Do
            End Do
            Fkl(k,l) = Fkl(k,l) + Ckl*F_kl
         End Do
      End Do
!
      Return
      End Subroutine Fck1

      Subroutine Fck2(AOInt,Dik,Fik,Cik,Djl,Fjl,Cjl)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dik(iBas,kBas), Fik(iBas,kBas),
     &       Djl(jBas,lBas), Fjl(jBas,lBas)
      Real*8 Cik, Cjl
      Integer i, j, k, l
      Real*8 F_jl, D_jl, Vijkl
!
      Do l = 1, lBas
         Do k = 1, kBas
            Do j = 1, jBas
               F_jl = Zero
               D_jl = Djl(j,l)*Cik
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)*ExFac
!
                  Fik(i,k) = Fik(i,k) + D_jl*Vijkl
                  F_jl = F_jl + Dik(i,k)*Vijkl
!
               End Do
               Fjl(j,l) = Fjl(j,l) + Cjl*F_jl
            End Do
         End Do
      End Do
!
      Return
      End Subroutine Fck2

      Subroutine Fck3(AOInt,
     &                Dij,Fij,Cij,Dkl,Fkl,Ckl,
     &                Dik,Fik,Cik,Djl,Fjl,Cjl)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dij(iBas,jBas), Fij(iBas,jBas),
     &       Dkl(kBas,lBas), Fkl(kBas,lBas),
     &       Dik(iBas,kBas), Fik(iBas,kBas),
     &       Djl(jBas,lBas), Fjl(jBas,lBas)
      Real*8 Cij, Ckl, Cik, Cjl
      Integer i, j, k, l
      Real*8 F_kl, D_kl, F_jl, D_jl, Vijkl
!
      Do l = 1, lBas
         Do k = 1, kBas
            F_kl = Zero
            D_kl = Dkl(k,l)*Cij
            Do j = 1, jBas
               F_jl = Zero
               D_jl = Djl(j,l)*Cik
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)
!
                  Fij(i,j) = Fij(i,j) +D_kl*Vijkl
                  F_kl = F_kl + Dij(i,j)*Vijkl
!
                  Fik(i,k) = Fik(i,k) + D_jl*Vijkl*ExFac
                  F_jl = F_jl + Dik(i,k)*Vijkl
!
               End Do
               Fjl(j,l) = Fjl(j,l) + Cjl*F_jl*ExFac
            End Do
            Fkl(k,l) = Fkl(k,l) + Ckl*F_kl
         End Do
      End Do
!
      Return
      End Subroutine Fck3

      Subroutine Fck4(AOInt,Dil,Fil,Cil,Djk,Fjk,Cjk)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dil(iBas,lBas), Fil(iBas,lBas),
     &       Djk(jBas,kBas), Fjk(jBas,kBas)
      Real*8 Cil, Cjk
      Integer i, j, k, l
      Real*8 F_jk, D_jk, Vijkl
!
!
      Do l = 1, lBas
         Do k = 1, kBas
            Do j = 1, jBas
               F_jk = Zero
               D_jk = Djk(j,k)*Cil
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)*ExFac
!
                  Fil(i,l) = Fil(i,l) + D_jk*Vijkl
                  F_jk = F_jk + Dil(i,l)*Vijkl
!
               End Do
               Fjk(j,k) = Fjk(j,k) + Cjk*F_jk
            End Do
         End Do
      End Do
!
      Return
      End Subroutine Fck4

      Subroutine Fck5(AOInt,
     &                Dij,Fij,Cij,Dkl,Fkl,Ckl,
     &                Dil,Fil,Cil,Djk,Fjk,Cjk)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dij(iBas,jBas), Fij(iBas,jBas),
     &       Dkl(kBas,lBas), Fkl(kBas,lBas),
     &       Dil(iBas,lBas), Fil(iBas,lBas),
     &       Djk(jBas,kBas), Fjk(jBas,kBas)
      Real*8 Cij, Ckl, Cil, Cjk
      Integer i, j, k, l
      Real*8 F_kl, D_kl, F_jk, D_jk, Vijkl
!
      Do l = 1, lBas
         Do k = 1, kBas
            F_kl = Zero
            D_kl = Dkl(k,l)*Cij
            Do j = 1, jBas
               F_jk = Zero
               D_jk = Djk(j,k)*Cil
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)
!
                  Fij(i,j) = Fij(i,j) +D_kl*Vijkl
                  F_kl = F_kl + Dij(i,j)*Vijkl
!
                  Fil(i,l) = Fil(i,l) + D_jk*Vijkl*ExFac
                  F_jk = F_jk + Dil(i,l)*Vijkl
!
               End Do
               Fjk(j,k) = Fjk(j,k) + Cjk*F_jk*ExFac
            End Do
            Fkl(k,l) = Fkl(k,l) + Ckl*F_kl
         End Do
      End Do
!
      Return
      End Subroutine Fck5

      Subroutine Fck6(AOInt,
     &                Dik,Fik,Cik,Djl,Fjl,Cjl,
     &                Dil,Fil,Cil,Djk,Fjk,Cjk)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dik(iBas,kBas), Fik(iBas,kBas),
     &       Djl(jBas,lBas), Fjl(jBas,lBas),
     &       Dil(iBas,lBas), Fil(iBas,lBas),
     &       Djk(jBas,kBas), Fjk(jBas,kBas)
      Real*8 Cik, Cjl, Cil, Cjk
      Real*8 F_jl, D_jl, F_jk, D_jk, Vijkl
      Integer i, j, k, l
!
      Do l = 1, lBas
         Do k = 1, kBas
            Do j = 1, jBas
               F_jl = Zero
               D_jl = Djl(j,l)*Cik
               F_jk = Zero
               D_jk = Djk(j,k)*Cil
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)
!
                  Fik(i,k) = Fik(i,k) + D_jl*Vijkl*ExFac
                  F_jl = F_jl + Dik(i,k)*Vijkl
!
                  Fil(i,l) = Fil(i,l) + D_jk*Vijkl*ExFac
                  F_jk = F_jk + Dil(i,l)*Vijkl
!
               End Do
               Fjl(j,l) = Fjl(j,l) + Cjl*F_jl*ExFac
               Fjk(j,k) = Fjk(j,k) + Cjk*F_jk*ExFac
            End Do
         End Do
      End Do
!
      Return
      End Subroutine Fck6

      Subroutine Fck7(AOInt,
     &                Dij,Fij,Cij,Dkl,Fkl,Ckl,
     &                Dik,Fik,Cik,Djl,Fjl,Cjl,
     &                Dil,Fil,Cil,Djk,Fjk,Cjk)
      Real*8 AOInt(iBas,jBas,kBas,lBas),
     &       Dij(iBas,jBas), Fij(iBas,jBas),
     &       Dkl(kBas,lBas), Fkl(kBas,lBas),
     &       Dik(iBas,kBas), Fik(iBas,kBas),
     &       Djl(jBas,lBas), Fjl(jBas,lBas),
     &       Dil(iBas,lBas), Fil(iBas,lBas),
     &       Djk(jBas,kBas), Fjk(jBas,kBas)
      Real*8 Cij, Ckl, Cik, Cjl, Cil, Cjk
      Integer i, j, k, l
      Real*8 F_kl, D_kl, F_jl, D_jl, F_jk, D_jk, Vijkl
!
      Do l = 1, lBas
         Do k = 1, kBas
            F_kl = Zero
            D_kl = Dkl(k,l)*Cij
            Do j = 1, jBas
               F_jl = Zero
               D_jl = Djl(j,l)*Cik
               F_jk = Zero
               D_jk = Djk(j,k)*Cil
               Do i = 1, iBas
                  Vijkl = AOInt(i,j,k,l)
!
                  Fij(i,j) = Fij(i,j) +D_kl*Vijkl
                  F_kl = F_kl + Dij(i,j)*Vijkl
!
                  Fik(i,k) = Fik(i,k) + D_jl*Vijkl*ExFac
                  F_jl = F_jl + Dik(i,k)*Vijkl
!
                  Fil(i,l) = Fil(i,l) + D_jk*Vijkl*ExFac
                  F_jk = F_jk + Dil(i,l)*Vijkl
!
               End Do
               Fjl(j,l) = Fjl(j,l) + Cjl*F_jl*ExFac
               Fjk(j,k) = Fjk(j,k) + Cjk*F_jk*ExFac
            End Do
            Fkl(k,l) = Fkl(k,l) + Ckl*F_kl
         End Do
      End Do
!
      Return
      End Subroutine Fck7

      End Subroutine FckAcc
