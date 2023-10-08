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
      Subroutine FckAcc_NoSym(iCmp, jCmp, kCmp, lCmp, Shijij,
     &                        iShell, nijkl, AOInt,FMat,DMat,nDens,
     &                        iAO,iAOst,iBas,jBas,kBas,lBas,ExFac)
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
!***********************************************************************
      use SOAO_Info, only: iAOtSO
      use Gateway_Info, only: CutInt
      use Constants, only: Zero, One, Four, Half
      Implicit None
      Integer nijkl, iCmp, jCmp, kCmp, lCmp, nDens
      Real*8 AOInt(nijkl,iCmp,jCmp,kCmp,lCmp), FMat(nDens),
     &       DMat(nDens)
      Logical Shij, Shkl, Shijij
      Integer iShell(4), iAO(4), iAOst(4)
      Integer iBas, jBas, kBas, lBas
      Real*8 ExFac

!     Local Arrays
      Integer i, j, iTri
      Real*8 Fac, Fac_C, Fac_E
      Integer in(4)
      Integer i1, i2, i3, i4, iSO, jSO, kSO, lSO
      Integer k, l, kl, jk, jl, ij, ik, il
      Real*8 D_kl, F_kl, D_jl, F_jl, D_jk, F_jk, AOijkl
!     Statement Function
!
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!
!
!     Write (*,*) DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,AOInt,1),
!    &            DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,AOInt,1,One  ,0)
!     If (iPrint.ge.49) Then
!        Write (*,*) ' FckAcc:AOIn',DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
!    &               AOInt,1,AOInt,1), DDot_(nijkl*iCmp*jCmp*kCmp*lCmp,
!    &               AOInt,1,One,0)
!     End If
!
!     Quadruple loop over elements of the basis functions angular
!     description. Loops are reduced to just produce unique SO integrals
!     Observe that we will walk through the memory in AOInt in a
!     sequential way.
!
      Shij = iShell(1).eq.iShell(2)
      Shkl = iShell(3).eq.iShell(4)
!
      Fac=One
      If (Shij) Fac=Fac*Half
      If (Shkl) Fac=Fac*Half
      If (Shijij) Fac=Fac*Half
      Fac_C=Four*Fac
      Fac_E=-Fac*ExFac
!
      Do 100 i1 = 1, iCmp
         in(1)=iAOtSO(iAO(1)+i1,0)+iAOst(1)
         Do 200 i2 = 1, jCmp
            in(2)=iAOtSO(iAO(2)+i2,0)+iAOst(2)
            Do 300 i3 = 1, kCmp
               in(3)=iAOtSO(iAO(3)+i3,0)+iAOst(3)
               Do 400 i4 = 1, lCmp
                  in(4)=iAOtSO(iAO(4)+i4,0)+iAOst(4)
!
!................ Get the starting SO indices.
!
                  iSO = in(1)
                  jSO = in(2)
                  kSO = in(3)
                  lSO = in(4)
!
                  nijkl=0
                  Do l = lSO, lSO+lBas-1
                     Do k = kSO, kSO+kBas-1
                        kl=iTri(k,l)
                        D_kl=DMat(kl)*Fac_C
                        F_kl=Zero
!
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
!
                              FMat(ij) = FMat(ij)+AOijkl*D_kl
                              F_kl     = F_kl    +AOijkl*DMat(ij)
!
                              FMat(ik) = FMat(ik)+AOijkl*D_jl
                              F_jl     = F_jl    +AOijkl*DMat(ik)
                              FMat(il) = FMat(il)+AOijkl*D_jk
                              F_jk     = F_jk    +AOijkl*DMat(il)
!
 10                           Continue
                           End Do
                           FMat(jl)=FMat(jl)+F_jl*Fac_E
                           FMat(jk)=FMat(jk)+F_jk*Fac_E
                        End Do
                        FMat(kl)=FMat(kl)+F_kl*Fac_C
                     End Do
                  End Do
!
 400           Continue
 300        Continue
 200     Continue
 100  Continue
!
      End Subroutine FckAcc_NoSym
