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
! Copyright (C) 1992,2007, Roland Lindh                                *
!***********************************************************************

subroutine PGet2_RI3(iCmp,jBas,kBas,lBas,iAO,iAOst,nijkl,PSO,nPSO,DSO,nDSO,ExFac,CoulFac,PMax,V_k,mV_k,ZpK,nSA,nAct)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified for 3-center RI gradients, March 2007           *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use SOAO_Info, only: iAOtSO
use pso_stuff, only: AOrb, lPSO, nnP, Thpkl
use Basis_Info, only: nBas, nBas_Aux
use Symmetry_Info, only: Mul, nIrrep
use RI_glob, only: BklK, CijK, CilK, CMOi, iAdrCVec, iOff_Ymnij, LuCVector, nAvec, nChOrb, nIJR, nYmnij, tbvec, Yij, Ymnij
use Constants, only: Zero, One, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), jBas, kBas, lBas, iAO(4), iAOst(4), nijkl, nPSO, nDSO, mV_k, nSA, nAct(0:7)
real(kind=wp), intent(out) :: PSO(nijkl,nPSO), PMax
real(kind=wp), intent(in) :: DSO(nDSO,nSA), ExFac, CoulFac, V_k(mV_k,nSA), Zpk(*)
integer(kind=iwp) :: i, i2, i3, i4, iAdr, ij, ik, il, imo, iMO1, iMO2, Indkl, iSO, iThpkl, iVec, j, j2, j23, j3, j4, jAOj, jC, &
                     jmo, jp, js, jSO, jSO_off, jSOj, jSym(0:7), k, kAct, kAOk, kmo, ks, kSO, kSOk, kSym(0:7), l, lAct, lAOl, &
                     lCVec, lda, lmo, lOper, ls, lSO, lSOl, lSym(0:7), MemSO2, mijkl, n2j, nCumnnP(0:7), nCumnnP2(0:7), nJ, njSym, &
                     nk, nkSym, nl, nlSym, ntmp
real(kind=wp) :: Cpu, Cpu1, Cpu2, ExFac_, Fac, temp, tmp, Wall, Wall1, Wall2
real(kind=wp), pointer :: Xki(:), Xli(:)
integer(kind=iwp), external :: iPntSO
real(kind=wp), external :: ddot_

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx(' In PGET_RI3:DSO ',[iD0Lbl],iComp,1,D0)
call RecPrt('V_K',' ',V_K,1,mV_K)
write(u6,*)
write(u6,*) 'Distribution of Ymnij'
do iSym=1,nIrrep
  if (nYmnij(iSym,1) > 0) then
    write(u6,*) 'iSym=',iSym
    do i=iOff_Ymnij(iSym,1)+1,iOff_Ymnij(iSym,1)+nYmnij(iSym,1)
      write(u6,*) 'Ymnij=',Ymnij(1)%A(i)
    end do
  end if
end do
write(u6,*) 'jbas,kbas,lbas=',jBas,kBas,lBas
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(Cpu1,Wall1)

Fac = Quart
lOper = 1
PMax = Zero
iSO = 1
PSO(:,:) = Zero

if (lPSO) then
  nCumnnP(0) = 0
  nBas_Aux(0) = nBas_Aux(0)-1
  do i=1,nIrrep-1
    nCumnnP(i) = nCumnnP(i-1)+nnP(i-1)*nBas_Aux(i-1)
  end do
  nBas_Aux(0) = nBas_Aux(0)+1
end if

! i2, j2, jBas: auxiliary basis
! i3, j3, kBas: valence basis
! i4, j4, lBas: valence basis
!
! Note when j2 is symmetric then we can have both Coulomb and
! exchange contributions, while for j2 asymmetric we will
! only have exchange contributions

MemSO2 = 0
do i2=1,iCmp(2)
  njSym = 0
  do j=0,nIrrep-1
    if (iAOtSO(iAO(2)+i2,j) > 0) then
      jSym(njSym) = j
      njSym = njSym+1
    end if
  end do
  do i3=1,iCmp(3)
    nkSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(3)+i3,j) > 0) then
        kSym(nkSym) = j
        nkSym = nkSym+1
      end if
    end do
    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do

      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        !nJ = nChOrb(j2,iSO)
        nJ = jBas

        if (lPSO) then
          ntmp = 0
          do j4=0,nIrrep-1
            j3 = Mul(j4+1,j2+1)-1
            if (j3 <= j4) nCumnnP2(j3) = ntmp
            if (j3 == j4) ntmp = ntmp+nTri_Elem(nAct(j3))
            if (j3 < j4) ntmp = ntmp+nAct(j3)*nAct(j4)
          end do
          do j4=0,nIrrep-1
            j3 = Mul(j4+1,j2+1)-1
            if (j3 > j4) nCumnnP2(j3) = nCumnnP2(j4)
          end do
        end if

        do ks=0,nkSym-1
          j3 = kSym(ks)
          j23 = Mul(j2+1,j3+1)-1
          nk = nYmnij(j3+1,1)
          kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)

          ! Pointers to the full list of the X_mu,i elements
          ! Note this list runs over all basis functions mu
          ! (kBas*iCmp(3)). Here we only want to pick up the
          ! subblock for a fixed iCmp(3) value.

          if ((nk < nChOrb(j3,iSO)) .and. (ExFac /= Zero) .and. (nk > 0)) then

            ! Offset to where the block starts (jbas,kbas,i3)

            lda = size(CMOi(1)%SB(j3+1)%A2,1)
            ik = 1+lda*(kSO-1)
            Xki(1:) => CMOi(1)%SB(j3+1)%A1(ik:)

            ! Loop over the auxiliary basis functions which has
            ! significant contributions to the k shell.

            imo = 1
            do k=1,nk
              kmo = Ymnij(1)%A(k+iOff_Ymnij(j3+1,1))

              call dcopy_(kBas,Xki(kmo:),nChOrb(j3,iSO),Yij(imo,1,1),nk)

              imo = imo+1
            end do
            ! Reset pointers
            Xki(1:nk*kBas) => Yij(1:nk*kBas,1,1)
            !call RecPrt('X(i,mu)C',' ',Xki,nk,kBas)
          else if ((ExFac /= Zero) .and. (nk > 0)) then
            lda = size(CMOi(1)%SB(j3+1)%A2,1)
            ik = 1+lda*(kSO-1)
            Xki(1:) => CMOi(1)%SB(j3+1)%A1(ik:)
            !call RecPrt('X(i,mu)R',' ',Xki,nk,kBas)
          else
            nullify(Xki)
          end if

          do ls=0,nlSym-1
            j4 = lSym(ls)
            if (j23 /= j4) cycle
            nl = nYmnij(j4+1,1)
            lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

            ! Pointers to the full list of the X_mu,i elements

            if ((nl < nChOrb(j4,iSO)) .and. (ExFac /= Zero) .and. (nl > 0)) then

              lda = size(CMOi(1)%SB(j4+1)%A2,1)
              il = 1+lda*(lSO-1)
              Xli(1:) => CMOi(1)%SB(j4+1)%A1(il:)
              imo = 1
              do l=1,nl
                lmo = Ymnij(1)%A(l+iOff_Ymnij(j4+1,1))

                call dcopy_(lBas,Xli(lmo:),nChOrb(j4,iSO),Yij(imo,2,1),nl)

                imo = imo+1
              end do
              ! Reset pointers
              Xli(1:nl*lBas) => Yij(1:nl*lBas,2,1)
              !call RecPrt('X(j,nu)C',' ',Xli,nl,lBas)
            else if ((ExFac /= Zero) .and. (nl > 0)) then
              lda = size(CMOi(1)%SB(j4+1)%A2,1)
              il = 1+lda*(lSO-1)
              Xli(1:) => CMOi(1)%SB(j4+1)%A1(il:)
              !call RecPrt('X(j,nu)R',' ',Xli,nl,lBas)
            else
              nullify(Xli)
            end if

            MemSO2 = MemSO2+1

            jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
            jSO_off = jSO-nBas(j2)

            ExFac_ = ExFac
            if (nJ*nk*nl == 0) ExFac_ = Zero
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Read a block of C_kl^J and transform it to AO basis.

            if (ExFac_ /= Zero) then

              ! Read C(i,j,J) for a fix i2 value

              lCVec = nIJR(j3+1,j4+1,iSO)*jBas
              iAdr = iAdrCVec(j2+1,j3+1,1)+nIJR(j3+1,j4+1,iSO)*(jSO_Off-1)
              call dDaFile(LuCVector(j2+1,1),2,Cijk,lCVec,iAdr)
              !call RecPrt('C(ij,K)',' ',CijK,nIJR(j3+1,j4+1,iSO),jBas)

              ! Extract only those C_kl^Js for which we deem k and l
              ! to belong to the shell-pair and to be of
              ! significance. Use temporary memory location at CilK.

              if (nk*nl < nChOrb(j3,iSO)*nChOrb(j4,iSO)) then
                ij = 1
                do j=1,nl
                  jmo = Ymnij(1)%A(j+iOff_Ymnij(j4+1,1))
                  do i=1,nk
                    imo = Ymnij(1)%A(i+iOff_Ymnij(j3+1,1))

                    jC = imo+nChOrb(j3,iSO)*(jmo-1)

                    ! For this particular ij combination pick the whole row.

                    call dcopy_(jBas,CijK(jC),nChOrb(j3,iSO)*nChOrb(j4,iSO),CilK(ij),nk*nl)
                    ij = ij+1
                  end do
                end do

                ! Copy back to original memory position.

                n2j = nk*nl*jBas
                CijK(1:n2j) = CilK(1:n2j)
              end if

              ! Transform according to Eq. 16 (step 4) and
              ! generate B_kl^J. This is a transformation from
              ! the MO basis, ij, to the AO basis mn.

              ! E(jK,m) = Sum_i C(i,jK)' * X(i,m)

              call dGEMM_('T','N',nl*jBas,kBas,nk,One,CijK,nk,Xki,nk,Zero,CilK,nl*jBas)

              ! B(Km,n) = Sum_j E(j, Km)' * X(j,n)

              call dGEMM_('T','N',jBas*kBas,lBas,nl,One,CilK,nl,Xli,nl,Zero,BklK,jBas*kBas)

            end if
            !                                                          *
            !***********************************************************
            !                                                          *
            ! Active term (CASSCF and SA-CASSCF)

            if (lPSO) then
              Thpkl(1:jBas*kBas*lBas) = Zero
              if (nAct(j3)*nAct(j4) /= 0) then
                do iVec=1,nAVec
                  iMO1 = 1
                  iMO2 = 1
                  if (iVec == 2) iMO2 = 2
                  if (iVec == 4) then
                    iMO1 = 2
                  end if

                  do jAOj=0,jBas-1
                    jSOj = jSO+jAOj-nBas(j2)
                    jp = nCumnnP(j2)+(jSOj-1)*nnP(j2)+nCumnnP2(j3)
                    do lAOl=0,lBas-1
                      lSOl = lSO+lAOl

                      if (j3 == j4) then
                        do kAct=1,nAct(j3)
                          ! Zpk(*,iVec)
                          tmp = ddot_(kAct,Zpk(jp+iTri(kAct,1)),1,AOrb(iMO1)%SB(j4+1)%A2(:,lSOl),1)
                          do lAct=kAct+1,nAct(j4)
                            tmp = tmp+Zpk(jp+iTri(lAct,kAct))*AOrb(iMO1)%SB(j4+1)%A2(lAct,lSOl)
                          end do
                          Cilk(kAct) = tmp
                        end do
                      else
                        if (j3 < j4) then
                          call dGeMV_('N',nAct(j3),nAct(j4),One,Zpk(jp+1),nAct(j3),AOrb(iMO1)%SB(j4+1)%A2(:,lSOl),1,Zero,CilK,1)
                        else
                          call dGeMV_('T',nAct(j4),nAct(j3),One,Zpk(jp+1),nAct(j4),AOrb(iMO1)%SB(j4+1)%A2(:,lSOl),1,Zero,CilK,1)
                        end if
                      end if

                      iThpkl = jAOj+lAOl*kBas*jBas+1
                      call dGeMV_('T',nAct(j3),kBas,One,AOrb(iMO2)%SB(j3+1)%A2(:,kSO),nAct(j3),Cilk,1,One,Thpkl(iThpkl),jBas)

                    end do
                  end do
                end do
              end if
            end if

            !                                                          *
            !***********************************************************
            !                                                          *
            if (ExFac_ /= Zero) then
            !                                                          *
            !***********************************************************
            !                                                          *

#             define _EXCHANGE_
              if (j3 /= j4) then
                if (lPSO) then
                  ! Exchange and active contributions
#                 define _ACTIVE_
#                 include "pget2_ri3.fh"
#                 undef _ACTIVE_
                else
                  ! Exchange contribution
#                 include "pget2_ri3.fh"
                end if
              else
#               define _COULOMB_
                if (lPSO) then
                  ! Coulomb, Exchange and active contributions
#                 define _ACTIVE_
#                 include "pget2_ri3.fh"
#                 undef _ACTIVE_
                else
                  ! Coulomb and Exchange contributions
#                 include "pget2_ri3.fh"
                end if
#               undef _COULOMB_
              end if
#             undef _EXCHANGE_
              !                                                        *
              !*********************************************************
              !                                                        *
            else if ((ExFac_ == Zero) .and. (j3 == j4)) then
              !                                                        *
              !*********************************************************
              !                                                        *

#             define _COULOMB_
              if (lPSO) then
                ! Coulomb and active contributions
#               define _ACTIVE_
#               include "pget2_ri3.fh"
#               undef _ACTIVE_
              else
                ! Coulomb only contribution
#               include "pget2_ri3.fh"
              end if
#             undef _COULOMB_
              !                                                        *
              !*********************************************************
              !                                                        *
            end if
            !                                                          *
            !***********************************************************
            !                                                          *

          end do
          nullify(Xki,Xli)
        end do
      end do

    end do
  end do
end do
if (nPSO /= MemSO2) then
  write(u6,*) ' PGET_RI3: nPSO /= MemSO2'
  write(u6,*) nPSO,MemSO2
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGET_RI3:PSO ',' ',PSO,nijkl,nPSO)
#endif

call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tbvec(1) = tbvec(1)+Cpu
tbvec(2) = tbvec(2)+Wall

return

end subroutine PGet2_RI3
