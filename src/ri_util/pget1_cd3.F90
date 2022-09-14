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

subroutine PGet1_CD3(PAO,ijkl,nPAO,iCmp,iAO,iAOst,iBas,jBas,kBas,lBas,kOp,DSO,DSO_Var,nDSO,ExFac,CoulFac,PMax,V_k,U_k,mV_k)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified for Cholesky 1-center gradients May 2007 by     *
!             R. Lindh                                                 *
!***********************************************************************

use Index_Functions, only: iTri
use pso_stuff, only: ij2K, iOff_ij2K
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use RI_glob, only: BklK, BMP2, CijK, CilK, CMOi, iAdrCVec, iMP2prpt, LuBVector, LuCVector, nChOrb, nIJR, tbvec
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), iBas, jBas, kBas, lBas, kOp(4), nDSO, mV_k
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: DSO(nDSO), DSO_Var(nDSO), ExFac, CoulFac, V_k(mV_k), U_k(mV_k)
integer(kind=iwp) :: i, i1, i2, i3, i4, iAdr, iAOi, ijVec, indexB, Indij, Indkl, iPAO, irc, iSO, iSOi, jAOj, jSO, jSOj, jSym, &
                     kAOk, kSO, kSOk, kSym, lAOl, lBVec, lSO, lSOl, lSym, nijkl, nKBas, nLBas, NumOrb
real(kind=wp) :: Cpu, Cpu1, Cpu2, Fac, Fac_ij, temp, tempJ, tempK, Wall, Wall1, Wall2
real(kind=wp), external :: Compute_B

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iComp = 1
call PrMtrx('DSO     ',[iD0Lbl],iComp,1,D0)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Quadruple loop over elements of the basis functions angular description.
! Observe that we will walk through the memory in PAO in a sequential way.
!
!Fac = Quart

call CWTime(Cpu1,Wall1)

Fac = Half
PMax = Zero
iPAO = 0
jSym = 1
kSym = 1
lSym = 1
NumOrb = nChOrb(kSym-1,1)

if ((ExFac /= Zero) .and. (NumOrb > 0) .and. (iMP2prpt /= 2)) then

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  do i1=1,iCmp(1)
    iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
    do iAOi=0,iBas-1
      iSOi = iSO+iAOi
      do i2=1,iCmp(2)
        jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj
          if (iSOi == jSOj) then
            Fac_ij = One
          else
            Fac_ij = Half
          end if
          Indij = iTri(iSOi,jSOj)
          ijVec = ij2K(iOff_ij2K(1)+Indij)

          if (ijVec /= 0) then
            iAdr = nIJR(kSym,lSym,1)*(ijVec-1)+iAdrCVec(jSym,kSym,1)
            call dDaFile(LuCVector(jSym,1),2,CijK,nIJR(kSym,lSym,1),iAdr)

            call dGEMM_('T','N',NumOrb,nKBas,NumOrb,One,CijK,NumOrb,CMOi(1)%SB(1)%A2(:,kSO),NumOrb,Zero,CilK,max(1,NumOrb))

            call dGEMM_('T','N',nKBas,nLBas,NumOrb,One,CilK,NumOrb,CMOi(1)%SB(1)%A2(:,lSO),NumOrb,Zero,BklK,max(1,nKBas))
          end if

          do i3=1,iCmp(3)
            kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
            do i4=1,iCmp(4)
              lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

              iPAO = i4+(i3-1)*iCmp(4)+(i2-1)*iCmp(4)*iCmp(3)+(i1-1)*iCmp(4)*iCmp(3)*iCmp(2)

              do lAOl=0,lBas-1
                lSOl = lSO+lAOl
                do kAOk=0,kBas-1
                  kSOk = kSO+kAOk
                  indexB = 1+(kAOk+(i3-1)*kBas)+(lAOl+(i4-1)*lBas)*nKBas
                  nijkl = iAOi+jAOj*iBas+kAOk*iBas*jBas+lAOl*iBas*jBas*kBas+1

                  Indkl = iTri(kSOk,lSOl)
                  temp = V_k(Indij)*DSO(Indkl)*CoulFac
                  if (ijVec /= 0) then
                    tempK = BklK(indexB)
                  else
                    tempK = Zero
                  end if

                  temp = temp-tempK*ExFac*Half*fac_ij
                  PMax = max(PMax,abs(Temp))
                  PAO(nijkl,iPAO) = Fac*temp

                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

else if ((iMP2prpt == 2) .and. (NumOrb > 0)) then

  nKBas = kBas*iCmp(3)
  nLBas = lBas*iCmp(4)

  kSO = iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)
  lSO = iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)

  do i1=1,iCmp(1)
    iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
    do iAOi=0,iBas-1
      iSOi = iSO+iAOi
      do i2=1,iCmp(2)
        jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj
          if (iSOi == jSOj) then
            Fac_ij = One
          else
            Fac_ij = Half
          end if
          Indij = iTri(iSOi,jSOj)
          ijVec = ij2K(iOff_ij2K(1)+Indij)
          if (ijVec /= 0) then
            iAdr = nIJR(kSym,lSym,1)*(ijVec-1)+iAdrCVec(jSym,kSym,1)
            call dDaFile(LuCVector(jSym,1),2,CijK,nIJR(kSym,lSym,1),iAdr)

            call dGEMM_('T','N',NumOrb,nKBas,NumOrb,One,CijK,NumOrb,CMOi(1)%SB(1)%A2(:,kSO),NumOrb,Zero,CilK,max(1,NumOrb))

            call dGEMM_('T','N',nKBas,nLBas,NumOrb,One,CilK,NumOrb,CMOi(1)%SB(1)%A2(:,lSO),NumOrb,Zero,BklK,max(1,nKBas))
            lBVec = nBas(0)*nBas(0)
            do i=1,2
              iAdr = 1+nBas(0)*nBas(0)*(ijVec-1)
              call dDaFile(LuBVector(i),2,Bmp2(:,i),lBVec,iAdr)
            end do

          end if
          do i3=1,iCmp(3)
            kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
            do i4=1,iCmp(4)
              lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

              iPAO = i4+(i3-1)*iCmp(4)+(i2-1)*iCmp(4)*iCmp(3)+(i1-1)*iCmp(4)*iCmp(3)*iCmp(2)

              do lAOl=0,lBas-1
                lSOl = lSO+lAOl
                do kAOk=0,kBas-1
                  kSOk = kSO+kAOk
                  indexB = 1+(kAOk+(i3-1)*kBas)+(lAOl+(i4-1)*lBas)*nKBas
                  nijkl = iAOi+jAOj*iBas+kAOk*iBas*jBas+lAOl*iBas*jBas*kBas+1

                  Indkl = iTri(kSOk,lSOl)
                  temp = V_k(Indij)*DSO(Indkl)*CoulFac

                  if (ijVec /= 0) then
                    tempK = BklK(indexB)
                  else
                    tempK = Zero
                  end if
                  temp = temp+U_k(indij)*DSO(indkl)*CoulFac
                  temp = temp+V_k(indij)*(DSO_Var(indkl)-DSO(indkl))*CoulFac
                  if (ijVec /= 0) then
                    tempJ = Compute_B(irc,kSOk,lSOl,0,nBas(0),2)
                    temp = temp+tempJ*CoulFac*fac_ij

                    tempK = tempK+Compute_B(irc,kSOk,lSOl,0,nBas(0),1)
                  end if
                  temp = temp-tempK*ExFac*Half*fac_ij
                  PMax = max(PMax,abs(Temp))
                  PAO(nijkl,iPAO) = Fac*temp

                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

else
  do i1=1,iCmp(1)
    do i2=1,iCmp(2)
      do i3=1,iCmp(3)
        do i4=1,iCmp(4)
          iSO = iAOtSO(iAO(1)+i1,kOp(1))+iAOst(1)
          jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
          kSO = iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
          lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

          iPAO = iPAO+1

          nijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl
            do kAOk=0,kBas-1
              kSOk = kSO+kAOk
              Indkl = iTri(kSOk,lSOl)
              do jAOj=0,jBas-1
                jSOj = jSO+jAOj
                do iAOi=0,iBas-1
                  iSOi = iSO+iAOi
                  nijkl = nijkl+1

                  ! V_k(ij)*D(kl)

                  Indij = iTri(iSOi,jSOj)

                  temp = V_k(Indij)*DSO(Indkl)*coulfac

                  PMax = max(PMax,abs(Temp))
                  PAO(nijkl,iPAO) = Fac*temp
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end if
if (iPAO /= nPAO) then
  write(u6,*) ' Error in PGet1_CD3!'
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1_CD3:PAO ',' ',PAO,ijkl,nPAO)
#endif

call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tbvec(1) = tbvec(1)+Cpu
tbvec(2) = tbvec(2)+Wall

return

end subroutine PGet1_CD3
