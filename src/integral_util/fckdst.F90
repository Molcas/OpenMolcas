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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine FckDst(TwoHam,nDens,Fij,iBas,jBas,iCmp,jCmp,ikop1,ikop2,Irrep,Shij,iAO1,iAO2,iAOst1,iAOst2,fact)

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: iChTbl, iOper, nIrrep
use SOAO_Info, only: iAOtSO, nSOInf
use Basis_Info, only: nBas
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDens, iBas, jBas, iCmp, jCmp, ikOp1, ikOp2, Irrep, iAO1, iAO2, iAOst1, iAOst2
real(kind=wp), intent(inout) :: TwoHam(nDens)
real(kind=wp), intent(in) :: Fij(0:iBas-1,0:jBas-1,iCmp,jCmp), Fact
logical(kind=iwp), intent(in) :: Shij
integer(kind=iwp) :: i1, i2, iAO, iAOi, iChO, iIR, iIrrep, ipF, ipFij, iPnt(0:7), iPntij, iSO, iSOi, iSOj, jAO, jAOj, jirr(0:7), &
                     jIrrep, jSO, jSOi, jSOj, l1, l2, NrOpr
real(kind=wp) :: Fac, x1, x2, x3, x4, xr

iChO = iOper(Irrep)
if (iChO == 0) then

  iPntij = 0
  do iIrrep=0,nIrrep-1
    iPnt(iIrrep) = iPntij
    iPntij = iPntij+nTri_Elem(nBas(iIrrep))
  end do

  ! Distribute contributions from the intermediate skeleton
  ! Fock matrix onto the symmetry adapted Fock matrix.

  iiR = NrOpr(ieor(ikOp1,ikOp2))
  do i1=1,iCmp
    do i2=1,jCmp
      do iIrrep=0,nIrrep-1
        if (max(iAO1+i1,iAO2+i2) > nSOInf) then
          write(u6,*) 'Fckdst: Max(iSO.jSO)>nSOInf (1)'
          call Abend()
        end if
        iSO = iAOtSO(iAO1+i1,iIrrep)+iAOst1
        jSO = iAOtSO(iAO2+i2,iIrrep)+iAOst2
        XR = real(iChTbl(iIrrep,iiR),kind=wp)

        if ((iSO < 0) .or. (jSO < 0)) cycle

        iPntij = iPnt(iIrrep)
        do jAOj=0,jBas-1
          do iAOi=0,iBas-1
            Fac = XR
            if (Shij .and. (i1 == i2) .and. (iAOi+iAOst1 == jAOj+iAOst2)) Fac = Two*XR
            jSOj = jSO+jAOj
            iSOi = iSO+iAOi
            ipFij = iPntij+iTri(iSOi,jSOj)
            TwoHam(ipFij) = TwoHam(ipFij)+Fact*Fac*Fij(iAOi,jAOj,i1,i2)
          end do
        end do

      end do ! iIrrep
    end do   ! i2
  end do     ! i1

else

  do iIrrep=0,nIrrep-1
    jIrr(iIrrep) = NrOpr(ieor(iOper(iIrrep),iChO))
  end do

  iPntij = 0
  iPnt(0:nIrrep-1) = -1
  do iIrrep=0,nIrrep-1
    jIrrep = jIrr(iIrrep)
    if (iIrrep > jIrrep) then
      iPnt(iIrrep) = iPntij
      iPntij = iPntij+nBas(jIrrep)*nBas(iIrrep)
    end if
  end do

  ! Distribute contributions from the intermediate skeleton
  ! Fock matrix onto the symmetry adapted Fock matrix.

  l1 = NrOpr(ikop1)
  l2 = NrOpr(ikop2)
  do i1=1,iCmp
    do i2=1,jCmp
      do iIrrep=0,nIrrep-1
        jIrrep = jIrr(iIrrep)
        if (iIrrep < jIrrep) cycle
        X1 = real(iChTbl(iIrrep,l1),kind=wp)
        X2 = real(iChTbl(jIrrep,l2),kind=wp)
        X3 = real(iChTbl(jIrrep,l1),kind=wp)
        X4 = real(iChTbl(iIrrep,l2),kind=wp)
        if (max(iAO1+i1,iAO2+i2) > nSOInf) then
          write(u6,*) 'Fckdst: Max(iSO.jSO)>nSOInf (2)'
          call Abend()
        end if
        iSOi = iAOtSO(iAO1+i1,iIrrep)+iAOst1
        jSOj = iAOtSO(iAO2+i2,jIrrep)+iAOst2
        jSOi = iAOtSO(iAO2+i2,iIrrep)+iAOst2
        iSOj = iAOtSO(iAO1+i1,jIrrep)+iAOst1
        if ((iSOi > -1) .and. (jSOj > -1)) then
          iPntij = iPnt(iIrrep)
          do jAO=0,jBas-1
            do iAO=0,iBas-1
              Fac = X1*X2
              ipF = iPntij+(jSOj+jAO-1)*nBas(iIrrep)+iSOi+iAO
              TwoHam(ipF) = TwoHam(ipF)+Fact*Fac*Fij(iAO,jAO,i1,i2)
            end do
          end do
        end if
        if ((jSOi > -1) .and. (iSOj > -1)) then
          iPntij = iPnt(iIrrep)
          do jAO=0,jBas-1
            do iAO=0,iBas-1
              Fac = X3*X4
              ipF = iPntij+nBas(iIrrep)*(iSOj+iAO-1)+jSOi+jAO
              TwoHam(ipF) = TwoHam(ipF)+Fact*Fac*Fij(iAO,jAO,i1,i2)
            end do
          end do
        end if
      end do
    end do
  end do

end if

return

end subroutine FckDst
