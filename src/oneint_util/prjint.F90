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
!               1993, Per Boussard                                     *
!***********************************************************************

subroutine PrjInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of ECP integrals.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
!             Physics, University of Stockholm, Sweden, October '93.   *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Real_Spherical, only: ipSph, RSph
use Symmetry_Info, only: iChTbl, nIrrep
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: ia, iaC, iAng, ib, iBk, iC, iCb, iCnt, iCnttp, iComp, iDCRT(0:7), iIC, iIrrep, ip, ipaC, ipCb, ipF1, ipf2, &
                     ipK1, ipK2, ipOff, ipP1, ipP2, ipTmp, ipZ1, ipZ2, ipZI1, ipZI2, iShll, lDCRT, llOper, LmbdT, mArr, mdc, nac, &
                     nBasisi, ncb, nDCRT, nExpi, nOp
real(kind=wp) :: C(3), Fact, Factor, TC(3), Xg
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif
integer(kind=iwp), external :: NrOpr

#include "macros.fh"
unused_var(P)
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(iChO)
unused_var(PtChrg)
unused_var(iAddPot)

#ifdef _DEBUGPRINT_
call RecPrt(' In PrjInt: A',' ',A,1,3)
call RecPrt(' In PrjInt: RB',' ',RB,1,3)
call RecPrt(' In PrjInt: Ccoor',' ',Ccoor,1,3)
write(u6,*) ' In PrjInt: la,lb=',' ',la,lb
#endif

rFinal(:,:,:,:) = Zero

llOper = lOper(1)
iComp = 1
mdc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%ECP) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

      call DCR(LmbdT,iStabM,nStabM,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,iDCRT,nDCRT)
      Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

      do lDCRT=0,nDCRT-1
        call OA(iDCRT(lDCRT),C,TC)
        do iAng=0,dbsc(iCnttp)%nPrj-1
          iShll = dbsc(iCnttp)%iPrj+iAng
          nExpi = Shells(iShll)%nExp
          nBasisi = Shells(iShll)%nBasis
          if ((nExpi == 0) .or. (nBasisi == 0)) cycle

#         ifdef _DEBUGPRINT_
          call RecPrt('Cff',' ',Shells(iShll)%pCff,nExpi,nBasisi)
#         endif
          ip = 1
          ipF1 = ip
          nac = nTri_Elem1(la)*nTri_Elem1(iAng)
          ip = ip+nAlpha*nExpi*nac
          ipP1 = ip
          ip = ip+3*nAlpha*nExpi
          ipZ1 = ip
          ip = ip+nAlpha*nExpi
          ipK1 = ip
          ip = ip+nAlpha*nExpi
          ipZI1 = ip
          ip = ip+nAlpha*nExpi
          if (ip-1 > nArr*nZeta) then
            call WarningMessage(2,'PrjInt: ip-1 > nArr*nZeta(1)')
            call Abend()
          end if
          mArr = (nArr*nZeta-(ip-1))/nZeta

          ! Calculate Effective center and exponent for <A|core>

          call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,Alpha,Shells(iShll)%Exp)
          call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

          ! Calculate Overlap <A|core>

          nHer = (la+iAng+2)/2
          call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,Array(ipZ1),Array(ipZI1),Array(ipK1),Array(ipP1),Array(ipF1), &
                      nAlpha*nExpi,iComp,la,iAng,A,TC,nHer,Array(ip),mArr,CCoor,nOrdOp)
          ip = ip-6*nAlpha*nExpi

          ipF2 = ip
          ncb = nTri_Elem1(iAng)*nTri_Elem1(lb)
          ip = ip+nExpi*nBeta*ncb
          ipP2 = ip
          ip = ip+3*nExpi*nBeta
          ipZ2 = ip
          ip = ip+nExpi*nBeta
          ipK2 = ip
          ip = ip+nExpi*nBeta
          ipZI2 = ip
          ip = ip+nExpi*nBeta
          if (ip-1 > nArr*nZeta) then
            call WarningMessage(2,'PrjInt: ip-1 > nArr*nZeta(2)')
            call Abend()
          end if
          mArr = (nArr*nZeta-(ip-1))/nZeta

          ! Calculate Effective center and exponent for <core|B>

          call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,Shells(iShll)%Exp,Beta)
          call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

          ! Calculate Overlap <core|B>

          nHer = (iAng+lb+2)/2
          call MltPrm(Shells(iShll)%Exp,nExpi,Beta,nBeta,Array(ipZ2),Array(ipZI2),Array(ipK2),Array(ipP2),Array(ipF2),nExpi*nBeta, &
                      iComp,iAng,lb,TC,RB,nHer,Array(ip),mArr,CCoor,nOrdOp)
          ip = ip-6*nExpi*nBeta
          ipTmp = ip
          ip = ip+max(nAlpha*nExpi*nac,nBeta*ncb*nBasisi)
          if (ip-1 > nArr*nZeta) then
            call WarningMessage(2,'PrjInt: ip-1 > nArr*nZeta(3)')
            call Abend()
          end if

          ! Calculate Contraction over components of the core
          ! orbitals of type <A|core>Bc<core|B> where we now have in
          ! Array(ipF1) the cartesian components of <A|core>, and
          ! similarily, in Array(ipF2), we have stored the cartesian
          ! components of <core|B>. Observe that the core orbitals
          ! orthonomal atomic orbitals. Hence, the transformation
          ! to the spherical harmonics has to be for normilized
          ! spherical harmonics.

          ! From the lefthandside overlap, form iKaC from ikac by
          ! 1) i,kac -> k,aci

          call DgeTMo(Array(ipF1),nAlpha,nAlpha,nExpi*nac,Array(ipTmp),nExpi*nac)

          ! 2) aciK =  k,aci * k,K

          call DGEMM_('T','N',nAlpha*nac,nBasisi,nExpi,One,Array(ipTmp),nExpi,Shells(iShll)%pCff,nExpi,Zero,Array(ipF1),nAlpha*nac)

          ! 3) Mult by shiftoperators aci,K -> Bk(K) * aci,K

          ipOff = ipF1-1
          do iBk=1,nBasisi
            Array(ipOff+1:ipOff+nAlpha*nac) = Shells(ishll)%Bk(iBk)*Array(ipOff+1:ipOff+nAlpha*nac)
            ipOff = ipOff+nAlpha*nac
          end do ! iBk

          ! 4) a,ciK -> ciKa

          call DgeTMo(Array(ipF1),nTri_Elem1(la),nTri_Elem1(la),nTri_Elem1(iAng)*nAlpha*nBasisi,Array(ipTmp), &
                      nTri_Elem1(iAng)*nAlpha*nBasisi)

          ! 5) iKa,C = c,iKa * c,C

          call DGEMM_('T','N',nAlpha*nBasisi*nTri_Elem1(la),(2*iAng+1),nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                      RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF1),nAlpha*nBasisi*nTri_Elem1(la))

          ! And (almost) the same thing for the righthand side, form
          ! KjCb from kjcb
          ! 1) jcb,K = k,jcb * k,K

          call DGEMM_('T','N',nBeta*ncb,nBasisi,nExpi,One,Array(ipF2),nExpi,Shells(iShll)%pCff,nExpi,Zero,Array(ipTmp),nBeta*ncb)

          ! 2)  j,cbK -> cbK,j

          call DgeTMo(Array(ipTmp),nBeta,nBeta,nBasisi*ncb,Array(ipF2),nBasisi*ncb)

          ! 3) bKj,C = c,bKj * c,C

          call DGEMM_('T','N',nTri_Elem1(lb)*nBasisi*nBeta,(2*iAng+1),nTri_Elem1(iAng),One,Array(ipF2),nTri_Elem1(iAng), &
                      RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipTmp),nTri_Elem1(lb)*nBasisi*nBeta)

          ! 4) b,KjC -> KjC,b

          call DgeTMo(Array(ipTmp),nTri_Elem1(lb),nTri_Elem1(lb),nBasisi*nBeta*(2*iAng+1),Array(ipF2),nBasisi*nBeta*(2*iAng+1))

          ! Next Contract (iKaC)*(KjCb) over K and C, producing ijab,
          ! by the following procedure:
          ! Loop over a and b
          !   Loop over C
          !     Contract iK(aC)*Kj(Cb), over K producing ij(aCb),
          !       accumulate to ij(ab)
          !   End loop C
          ! End Loop b and a

          do ib=1,nTri_Elem1(lb)
            do ia=1,nTri_Elem1(la)

              do iC=1,(2*iAng+1)
                iaC = (iC-1)*nTri_Elem1(la)+ia
                ipaC = (iaC-1)*nAlpha*nBasisi+ipF1
                iCb = (ib-1)*(2*iAng+1)+iC
                ipCb = (iCb-1)*nBasisi*nBeta+ipF2

                iIC = 0
                do iIrrep=0,nIrrep-1
                  if (.not. btest(llOper,iIrrep)) cycle
                  iIC = iIC+1
                  nOp = NrOpr(iDCRT(lDCRT))
                  Xg = real(iChTbl(iIrrep,nOp),kind=wp)
                  Factor = Xg*Fact
                  call DGEMM_('N','N',nAlpha,nBeta,nBasisi,Factor,Array(ipaC),nAlpha,Array(ipCb),nBasisi,One,rFinal(1,ia,ib,iIC), &
                              nAlpha)
                end do ! iIrrep

              end do ! iC
            end do   ! ia
          end do     ! ib
        end do ! iAng

      end do ! lDCRT
    end do ! iCnt
  end if
  mdc = mdc+dbsc(iCnttp)%nCntr
end do ! iCnttp

#ifdef _DEBUGPRINT_
write(u6,*) ' Result in PrjInt'
do ia=1,(la+1)*(la+2)/2
  do ib=1,(lb+1)*(lb+2)/2
    write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
    call RecPrt(Label,' ',rFinal(:,ia,ib,1),nAlpha,nBeta)
  end do
end do
#endif

return

end subroutine PrjInt
