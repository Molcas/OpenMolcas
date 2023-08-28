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
! Copyright (C) 1994, Roland Lindh                                     *
!               1994, Luis Seijo                                       *
!***********************************************************************

subroutine SROInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of MP integrals.          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Luis Seijo, Dept. of Applied Phys-  *
!             ical Chemistry, the Free University of Madrid, Spain,    *
!             September '94.                                           *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Real_Spherical, only: ipSph, RSph
use Symmetry_Info, only: iChTbl, nIrrep
use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, iaC, iAng, ib, iC, iCb, iCnt, iCnttp, iComp, iDCRT(0:7), iIC, iIrrep, ip, ipaC, ipC, ipCb, ipF1, ipF2, &
                     ipK1, ipK2, ipP1, ipP2, iPrint, ipTmp, ipZ1, ipZ2, ipZI1, ipZI2, iRout, iShll, l, lDCRT, llOper, LmbdT, mArr, &
                     mdc, nac, ncb, nDCRT, nExpi, nOp
real(kind=wp) :: C(3), Fact, Factor, TC(3), Xg
character(len=80) :: Label
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(iChO)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 191
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In SROInt: A',' ',A,1,3)
  call RecPrt(' In SROInt: RB',' ',RB,1,3)
  call RecPrt(' In SROInt: CoorO',' ',CoorO,1,3)
  call RecPrt(' In SROInt: P',' ',P,nZeta,3)
  write(u6,*) ' In SROInt: la,lb=',' ',la,lb
end if

rFinal(:,:,:,:) = Zero

llOper = lOper(1)
iComp = 1
mdc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%ECP .and. (dbsc(iCnttp)%nSRO > 0)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

      call DCR(LmbdT,iStabM,nStabM,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,iDCRT,nDCRT)
      Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

      do lDCRT=0,nDCRT-1
        call OA(iDCRT(lDCRT),C,TC)
        do iAng=0,dbsc(iCnttp)%nSRO-1
          iShll = dbsc(iCnttp)%iSRO+iAng
          nExpi = Shells(iShll)%nExp
          if (nExpi == 0) cycle

          ip = 1
          ipC = ip
          ip = ip+nExpi**2

          if (iPrint >= 49) call RecPrt(' The Akl matrix',' ',Shells(iShll)%Akl(:,:,1),nExpi,nExpi)
          Array(ipC:ipC+nExpi**2-1) = pack(Shells(iShll)%Akl(:,:,1),.true.)
          if (EQ(A,RB) .and. EQ(A,TC) .and. dbsc(iCnttp)%NoPair) then
            if (iPrint >= 49) call RecPrt(' The Adl matrix',' ',Shells(iShll)%Akl(:,:,2),nExpi,nExpi)
            Array(ipC:ipC+nExpi**2-1) = Array(ipC:ipC+nExpi**2-1)+pack(Shells(iShll)%Akl(:,:,2),.true.)
          end if

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
            call WarningMessage(2,'SROInt: ip-1 > nArr*nZeta(1)')
            write(u6,*) ' nArr, nZeta=',nArr,nZeta
            write(u6,*) ' nac, nAlpha=',nac,nAlpha
            write(u6,*) ' nExpi=',nExpi
            call Abend()
          end if
          mArr = (nArr*nZeta-(ip-1))/nZeta

          ! Calculate Effective center and exponent for <A|alm>

          call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,Alpha,Shells(iShll)%Exp)
          call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

          ! Calculate Overlap <A|alm>

          nHer = (la+iAng+2)/2
          call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,Array(ipZ1),Array(ipZI1),Array(ipK1),Array(ipP1),Array(ipF1), &
                      nAlpha*nExpi,iComp,la,iAng,A,TC,nHer,Array(ip),mArr,CoorO,nOrdOp)
          if (iPrint >= 99) call RecPrt('<a|srbs>',' ',Array(ipF1),nAlpha*nExpi,nac)
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
            call WarningMessage(2,'SROInt: ip-1 > nArr*nZeta(2)')
            call Abend()
          end if
          mArr = (nArr*nZeta-(ip-1))/nZeta

          ! Calculate Effective center and exponent for <blm|B>

          call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,Shells(iShll)%Exp,Beta)
          call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

          ! Calculate Overlap <blm|B>

          nHer = (iAng+lb+2)/2
          call MltPrm(Shells(iShll)%Exp,nExpi,Beta,nBeta,Array(ipZ2),Array(ipZI2),Array(ipK2),Array(ipP2),Array(ipF2),nExpi*nBeta, &
                      iComp,iAng,lb,TC,RB,nHer,Array(ip),mArr,CoorO,nOrdOp)
          if (iPrint >= 99) call RecPrt('<srbs|b>',' ',Array(ipF2),nExpi*nBeta,ncb)
          ip = ip-6*nExpi*nBeta
          ipTmp = ip
          ip = ip+max(nAlpha*nExpi*nac,nExpi*nBeta*ncb)
          if (ip-1 > nArr*nZeta) then
            call WarningMessage(2,'SROInt: ip-1 > nArr*nZeta(3)')
            call Abend()
          end if

          ! Calculate Contraction over the spectral resolvent basis
          ! set of the type <A|alm>A(l;ab)<blm|B> where we now have in
          ! Array(ipF1) the cartesian components of <A|alm>, and
          ! similarily, in Array(ipF2), we have stored the cartesian
          ! components of <alm|B>. Observe that as opposed to the
          ! projection operator that this contraction is done in the
          ! primitive basis.

          ! From the lefthandside overlap, form ikaC from ikac by
          ! 1) ika,c -> c,ika

          call DgeTMo(Array(ipF1),nAlpha*nExpi*nTri_Elem1(la),nAlpha*nExpi*nTri_Elem1(la),nTri_Elem1(iAng),Array(ipTmp), &
                      nTri_Elem1(iAng))

          ! 2) ika,C = c,ika * c,C

          call DGEMM_('T','N',nAlpha*nExpi*nTri_Elem1(la),(2*iAng+1),nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                      RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF1),nAlpha*nExpi*nTri_Elem1(la))
          if (iPrint >= 99) call RecPrt('<A|srbs>',' ',Array(ipF1),nAlpha*nExpi,nTri_Elem1(la)*(2*iAng+1))

          ! And (almost) the same thing for the righthand side, form
          ! kjCb from kjcb
          ! 1) kj,cb -> cb,kj

          call DgeTMo(Array(ipF2),nBeta*nExpi,nBeta*nExpi,nTri_Elem1(iAng)*nTri_Elem1(lb),Array(ipTmp), &
                      nTri_Elem1(iAng)*nTri_Elem1(lb))

          ! 2) bkj,C = c,bkj * c,C

          call DGEMM_('T','N',nTri_Elem1(lb)*nExpi*nBeta,(2*iAng+1),nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                      RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF2),nTri_Elem1(lb)*nExpi*nBeta)

          ! 3) b,kjC -> kjC,b

          call DgeTMo(Array(ipF2),nTri_Elem1(lb),nTri_Elem1(lb),nExpi*nBeta*(2*iAng+1),Array(ipTmp),nExpi*nBeta*(2*iAng+1))
          l = nExpi*nBeta*(2*iAng+1)*nTri_Elem1(lb)
          Array(ipF2:ipF2+l-1) = Array(ipTmp:ipTmp+l-1)
          if (iPrint >= 99) call RecPrt('<srbs|B>',' ',Array(ipF2),nExpi*nBeta,(2*iAng+1)*nTri_Elem1(lb))

          ! Next Contract (ikaC)*(klC)*(ljCb) over k,l and C,
          ! producing ijab,
          ! by the following procedure:
          ! Loop over a and b
          !   Loop over C
          !     Contract ik(aC)*kl(C), over k producing il(aC),
          !     Contract il(aC)*lj(Cb), over l producing ij(aCb)
          !       accumulate to ij(ab)
          !   End loop C
          ! End Loop b and a

          do ib=1,nTri_Elem1(lb)
            do ia=1,nTri_Elem1(la)
              if (iPrint >= 99) write(u6,*) ' ia,ib=',ia,ib

              do iC=1,(2*iAng+1)
                if (iPrint >= 99) write(u6,*) ' iC,=',iC
                iaC = (iC-1)*nTri_Elem1(la)+ia
                ipaC = (iaC-1)*nAlpha*nExpi+ipF1
                iCb = (ib-1)*(2*iAng+1)+iC
                ipCb = (iCb-1)*nExpi*nBeta+ipF2

                iIC = 0
                if (iPrint >= 99) then
                  call RecPrt('<ia|iC>',' ',Array(ipaC),nAlpha,nExpi)
                  call RecPrt('<iC|ib>',' ',Array(ipCb),nExpi,nBeta)
                end if
                do iIrrep=0,nIrrep-1
                  if (.not. btest(llOper,iIrrep)) cycle
                  if (iPrint >= 99) write(u6,*) ' iIC=',iIC
                  iIC = iIC+1
                  nOp = NrOpr(iDCRT(lDCRT))
                  Xg = real(iChTbl(iIrrep,nOp),kind=wp)
                  Factor = Xg*Fact
                  call DGEMM_('N','N',nAlpha,nExpi,nExpi,One,Array(ipaC),nAlpha,Array(ipC),nExpi,Zero,Array(ipTmp),nAlpha)
                  call DGEMM_('N','N',nAlpha,nBeta,nExpi,Factor,Array(ipTmp),nAlpha,Array(ipCb),nExpi,One,rFinal(1,ia,ib,iIC), &
                              nAlpha)
                end do

              end do
            end do
          end do

        end do
      end do
    end do
  end if
  mdc = mdc+dbsc(iCnttp)%nCntr
end do

if (iPrint >= 99) then
  write(u6,*) ' Result in SROInt'
  do ia=1,(la+1)*(la+2)/2
    do ib=1,(lb+1)*(lb+2)/2
      write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
      call RecPrt(Label,' ',rFinal(:,ia,ib,1),nAlpha,nBeta)
    end do
  end do
end if

return

end subroutine SROInt
