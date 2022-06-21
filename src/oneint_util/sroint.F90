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

use Basis_Info
use Center_Info
use Real_Spherical
use Symmetry_Info, only: nIrrep, iChTbl

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
real*8 C(3), TC(3)
integer iDCRT(0:7), iTwoj(0:7)
character*80 Label
logical EQ
data iTwoj/1,2,4,8,16,32,64,128/
! Statement function for Cartesian index
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 191
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In SROInt: A',' ',A,1,3)
  call RecPrt(' In SROInt: RB',' ',RB,1,3)
  call RecPrt(' In SROInt: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In SROInt: P',' ',P,nZeta,3)
  write(6,*) ' In SROInt: la,lb=',' ',la,lb
end if

call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

llOper = lOper(1)
iComp = 1
mdc = 0
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%ECP) .or. (dbsc(iCnttp)%nSRO <= 0)) then
    mdc = mdc+dbsc(iCnttp)%nCntr
    cycle
  end if
  do iCnt=1,dbsc(iCnttp)%nCntr
    C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,iDCRT,nDCRT)
    Fact = dble(nStabM)/dble(LmbdT)

    do lDCRT=0,nDCRT-1
      call OA(iDCRT(lDCRT),C,TC)
      do iAng=0,dbsc(iCnttp)%nSRO-1
        iShll = dbsc(iCnttp)%iSRO+iAng
        nExpi = Shells(iShll)%nExp
        if (nExpi == 0) cycle

        ip = 1
        ipC = ip
        ip = ip+nExpi**2

        if (iPrint >= 49) call RecPrt(' The Akl matrix',' ',Shells(iShll)%Akl(1,1,1),nExpi,nExpi)
        call dcopy_(nExpi**2,Shells(iShll)%Akl(1,1,1),1,Array(ipC),1)
        if (EQ(A,RB) .and. EQ(A,TC) .and. dbsc(iCnttp)%NoPair) then
          if (iPrint >= 49) call RecPrt(' The Adl matrix',' ',Shells(iShll)%Akl(1,1,2),nExpi,nExpi)
          call DaXpY_(nExpi**2,One,Shells(iShll)%Akl(1,1,2),1,Array(ipC),1)
        end if

        ipF1 = ip
        nac = nElem(la)*nElem(iAng)
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
          write(6,*) ' nArr, nZeta=',nArr,nZeta
          write(6,*) ' nac, nAlpha=',nac,nAlpha
          write(6,*) ' nExpi=',nExpi
          call Abend()
        end if
        mArr = (nArr*nZeta-(ip-1))/nZeta

        ! Calculate Effective center and exponent for <A|alm>

        call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,Alpha,Shells(iShll)%Exp)
        call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

        ! Calculate Overlap <A|alm>

        nHer = (la+iAng+2)/2
        call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,Array(ipZ1),Array(ipZI1),Array(ipK1),Array(ipP1),Array(ipF1), &
                    nAlpha*nExpi,iComp,la,iAng,A,TC,nHer,Array(ip),mArr,CCoor,nOrdOp)
        if (iPrint >= 99) call RecPrt('<a|srbs>',' ',Array(ipF1),nAlpha*nExpi,nac)
        ip = ip-6*nAlpha*nExpi

        ipF2 = ip
        ncb = nElem(iAng)*nElem(lb)
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
                    iComp,iAng,lb,TC,RB,nHer,Array(ip),mArr,CCoor,nOrdOp)
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

        call DgeTMo(Array(ipF1),nAlpha*nExpi*nElem(la),nAlpha*nExpi*nElem(la),nElem(iAng),Array(ipTmp),nElem(iAng))

        ! 2) ika,C = c,ika * c,C

        call DGEMM_('T','N',nAlpha*nExpi*nElem(la),(2*iAng+1),nElem(iAng),1.0d0,Array(ipTmp),nElem(iAng),RSph(ipSph(iAng)), &
                    nElem(iAng),0.0d0,Array(ipF1),nAlpha*nExpi*nElem(la))
        if (iPrint >= 99) call RecPrt('<A|srbs>',' ',Array(ipF1),nAlpha*nExpi,nElem(la)*(2*iAng+1))

        ! And (almost) the same thing for the righthand side, form
        ! kjCb from kjcb
        ! 1) kj,cb -> cb,kj

        call DgeTMo(Array(ipF2),nBeta*nExpi,nBeta*nExpi,nElem(iAng)*nElem(lb),Array(ipTmp),nElem(iAng)*nElem(lb))

        ! 2) bkj,C = c,bkj * c,C

        call DGEMM_('T','N',nElem(lb)*nExpi*nBeta,(2*iAng+1),nElem(iAng),1.0d0,Array(ipTmp),nElem(iAng),RSph(ipSph(iAng)), &
                    nElem(iAng),0.0d0,Array(ipF2),nElem(lb)*nExpi*nBeta)

        ! 3) b,kjC -> kjC,b

        call DgeTMo(Array(ipF2),nElem(lb),nElem(lb),nExpi*nBeta*(2*iAng+1),Array(ipTmp),nExpi*nBeta*(2*iAng+1))
        call dcopy_(nExpi*nBeta*(2*iAng+1)*nElem(lb),Array(ipTmp),1,Array(ipF2),1)
        if (iPrint >= 99) call RecPrt('<srbs|B>',' ',Array(ipF2),nExpi*nBeta,(2*iAng+1)*nElem(lb))

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

        do ib=1,nElem(lb)
          do ia=1,nElem(la)
            if (iPrint >= 99) write(6,*) ' ia,ib=',ia,ib

            do iC=1,(2*iAng+1)
              if (iPrint >= 99) write(6,*) ' iC,=',iC
              iaC = (iC-1)*nElem(la)+ia
              ipaC = (iaC-1)*nAlpha*nExpi+ipF1
              iCb = (ib-1)*(2*iAng+1)+iC
              ipCb = (iCb-1)*nExpi*nBeta+ipF2

              iIC = 0
              if (iPrint >= 99) then
                call RecPrt('<ia|iC>',' ',Array(ipaC),nAlpha,nExpi)
                call RecPrt('<iC|ib>',' ',Array(ipCb),nExpi,nBeta)
              end if
              do iIrrep=0,nIrrep-1
                if (iand(llOper,iTwoj(iIrrep)) == 0) cycle
                if (iPrint >= 99) write(6,*) ' iIC=',iIC
                iIC = iIC+1
                nOp = NrOpr(iDCRT(lDCRT))
                Xg = dble(iChTbl(iIrrep,nOp))
                Factor = Xg*Fact
                call DGEMM_('N','N',nAlpha,nExpi,nExpi,One,Array(ipaC),nAlpha,Array(ipC),nExpi,Zero,Array(ipTmp),nAlpha)
                call DGEMM_('N','N',nAlpha,nBeta,nExpi,Factor,Array(ipTmp),nAlpha,Array(ipCb),nExpi,One,final(1,ia,ib,iIC),nAlpha)
              end do

            end do
          end do
        end do

      end do
    end do
  end do
  mdc = mdc+dbsc(iCnttp)%nCntr
end do

if (iPrint >= 99) then
  write(6,*) ' Result in SROInt'
  do ia=1,(la+1)*(la+2)/2
    do ib=1,(lb+1)*(lb+2)/2
      write(Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
      call RecPrt(Label,' ',final(1,ia,ib,1),nAlpha,nBeta)
    end do
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Zeta)
  call Unused_real_array(ZInv)
  call Unused_real_array(rKappa)
  call Unused_integer(nRys)
  call Unused_integer_array(iChO)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine SROInt
