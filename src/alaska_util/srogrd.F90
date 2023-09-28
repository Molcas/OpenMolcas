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
! Copyright (C) 1994,1995, Roland Lindh                                *
!               1994, Luis Seijo                                       *
!***********************************************************************

subroutine SROGrd( &
#                 define _CALLING_
#                 include "grd_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of MP integrals.          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, and Luis Seijo, Dept. of Applied Phys-  *
!             ical Chemistry, the Free University of Madrid, Spain,    *
!             September '94.                                           *
!                                                                      *
!             Modified to gradients April '95. R. Lindh                *
!***********************************************************************

use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Real_Spherical, only: ipSph, RSph
use Symmetry_Info, only: iOper
use Index_Functions, only: nTri_Elem1
use Disp, only: Dirct, IndDsp
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
integer(kind=iwp) :: i, ia, iaC, iAng, ib, iC, iCar, iCb, iCent, iCmp, iDCRT(0:7), iGamma, iIrrep, ip, ipA, ipaC, ipAxyz, ipB, &
                     ipBxyz, ipC, ipCb, ipCxyz, ipF1, ipF1a, ipF2, ipF2a, ipK1, ipK2, ipP1, ipP2, ipQ1, iPrint, ipRxyz, ipTmp, &
                     ipZ1, ipZ2, ipZI1, ipZI2, iRout, iShll, iStrt, iuvwx(4), iVec, j, JndGrd(3,4), kCnt, kCnttp, kdc, ld, lDCRT, &
                     LmbdT, lOp(4), mGrad, mVec, mVecAC, mVecCB, nac, ncb, nDAO, nDCRT, nDisp, nExpi, n_Her, nVecAC, nVecCB
real(kind=wp) :: C(3), Fact, TC(3)
character(len=80) :: Label
logical(kind=iwp) :: ABeq(3), EQ, JfGrad(3,4)
real(kind=wp), external :: DNrm2_
logical(kind=iwp), external :: TF
#include "print.fh"

#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(nHer)

iRout = 191
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In SROGrd: A',' ',A,1,3)
  call RecPrt(' In SROGrd: RB',' ',RB,1,3)
  call RecPrt(' In SROGrd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In SROGrd: P',' ',P,nZeta,3)
  write(u6,*) ' In SROGrd: la,lb=',' ',la,lb
end if

nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
iIrrep = 0
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = iOper(kOp(1))
lOp(2) = iOper(kOp(2))

kdc = 0
do kCnttp=1,nCnttp
  if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
  if (.not. dbsc(kCnttp)%ECP) cycle
  if (dbsc(kCnttp)%nSRO <= 0) cycle
  do kCnt=1,dbsc(kCnttp)%nCntr
    C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    call DCR(LmbdT,iStabM,nStabM,dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)

    iuvwx(3) = dc(kdc+kCnt)%nStab
    iuvwx(4) = dc(kdc+kCnt)%nStab
    JndGrd(:,1:2) = IndGrd(:,:)
    do i=1,3
      do j=1,2
        JfGrad(i,j) = IfGrad(i,j)
      end do
    end do

    nDisp = IndDsp(kdc+kCnt,iIrrep)
    do iCar=0,2
      JfGrad(iCar+1,3) = .false.
      iCmp = 2**iCar
      if (TF(kdc+kCnt,iIrrep,iCmp) .and. (.not. dbsc(kCnttp)%pChrg)) then
        nDisp = nDisp+1
        if (Dirct(nDisp)) then
          JndGrd(iCar+1,1) = abs(JndGrd(iCar+1,1))
          JndGrd(iCar+1,2) = abs(JndGrd(iCar+1,2))
          JndGrd(iCar+1,3) = -nDisp
          JfGrad(iCar+1,1) = .true.
          JfGrad(iCar+1,2) = .true.
        else
          JndGrd(iCar+1,3) = 0
        end if
      else
        JndGrd(iCar+1,3) = 0
      end if
    end do
    JndGrd(:,4) = 0
    JfGrad(1,4) = .false.
    JfGrad(2,4) = .false.
    JfGrad(3,4) = .false.
    mGrad = 0
    do iCar=1,3
      do i=1,2
        if (JfGrad(iCar,i)) mGrad = mGrad+1
      end do
    end do
    if (mGrad == 0) cycle

    do lDCRT=0,nDCRT-1
      lOp(3) = iDCRT(lDCRT)
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      if (EQ(A,RB) .and. EQ(A,TC)) cycle
      do iAng=0,dbsc(kCnttp)%nSRO-1
        iShll = dbsc(kCnttp)%iSRO+iAng
        nExpi = Shells(iShll)%nExp
        if (nExpi == 0) cycle

        ip = 1
        ipC = ip
        ip = ip+nExpi**2

        if (iPrint >= 49) call RecPrt(' The Akl matrix',' ',Shells(iShll)%Akl(1,1,1),nExpi,nExpi)
        call dcopy_(nExpi**2,Shells(iShll)%Akl(1,1,1),1,Array(ipC),1)
        if (EQ(A,RB) .and. EQ(A,TC) .and. dbsc(kCnttp)%NoPair) then
          call DaXpY_(nExpi**2,One,Shells(iShll)%Akl(1,1,2),1,Array(ipC),1)
          if (iPrint >= 49) call RecPrt(' The Adl matrix',' ',Shells(iShll)%Akl(1,1,2),nExpi,nExpi)
        end if

        ipF1 = ip
        nac = nTri_Elem1(la)*nTri_Elem1(iAng)*4
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
          write(u6,*) '  ip-1 > nArr*nZeta(1) in SROGrd'
          write(u6,*) ' nArr, nZeta=',nArr,nZeta
          write(u6,*) ' nac, nAlpha=',nac,nAlpha
          write(u6,*) ' nExpi=',nExpi
          call Abend()
        end if

        ! Calculate Effective center and exponent for <A|alm>

        call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,Alpha,Shells(iShll)%Exp)
        call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

        ! Calculate Overlap <A|core> and derivative <A'|core>

        n_Her = ((la+1)+iAng+2)/2
        ipAxyz = ip
        ip = ip+nAlpha*nExpi*3*n_Her*(la+2)
        ipCxyz = ip
        ip = ip+nAlpha*nExpi*3*n_Her*(iAng+1)
        ipRxyz = ip
        ip = ip+nAlpha*nExpi*3*n_Her*(nOrdOp+1)
        ipQ1 = ip
        ip = ip+nAlpha*nExpi*3*(la+2)*(iAng+1)*(nOrdOp+1)
        ipA = ip
        ip = ip+nAlpha*nExpi
        if (ip-1 > nArr*nZeta) then
          write(u6,*) '  ip-1 > nArr*nZeta(1b) in SROGrd'
          call Abend()
        end if
        ABeq(1) = A(1) == TC(1)
        ABeq(2) = A(2) == TC(2)
        ABeq(3) = A(3) == TC(3)
        call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,A,Array(ipAxyz),la+1,HerR(iHerR(n_Her)),n_Her,ABeq)
        call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,TC,Array(ipCxyz),iAng,HerR(iHerR(n_Her)),n_Her,ABeq)
        ABeq(1) = .false.
        ABeq(2) = .false.
        ABeq(3) = .false.
        call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(n_Her)),n_Her,ABeq)
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipAxyz)=',DNrm2_(nAlpha*nExpi*3*n_Her*(la+2),Array(ipAxyz),1)
          write(u6,*) ' Array(ipCxyz)=',DNrm2_(nAlpha*nExpi*3*n_Her*(iAng+1),Array(ipCxyz),1)
          write(u6,*) ' Array(ipRxyz)=',DNrm2_(nAlpha*nExpi*3*n_Her*(nOrdOp+1),Array(ipRxyz),1)
        end if
        call Assmbl(Array(ipQ1),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipCxyz),iAng,nAlpha*nExpi,HerW(iHerW(n_Her)),n_Her)
        iStrt = ipA
        do iGamma=1,nExpi
          call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
          iStrt = iStrt+nAlpha
        end do
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*nExpi,Array(ipA),1)
        end if
        call rKappa_Zeta(Array(ipK1),Array(ipZ1),nExpi*nAlpha)
        ld = 1
        call CmbnAC(Array(ipQ1),nAlpha*nExpi,la,iAng,Array(ipK1),Array(ipF1),Array(ipA),JfGrad(1,1),ld,nVecAC)
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipQ1)=',DNrm2_(nAlpha*nExpi*3*(la+2)*(iAng+1)*(nOrdOp+1),Array(ipQ1),1)
          write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*nExpi,Array(ipA),1)
        end if
        ip = ip-nAlpha*nExpi*(6+3*n_Her*(la+2)+3*n_Her*(iAng+1)+3*n_Her*(nOrdOp+1)+3*(la+2)*(iAng+1)*(nOrdOp+1)+1)

        ipF2 = ip
        ncb = nTri_Elem1(iAng)*nTri_Elem1(lb)*4
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
          write(u6,*) '  ip-1 > nArr*nZeta(2) in SROGrd'
          call Abend()
        end if

        ! Calculate Effective center and exponent for <blm|B>

        call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,Shells(iShll)%Exp,Beta)
        call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

        ! Calculate Overlap <core|B> and <core|B'>

        n_Her = ((iAng+1)+lb+2)/2
        ipCxyz = ip
        ip = ip+nBeta*nExpi*3*n_Her*(iAng+1)
        ipBxyz = ip
        ip = ip+nBeta*nExpi*3*n_Her*(lb+2)
        ipRxyz = ip
        ip = ip+nBeta*nExpi*3*n_Her*(nOrdOp+1)
        ipQ1 = ip
        ip = ip+nBeta*nExpi*3*(iAng+1)*(lb+2)*(nOrdOp+1)
        ipB = ip
        ip = ip+nBeta*nExpi
        if (ip-1 > nArr*nZeta) then
          write(u6,*) '  ip-1 > nArr*nZeta(2b) in SROGrd'
          call Abend()
        end if
        ABeq(1) = TC(1) == RB(1)
        ABeq(2) = TC(2) == RB(2)
        ABeq(3) = TC(3) == RB(3)
        call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,TC,Array(ipCxyz),iAng,HerR(iHerR(n_Her)),n_Her,ABeq)
        call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(n_Her)),n_Her,ABeq)
        ABeq(1) = .false.
        ABeq(2) = .false.
        ABeq(3) = .false.
        call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(n_Her)),n_Her,ABeq)
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipCxyz)=',DNrm2_(nBeta*nExpi*3*n_Her*(iAng+1),Array(ipCxyz),1)
          write(u6,*) ' Array(ipBxyz)=',DNrm2_(nBeta*nExpi*3*n_Her*(lb+2),Array(ipBxyz),1)
          write(u6,*) ' Array(ipRxyz)=',DNrm2_(nBeta*nExpi*3*n_Her*(nOrdOp+1),Array(ipRxyz),1)
        end if
        call Assmbl(Array(ipQ1),Array(ipCxyz),iAng,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,nExpi*nBeta,HerW(iHerW(n_Her)),n_Her)
        iStrt = ipB
        do iGamma=1,nExpi
          call dcopy_(nBeta,Beta,1,Array(iStrt),nExpi)
          iStrt = iStrt+1
        end do
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipB)=',DNrm2_(nExpi*nBeta,Array(ipB),1)
        end if
        call rKappa_Zeta(Array(ipK2),Array(ipZ2),nExpi*nBeta)
        ld = 1
        call CmbnCB(Array(ipQ1),nExpi*nBeta,iAng,lb,Array(ipK2),Array(ipF2),Array(ipB),JfGrad(1,2),ld,nVecCB)
        if (iPrint >= 49) then
          write(u6,*) ' Array(ipQ1)=',DNrm2_(nExpi*nBeta*3*(la+2)*(iAng+1)*(nOrdOp+1),Array(ipQ1),1)
          write(u6,*) ' Array(ipB)=',DNrm2_(nExpi*nBeta,Array(ipB),1)
        end if
        ip = ip-nBeta*nExpi*(6+3*n_Her*(lb+2)+3*n_Her*(iAng+1)+3*n_Her*(nOrdOp+1)+3*(lb+2)*(iAng+1)*(nOrdOp+1)+1)
        nac = nTri_Elem1(la)*nTri_Elem1(iAng)*nVecAC
        ncb = nTri_Elem1(iAng)*nTri_Elem1(lb)*nVecCB
        ipTmp = ip
        ip = ip+max(nAlpha*nExpi*nac,nExpi*nBeta*ncb)
        if (ip-1 > nArr*nZeta) then
          write(u6,*) '  ip-1 > nArr*nZeta(3) in SROGrd'
          call Abend()
        end if
        nac = nTri_Elem1(la)*nTri_Elem1(iAng)
        ncb = nTri_Elem1(iAng)*nTri_Elem1(lb)

        ! Calculate Contraction over the spectral resolvent basis
        ! set of the type <A|alm>A(l;ab)<blm|B> where we now have in
        ! Array(ipF1) the cartesian components of <A|alm>, and
        ! similarily, in Array(ipF2), we have stored the cartesian
        ! components of <alm|B>. Observe that as opposed to the
        ! projection operator that this contraction is done in the
        ! primitive basis.

        ! From the lefthandside overlap, form ikaCx from ikacx by
        ! 1) ika,cx -> cx,ika

        call DgeTMo(Array(ipF1),nAlpha*nExpi*nTri_Elem1(la),nAlpha*nExpi*nTri_Elem1(la),nTri_Elem1(iAng)*nVecAC,Array(ipTmp), &
                    nTri_Elem1(iAng)*nVecAC)

        ! 2) xika,C = c,xika * c,C

        call DGEMM_('T','N',nVecAC*nAlpha*nExpi*nTri_Elem1(la),(2*iAng+1),nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                    RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF1),nVecAC*nAlpha*nExpi*nTri_Elem1(la))

        ! 3) x,ikaC -> ikaC,x

        call DGetMo(Array(ipF1),nVecAC,nVecAC,nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1),Array(ipTmp), &
                    nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1))
        call dcopy_(nVecAC*nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1),Array(ipTmp),1,Array(ipF1),1)

        ! And (almost) the same thing for the righthand side, form
        ! kjCbx from kjcbx
        ! ) kj,cbx -> cbx,kj

        call DgeTMo(Array(ipF2),nBeta*nExpi,nBeta*nExpi,ncb*nVecCB,Array(ipTmp),ncb*nVecCB)

        ! 2) bxkj,C = c,bxkj * c,C

        call DGEMM_('T','N',nTri_Elem1(lb)*nVecCB*nExpi*nBeta,(2*iAng+1),nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                    RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF2),nTri_Elem1(lb)*nVecCB*nExpi*nBeta)

        ! 3) bx,kjC -> kjC,bx

        call DgeTMo(Array(ipF2),nTri_Elem1(lb)*nVecCB,nTri_Elem1(lb)*nVecCB,nExpi*nBeta*(2*iAng+1),Array(ipTmp), &
                    nExpi*nBeta*(2*iAng+1))
        call dcopy_(nExpi*nBeta*(2*iAng+1)*nTri_Elem1(lb)*nVecCB,Array(ipTmp),1,Array(ipF2),1)

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

        rFinal(:,:,:,1,:) = Zero

        mVec = 0
        mVecAC = 1
        mVecCB = 1
        do iCar=1,3
          do iCent=1,2
            if (JfGrad(iCar,iCent)) then
              mVec = mVec+1
              if (iCent == 1) then
                mVecAC = mVecAC+1
                ipF1a = ipF1+(mVecAC-1)*nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1)
                ipF2a = ipF2
              else
                ipF1a = ipF1
                mVecCB = mVecCB+1
                ipF2a = ipF2+(mVecCB-1)*nExpi*nBeta*(2*iAng+1)*nTri_Elem1(lb)
              end if

              do ib=1,nTri_Elem1(lb)
                do ia=1,nTri_Elem1(la)
                  if (iPrint >= 99) write(u6,*) ' ia,ib=',ia,ib

                  do iC=1,(2*iAng+1)
                    if (iPrint >= 99) write(u6,*) ' iC,=',iC
                    iaC = (iC-1)*nTri_Elem1(la)+ia
                    ipaC = (iaC-1)*nAlpha*nExpi+ipF1a
                    iCb = (ib-1)*(2*iAng+1)+iC
                    ipCb = (iCb-1)*nExpi*nBeta+ipF2a

                    if (iPrint >= 99) then
                      call RecPrt('<ia|iC>',' ',Array(ipaC),nAlpha,nExpi)
                      call RecPrt('<iC|ib>',' ',Array(ipCb),nExpi,nBeta)
                    end if

                    call DGEMM_('N','N',nAlpha,nExpi,nExpi,One,Array(ipaC),nAlpha,Array(ipC),nExpi,Zero,Array(ipTmp),nAlpha)
                    call DGEMM_('N','N',nAlpha,nBeta,nExpi,Fact,Array(ipTmp),nAlpha,Array(ipCb),nExpi,One,rFinal(:,ia,ib,1,mVec), &
                                nAlpha)

                  end do
                end do
              end do

            end if
          end do
        end do

        if (iPrint >= 49) then
          do iVec=1,mVec
            write(u6,*) iVec,sqrt(DNrm2_(nZeta*nTri_Elem1(la)*nTri_Elem1(lb),rFinal(:,:,:,1,iVec),1))
          end do
        end if
        if (iPrint >= 99) then
          write(u6,*) ' Result in SROGrd'
          do ia=1,nTri_Elem1(la)
            do ib=1,nTri_Elem1(lb)
              do iVec=1,mVec
                write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
                call RecPrt(Label,' ',rFinal(:,ia,ib,1,iVec),nAlpha,nBeta)
              end do
            end do
          end do
        end if

        ! Distribute contributions to the gradient

        call Distg1X(rFinal,DAO,nZeta,nDAO,mVec,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

      end do
    end do
  end do
end do

return

end subroutine SROGrd
