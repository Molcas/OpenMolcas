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
! Copyright (C) 1991,2001, Roland Lindh                                *
!***********************************************************************

subroutine PCMInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of nuclear attraction     *
!         integrals.                                                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden, January '91                             *
!                                                                      *
!             Modified to PCM-integrals, by RL June '01, Napoli, Italy.*
!***********************************************************************

use PCM_arrays, only: C_Tessera, nTiles, q_Tessera
use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Rys_interfaces, only: cff2d_kernel, modu2_kernel, rys2d_kernel, tval_kernel
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Symmetry_Info, only: ChOper
use Definitions, only: u6
#endif

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: iAnga(4), iDCRT(0:7), ipIn, iTile, jStab_(0:0), lDCRT, LmbdT, mabMax, mabMin, nDCRT, nFlop, nMem, nOp, NrOpr, &
                     nStab_, nT
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ii
#endif
real(kind=wp) :: C(3), Coora(3,4), CoorAC(3,2), Coori(3,4), Fact, qTessera, TC(3)
logical(kind=iwp) :: NoSpecial
procedure(cff2d_kernel) :: XCff2D
procedure(modu2_kernel) :: Fake
procedure(rys2d_kernel) :: XRys2D
procedure(tval_kernel) :: TNAI
logical(kind=iwp), external :: EQ

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(nHer)
unused_var(CoorO)
unused_var(nOrdOp)
unused_var(PtChrg)
unused_var(iAddPot)

rFinal(:,:,:,:) = Zero

iAnga(1) = la
iAnga(2) = lb
iAnga(3) = 0
iAnga(4) = 0
Coora(:,1) = A(:)
Coora(:,2) = RB(:)
Coori(:,1:2) = Coora(:,1:2)
mabMin = nTri3_Elem1(max(la,lb)-1)
if (EQ(A,RB)) mabMin = nTri3_Elem1(la+lb-1)
mabMax = nTri3_Elem1(la+lb)-1

! Compute FLOP's and size of work array which Hrr will use.

call mHrr(la,lb,nFLOP,nMem)

! Find center to accumulate angular momentum on. (HRR)

if (la >= lb) then
  CoorAC(:,1) = A(:)
else
  CoorAC(:,1) = RB(:)
end if

! The coordinates of the individual tiles are not stabilized by any
! operator but the unit operator.

nStab_ = 1
jStab_ = 0

! Loop over tiles.

do iTile=1,nTiles
  QTessera = Q_Tessera(iTile)
  C(:) = C_Tessera(:,iTile)
# ifdef _DEBUGPRINT_
  call RecPrt('C',' ',C,1,3)
# endif

  ! Find the DCR for M and S

  call DCR(LmbdT,iStabM,nStabM,jStab_,nStab_,iDCRT,nDCRT)
  Fact = One/real(LmbdT,kind=wp)

# ifdef _DEBUGPRINT_
  write(u6,*) ' m      =',nStabM
  write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
  write(u6,*) ' s      =',nStab_
  write(u6,'(9A)') '(S)=',ChOper(jStab_)
  write(u6,*) ' LambdaT=',LmbdT
  write(u6,*) ' t      =',nDCRT
  write(u6,'(9A)') '(T)=',(ChOper(iDCRT(ii)),ii=0,nDCRT-1)
# endif

  do lDCRT=0,nDCRT-1
    call OA(iDCRT(lDCRT),C,TC)
    CoorAC(:,2) = TC(:)
    Coori(:,3) = TC(:)
    Coori(:,4) = TC(:)
    Coora(:,3) = TC(:)
    Coora(:,4) = TC(:)

    ! Compute integrals with the Rys quadrature.

    nT = nZeta
    NoSpecial = .true.
    call Rys(iAnga,nT,Zeta,ZInv,nZeta,[One],[One],1,P,nZeta,TC,1,rKappa,[One],Coori,Coora,CoorAC,mabmin,mabmax,0,0,Array, &
             nArr*nZeta,TNAI,Fake,XCff2D,XRys2D,NoSpecial)

    ! Use the HRR to compute the required primitive integrals.

    call HRR(la,lb,A,RB,Array,nZeta,nMem,ipIn)

    ! Accumulate contributions to the symmetry adapted operator

    nOp = NrOpr(iDCRT(lDCRT))
    call SymAdO(Array(ipIn),nZeta,la,lb,nComp,rFinal,nIC,nOp,lOper,iChO,-Fact*QTessera)
#   ifdef _DEBUGPRINT_
    write(u6,*) Fact*QTessera
    call RecPrt('PCMInt: Array(ipIn)',' ',Array(ipIn),nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nComp)
    call RecPrt('PCMInt: rFinal',' ',rFinal,nZeta,nTri_Elem1(la)*nTri_Elem1(lb)*nIC)
#   endif

  end do
end do

end subroutine PCMInt
