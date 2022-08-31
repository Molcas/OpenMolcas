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

subroutine MkFck(iAnga,iCmp,Shijij,iShll,iShell,iBasi,jBasj,kBask,lBasl,iAO,iAOst,nOp,jOp,Dij,mDij,nDij,ij1,ij2,ij3,ij4,Dkl,mDkl, &
                 nDkl,kl1,kl2,kl3,kl4,Dik,mDik,nDik,ik1,ik2,ik3,ik4,Dil,mDil,nDil,il1,il2,il3,il4,Djk,mDjk,nDjk,jk1,jk2,jk3,jk4, &
                 Djl,mDjl,nDjl,jl1,jl2,jl3,jl4,AOInt,nAO,TwoHam,nFock,Scrtch2,nS2,FckTmp,nFT,pert,iuvwx,iCent,iCar,indgrd,ipDisp)
!***********************************************************************
!                                                                      *
! Object: Driver for the generation of the two electron contribution   *
!         to the Fock Matrix directly from the two electron integrals. *
!                                                                      *
!     Author:  Anders Bernhardsson 1995                                *
!***********************************************************************

use Symmetry_Info, only: nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), iCmp(4), iShll(4), iShell(4), iBasi, jBasj, kBask, lBasl, iAO(4), iAOst(4), nOp(4), &
                                 jOp(6), mDij, nDij, ij1, ij2, ij3, ij4, mDkl, nDkl, kl1, kl2, kl3, kl4, mDik, nDik, ik1, ik2, &
                                 ik3, ik4, mDil, nDil, il1, il2, il3, il4, mDjk, nDjk, jk1, jk2, jk3, jk4, mDjl, nDjl, jl1, jl2, &
                                 jl3, jl4, nAO, nFock, nS2, nFT, iuvwx, iCent, iCar, indgrd(3,4,0:7), ipdisp(*)
logical(kind=iwp), intent(in) :: Shijij, pert(0:7)
real(kind=wp), intent(in) :: Dij(mDij,nDij), Dkl(mDkl,nDkl), Dik(mDik,nDik), Dil(mDil,nDil), Djk(mDjk,nDjk), Djl(mDjl,nDjl), &
                             AOInt(nAO)
real(kind=wp), intent(inout) :: TwoHam(nFock)
real(kind=wp), intent(out) :: Scrtch2(nS2), FckTmp(nFT)
integer(kind=iwp) :: nijkl
real(kind=wp) :: Fact

! Just to make a nice interface

!iRout = 12
!iPrint = nPrint(iRout)
nijkl = iBasi*jBasj*kBask*lBasl

! Accumulate contributions directly to the symmetry adapted Fock matrix.

Fact = real(iuvwx,kind=wp)/real(nIrrep,kind=wp)

call FckAcc_mck(iAnga,iCmp(1),iCmp(2),iCmp(3),iCmp(4),Shijij,iShll,iShell,nOp,nijkl,AOInt,TwoHam,nFock,Scrtch2,nS2,iAO,iAOst, &
                iBasi,jBasj,kBask,lBasl,Dij(1,jOp(1)),ij1,ij2,ij3,ij4,Dkl(1,jOp(2)),kl1,kl2,kl3,kl4,Dik(1,jOp(3)),ik1,ik2,ik3,ik4, &
                Dil(1,jOp(4)),il1,il2,il3,il4,Djk(1,jOp(5)),jk1,jk2,jk3,jk4,Djl(1,jOp(6)),jl1,jl2,jl3,jl4,FckTmp,nFT,fact,iCar, &
                iCent,pert,indgrd,ipdisp)

return

end subroutine MkFck
