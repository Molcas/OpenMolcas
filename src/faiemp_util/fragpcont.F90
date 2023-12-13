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
! Copyright (C) Ben Swerts                                             *
!               2016, Liviu Ungur                                      *
!***********************************************************************

subroutine FragpCont(F1,mi,mK,ma,mC,F2,mL,mj,mD,mb,W,rFinal,Factor)
!***********************************************************************
!                                                                      *
! Object: Specialized contraction of 3 4D matrices.                    *
!         It is very inefficient and should be rewritten using BLAS    *
!                                                                      *
! Called from: FragPInt and FragPGrd                                   *
!                                                                      *
!     Author: Ben Swerts                                               *
!   Modified: Liviu Ungur                                              *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mi, mK, ma, mC, mL, mj, mD, mb
real(kind=wp), intent(in) :: F1(mi,mK,ma,mC), F2(mL,mj,mD,mb), W(mK,mC,mL,mD), Factor
real(kind=wp), intent(inout) :: rFinal(mi,mj,ma,mb)
! local variables
integer(kind=iwp) :: ib, ia, ij, ii, iC, iK, iD, iL, j1, j2, j12
real(kind=wp) :: xt, xt2, WW1(mK*mC*mL*mD), F12(mK*mC*mL*mD)
logical(kind=iwp) :: DBG
real(kind=wp), external :: ddot_ !, dnrm2_

DBG = .false.

if (DBG) write(u6,'(A)') 'Enter FragpCont'
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mi=',mi
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mK=',mK
if (DBG) write(u6,'(A,i5)') 'FragpCont:  ma=',ma
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mC=',mC
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mL=',mL
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mj=',mj
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mD=',mD
if (DBG) write(u6,'(A,i5)') 'FragpCont:  mb=',mb
if (DBG) write(u6,'(A,F20.14)') 'FragpCont:  Factor=',Factor

if (DBG) write(u6,'(A,i5)') 'FragpCont: F1(mi,mK,ma,mC):'
if (DBG) then
  do ii=1,mi
    do iK=1,mK
      do ia=1,ma
        do iC=1,mC
          write(u6,'(4(a,i3),a,es22.14)') 'FragpCont: F1(',ii,',',iK,',',ia,',',iC,')=',F1(ii,iK,ia,iC)
        end do
      end do
    end do
  end do

  write(u6,'(A,i5)') 'FragpCont: F2(mL,mj,mD,mb):'
  do iL=1,mL
    do ij=1,mj
      do iD=1,mD
        do ib=1,mb
          write(u6,'(4(a,i3),a,es22.14)') 'FragpCont: F2(',iL,',',ij,',',iD,',',ib,')=',F2(iL,ij,iD,ib)
        end do
      end do
    end do
  end do

  write(u6,'(A,i5)') 'FragpCont: W(mK,mC,mL,mD):'
  do iK=1,mK
    do iC=1,mC
      do iL=1,mL
        do iD=1,mD
          write(u6,'(4(a,i3),a,es22.14)') 'FragpCont:  W(',iK,',',iC,',',iL,',',iD,')=',W(iK,iC,iL,iD)
        end do
      end do
    end do
  end do
end if !DBG

!write(u6,*) 'norm of FINAL in fragpcont',dnrm2_(mi*mj*ma*mb,rFinal(1:mi,1:mj,1:ma,1:mb),1)
!T = F1*W
!call DGEMM_('T','N',iBas*nElem(la)*nAlpha,iSize,nElem(iAng),One,Array(ipTmp),nElem(iAng),RSph(ipSph(iAng)),nElem(iAng),Zero, &
!            Array(ipF1),nAlpha*iBas*nElem(la))
!rFinal = T*F2
!call DGEMM_('T','N',iBas*nElem(la)*nAlpha,iSize,nElem(iAng),One,Array(ipTmp),nElem(iAng),RSph(ipSph(iAng)),nElem(iAng),Zero, &
!            Array(ipF1),nAlpha*iBas*nElem(la))
!write(u6,*) 'factor=',factor
do ib=1,mb
  do ia=1,ma
    do ij=1,mj
      do ii=1,mi
        xt = Zero
        xt2 = Zero
        WW1 = Zero
        F12 = Zero

        j1 = 0
        j12 = 0
        do iC=1,mC
          do iK=1,mK
            j1 = j1+1

            j2 = 0
            do iD=1,mD
              do iL=1,mL
                j2 = j2+1
                j12 = j12+1

                WW1(j12) = W(iK,iC,iL,iD)
                F12(j12) = F1(ii,iK,ia,iC)*F2(iL,ij,iD,ib)

                xt = xt+F1(ii,iK,ia,iC)*W(iK,iC,iL,iD)*F2(iL,ij,iD,ib)
              end do
            end do
          end do
        end do

        xt2 = xt2+ddot_(mC*mK*mD*mL,F12,1,WW1,1)

        !write(u6,*) 'factor, xt, xt2, diff,:',factor,xt,xt2,xt-xt2

        ! TMP
        !call dgemm_('N',m,n,One,F1(ii,:,ia,:),lda,W(:,:,iL,iD),incx,beta,TMP,incy)
        !
        !call daxpy_(n,a,x,incx,y,incy)

        rFinal(ii,ij,ia,ib) = rFinal(ii,ij,ia,ib)+Factor*xt

      end do !ii
    end do !ij
  end do !ia
end do !ib

return

end subroutine FragpCont
