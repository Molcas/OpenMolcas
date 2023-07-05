!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ProjSym(nCent,Ind,A,iDCRs,B,dB,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,iq,rMult)

use Slapaf_Info, only: jStab, nStab, Smmtrc

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "warnings.h"
#include "real.fh"
real*8 Tx(3,MxAtom), A(3,nCent), B(3,nCent), ATemp(3), dB(3,nCent,3,nCent), BM(nB_Tot), dBM(ndB_Tot)
integer Ind(nCent), iDCRs(nCent), iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
logical Proc_dB

#ifdef _DEBUGPRINT_
call RecPrt('B',' ',B,3,nCent)
call RecPrt('dB',' ',dB,3*nCent,3*nCent)
write(6,*) iDCRs
#endif

! Set up the T-matrix

! Project away nonsymmetric displacements

call dcopy_(3*nCent,[One],0,Tx,1)
do i=1,nCent
  call NonSym(nStab(Ind(i)),jStab(0,Ind(i)),A(1,i),Tx(1,i))

  ! Rotate vector back to the unique center

  call OA(iDCRS(i),Tx(1:3,i),ATemp)
  Tx(:,i) = ATemp(:)
end do

! Create BqR

nq = 0
do i=1,nCent
  do ixyz=1,3
    if (Smmtrc(ixyz,Ind(i))) then
      iDim = 0
      do j=1,Ind(i)
        jxyz_Max = 3
        if (j == Ind(i)) jxyz_Max = ixyz
        do jxyz=1,jxyz_Max
          if (Smmtrc(jxyz,j)) iDim = iDim+1
        end do
      end do
      mB_Tot = mB_Tot+1
      nq = nq+1
      BM(mB_Tot) = Tx(ixyz,i)*B(ixyz,i)
      iBM(mB_Tot) = iDim
    end if
  end do
end do
nqB(iq) = nq

! Create dBqR

if (Proc_dB) then
  do i=1,nCent
    do ixyz=1,3
      if (Smmtrc(ixyz,Ind(i))) then
        iDim = 0
        do j=1,Ind(i)
          jxyz_Max = 3
          if (j == Ind(i)) jxyz_Max = ixyz
          do jxyz=1,jxyz_Max
            if (Smmtrc(jxyz,j)) iDim = iDim+1
          end do
        end do

        do j=1,nCent
          do jxyz=1,3
            if (Smmtrc(jxyz,Ind(j))) then

              jDim = 0
              do k=1,Ind(j)
                kxyz_Max = 3
                if (k == Ind(j)) kxyz_Max = jxyz
                do kxyz=1,kxyz_Max
                  if (Smmtrc(kxyz,k)) jDim = jDim+1
                end do
              end do

              mdB_Tot = mdB_Tot+1
              dBM(mdB_Tot) = rMult* Tx(ixyz,i)*dB(ixyz,i,jxyz,j)*Tx(jxyz,j)
              idBM(1,mdB_Tot) = iDim
              idBM(2,mdB_Tot) = jDim

            end if
          end do
        end do

      end if
    end do
  end do
end if

return

end subroutine ProjSym
