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

subroutine Mod_P2(P2mo,nP2Act,D1mo,nD1mo,DS1mo,ExFac,nDet)

use nq_Info, only: iOff_Ash, mIrrep, nAsh
use Constants, only: One, Two, Quart
use Definitions, only: wp, u6

implicit none
integer nP2Act, nD1mo, nDet
real*8 P2mo(nP2Act), D1mo(nD1mo), DS1mo(nD1mo), ExFac
integer iOff_, iIrrep, jIrrep, kIrrep, ijIrrep, ijkIrrep, k_, k, l_, l, kl, i, i_, j, j_, il, ik, ij, ijkl, jk, jl
real*8 P2Act, Fact
! Statement function
integer iTri
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
iOff_ = 0
do iIrrep=0,mIrrep-1
  iOff_Ash(iIrrep) = iOff_
  iOff_ = iOff_+nAsh(iIrrep)
end do
!
!***********************************************************************
!
! active space Ptvxy

if (nDet == 1) then

  P2Act = real(nP2Act,kind=wp)
  call Put_Temp('nP2Act  ',[P2Act],1)
  call Put_Temp('P2_RAW  ',P2mo,nP2Act)

  do iIrrep=0,mIrrep-1
    do jIrrep=0,mIrrep-1
      ijIrrep = ieor(iIrrep,jIrrep)
      do kIrrep=0,mIrrep-1
        ijkIrrep = ieor(ijIrrep,kIrrep)

        do k_=1,nASh(kIrrep)
          k = iOff_Ash(kIrrep)+k_
          do l_=1,nASh(ijkIrrep)
            l = iOff_Ash(ijkIrrep)+l_
            if (l > k) Go To 100
            kl = iTri(k,l)
            do i_=1,nASh(iIrrep)
              i = iOff_Ash(iIrrep)+i_
              il = iTri(i,l)
              ik = iTri(i,k)
              do j_=1,nASh(jIrrep)
                j = iOff_Ash(jIrrep)+j_
                if (j > i) Go To 200
                ij = iTri(i,j)
                if (kl > ij) Go To 200
                ijkl = iTri(ij,kl)
                jk = iTri(j,k)
                jl = iTri(j,l)

                Fact = One
                if ((iIrrep == jIrrep) .and. (k == l)) Fact = Two
                P2mo(ijkl) = Fact*P2mo(ijkl)

                if (iIrrep == ijkIrrep) P2mo(ijkl) = P2mo(ijkl)+(One-ExFac)*Quart*(D1mo(jk)*D1mo(il)+DS1mo(jk)*DS1mo(il))

                if (iIrrep == kIrrep) P2mo(ijkl) = P2mo(ijkl)+(One-ExFac)*Quart*(D1mo(jl)*D1mo(ik)+DS1mo(jl)*DS1mo(ik))

                P2mo(ijkl) = P2mo(ijkl)/Fact

200             continue
              end do
            end do
100         continue
          end do
        end do

      end do
    end do
  end do
  call Put_Temp('P2_KS   ',P2mo,nP2Act)
else
  write(u6,*) ' Not implemented yet!!! nDet=',nDet
  call Abend()
end if

return

end subroutine Mod_P2
