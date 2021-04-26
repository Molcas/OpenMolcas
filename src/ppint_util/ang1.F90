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

subroutine ang1(ang,dfac,nanb,lalb,mamb,lamu,lmf,lml,lmx,lmy,lmz,ltot1,xab,yab,zab,xk,yk,zk,zlm)
! compute type 1 angular integrals

implicit real*8(a-h,o-z)
parameter(a0=0.0d0,a1=1.0d0)
dimension ang(ltot1,*), dfac(*), lmf(*), lml(*), lmx(*), lmy(*), lmz(*), xab(*), yab(*), zab(*), zlm(*)

call wzero(ltot1*lamu,ang,1)
do n=1,nanb
  if (xab(n) == a0) go to 96
  do l=1,lalb
    if (yab(l) == a0) go to 94
    do m=1,mamb
      if (zab(m) == a0) go to 92
      nlm = ((n-2)+l)+m
      lamlo = mod(nlm-1,2)+1
      lamhi = min(nlm,lamu)
      if (lamlo > lamhi) go to 92
      do lam=lamlo,lamhi,2
        l2 = lam+lam-1
        angt = a0
        loc = (lam-1)**2
        do mu1=1,l2
          istart = lmf(loc+mu1)
          if ((mod(n,2) == mod(lmx(istart),2)) .or. (mod(l,2) == mod(lmy(istart),2)) .or. (mod(m,2) == mod(lmz(istart),2))) go to 80
          pre = a0
          aint = a0
          iend = lml(loc+mu1)
          do i=istart,iend
            indx = lmx(i)
            indy = lmy(i)
            indz = lmz(i)
            if (indx == 0) then
              xkp = a1
            else
              xkp = xk**indx
            end if
            if (indy == 0) then
              ykp = a1
            else
              ykp = yk**indy
            end if
            if (indz == 0) then
              zkp = a1
            else
              zkp = zk**indz
            end if
            pre = pre+zlm(i)*xkp*ykp*zkp
            aint = aint+zlm(i)*dfac(n+indx)*dfac(l+indy)*dfac(m+indz)/dfac((n+indx)+(l+indy)+(m+indz))
          end do
          angt = angt+pre*aint
          80        continue
        end do
        ang(nlm,lam) = ang(nlm,lam)+((xab(n)*yab(l))*zab(m))*angt
      end do
      92    continue
    end do
    94  continue
  end do
  96 continue
end do

return

end subroutine ang1
