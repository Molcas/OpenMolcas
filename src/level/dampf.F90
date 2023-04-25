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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine dampF(r,RHOAB,NCMM,MMLR,IVSR,IDSTT,DM)

implicit none
integer NCMM, MMLR(NCMM), IVSR, IDSTT, IDFF, FIRST, m, MM
real*8 r, RHOAB, bTT(-2:2), cDS(-4:0), bDS(-4:0), br, XP, YP, DM(NCMM), bpm(20,-2:0), cpm(20,-2:0), ZK
data bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
data bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
data cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
data FIRST/1/
save FIRST, bpm, cpm

!write(6,*) 'Made it inside of dampF! IVSR=',IVSR
if (NCMM > 4) then
  write(6,*) 'IDSTT=',IDSTT
end if
if (FIRST == 1) then
  do m=1,20
    do IDFF=-2,0
      !bpm(m,IDFF) = bDS(IDFF)/dfloat(m)
      bpm(m,IDFF) = bDS(IDFF)/dble(m)
      !cpm(m,IDFF) = cDS(IDFF)/dsqrt(dfloat(m))
      cpm(m,IDFF) = cDS(IDFF)/dsqrt(dble(m))
    end do
  end do
  FIRST = 0
end if
br = RHOAB*r
!write(6,*) 'NCMM=',NCMM
do m=1,NCMM
  MM = MMLR(m)
  XP = dexp(-(bpm(MM,IVSR)+cpm(MM,IVSR)*br)*br)
  YP = 1.d0-XP
  ZK = MM-1.d0
  DM(m) = YP**(MM-1)
  !... Actually ...  DM(m)= YP**(MM + IVSR/2)  :  set it up this way to
  !   avoid taking exponential of a logarithm for fractional powers (slow)
  if (IVSR == -4) then
    ZK = ZK-1.d0
    DM(m) = DM(m)/YP
  end if
  if (IVSR == -3) then
    ZK = ZK-0.5d0
    DM(m) = DM(m)/dsqrt(YP)
  end if
  if (IVSR == -1) then
    ZK = ZK+0.5d0
    DM(m) = DM(m)*dsqrt(YP)
  end if
  if (IVSR == 0) then
    ZK = MM
    DM(m) = DM(m)*YP
  end if
  if (IVSR == -9) then
  end if
end do
br = bTT(1) !Make sure that it's "referenced" in subroutine too!

return

!600 format(/,' *** ERROR ***  For  IDSTT=',i3,'   IVSR=',i3,' no damping function is defined')
!602 format( /,' ***ERROR ***  RHOAB=', F7.4,'  yields an invalid Damping Function definition')

end subroutine dampF
