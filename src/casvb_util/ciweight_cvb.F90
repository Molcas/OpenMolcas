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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

!***********************************************************************
!*                                                                     *
!*  CIWEIGHT  := Chirgwin-Couson/Lowdin/inverse-overlap weights of     *
!*            := full CASSCF vector and residual.                      *
!*                                                                     *
!***********************************************************************
subroutine ciweight_cvb(civec,civbs,civb,citmp,vec5,orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: civec(*), civbs(*), civb(*), citmp(*), vec5(*), orbs(norb,norb), sorbs(norb,norb), orbinv(norb,norb), &
                 owrk(norb,norb), gjorb(*), gjorb2(*), gjorb3(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: icitmp, icivb, icivbs, icivec, ion, ionmax, ionmin, iretval1, iretval2, ivec5, k1, k10, k11, k12, k13, k14, &
                     k15, k16, k17, k18, k19, k2, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k3, k30, k31, k32, k33, k34, &
                     k4, k5, k6, k7, k8, k9, mxasg, mxdetcas, mxrem, mxsng, ncnfcas
integer(kind=iwp), external :: mstacki_cvb, mstackr_cvb

ionmin = max(nel-norb,0)
ionmax = nbet
mxrem = norb-ionmin
mxsng = nel-2*ionmin
mxasg = nalf-ionmin
call icomb_cvb(mxsng,mxasg,mxdetcas)
! Work out number of configurations in CASSCF vector:
ncnfcas = 0
do ion=ionmin,ionmax
  call icomb_cvb(norb,ion,iretval1)
  call icomb_cvb(norb-ion,nel-2*ion,iretval2)
  ncnfcas = ncnfcas+iretval1*iretval2
end do

k1 = mstackr_cvb(ionmax-ionmin+1)
k2 = mstackr_cvb(ionmax-ionmin+1)
k3 = mstackr_cvb(ionmax-ionmin+1)
k4 = mstackr_cvb(ionmax-ionmin+1)
k5 = mstackr_cvb(ionmax-ionmin+1)
k6 = mstackr_cvb(ionmax-ionmin+1)
k7 = mstacki_cvb(norb+1)
k8 = mstacki_cvb(norb+1)
k9 = mstacki_cvb((norb+1)*(nalf+1))
k10 = mstacki_cvb((norb+1)*(nbet+1))
k11 = mstacki_cvb(norb)
k12 = mstacki_cvb(norb)
k13 = mstacki_cvb(norb+1)
k14 = mstacki_cvb(norb+1)
k15 = mstacki_cvb(norb+1)
k16 = mstacki_cvb((norb+1)*(ionmax+1))
k17 = mstacki_cvb(norb)
k18 = mstacki_cvb(norb)
k19 = mstacki_cvb(norb+1)
k20 = mstacki_cvb(norb+1)
k21 = mstacki_cvb(norb+1)
k22 = mstacki_cvb((mxrem+1)*(mxsng+1))
k23 = mstacki_cvb(norb)
k24 = mstacki_cvb(norb)
k25 = mstacki_cvb(norb+1)
k26 = mstacki_cvb(norb+1)
k27 = mstacki_cvb(norb+1)
k28 = mstacki_cvb((mxsng+1)*(mxasg+1))
k29 = mstacki_cvb(norb)
k30 = mstacki_cvb(norb)
if (mod(iciweights,8) > 3) then
  k31 = mstackr_cvb(ncnfcas)
  k32 = mstackr_cvb(ncnfcas)
  k33 = mstacki_cvb(mxdetcas)
  k34 = mstacki_cvb(mxdetcas)
else
  k31 = 0
  k32 = 0
  k33 = 0
  k34 = 0
end if
icivec = nint(civec(1))
icivbs = nint(civbs(1))
icivb = nint(civb(1))
icitmp = nint(citmp(1))
ivec5 = nint(vec5(1))
call ciweight2_cvb(civec,civbs,civb,citmp,vec5,orbs,sorbs,orbinv,owrk,gjorb,gjorb2,gjorb3,work(iaddr_ci(icitmp)), &
                   work(iaddr_ci(icivbs)),work(iaddr_ci(icivec)),work(iaddr_ci(icivb)),work(iaddr_ci(ivec5)),work(k1),work(k2), &
                   work(k3),work(k4),work(k5),work(k6),iwork(k7),iwork(k8),iwork(k9),iwork(k10),iwork(k11),iwork(k12),iwork(k13), &
                   iwork(k14),iwork(k15),iwork(k16),iwork(k17),iwork(k18),iwork(k19),iwork(k20),iwork(k21),iwork(k22),iwork(k23), &
                   iwork(k24),iwork(k25),iwork(k26),iwork(k27),iwork(k28),iwork(k29),iwork(k30),work(k31),work(k32),iwork(k33), &
                   iwork(k34),ionmin,ionmax,mxrem,mxsng,mxasg,ncnfcas,mxdetcas)
call mfreer_cvb(k1)

return

end subroutine ciweight_cvb
