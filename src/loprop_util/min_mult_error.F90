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

subroutine Min_Mult_Error(EC,A,B,Ci,Cj,rMP,xrMP,xxrMP,xnrMP,lMax,nij,nElem,iAtom,jAtom,nAtoms,nPert,C_o_C,Scratch_New,Scratch_Org, &
                          iPlot,T_Values,iWarnings,Num_Warnings)

implicit real*8(A-H,O-Z)
#include "real.fh"
dimension EC(3,nij), C_o_C(3), A(3,nij), B(3,nij), Ci(3), Cj(3), R_ij(3)
dimension rMP(nij,0:nElem-1,0:nPert-1), xnrMP(nij,nElem)
dimension xrMP(nij,nElem), xxrMP(nij,nElem), T_Values(nij)
dimension Scratch_Org(nij*(2+lMax+1)), Scratch_New(nij*(2+lMax+1))
dimension iWarnings(nij)
parameter(Error_Threshold=1.0D-12,Delta_Threshold=1.0D-12)
external Error_for_t, Golden

ij = iAtom*(iAtom-1)/2+jAtom
iPrint_Errors = 0
l = lMax-1
do i=1,3
  R_ij(i) = Cj(i)-Ci(i)
end do

t_min = Zero
t_max = Zero
do i=1,3
  if (R_ij(i) /= 0.0d0) then
    t_min = (Ci(i)-EC(i,ij))/R_ij(i)
    t_max = (Cj(i)-EC(i,ij))/R_ij(i)
  end if
end do
Delta_Orig = 0.1d0
Delta = Delta_Orig

! Check the range of possible t-values for minima

num_min = 0
iSlope = 0
if (iPlot == 1) then
  write(6,*)
  write(6,*) 'iAtom, jAtom = ',iAtom,jAtom
end if
t_temp = t_min
Error_Best = -1.0d0
Error_Old = Zero
t_best = 0.0d0
50 continue
Error = Error_for_t(t_temp,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org, &
                    iPrint_Errors)
if (iPlot == 1) then
  write(6,'(1X,A,F5.2,F16.12)') 't, Error = ',t_temp,Error
  call xFlush(6)
end if
Delta_Error = Error-Error_Old
Error_Old = Error
iSlope_Old = iSlope
if (abs(Delta_Error) < Error_Threshold) then
  iSlope = 0
else if (Delta_Error < Zero) then
  iSlope = -1
else
  iSlope = 1
end if
if ((iSlope_Old < 0) .and. (iSlope >= 0)) then
  num_min = num_min+1
end if
if ((Error < Error_Best) .or. (Error_Best < Zero)) then
  Error_Best = Error
  t_best = t_temp
end if
t_temp = t_temp+Delta/10.0d0
if (t_temp <= t_max+Delta/100.0d0) goto 50

! Any warnings from scan?

if (num_min > 1) then
  iWarnings(ij) = 1
  Num_Warnings = Num_Warnings+1
end if

! Find minima with Golden Section Search
!
! (It's assumed that the error function is either constant or has at
! least one minima somewhere, although not necessarily within the
! allowed range).

ax = t_best
bx = t_best+Delta
! First make an initial bracketing of the minima
call mnBrak(ax,bx,cx,fa,fb,fc,Error_for_t,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New, &
            Scratch_Org,iPrint_Errors)

if (abs(fa-fc) < Error_Threshold) then
  iWarnings(ij) = 4
  Num_Warnings = Num_Warnings+1
  t_final = Zero
else
  Error = Golden(ax,bx,cx,Error_for_t,Delta_Threshold,Error_Threshold,t_final,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax, &
                 nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
end if

! Check that the minima is within the allowed range

if (t_final > t_max) then
  t_final = t_max
  iWarnings(ij) = 2
  Num_Warnings = Num_Warnings+1
else if (t_final < t_min) then
  t_final = t_min
  iWarnings(ij) = 2
  Num_Warnings = Num_Warnings+1
end if

T_values(ij) = t_final
B(1,ij) = EC(1,ij)+t_final*R_ij(1)
B(2,ij) = EC(2,ij)+t_final*R_ij(2)
B(3,ij) = EC(3,ij)+t_final*R_ij(3)

return

end subroutine Min_Mult_Error
