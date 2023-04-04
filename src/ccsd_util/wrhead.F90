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

subroutine wrhead()
! this routine writes head of the output file
! parameters transported through cmm common

use ccsd_global, only: ccconv, cycext, eps, filerst, firstext, fullprint, iokey, ispin, keyrst, keysa, lsym, maxiter, mchntyp, &
                       mhkey, mmul, noa, nob, noop, norb, nsym, nva, nvb, shifto, shiftv, slim, title, typden, yesext
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i, nhelp, nhelp1, nhelp2

!1 write header and title

!call HelloPN()
if (noop == 1) then
  write(u6,'(6X,A)') ' No Operation is required'
  write(u6,'(6X,A)') ' Happy Landing'
  call Finish(0)
end if

if (title(1:1) /= ' ') then
  write(u6,*)
  write(u6,'(6X,112A1)') ('*',i=1,112)
  write(u6,'(6X,A1,19X,A72,19X,A1)') '*',title,'*'
  write(u6,'(6X,112A1)') ('*',i=1,112)
end if

!2 write occupations and matrix multiplication tablaux

write(u6,*)
write(u6,*)
write(u6,'(6X,A)') 'Wave function specifications:'
write(u6,'(6X,A)') '-----------------------------'
write(u6,*)
write(u6,'(6X,A,T45,I6)') 'Spin mutiplicity',ispin
write(u6,'(6X,A,T45,I6)') 'State symmetry',lsym
write(u6,*)
write(u6,'(6X,A)') 'Orbital specifications:'
write(u6,'(6X,A)') '-----------------------'
write(u6,*)
write(u6,'(6X,A,T47,8I4)') 'Symmetry species',(nhelp,nhelp=1,nsym)
write(u6,'(6X,A,T47,8I4)') 'Total no. of orbitals',(norb(nhelp),nhelp=1,nsym)
write(u6,'(6X,A,T47,8I4)') 'No. of occupied orbitals with alpha spin',(noa(nhelp),nhelp=1,nsym)
write(u6,'(6X,A,T47,8I4)') 'No. of occupied orbitals with beta spin',(nob(nhelp),nhelp=1,nsym)
write(u6,'(6X,A,T47,8I4)') 'No. of virtual orbitals with alpha spin',(nva(nhelp),nhelp=1,nsym)
write(u6,'(6X,A,T47,8I4)') 'No. of virtual orbitals with beta spin',(nvb(nhelp),nhelp=1,nsym)
write(u6,*)
if (fullprint > 1) then
  write(u6,*)
  write(u6,16) nsym
  write(u6,17)
  write(u6,*)
  write(u6,18) (nhelp,nhelp=1,nsym)
  write(u6,*)
  do nhelp1=1,nsym
    write(u6,20) nhelp1,(mmul(nhelp1,nhelp),nhelp=1,nsym)
  end do
  write(u6,*)
end if
write(u6,*)
write(u6,'(6X,A)') 'Methods and options:'
write(u6,'(6X,A)') '--------------------'
write(u6,*)
write(u6,'(6X,A,T45,I3)') 'Max no. of iterations',maxiter
if (typden == 0) then
  write(u6,'(6X,A,T45,A)') 'Type of denominators','diagonal Fock matrix elements'
else if (typden == 1) then
  write(u6,'(6X,A,T45,A)') 'Type of denominators','spin averaged diagonal Fock matrix elements'
else
  write(u6,'(6X,A,T45,A)') 'Type of denominators','orbital energies'
end if
write(u6,'(6X,A,T45,F22.14)') 'energy convergence criterium',ccconv
!if (yesext /= 0) then
!end if

!6 write orbital energies per symmetry

write(u6,*)
if (fullprint > 0) then
  write(u6,51)
  write(u6,*)

  nhelp2 = 1
  do nhelp=1,nsym
    write(u6,54) nhelp
    do nhelp1=1,norb(nhelp)
      write(u6,55) nhelp1,eps(nhelp2)
      nhelp2 = nhelp2+1
    end do
  end do
  write(u6,*)
end if

!6 extrapolation parameters

if (yesext == 0) then
  write(u6,61)
else
  write(u6,62)
  write(u6,63) firstext
  write(u6,64) cycext
end if
write(u6,*)

!7  convergence criterion

!FUE write(u6,71) ccconv

!8 spin adaptation

if (keysa == 0) then
  write(u6,81)
else if (keysa == 1) then
  write(u6,82)
else if (keysa == 2) then
  write(u6,83)
else if (keysa == 3) then
  write(u6,84)
else if (keysa == 4) then
  write(u6,85)
end if

!9 restart status
if (keyrst == 0) then
  write(u6,91)
else if (keyrst == 1) then
  write(u6,92) filerst
else if (keyrst == 2) then
  write(u6,93) filerst
end if

!10 write matrix multiplication performance

if (mchntyp == 1) then
  write(u6,101)
else if (mchntyp == 2) then
  write(u6,102) slim
end if

!11 write denominator shifts
write(u6,111) shifto
write(u6,112) shiftv

!12 write workspace

!if (maxspace == 0) then
!  write(u6,121)
!else
!  write(u6,122) maxspace
!end if

!13 level of printing

if (fullprint == 0) then
  write(u6,131)
else if (fullprint == 1) then
  write(u6,132)
else if (fullprint == 2) then
  write(u6,133)
else if (fullprint == 3) then
  write(u6,134)
end if

!14 I/O handling

if (iokey == 1) then
  write(u6,141)
else
  write(u6,142)
end if

!15 Matrix handling

if (mhkey == 1) then
  write(u6,151)
else
  write(u6,152)
end if

!* can be added

write(u6,*)
write(u6,*)

return

16 format(' NUMBER OF IRREPS             :',i3)
17 format(' MATRIX MULTIPLICATION TABLE')
18 format(' IRREP #',9x,8(i3,2x))
20 format(' IRREP #',i3,6x,8(i3,2x))
51 format(' LIST OF ORBITAL ENERGIES')
54 format(' IRREDUCIBLE REPRESENTATION NO:',i2)
55 format(' ORBITAL NO:',i3,5x,f16.10)
61 format(' DIIS EXTRAPOLATION USED     : NO')
62 format(' DIIS EXTRAPOLATION USED     : YES')
63 format(' FIRST ITERATION OF EXT.     :',i3)
64 format(' EXTRAPOLATION CYCLE         :',i3)
!FUE 71 format (' ENERGY CONVERGENCE CRITER.  : ',f10.5)
81 format(' SPIN ADAPTATION             : NONE ')
82 format(' SPIN ADAPTATION             : T2 DDVV ')
83 format(' SPIN ADAPTATION             : T2 DDVV + T1 DV ')
84 format(' SPIN ADAPTATION             : T1 AND T2 FULL ')
85 format(' SPIN ADAPTATION             : T2 FULL WITHOUT SDVS')
91 format(' RESTART STATUS              : NONE ')
92 format(' RST. INF. WILL BE SAVED IN  : ',a6)
93 format(' RST. INF. WILL BE LOAD FROM : ',a6)
101 format(' PREFERENCE MATRIX MULT.     : NORMAL')
102 format(' PREFERENCE MATRIX MULT.     : TRANSP ; LIMIT =',f12.5)
111 format(' DENOMINATOR SHIFT FOR OCC.  : ',f12.5)
112 format(' DENOMINATOR SHIFT FOR VIRT. : ',f12.5)
!121 format (' MAXIMAL ALLOWED WORK SPACE  : UNLIMITED')
!122 format (' MAXIMAL ALLOWED WORK SPACE  : ',i10)
131 format(' LEVEL OF OUTPUT PRINTING    : MINIMAL')
132 format(' LEVEL OF OUTPUT PRINTING    : MEDIUM')
133 format(' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
134 format(' LEVEL OF OUTPUT PRINTING    : DEBUG')
141 format(' INPUT/OUTPUT HANDLING       : Standard SQ ')
142 format(' INPUT/OUTPUT HANDLING       : Molcas4  DA ')
151 format(' MATRIX OPERATIONS           : ESSL        ')
152 format(' MATRIX OPERATIONS           : Fortran code')

end subroutine wrhead
