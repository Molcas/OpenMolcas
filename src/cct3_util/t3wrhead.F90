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

subroutine t3wrhead()
! this routine writes head of the output file
! parameters transported through cmm common

use Definitions, only: iwp, u6

implicit none
#include "t31.fh"
integer(kind=iwp) :: nhelp, nhelp1, nhelp2

!1 write title

if (noop == 1) then
  write(u6,'(6X,A)') ' No operation is required'
  write(u6,'(6X,A)') ' Happy Landing'
  call Finish(0)
end if

!GG do nhelp=1,ntit
!GG   write(u6,9) title(nhelp)
!GG end do
!GG 9 format(A72)
!GG write(u6,*)

!2 write occupations and matrix multiplication tablaux

write(u6,11) (norb(nhelp),nhelp=1,nsym)
11 format(' NORB ',8(i3,2x))
write(u6,12) (noa(nhelp),nhelp=1,nsym)
12 format(' NOA  ',8(i3,2x))
write(u6,13) (nob(nhelp),nhelp=1,nsym)
13 format(' NOB  ',8(i3,2x))
write(u6,14) (nva(nhelp),nhelp=1,nsym)
14 format(' NVA  ',8(i3,2x))
write(u6,15) (nvb(nhelp),nhelp=1,nsym)
15 format(' NVB  ',8(i3,2x))

if (fullprint > 1) then
  write(u6,*)
  write(u6,*)
  write(u6,17)
17 format(' MATRIX MULTIPLICATION TABLE')
  write(u6,*)
  write(u6,18) (nhelp,nhelp=1,nsym)
18 format(' IRREP #',9x,8(i3,2x))
  write(u6,*)
  do nhelp1=1,nsym
    write(u6,20) nhelp1,(mmul(nhelp1,nhelp),nhelp=1,nsym)
  end do
20 format(' IRREP #',i3,6x,8(i3,2x))
  write(u6,*)
end if

write(u6,*)
write(u6,16) nsym
16 format(' NUMBER OF IRREPS             :',i3)
write(u6,161) ispin
161 format(' MULTIPLICITY                 :',i3)
write(u6,162) lsym
162 format(' OVERALL SYMMETRY STATE       :',i3)

!4 write type of triples

if (typt3 == 0) then
  write(u6,31)
31 format(' METHOD                       : CCSD')
else if (typt3 == 1) then
  write(u6,32)
32 format(' METHOD                       : CCSD+T(CCSD)')
else if (typt3 == 2) then
  write(u6,33)
33 format(' METHOD                       : CCSD+T(CCSD)+<T3(1)WT1> = CCSD(T)')
else if (typt3 == 3) then
  write(u6,34)
34 format(' METHOD                       : CCSD+T(CCSD)+<T3(1)WT1>+<T3(1)UT2> = CCSD(T)')
end if

!5 write type of fok division

if (typden == 0) then
  write(u6,41)
else if (typden == 1) then
  write(u6,42)
else
  write(u6,43)
end if
41 format(' TYPE OF DENOMINATOR          : DIAGONAL')
43 format(' TYPE OF DENOMINATOR          : ORBITAL ENERGIES')
42 format(' TYPE OF DENOMINATOR          : (FAA+FBB)/2')

!6 write orbital energies per symmetry

if (fullprint > 0) then
  write(u6,*)
  write(u6,51)
51 format(' LIST OF ORBITAL ENERGIES')
  write(u6,*)

  nhelp2 = 1
  do nhelp=1,nsym
    write(u6,54) nhelp
54  format(' IRREDUCIBLE REPRESENTATION NO:',i2)
    do nhelp1=1,norb(nhelp)
      write(u6,55) nhelp1,eps(nhelp2)
55    format(' ORBITAL NO:',i3,5x,f16.10)
      nhelp2 = nhelp2+1
    end do
  end do
  write(u6,*)
end if

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
81 format(' SPIN ADAPTATION             : NONE ')
82 format(' SPIN ADAPTATION             : T2 DDVV ')
83 format(' SPIN ADAPTATION             : T2 DDVV + T1 DV ')
84 format(' SPIN ADAPTATION             : T1 AND T2 FULL ')
85 format(' SPIN ADAPTATION             : T2 FULL WITHOUT SDVS')

!9 type file whre CCSD results are
write(u6,91) filerst
91 format(' CCSD RESULTS LOAD FROM FILE : ',a6)

!10 write matrix multiplication performance

if (mchntyp == 1) then
  write(u6,101)
else if (mchntyp == 2) then
  write(u6,102) slim
end if
101 format(' PREFERENCE MATRIX MULT.     : NORMAL')
102 format(' PREFERENCE MATRIX MULT.     : TRANSP ; LIMIT =',d12.5)

!11 write denominator shifts
write(u6,111) shifto
write(u6,112) shiftv
111 format(' DENOMINATOR SHIFT FOR OCC.  : ',d12.5)
112 format(' DENOMINATOR SHIFT FOR VIRT. : ',d12.5)

!12 write workspace

!if (maxspace == 0) then
!  write(u6,121)
!else
!  write(u6,122) maxspace
!end if
!121 format(' MAXIMAL ALLOWED WORK SPACE  : Unlimited')
!122 format(' MAXIMAL ALLOWED WORK SPACE  : ',i10)

!13 level of printing

if (fullprint == 0) then
  write(u6,131)
131 format(' LEVEL OF OUTPUT PRINTING    : MINIMAL')
else if (fullprint == 1) then
  write(u6,132)
132 format(' LEVEL OF OUTPUT PRINTING    : MEDIUM')
else
  write(u6,133)
133 format(' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
end if
!14 I/O handling

if (iokey == 1) then
  write(u6,141)
141 format(' INPUT/OUTPUT HANDLING       : Standard SQ ')
else
  write(u6,142)
142 format(' INPUT/OUTPUT HANDLING       : Molcas4  DA ')
end if
!
if (mhkey == 1) then
  write(u6,151)
151 format(' MATRIX HANDLING             : ESSL        ')
else
  write(u6,152)
152 format(' MATRIX HANDLING             : Fortran code')
end if

!20 Print parameters of IJ cycle segmentation

if (ijsegkey == 0) then
  write(u6,201)
201 format(' IJ CYCLE SEGMENTED          : NO          ')
else
  write(u6,202)
  write(u6,203) symimin
  write(u6,204) imin
  write(u6,205) symjmin
  write(u6,206) jmin
  write(u6,207) symimax
  write(u6,208) imax
  write(u6,209) symjmax
  write(u6,210) jmax
  write(u6,211)
202 format(' IJ CYCLE SEGMENTED          : YES         ')
203 format(' SYMI minimal                : ',i4)
204 format('    I minimal                : ',i4)
205 format(' SYMJ minimal                : ',i4)
206 format('    J minimal                : ',i4)
207 format(' SYMI maximal                : ',i4)
208 format('    I maximal                : ',i4)
209 format(' SYMJ maximal                : ',i4)
210 format('    J maximal                : ',i4)
211 format(' Be very careful in using IJSEgmentation technique')
end if

!* can be added

write(u6,*)
write(u6,*)

return

end subroutine t3wrhead
