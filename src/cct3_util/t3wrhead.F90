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

use CCT3_global, only: eps, filerst, fullprint, ijsegkey, imax, imin, iokey, ispin, jmax, jmin, keysa, lsym, mchntyp, mhkey, mmul, &
                       noa, nob, noop, norb, nsym, nva, nvb, slim, shifto, shiftv, symimax, symimin, symjmax, symjmin, typden, typt3
use Definitions, only: iwp, u6

implicit none
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
!GG write(u6,*)

!2 write occupations and matrix multiplication tablaux

write(u6,11) (norb(nhelp),nhelp=1,nsym)
write(u6,12) (noa(nhelp),nhelp=1,nsym)
write(u6,13) (nob(nhelp),nhelp=1,nsym)
write(u6,14) (nva(nhelp),nhelp=1,nsym)
write(u6,15) (nvb(nhelp),nhelp=1,nsym)

if (fullprint > 1) then
  write(u6,*)
  write(u6,*)
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
write(u6,16) nsym
write(u6,161) ispin
write(u6,162) lsym

!4 write type of triples

if (typt3 == 0) then
  write(u6,31)
else if (typt3 == 1) then
  write(u6,32)
else if (typt3 == 2) then
  write(u6,33)
else if (typt3 == 3) then
  write(u6,34)
end if

!5 write type of fok division

if (typden == 0) then
  write(u6,41)
else if (typden == 1) then
  write(u6,42)
else
  write(u6,43)
end if

!6 write orbital energies per symmetry

if (fullprint > 0) then
  write(u6,*)
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

!9 type file whre CCSD results are
write(u6,91) filerst

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
else
  write(u6,133)
end if
!14 I/O handling

if (iokey == 1) then
  write(u6,141)
else
  write(u6,142)
end if

if (mhkey == 1) then
  write(u6,151)
else
  write(u6,152)
end if

!20 Print parameters of IJ cycle segmentation

if (ijsegkey == 0) then
  write(u6,201)
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
end if

!* can be added

write(u6,*)
write(u6,*)

return

!GG 9 format(A72)
11 format(' NORB ',8(i3,2x))
12 format(' NOA  ',8(i3,2x))
13 format(' NOB  ',8(i3,2x))
14 format(' NVA  ',8(i3,2x))
15 format(' NVB  ',8(i3,2x))
16 format(' NUMBER OF IRREPS             :',i3)
17 format(' MATRIX MULTIPLICATION TABLE')
18 format(' IRREP #',9x,8(i3,2x))
20 format(' IRREP #',i3,6x,8(i3,2x))
31 format(' METHOD                       : CCSD')
32 format(' METHOD                       : CCSD+T(CCSD)')
33 format(' METHOD                       : CCSD+T(CCSD)+<T3(1)WT1> = CCSD(T)')
34 format(' METHOD                       : CCSD+T(CCSD)+<T3(1)WT1>+<T3(1)UT2> = CCSD(T)')
41 format(' TYPE OF DENOMINATOR          : DIAGONAL')
42 format(' TYPE OF DENOMINATOR          : (FAA+FBB)/2')
43 format(' TYPE OF DENOMINATOR          : ORBITAL ENERGIES')
51 format(' LIST OF ORBITAL ENERGIES')
54 format(' IRREDUCIBLE REPRESENTATION NO:',i2)
55 format(' ORBITAL NO:',i3,5x,f16.10)
81 format(' SPIN ADAPTATION             : NONE ')
82 format(' SPIN ADAPTATION             : T2 DDVV ')
83 format(' SPIN ADAPTATION             : T2 DDVV + T1 DV ')
84 format(' SPIN ADAPTATION             : T1 AND T2 FULL ')
85 format(' SPIN ADAPTATION             : T2 FULL WITHOUT SDVS')
91 format(' CCSD RESULTS LOAD FROM FILE : ',a6)
101 format(' PREFERENCE MATRIX MULT.     : NORMAL')
102 format(' PREFERENCE MATRIX MULT.     : TRANSP ; LIMIT =',d12.5)
111 format(' DENOMINATOR SHIFT FOR OCC.  : ',d12.5)
112 format(' DENOMINATOR SHIFT FOR VIRT. : ',d12.5)
!121 format(' MAXIMAL ALLOWED WORK SPACE  : Unlimited')
!122 format(' MAXIMAL ALLOWED WORK SPACE  : ',i10)
131 format(' LEVEL OF OUTPUT PRINTING    : MINIMAL')
132 format(' LEVEL OF OUTPUT PRINTING    : MEDIUM')
133 format(' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
141 format(' INPUT/OUTPUT HANDLING       : Standard SQ ')
142 format(' INPUT/OUTPUT HANDLING       : Molcas4  DA ')
151 format(' MATRIX HANDLING             : ESSL        ')
152 format(' MATRIX HANDLING             : Fortran code')
161 format(' MULTIPLICITY                 :',i3)
162 format(' OVERALL SYMMETRY STATE       :',i3)
201 format(' IJ CYCLE SEGMENTED          : NO          ')
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

end subroutine t3wrhead
