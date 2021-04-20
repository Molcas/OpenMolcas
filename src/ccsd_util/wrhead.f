************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine wrhead
c
c     this routine write head of ther output file
c     parameters transported through cmm common
c
#include "ccsd1.fh"
c
c     help variables
c
       integer nhelp,nhelp1,nhelp2
c
c1    write header and title
c
*      Call HelloPN
      if (noop.eq.1) then
      write(6,'(6X,A)') ' No Operation is required'
      write(6,'(6X,A)') ' Happy Landing'
      Call Finish(0)
      end if
c
      if ( title(1:1).ne.' ' ) then
          write(6,*)
          write(6,'(6X,112A1)') ('*',i=1,112)
          write(6,'(6X,A1,19X,A72,19X,A1)') '*',title,'*'
          write(6,'(6X,112A1)') ('*',i=1,112)
      end if
c
c2    write occupations and matrix multiplication tablaux
c
       Write(6,*)
       Write(6,*)
       Write(6,'(6X,A)')'Wave function specifications:'
       Write(6,'(6X,A)')'-----------------------------'
       Write(6,*)
       Write(6,'(6X,A,T45,I6)') 'Spin mutiplicity',ispin
       Write(6,'(6X,A,T45,I6)')   'State symmetry',lsym
       Write(6,*)
       Write(6,'(6X,A)')'Orbital specifications:'
       Write(6,'(6X,A)')'-----------------------'
       Write(6,*)
       Write(6,'(6X,A,T47,8I4)')
     & 'Symmetry species',
     & (nhelp,nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'Total no. of orbitals',
     & (norb(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of occupied orbitals with alpha spin',
     & (noa(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of occupied orbitals with beta spin',
     & (nob(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of virtual orbitals with alpha spin',
     & (nva(nhelp),nhelp=1,nsym)
       Write(6,'(6X,A,T47,8I4)')
     & 'No. of virtual orbitals with beta spin',
     & (nvb(nhelp),nhelp=1,nsym)
       Write(6,*)
           if (fullprint.gt.1) then
           write(6,*)
           write(6,16) nsym
 16        format (' NUMBER OF IRREPS             :',i3)
           write(6,17)
17         format (' MATRIX MULTIPLICATION TABLE')
           write(6,*)
           write(6,18) (nhelp,nhelp=1,nsym)
18         format (' IRREP #',9x,8(i3,2x))
           write(6,*)
           do 19 nhelp1=1,nsym
           write(6,20) nhelp1,(mmul(nhelp1,nhelp),nhelp=1,nsym)
19         continue
20         format (' IRREP #',i3,6x,8(i3,2x))
           write(6,*)
           end if
       Write(6,*)
       Write(6,'(6X,A)')'Methods and options:'
       Write(6,'(6X,A)')'--------------------'
       Write(6,*)
       Write(6,'(6X,A,T45,I3)') 'Max no. of iterations',maxiter
       if (typden.eq.0) then
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'diagonal Fock matrix elements'
       else if (typden.eq.1) then
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'spin averaged diagonal Fock matrix elements'
       else
         Write(6,'(6X,A,T45,A)') 'Type of denominators',
     &   'orbital energies'
       end if
       Write(6,'(6X,A,T45,F22.14)') 'energy convergence criterium',
     & ccconv
       if (yesext.ne.0) then
       end if
c
c6    write orbital energies per symmetry
c
       write(6,*)
       if (fullprint.gt.0) then
       write(6,51)
 51     format (' LIST OF ORBITAL ENERGIES')
       write(6,*)
c
       nhelp2=1
       do nhelp=1,nsym
       write(6,54) nhelp
 54     format (' IRREDUCIBLE REPRESENTATION NO:',i2)
       do nhelp1=1,norb(nhelp)
       write(6,55) nhelp1,eps(nhelp2)
 55     format (' ORBITAL NO:',i3,5x,f16.10)
       nhelp2=nhelp2+1
       end do
       end do
       write(6,*)
       end if
c
c6    extrapolation parameters
c
       if (yesext.eq.0) then
       write(6,61)
 61     format (' DIIS EXTRAPOLATION USED     : NO')
       else
        write(6,62)
 62      format (' DIIS EXTRAPOLATION USED     : YES')
        write(6,63) firstext
 63      format (' FIRST ITERATION OF EXT.     :',i3)
        write(6,64) cycext
 64      format (' EXTRAPOLATION CYCLE         :',i3)
       end if
       write(6,*)
c
c7    convergence criterion
c
CFUE       write(6,71) ccconv
CFUE 71     format (' ENERGY CONVERGENCE CRITER.  : ',d10.5)
c
c8    spin adaptation
c
       if (keysa.eq.0) then
       write(6,81)
       else if (keysa.eq.1) then
       write(6,82)
       else if (keysa.eq.2) then
       write(6,83)
       else if (keysa.eq.3) then
       write(6,84)
       else if (keysa.eq.4) then
       write(6,85)
       end if
 81     format (' SPIN ADAPTATION             : NONE ')
 82     format (' SPIN ADAPTATION             : T2 DDVV ')
 83     format (' SPIN ADAPTATION             : T2 DDVV + T1 DV ')
 84     format (' SPIN ADAPTATION             : T1 AND T2 FULL ')
 85     format (' SPIN ADAPTATION             : T2 FULL WITHOUT SDVS')
c
c9    restart status
       if (keyrst.eq.0) then
       write(6,91)
       else if (keyrst.eq.1) then
       write(6,92) filerst
       else if (keyrst.eq.2) then
       write(6,93) filerst
       end if
 91     format (' RESTART STATUS              : NONE ')
 92     format (' RST. INF. WILL BE SAVED IN  : ',a6)
 93     format (' RST. INF. WILL BE LOAD FROM : ',a6)
c
c10   write matrix multiplication performance
c
       if (mchntyp.eq.1) then
       write(6,101)
       else if (mchntyp.eq.2) then
       write(6,102) slim
       end if
 101    format (' PREFERENCE MATRIX MULT.     : NORMAL')
 102    format (' PREFERENCE MATRIX MULT.     : TRANSP ; LIMIT =',d12.5)
c
c11   write denominator shifts
       write(6,111) shifto
       write(6,112) shiftv
 111    format (' DENOMINATOR SHIFT FOR OCC.  : ',d12.5)
 112    format (' DENOMINATOR SHIFT FOR VIRT. : ',d12.5)
c
c12   write workspace
c
c      if (maxspace.eq.0) then
c      write(6,121)
c      else
c      write(6,122) maxspace
c      end if
c121    format (' MAXIMAL ALLOWED WORK SPACE  : UNLIMITED')
c122    format (' MAXIMAL ALLOWED WORK SPACE  : ',i10)
c
c13   level of printing
c
       if (fullprint.eq.0) then
       write(6,131)
 131    format (' LEVEL OF OUTPUT PRINTING    : MINIMAL')
       else if (fullprint.eq.1) then
       write(6,132)
 132    format (' LEVEL OF OUTPUT PRINTING    : MEDIUM')
       else if (fullprint.eq.2) then
       write(6,133)
 133    format (' LEVEL OF OUTPUT PRINTING    : MAXIMAL')
       else if (fullprint.eq.3) then
       write(6,134)
 134    format (' LEVEL OF OUTPUT PRINTING    : DEBUG')
       end if
c
c14   I/O handling
c
       if (iokey.eq.1) then
       write(6,141)
 141    format (' INPUT/OUTPUT HANDLING       : Standard SQ ')
       else
       write(6,142)
 142    format (' INPUT/OUTPUT HANDLING       : Molcas4  DA ')
       end if
c
c15    Matrix handling
c
       if (mhkey.eq.1) then
       write(6,151)
 151    format (' MATRIX OPERATIONS           : ESSL        ')
       else
       write(6,152)
 152    format (' MATRIX OPERATIONS           : Fortran code')
       end if
c
c*    can be added
c
       write(6,*)
       write(6,*)
c
       return
       end
