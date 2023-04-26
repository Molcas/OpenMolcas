************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Anders Bernhardsson                              *
************************************************************************
      Subroutine RHS_PT2(rkappa,CLag,SLag)

        Implicit Real*8(a-h,o-z)
        Real*8 rKappa(*),CLag(*),SLag(*)

! The RHS array for CASPT2 has been already calculated in the
! CASPT2 module, so here we only need to read it from file

#include "Pointers.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "Files_mclr.fh"

!     Read in a and b part of effective gradient from CASPT2

      nOLag = 0
      nCLag = 0
      DO i = 1, nSym
        nOLag = nOLag + nOrb(i)*nOrb(i)
        nCLag = nCLag + nRoots*nCSF(i)
      END DO
      nSLag = nRoots*nRoots

      Do i = 1, nCLag
        Read (LUPT2,*,END=200) CLag(i)
      End Do
      Do i = 1, nOLag
        Read (LUPT2,*,END=200) tmp ! rKappa(i)
        rKappa(i) = rKappa(i) + tmp
      End Do
      Do i = 1, nSLag
        Read (LUPT2,*,END=200) SLag(i)
      End Do

      return

  200 continue
      write(6,*)
      write(6,'(1x,"The file which has to be written in CASPT2 module ",
     *            "does not exist in RHS_PT2.")')
      write(6,'(1x,"For single-point gradient calculation, you need ",
     *            "GRAD or GRDT keyword in &CASPT2.")')
      write(6,'(1x,"For geometry optimization, you do not need ",
     *            "anything, so something is wrong with the code.")')
      write(6,*)
      call abend()

      End