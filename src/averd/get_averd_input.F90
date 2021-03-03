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
!---------------------------------------------------------------------
! Read input for averd.
!
! Title         -Title
! Wset          -Set of weights. Dimension is Nset.
! iPrint        -How much print.(1=minimal,2=print average orbitals,
!                5=print all orbitals,99=wacko!)
! Nset          -Number of orbitals to read
! DensityBased  -Is the procedure density or orbital based
! ThrOcc        -Print which orbitals have occupation number below
!                this threshold.
!-----------------------------------------------------------------------

subroutine Get_Averd_input(Title,Wset,iPrint,Nset,DensityBased,ThrOcc)

implicit real*8(a-h,o-z)

#include "mxave.fh"
#include "warnings.fh"

character*72 Title
dimension Wset(MxSets)
logical DensityBased

character*180 Key
character*4 Kword
character*180 Get_Ln
integer iCLast
external Get_Ln, iCLast

!-- Call subroutines that handle the input.

LuRd = 21
call SpoolInp(LuRd)
rewind(LuRd)
call RdNLst(LuRd,'AVERD')

!-- Label 1000 is the top.

1000 continue

!-- Get_Ln read the keyword and skips line starting with *
!   or is empty.

Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)

!-- The keywords...

if (Kword(1:4) == 'WSET') Go To 101
if (Kword(1:4) == 'PRIN') Go To 102
if (Kword(1:4) == 'TITL') Go To 103
if (Kword(1:4) == 'ORBI') Go To 104
if (Kword(1:4) == 'OCCU') Go To 105
if (Kword(1:4) == 'END ') Go To 9999

!-- ...and what happens if something else is encountered.

iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(6,*) ' '
write(6,'(1X,A,A)') Kword(1:Last),' is not a valid keyword!'
write(6,*) ' ERROR!'
call Quit(_RC_INPUT_ERROR_)

!-- Read weights.

101 continue
Key = Get_Ln(LuRd)
call Get_I1(1,Nset)
Key = Get_Ln(LuRd)
call Get_F(1,Wset,Nset)
Go To 1000

!-- How much print?

102 continue
Key = Get_Ln(LuRd)
call Get_I1(1,iPrint)
Go To 1000

!-- Title

103 continue
Key = Get_Ln(LuRd)
Title = Key(1:len(Title))
Go To 1000

!-- Should it be density based, or orbital based.

104 continue
DensityBased = .false.
Go To 1000

!-- I want to be told which orbitals have occ.num. below threshold.

105 continue
Key = Get_Ln(LuRd)
call Get_F1(1,ThrOcc)
Go To 1000

!-- A most Graceful Exit.

9999 continue

return

end subroutine Get_Averd_input
