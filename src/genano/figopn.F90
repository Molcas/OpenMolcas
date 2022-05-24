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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine prints the prologue part to the postscript figure of    *
! occupation numbers.                                                  *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine FigOpn(lu)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lu

write(lu,'(a)') '%!PS'
write(lu,'(a)') '%%BoundingBox: 70 95 550 735'
write(lu,'(a)') '/TF /Helvetica findfont 18 scalefont def'
write(lu,'(a)') '/SF /Symbol findfont 18 scalefont def'
write(lu,'(a)') '/HF /Helvetica findfont 24 scalefont def'
write(lu,'(a)') '/Y {0.5 sub neg 75 mul} def'
write(lu,'(a)') '/X {75 mul} def'
write(lu,'(a)') '/Draw {'
write(lu,'(a)') '  newpath'
write(lu,'(a)') '  2.0 setlinewidth'
write(lu,'(a)') '  moveto'
write(lu,'(a)') '  0.5 X 0 rlineto'
write(lu,'(a)') '  currentpoint /y exch def /x exch def'
write(lu,'(a)') '  stroke'
write(lu,'(a)') '  x y moveto'
write(lu,'(a)') '} def'
write(lu,'(a)') '/Label {'
write(lu,'(a)') '  TF setfont'
write(lu,'(a)') '  6 -6 rmoveto show'
write(lu,'(a)') '} def'
write(lu,'(a)') '/Connect {'
write(lu,'(a)') '  [5 5] 0 setdash'
write(lu,'(a)') '  newpath'
write(lu,'(a)') '  0.5 setlinewidth'
write(lu,'(a)') '  moveto'
write(lu,'(a)') '  exch 0.5 X add exch'
write(lu,'(a)') '  lineto'
write(lu,'(a)') '  stroke'
write(lu,'(a)') '  [] 0 setdash'
write(lu,'(a)') '} def'
write(lu,'(a)') '/CenterLine {'
!***********************************************************************
write(lu,'(a)') '  dup stringwidth pop 0.5 mul neg 0 rmoveto show'
write(lu,'(a)') '} def'
write(lu,'(a)') '/RightLine {'
write(lu,'(a)') '  dup stringwidth pop neg 0 rmoveto show'
write(lu,'(a)') '} def'
write(lu,'(a)') '%.ManualFeed'
write(lu,'(a)') 'gsave'
write(lu,'(a)') '144 144 translate'
write(lu,'(a)') '-0.30 X -6.6 Y moveto'
write(lu,'(a)') 'TF setfont (lg\050) CenterLine'
write(lu,'(a)') 'SF setfont (h) show'
write(lu,'(a)') 'TF setfont (\051) show'
write(lu,'(a)') '%---'
write(lu,'(a)') 'TF setfont'
write(lu,'(a)') '0.25 X 1.0 Y moveto (s) CenterLine'
write(lu,'(a)') '1.25 X 1.0 Y moveto (p) CenterLine'
write(lu,'(a)') '2.25 X 1.0 Y moveto (d) CenterLine'
write(lu,'(a)') '3.25 X 1.0 Y moveto (f) CenterLine'
write(lu,'(a)') '4.25 X 1.0 Y moveto (g) CenterLine'
write(lu,'(a)') '%---'
write(lu,'(a)') '3.0 setlinewidth'
write(lu,'(a)') 'newpath'
!***********************************************************************
write(lu,'(a)') '-0.25 X  0.5 Y moveto -0.25 X -6.5 Y lineto'
write(lu,'(a)') '-0.40 X  0.0 Y moveto -0.25 X  0.0 Y lineto'
write(lu,'(a)') '-0.40 X -1.0 Y moveto -0.25 X -1.0 Y lineto'
write(lu,'(a)') '-0.40 X -2.0 Y moveto -0.25 X -2.0 Y lineto'
write(lu,'(a)') '-0.40 X -3.0 Y moveto -0.25 X -3.0 Y lineto'
write(lu,'(a)') '-0.40 X -4.0 Y moveto -0.25 X -4.0 Y lineto'
write(lu,'(a)') '-0.40 X -5.0 Y moveto -0.25 X -5.0 Y lineto'
write(lu,'(a)') '-0.40 X -6.0 Y moveto -0.25 X -6.0 Y lineto'
write(lu,'(a)') 'TF setfont'
write(lu,'(a)') '-0.45 X  0.0 Y 6 sub moveto (0.0) RightLine'
write(lu,'(a)') '-0.45 X -1.0 Y 6 sub moveto (\261 1.0) RightLine'
write(lu,'(a)') '-0.45 X -2.0 Y 6 sub moveto (\261 2.0) RightLine'
write(lu,'(a)') '-0.45 X -3.0 Y 6 sub moveto (\261 3.0) RightLine'
write(lu,'(a)') '-0.45 X -4.0 Y 6 sub moveto (\261 4.0) RightLine'
write(lu,'(a)') '-0.45 X -5.0 Y 6 sub moveto (\261 5.0) RightLine'
write(lu,'(a)') '-0.45 X -6.0 Y 6 sub moveto (\261 6.0) RightLine'
write(lu,'(a)') 'stroke'

return

end subroutine FigOpn
