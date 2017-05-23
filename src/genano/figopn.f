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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine prints the prologue part to the postscript figure of    *
* occupation numbers.                                                  *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine FigOpn(lu)
      Write(lu,'(a)') '%!PS'
      Write(lu,'(a)') '%%BoundingBox: 70 95 550 735'
      Write(lu,'(a)') '/TF /Helvetica findfont 18 scalefont def'
      Write(lu,'(a)') '/SF /Symbol findfont 18 scalefont def'
      Write(lu,'(a)') '/HF /Helvetica findfont 24 scalefont def'
      Write(lu,'(a)') '/Y {0.5 sub neg 75 mul} def'
      Write(lu,'(a)') '/X {75 mul} def'
      Write(lu,'(a)') '/Draw {'
      Write(lu,'(a)') '  newpath'
      Write(lu,'(a)') '  2.0 setlinewidth'
      Write(lu,'(a)') '  moveto'
      Write(lu,'(a)') '  0.5 X 0 rlineto'
      Write(lu,'(a)') '  currentpoint /y exch def /x exch def'
      Write(lu,'(a)') '  stroke'
      Write(lu,'(a)') '  x y moveto'
      Write(lu,'(a)') '} def'
      Write(lu,'(a)') '/Label {'
      Write(lu,'(a)') '  TF setfont'
      Write(lu,'(a)') '  6 -6 rmoveto show'
      Write(lu,'(a)') '} def'
      Write(lu,'(a)') '/Connect {'
      Write(lu,'(a)') '  [5 5] 0 setdash'
      Write(lu,'(a)') '  newpath'
      Write(lu,'(a)') '  0.5 setlinewidth'
      Write(lu,'(a)') '  moveto'
      Write(lu,'(a)') '  exch 0.5 X add exch'
      Write(lu,'(a)') '  lineto'
      Write(lu,'(a)') '  stroke'
      Write(lu,'(a)') '  [] 0 setdash'
      Write(lu,'(a)') '} def'
      Write(lu,'(a)') '/CenterLine {'
************************************************************************
      Write(lu,'(a)') '  dup stringwidth pop 0.5 mul neg 0 rmoveto show'
      Write(lu,'(a)') '} def'
      Write(lu,'(a)') '/RightLine {'
      Write(lu,'(a)') '  dup stringwidth pop neg 0 rmoveto show'
      Write(lu,'(a)') '} def'
      Write(lu,'(a)') '%.ManualFeed'
      Write(lu,'(a)') 'gsave'
      Write(lu,'(a)') '144 144 translate'
      Write(lu,'(a)') '-0.30 X -6.6 Y moveto'
      Write(lu,'(a)') 'TF setfont (lg\050) CenterLine'
      Write(lu,'(a)') 'SF setfont (h) show'
      Write(lu,'(a)') 'TF setfont (\051) show'
      Write(lu,'(a)') '%---'
      Write(lu,'(a)') 'TF setfont'
      Write(lu,'(a)') '0.25 X 1.0 Y moveto (s) CenterLine'
      Write(lu,'(a)') '1.25 X 1.0 Y moveto (p) CenterLine'
      Write(lu,'(a)') '2.25 X 1.0 Y moveto (d) CenterLine'
      Write(lu,'(a)') '3.25 X 1.0 Y moveto (f) CenterLine'
      Write(lu,'(a)') '4.25 X 1.0 Y moveto (g) CenterLine'
      Write(lu,'(a)') '%---'
      Write(lu,'(a)') '3.0 setlinewidth'
      Write(lu,'(a)') 'newpath'
************************************************************************
      Write(lu,'(a)') '-0.25 X  0.5 Y moveto -0.25 X -6.5 Y lineto'
      Write(lu,'(a)') '-0.40 X  0.0 Y moveto -0.25 X  0.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -1.0 Y moveto -0.25 X -1.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -2.0 Y moveto -0.25 X -2.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -3.0 Y moveto -0.25 X -3.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -4.0 Y moveto -0.25 X -4.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -5.0 Y moveto -0.25 X -5.0 Y lineto'
      Write(lu,'(a)') '-0.40 X -6.0 Y moveto -0.25 X -6.0 Y lineto'
      Write(lu,'(a)') 'TF setfont'
      Write(lu,'(a)') '-0.45 X  0.0 Y 6 sub moveto (0.0) RightLine'
      Write(lu,'(a)') '-0.45 X -1.0 Y 6 sub moveto (\261 1.0) RightLine'
      Write(lu,'(a)') '-0.45 X -2.0 Y 6 sub moveto (\261 2.0) RightLine'
      Write(lu,'(a)') '-0.45 X -3.0 Y 6 sub moveto (\261 3.0) RightLine'
      Write(lu,'(a)') '-0.45 X -4.0 Y 6 sub moveto (\261 4.0) RightLine'
      Write(lu,'(a)') '-0.45 X -5.0 Y 6 sub moveto (\261 5.0) RightLine'
      Write(lu,'(a)') '-0.45 X -6.0 Y 6 sub moveto (\261 6.0) RightLine'
      Write(lu,'(a)') 'stroke'
      Return
      End
