"***********************************************************************
" This file is part of OpenMolcas.                                     *
"                                                                      *
" OpenMolcas is free software; you can redistribute it and/or modify   *
" it under the terms of the GNU Lesser General Public License, v. 2.1. *
" OpenMolcas is distributed in the hope that it will be useful, but it *
" is provided "as is" and without any express or implied warranties.   *
" For more details see the full text of the license in the file        *
" LICENSE or in <http://www.gnu.org/licenses/>.                        *
"***********************************************************************
if did_filetype() " filetype already set...
  finish          " ...don't do these checks
endif
for l in getline(1,10)
  if (l =~ '^\s*This \(copy\|run\) of MOLCAS')
    setfiletype molcasoutput
    break
  endif
endfor
