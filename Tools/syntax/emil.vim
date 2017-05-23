" Vim syntax file
"***********************************************************************
" This file is part of OpenMolcas.                                     *
"                                                                      *
" OpenMolcas is free software; you can redistribute it and/or modify   *
" it under the terms of the GNU Lesser General Public License, v. 2.1. *
" OpenMolcas is distributed in the hope that it will be useful, but it *
" is provided "as is" and without any express or implied warranties.   *
" For more details see the full text of the license in the file        *
" LICENSE or in <http://www.gnu.org/licenses/>.                        *
"                                                                      *
" Copyright (C) 2013, Steven Vancoillie                                *
"               2017, Ignacio Fdez. Galván                             *
"***********************************************************************
" Language:	Molcas input (*.input)
" Version:	8.1
" Last Change:	2017 January 11
" Maintainer:	Ignacio Fdez. Galván <Ignacio.Fernandez@kemi.uu.se>
"
" shell variables, these are substituted by emil
syn match molcasVariable "\$\w\+"
" emil commands, starts with 1 or more '>'
syn match molcasEmil "^[ \t]*>\+[ \t]*\w\+"
" module, word starting with '&'
syn match molcasModule "^[ \t]*&\w*"
" comments: lines starting with '*', from '//' to the end of line, and regions between '/*' and '*/'
syn match molcasComment excludenl "^[ \t]*\*.*$"
syn match molcasComment excludenl "\/\/.*$"
syn region molcasComment start="/\*" end="\*/"
" numbers, taken from fortran syntax
syn match molcasNumber	display "\(\<\|[-+]\)\d\+\(\.\d\+\)\=\([dDeE][-+]\=\d\+\)\="

" set the colors for the syntax elements
hi def link molcasVariable Identifier
hi def link molcasEmil Statement
hi def link molcasModule Special
hi def link molcasComment Comment
hi def link molcasNumber Number
