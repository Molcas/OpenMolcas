" Vim syntax file
" Language:	Molcas input (*.input)
" Version:	8.1
" Last Change:	2022 November 1
" Maintainer:	Ignacio Fdez. Galván <Ignacio.Fernandez@kemi.uu.se>
"
" shell variables, these are substituted by emil
syn match molcasVariable "\$\w\+"
" emil commands, starts with 1 or more '>'
syn match molcasEmil "^[ \t]*>\+[ \t]*\(\(do\|end\|go\)\c[ \t]\+\)\?\w\+"
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
