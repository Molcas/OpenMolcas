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
" Copyright (C) 2013-2016,2020, Ignacio Fdez. Galván                   *
"***********************************************************************
" Language:     Molcas output (*.log)
" Version:      7.9
" Last Change:  2020 March 4
" Maintainer:   Ignacio Fdez. Galván <Ignacio.Fernandez@kemi.uu.se>

" Folds can be open and closed with "zo" and "zc"
" Or all folds at once with "zR" and "zM"
" Disable automatic folding with ":set nofoldenable" (enable with ":set foldenable")
" This assumes ":set foldmethod=syntax"

" Header (license, logo, version...)
syn region molcasHeader
      \ start="^\s*This \(copy\|run\) of MOLCAS"
      \ end="^\n\n"
      \ fold

" Make folds for each module (except auto)
syn region molcasFoldModule
      \ start="^--- Start Module:\(\s*auto\)\@!.*$"
      \ end="^--- Stop Module:\(\s*auto\)\@!.*$"
      \ fold
      \ contains=molcasFoldBlock,molcasResult,molcasMessage
      \ keepend

" Make folds for each block of data
syn region molcasFoldBlock
      \ start="^++ .*$"
      \ end="^--\( .*\)\?$"
      \ end="^--- Stop Module:\(\s*auto\)\@!.*$"
      \ fold
      \ contains=molcasFoldBlock,molcasResult,molcasMessage

" Info about links, etc
syn match molcasInfo "^\*\*\*.*\(\n    .*\)*$" contains=molcasFile
syn match molcasFile "\(\*\*\* files:.*\)\@<=\<[^ ]*\>"
syn match molcasFile "\(          .*\)\@<=\<[^ ]*\>" contained
syn match molcasFile "\(\*\*\* symbolic link created:.*\)\@<=\<[^ ]*\>"

" Timing info
syn match molcasTiming "^--- Module\s\+.*\s\+spent\s\+.*$"
syn match molcasTiming "^\s*Timing:.*"

" Identify important results
syn match molcasResult "^::.*$"

" Emil commands in the output
syn match molcasEmil "^\s*>\+\s*.*$"

" Messages and warnings
"syn match molcasMessage "\(^ ###.*\n\)\+"

" Errors
syn match molcasError "^Non-zero return code.*"
syn match molcasError "^ *\.\+\n\(^ *\.\.\..*\.\n\)*"

" Return code
syn match molcasReturn "^\.#\+\.\n\(^\.#.*#\.\n\)*\.#\+\.\n"

" Set colors
hi link molcasHeader PreProc
hi link molcasInfo Comment
hi link molcasFile Special
hi link molcasTiming Identifier
hi link molcasResult Number
hi link molcasEmil Statement
hi link molcasMessage ToDo
hi link molcasError Error
hi link molcasReturn Identifier
" Large parts of the output colored only if debugging
if exists("molcas_output_debug")
  hi link molcasFoldModule Statement
  hi link molcasFoldBlock Type
endif

syn sync fromstart

" Keep the first line as title of the fold
set foldtext=v:folddashes.substitute(getline(v:foldstart),'','','g') 
