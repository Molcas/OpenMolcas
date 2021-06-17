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
! Copyright (C) 2000, Per-Olof Widmark                                 *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

function ixNuclearChargeFromSymbol(Symbol,Rc,Opt)
!***********************************************************************
!                                                                      *
! This function returns the nuclear charge of an atom based on the     *
! chemical symbol.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
! Written: Feb. 2000                                                   *
! Modified: March 2017, Ignacio Fdez. Galvan (use periodic_table.fh)   *
!                                                                      *
!***********************************************************************

implicit none
#include "proputil.fh"
integer ixNuclearChargeFromSymbol
!----------------------------------------------------------------------*
! Parameters.                                                          *
!----------------------------------------------------------------------*
integer StopOnError
parameter(StopOnError=_OPT_STOP_ON_ERROR_)
!----------------------------------------------------------------------*
! Dummy parameters.                                                    *
!----------------------------------------------------------------------*
character*(*) Symbol
integer Rc
integer Opt
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
#include "periodic_table.fh"
character*2 Sym1, Sym2
integer Index
integer i
!----------------------------------------------------------------------*
! External references.                                                 *
!----------------------------------------------------------------------*
external UpCase

!----------------------------------------------------------------------*
! Locate symbol in table.                                              *
!----------------------------------------------------------------------*
Index = 0
Sym1 = adjustl(Symbol)
call UpCase(Sym1)
do i=1,Num_Elem
  Sym2 = adjustl(PTab(i))
  call UpCase(Sym2)
  if (Sym1 == Sym2) Index = i
end do
!----------------------------------------------------------------------*
! Are we successful.                                                   *
!----------------------------------------------------------------------*
if (Index == 0) then
  write(6,'(a)') '***'
  write(6,'(a)') '*** NuclearChargeBySymbol: error'
  write(6,'(2a)') '***    unknown atom: ',Symbol
  write(6,'(a)') '***'
  if (iand(Opt,StopOnError) /= 0) call Quit_OnUserError()
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
ixNuclearChargeFromSymbol = Index

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(Rc)

end function ixNuclearChargeFromSymbol
