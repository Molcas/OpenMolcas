.. _TUT\:sec\:environment:

Environment and EMIL Commands
=============================

The following are basic and most common commands for the |molcas| environment variables and input language (EMIL):

.. class:: variablelist

:variable:`MOLCAS`
  |molcas| home directory.

:variable:`MOLCAS_MEM`
  Memory definition in Mb. Default 1024.

:variable:`MOLCAS_PRINT`
  Printing level: 2 Normal, 3 Verbose

:variable:`MOLCAS_PROJECT`
  Name used for the project/files.

:variable:`MOLCAS_WORKDIR`
  Scratch directory for intermediate files.

.. class:: commandlist

:command:`>> Do While`
  Start of a loop in an input file for geometry optimization with conditional termination.

:command:`>> Foreach`
  Start of a loop in an input file over a number of items.

:command:`>> EndDo`
  End of a loop in an input file.

:command:`>> If ( condition )`
  Start of If block.

:command:`>> EndIf`
  End of If block.

:command:`>> Label Mark`
  Setting the label "Mark" in the input.

:command:`>> Goto Mark`
  Forward jump to the label "Mark" skipping that part of the input.
