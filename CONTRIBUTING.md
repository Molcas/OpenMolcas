Attribution
===========

If you add a new file, please include the OpenMolcas header and copyright
information, for example:

```
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
* Copyright (C) 2014,2016, John Doe                                    *
*               2017, Jane Doe                                         *
************************************************************************
```

* Use full names when possible
* Use a single line for each author, and a single author on each line
* Include all years where modifications were done, use ranges only if there
  were modifications in all the intermediate years.

If you modify an existing file, and the modification is significant enough, add
your name to the header in the same format, or the year if your name is already
there. What exactly is "significant enough" is rather subjective; changing a
format to print more decimal figures is probably not, rearranging the loop
structure can be. It is up to you to decide whether or not your modifications
are significant enough.

If you copy a file, copy the header too, and then do as above if you modify it.
Do not remove any name unless you are absolutely sure that all contributions by
that person have been removed from the file.

You can include more detailed information about yourself, the origin of any
external code, the nature of the modifications, etc. below the header.

Finally, add any new names to the list in `CONTRIBUTORS.md`. Note that in order
to maintain the formatting with online tools (see
[here](https://gitlab.com/Molcas/OpenMolcas/blob/master/CONTRIBUTORS.md)), each
name is a separate line that *ends with two spaces*. The list is sorted
alphabetically by last name (using "international" sorting order, where Á = A,
Ö = O etc.).


Commit messages
===============

When writing commit messages with `git`, try to be informative but brief. The
best is to follow a simple format (see
[here](https://chris.beams.io/posts/git-commit/)):

1. One line with a short title.
2. A blank line.
3. Some explanation of what was modified and why it was needed or useful.

For example:

```
Start work on new_shiny_feature

The old_dull_feature has some limitations, this new_shiny_feature
will solve them. For the moment just add a dummy keyword.
```

In many tools and on the [online
viewer](https://gitlab.com/Molcas/OpenMolcas/commits/master) only the first
line is displayed by default and the rest of the message is available if
desired.
