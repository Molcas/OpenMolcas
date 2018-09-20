# -*- coding: utf-8 -*-
#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2015,2017, Ignacio Fdez. Galván                        *
#***********************************************************************

# Set global citation order
citeorder = {}

# Force updating references if anything else changes
def update_bib(app, env, added, changed, removed):
  updated = added.union(changed).union(removed)
  if (app.config.ref_file and (len(updated) > 0)):
    if (app.config.ref_file not in updated):
      return [app.config.ref_file]
  return []

# Put the references file at the end of the list
def bib_at_end(app, env, docnames):
  if (app.config.ref_file):
    try:
      docnames.remove(app.config.ref_file)
      docnames.append(app.config.ref_file)
    except ValueError:
      pass

# Build "citeorder" according to the citation order,
# following the document order in sphinx.
# This has to be done at source-read, because the
# bibliography is sorted before the doctree-read stage
def sort_citations(app, docname, source):
  if (docname == app.config.ref_file):
    env = app.builder.env
    cited = env.bibtex_cache._cited
    rel = env.collect_relations()
    # each element is [up, prev, next],
    # so start with the one with prev=None, and follow next
    doc = [x for x in rel if rel[x][1] is None][0]
    i = 0
    while (doc is not None):
      for cite in cited[doc]:
        if (cite not in citeorder):
          citeorder[cite] = i
          i = i+1
      doc = rel[doc][2]
 
def setup(app):
  app.connect('env-get-outdated', update_bib)
  app.connect('env-before-read-docs', bib_at_end)
  app.connect('source-read', sort_citations)
  app.add_config_value('ref_file', 'references', 'html')


### Bibliography Styles ###

from pybtex.style.sorting import BaseSortingStyle
from pybtex.style.formatting.unsrt import dashify, Style as UnsrtStyle
from pybtex.style.formatting import toplevel
from pybtex.style.template import join, words, field, optional, first_of, names, sentence, tag, optional_field, href
from pybtex.style.template import first_of, optional, join, field
from pybtex.plugin import register_plugin
from pybtex.richtext import Text, Symbol

# Trivial sorting style that just uses "citeorder" as key
class CiteStyle(BaseSortingStyle):
  def sorting_key(self, entry):
    try:
      return citeorder[entry.key]
    except:
      return 9999999

# Formatting style
class MolcasStyle(UnsrtStyle):

  default_sorting_style = 'cite'

  date = words [optional_field('month'), field('year')]

  def format_names(self, role, as_sentence=True):
    formatted_names = names(role, sep=', ', sep2 = ', ', last_sep=', ')
    if as_sentence:
      return sentence(capfirst=False) [formatted_names]
    else:
      return formatted_names

  def format_article(self, e):
    pages = first_of [
      # article id with total pages
      optional [
        join [
          field('articleid'),
          optional [
            '(',join[u'1–', optional_field('pagetotal')],')'
          ]
        ]
      ],
      # pages only
      field('pages', apply_func=dashify),
    ]
    volume_and_pages = first_of [
      # volume and pages, with optional issue number
      optional [
        join [
          tag('strong') [field('volume')],
          optional['[', field('number'),']'],
          ' (', field('year'), ')',
          ' ', pages
        ],
      ],
      # pages only
      words ['pages', pages],
    ]
    template = toplevel [
      self.format_names('author'),
      #self.format_title(e, 'title'),
      sentence(capfirst=False) [
        tag('em') [field('journal')],
        optional[ volume_and_pages ],
      ],
      sentence(capfirst=False) [
        optional_field('note')
      ],
      self.format_web_refs(e),
    ]
    return template.format_data(e)

  def format_inbook(self, e):
    template = toplevel [
      sentence [self.format_names('author')],
      sentence [
        #self.format_btitle(e, 'title'),
        self.format_chapter_and_pages(e),
      ],
      self.format_volume_and_series(e),
      sentence [
        field('publisher'),
        optional_field('address'),
        optional [
          words [field('edition'), 'edition']
        ],
        self.date,
        optional_field('note'),
      ],
      self.format_web_refs(e),
    ]
    return template.format_data(e)

  def format_incollection(self, e):
    template = toplevel [
      sentence [self.format_names('author')],
      #self.format_title(e, 'title'),
      words [
        'In',
        sentence(capfirst=False) [
          optional[ self.format_editor(e, as_sentence=False) ],
          self.format_btitle(e, 'booktitle', as_sentence=False),
          self.format_volume_and_series(e, as_sentence=False),
          self.format_chapter_and_pages(e),
        ],
      ],
      sentence [
        optional_field('publisher'),
        optional_field('address'),
        self.format_edition(e),
        self.date,
      ],
      self.format_web_refs(e),
    ]
    return template.format_data(e)

  def format_doi(self, e):
    return href [
      join [
        'https://doi.org/',
        field('doi')
      ],
      join [
        'doi:',
        field('doi')
      ]
    ]

register_plugin('pybtex.style.sorting', 'cite', CiteStyle)
register_plugin('pybtex.style.formatting', 'molcas', MolcasStyle)
