r"""
I/O module for signatures (:mod: `xsig.io`)
===========================================

.. currentmodule:: xsig.io

Gene expression signatures are stored in various databases like
MSigDB_ and GeneSigDB_ and in formats that are more or less alike. This
module provides access to several such resources and allows fetching and
local storing of signatures.

References
----------
.. MSigDB http://www.broadinstitute.org/gsea/msigdb
.. GeneSigDB http://compbio.dfci.harvard.edu/genesigdb/
"""
from __future__ import absolute_import, division, print_function

__author__ = 'Vlad Popovici'
__version__ = 0.1

import xml.etree.ElementTree as ET
from .signature import GeneSignature
from .annot import Annotation


def load_xml_msigdb(fname):
    r"""
    Loads a whole database from an XML file compliant with MSigDB_ specifications.
    The user should download such a file from MSigDB_ site - registration required.

    Note: MSigDB groups the signatures into seven collections (currently), each
    with its specific domain. Here, however, all the signatures will be returned
    as a dictionary where the keys begin with 'cx_' (the collection ID from
    MSigDB_).

    To use a particular signature from the collection, you need to know its ID.
    For this, search for the signature on the MSigDB_ web site and use the
    collection and systematic name to compose the signature::

        collection_systematicname

    For example, NUCLEOPLASM signature has a collection name C5 and a systematic
    name M12345, thus its ID will
    be::

        C5_M12345

    This function will also filter out all signatures that are not related to
    human genome (i.e. ORGANISM must be 'Homo sapiens').

    Parameters
    ----------
    fname: file name, string
        XML file name with MSigDB data

    Return
    ------
        .. a dictionary with the gene signatures
        .. a list of discarded signatures

    References
    ----------

    .. _MSigDB: http://www.broadinstitute.org/gsea/msigdb
    """

    xml_file = ET.parse(fname)
    xml_root = xml_file.getroot()

    db = {}
    discarded = []

    for sig in xml_root.findall('GENESET'):
        if sig.attrib['ORGANISM'] != 'Homo sapiens':
            discarded.append(sig.attrib['SYSTEMATIC_NAME'])
            continue
        gid = sig.attrib['MEMBERS_EZID'].split(',')
        url = 'http://www.broadinstitute.ors/gsea/msigdb/geneset_page.jsp?geneSetName=' + \
              sig.attrib['STANDARD_NAME']
        key = sig.attrib['CATEGORY_CODE'] + '_' + sig.attrib['SYSTEMATIC_NAME']
        db[key] = GeneSignature(gid, info=url)

    return db, discarded


# end load_xml_msigdb



def load_gmt_genesigdb(fname, annot):
    r"""
    Loads a whole database from a file formatted according the MSigDB's GMT_ specifications.
    There are a number of transformations applied to the data, as follows:
    -the signature ID from the original table is prepended an "GS_"
    -the signatures not referring to Homo sapiens are filtered out
    -the gene symbols are converted to EntrezID using the annotation passed as parameter

    Parameters
    ----------
    fname: file name, string
        GMT file with all signatures
    annot: annotation object, xsig.Annotation
        An annotation object used for converting gene symbols to EntrezIDs.

    Return
    ------
        .. a dictionary with the gene signatures
        .. a list of discarded signatures

    References
    ----------
    .. _GMT: http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
    """

    f = open(fname, 'r')

    db = {}
    discarded = []

    for line in f:
        # parse the line: signature ID, organism, list of symbols
        tokens = line.split('\t')
        if not tokens[1].startswith('Human'):
            discarded.append(tokens[0])
            continue
        key = 'GS_' + tokens[0]
        url = 'http://compbio.dfci.harvard.edu/genesigdb/publicationSearch.jsp?searchQuery=' + \
              tokens[0]
        # last name in the list may end with '\r\n'...
        tokens[-1] = tokens[1].strip()
        gid = []
        try:
            gid = [annot.id_from_symbol(g) for g in tokens[2:]]
        except KeyError:
            # if some symbols cannot be matched, the signature is discarded
            discarded.append(tokens[0])
            continue

        db[key] = GeneSignature(gid, info=url)

    return db, discarded
