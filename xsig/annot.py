r"""
Genome annotation table (:mod:`xsig.annot`)
===========================================

.. currentmodule:: xsig.annot

An annotation table provides a fast access to a limited view of a genome
annotation. It stores only limited information about the genes:

- EntrezID
- official symbol
- alternate symbol
- chromosomal position
- synonyms
- gene title (short description)

"""
from __future__ import absolute_import, division, print_function

__author__ = 'Vlad Popovici'
__version__ = 0.2

import pandas as pd

class Annotation(object):
    r"""
    A minimalistic class for storing genome annotations.
    """

    def __init__(self, v=0):
        self.__version = v
        self.__data = None
        self.__gene_to_id = None


    @property
    def version(self):
        return self.__version


    @version.setter
    def version(self, v):
        self.__version = str(v)


    def read_from_table(self, fname, ver):
        r"""
        Reads data from a tab-delimited file.

        Parameters
        ----------
        fname: file name, str
            The file containing the annotation.
        ver: annotation version, str
            A version identifier (usually the date when the annotation was generated).
        """

        self.__data = pd.read_table(fname, header=None,
                                    names=['EntrezID','Symbol', 'Symbol_alt', 'Chr_pos', 'Synonyms', 'Title'],
                                    na_values=['-'], dtype=str)
        self.__data.index = self.__data.EntrezID    # keep it as unicode
        self.version = ver
        self.__gene_to_id = dict(zip(self.__data.Symbol, self.__data.EntrezID))  # reversed indexing


    def save(self, fname):
        r"""
        Saves the annotation in a HDF5 file, under a key of the form "annot-VERSION".

        Parameters
        ----------
        fname: file name, str
            The name of the file to store the annotation. Note that '.h5' extension is automatically
            added.
        """
        self.__data.to_hdf(fname, 'annot_'+self.__version, mode='a', format='fixed', complib='zlib', complevel=9)


    def load(self, fname, ver):
        r"""
        Loads an annotation from a HDF5 file, from the key "annot-VERSION", where the version is
        specified by the user.

        Parameters
        ----------
        fname: file name, str
            The name of the HDF5 file.
        ver: version, str
            Annotation version.
        """
        self.version = ver
        self.__data = pd.read_hdf(fname, 'annot-'+ver)
        self.__gene_to_id = dict(zip(self.__data.Symbol, self.__data.EntrezID))  # reversed indexing


    def id_from_symbol(self, gene):
        return self.__gene_to_id[gene]


    def gene_from_id(self, id):
        return self.__data.Symbol[id]


    def info_id(self, id):
        return self.__data.ix[id]


    def info_gene(self, gene):
        return self.__data.ix[self.id_from_symbol(gene)]


    def id_to_gene_map(self):
        return dict(self.__data.index, self.__data.Symbol)


    def gene_to_id_map(self):
        return self.__gene_to_id

    ## end Annotation class