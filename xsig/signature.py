r"""
Gene signature (:mod:`xsig.signature`)
======================================

.. currentmodule:: xsig.signature

A gene expression signature is a set of genes that are supposed to have a
correlated impact on the phenotype. The genes in the set can be weighted
and an "activity score" can be computed for the gene set. Examples of
such gene sets/signatures are the pathways or predictive/prognostic gene
signatures.

This module assumes that gene expression data is represented in terms of
Pandas DataFrame objects, with genes (or probesets, probes, etc) given by
columns and samples by rows.
"""

from __future__ import absolute_import, division, print_function

__author__ = 'Vlad Popovici'
__version__ = 0.1

import numpy as np
import pandas as pd

from pandas.core.series import Series
from pandas.core.frame import DataFrame


class GeneSignature(object):
    r"""
    A gene signature consists of a list of genes and, optionally, a list
    of coefficients which are combined to produce a score. This score can
    be

    - a gene set enrichment score (like in [#GSEA]_ methods)
    - a prognostic score (a continuous value that can be linked to survival)
    - a predictive score (usually a binary value predicting a class
    membership)
    - etc


    References
    ----------
    .. [#GSEA] Subramanian A et al. Gene set enrichment analysis: a knowledge-based
    approach for interpreting genome-wide expression profiles. PNAS 102(43):15545-50,
    2005. PMID: 16199517
    """

    def __init__(self, gid, w=None, b=None, info=''):
        r"""
        Initialize a basic gene signature with a set of genes and optional weights
        and bias (threshold).

        Parameters
        ----------
        gid: set of genes, list-like
            A list of gene EntrezIDs.
        w: weights, array
            A vector of weights. No normalization is performed. If not specified, it
            defaults to a vector of 1s.
        b: bias, float
            A bias (threshold) to be used for binarizing the score. It defaults to 0.
        info: about signature, string
            Some additional information about the signature - e.g. an URL etc.

        Internally, the gene set is stored as a pandas.core.series.Series object.
        """

        if b is None:
            self._bias = 0.0
        else:
            self._bias = b

        if w is None:
            self._gset = Series(dict(zip(gid, np.ones(len(gid), dtype=np.float64))), copy=True)
        else:
            if isinstance(w, list):
                if len(w) == 1:
                    _w = np.ndarray((len(gid), ), dtype=np.float64)
                    _w.fill(w[0])
                    self._gset = Series(dict(zip(gid, _w)), copy=True)
                elif len(w) == len(gid):
                    self._gset = Series(dict(zip(gid, w)), dtype=np.float64, copy=True)
                else:
                    raise ValueError('length of weights does not match the length of the gene list')
            else:
                # assume w is a constant - try to initialize the series:
                _w = np.ndarray((len(gid), ), dtype=np.float64)
                _w.fill(w)
                self._gset = Series(dict(zip(gid, _w)), copy=True)

        self.info = info

        return


    @property
    def gset(self):
        return self._gset


    @property
    def bias(self):
        return self._bias


    def genes(self):
        return self._gset.index.tolist()


    def remap(self, old_to_new):
        r"""
        Updates the symbols/gene IDs from the signature.

        Parameter
        ---------
        old_to_new: IDs map, dict
            The mapping from old symbols to new symbols, given as a dictionary,
            {old: new,...}

        Returns
        -------
            A list of symbols that could not be remapped, or None if all symbols
            could be updated.
        """

## end class GeneSignature


class GSS:
    r"""
    Gene Set/Signatures Scores: a namespace for various methods of computing scores
    for gene sets.
    """

    @staticmethod
    def average(gs, X, weighted=False):
        r"""
        Simple averaging of the genes in the set.
        """
        if not isinstance(X, pd.core.frame.DataFrame) and not isinstance(X, pd.core.series.Series):
            raise ValueError('X must be a Series or DataFrame instance')

        if weighted:
            s = gs.gset.sum()
            return X[gs.genes()].dot(gs.gset) / s
        else:
            if isinstance(X, pd.core.frame.DataFrame):
                # a set of samples, given in an expression matrix
                return X[gs.genes()].mean(axis=1)
            else:
                return X[gs.genes()].mean()

        return None
