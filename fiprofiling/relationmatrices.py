#!/usr/bin/env python3
import warnings

import numpy as np
import pandas as pd

from fiprofiling import fingerprints

class RM:
    def __init__(self, dataframe=pd.DataFrame(), *args,
                 num_datapoints=None, name=None):
        self.df = dataframe
        self.num_datapoints = num_datapoints
        if name:
            self.name = name
        else:
            self.name = "Unnamed " + type(self).__name__

    @staticmethod
    def csv2df__count(csv):
        df = pd.read_csv(csv, comment='#')
        df.set_index('Unnamed: 0', inplace=True)
        df.columns = [str(x) for x in df.columns]
        df.index = [str(x) for x in df.index]
        datapoint_count = None
        # try to find the integer as comment in the file.
        # TODO: come up with a more sensible way of handling metadata in csv
        csv.seek(0)
        for line in csv:
            if line.startswith('#'):
                datapoint_count = int(line.lstrip('#'))
                break
        return (df, datapoint_count)

    @classmethod
    def from_csv(cls, csv, name=None):
        df, num_datapoints = cls.csv2df__count(csv)
        return cls(df, num_datapoints=num_datapoints, name=name)

    @classmethod
    def from_fingerprints(cls, fps):
        raise NotImplementedError

    @classmethod
    def from_dataframe(cls, df, *args, num_datapoints=None, name=None):
        return cls(df, num_datapoints=num_datapoints, name=name)

    @classmethod
    def _otheRM_typecheck(cls, otheRM):
        if not isinstance(otheRM, cls):
            raise TypeError(
                "The provided other relation matrix needs to be of {0} type, not {1}".format(
                    [type(cls), type(otheRM)]))
        return True

    def is_square(self):
        x, y = self.df.shape
        return x == y

    def offdiag_values(self, *args, skipna=True):
        """Yields consecutively all off-diagonal values of the relation matrix.
        If skipna=True (default), all NaN values are skipped"""
        x, y = self.df.shape
        if not skipna:
            for i in range(x):
                for j in range(y):
                    if i != j:
                        yield self.df.iloc[i, j]
        else:
            for i in range(x):
                for j in range(y):
                    if i != j:
                        value = self.df.iloc[i, j]
                        if np.isfinite(value):
                            yield value

    def mean(self, *args, offdiag_only=True):
        if offdiag_only:
            offdiags = list(self.offdiag_values())
            return float(sum(offdiags)) / len(offdiags)
        else:
            return float(self.df.sum().sum()/self.df.count().sum())

    def triangular(self, strict=False):
        offset = int(strict) # strict ~ 1, non-strict ~ 0
        mask = np.triu(np.ones(self.df.shape), k=offset).astype(np.bool)
        return self.__class__(self.df.where(mask),
                              num_datapoints=self.num_datapoints,
                              name=self.name)

    def difference(self, otheRM):
        self._otheRM_typecheck(otheRM)
        diffname = self.name + ' - ' + otheRM.name
        difference_df = self.df - otheRM.df
        return self.__class__(difference_df,
                              num_datapoints=min((self.num_datapoints,
                                                  otheRM.num_datapoints)),
                              name=diffname)

    def distance(self, otheRM, *args, imputation='none', type='ratio',
                    minval=-1, maxval=1):
        """Returns a distance value between this relational matrix and the other
        provided as an argument"""
        diffmatrix = self.difference(otheRM)
        if imputation and imputation.lower() != 'none':
            if imputation == 'mean':
                fillval = diffmatrix.mean(offdiag_only=False)
                diffmatrix.df.fillna(fillval, inplace=True)
            elif imputation == 'zero':
                fillval = 0
                diffmatrix.df.fillna(fillval, inplace=True)
            else:
                raise NotImplementedError(
                    "Imputation method '{0}' is not implemented.".format(imputation))
        if type == 'euclidean':
            dist = np.sqrt(np.square(diffmatrix.df).values.sum())
        elif type == 'ratio':
            #Can't use straight mean of absolute values + skipna=True, bacause
            #the difference between NPMI values can be as much as 2 per relation
            diffsum = diffmatrix.df.abs().sum().sum() # sum of all differences
            filledvals = diffmatrix.df.count().sum() # sum of non-nan values
            # sum(|diffrm|) / (sum(non-NaN values)*max_possible_difference)
            dist = diffsum / (filledvals * (maxval - minval))
        return dist

    def is_equal(self, otheRM):
        self._otheRM_typecheck(otheRM)
        return self.df.equals(otheRM.df) \
            and self.num_datapoints == otheRM.num_datapoints

    def is_symmetric(self):
        return (self.df.transpose() == self.df).all().all()

    def to_csv(self, csvfile):
        self.df.to_csv(csvfile)
        csvfile.write("#{0}".format(self.num_datapoints))


class CORM(RM):
    @classmethod
    def from_fingerprints(cls, fps, *args, fpformat=None, name=None, labels=None):
        if not fpformat:
            # convert straight to numpy array
            arrayable_fps = (np.asarray(fp) for fp in fps)
        elif fpformat == 'bintext':
            arrayable_fps = (fingerprints.fp2nparray(fingerprints.bintext2fp(fp)) for fp in fps)
        elif fpformat == 'hextext':
            arrayable_fps = (fingerprints.fp2nparray(fingerprints.hextext2fp(fp)) for fp in fps)
        else:
            raise ValueError("Unknown fpformat option", fpformat)
        first_array = next(arrayable_fps)
        if not labels:
            labels = [str(i) for i in range(len(first_array))]
        coorm = cls(pd.DataFrame(np.outer(first_array, first_array),
                                 columns=labels,
                                 index=labels),
                    num_datapoints=1,
                    name=name)
        for arrayable_fp in arrayable_fps:
            coorm.add_arrayable_fingerprint(arrayable_fp)
        return coorm

    def add_fingerprint(self, fp, *args, fpformat=None):
        if not fp:
            warnings.warn(RuntimeWarning("Encountered empty-looking fingerprint string. Skipping:", fp))
        elif not fpformat:
            self.add_arrayable_fingerprint(np.asarray(fp))
        elif fpformat == 'bintext':
            self.add_arrayable_fingerprint(fingerprints.fp2nparray(fingerprints.bintext2fp(fp)))
        elif fpformat == 'hextext':
            self.add_arrayable_fingerprint(fingerprints.fp2nparray(fingerprints.hextext2fp(fp)))
        else:
            raise ValueError("Unknown fpformat option", fpformat)

    def add_arrayable_fingerprint(self, afp):
        outer_product = np.outer(afp, afp)
        self.df += outer_product.astype(np.int)
        self.num_datapoints += 1

    def __add__(self, other):
        datapoints = self.num_datapoints + other.num_datapoints
        return self.from_dataframe(self.df + other.df,
                                   num_datapoints=datapoints,
                                   name=self.name)


class COPRM(RM):
    # PROBRM = COORM/sum(all datapoints)
    @classmethod
    def from_CORM(cls, coorm):
        probability_df = coorm.df.divide(coorm.num_datapoints)
        return cls(probability_df, num_datapoints=coorm.num_datapoints,
                   name=coorm.name)


class IndependentCOPRM(RM):
    """Expected probability of occurrence with the assumption that
    all observed variables are independent."""
    @classmethod
    def from_CORM(cls, coorm):
        return cls.from_COPRM(COPRM.from_CORM(coorm))

    @classmethod
    def from_COPRM(cls, probrm):
        iprm_df = cls.probrm_df2iprm_df(probrm.df)
        return cls(iprm_df, num_datapoints=probrm.num_datapoints,
                   name=probrm.name)

    @staticmethod
    def probrm_df2iprm_df(probrm_df):
        diagonal = np.diag(probrm_df)
        # get all statistical combinations of the individual probabilities
        iprm_df = pd.DataFrame(np.outer(diagonal, diagonal),
                               columns=probrm_df.columns.values,
                               index=probrm_df.index.values)
        # apply back the diagonal, because P(A AND A) = P(A), not P(A)*P(A)
        for i, diagvalue in enumerate(diagonal):
            iprm_df.iloc[i, i] = diagvalue
        return iprm_df


class PMIRM(RM):
    # PMI = log2(p(a, b) / (p(a)p(b)))
    # PMI = log2(COPRM / IndependentCOPRM)
    @staticmethod
    def probability_df2pmi_df(prob_df):
        indep_prob_df = IndependentCOPRM.probrm_df2iprm_df(prob_df)
        # replace zeros for NaNs to avoid subsequent division by zero problems
        indep_prob_df.replace(0, np.nan, inplace=True)
        ratio_df = prob_df / indep_prob_df
        # replace zeros for NaNs to avoid log(0) = -inf values in result
        ratio_df.replace(0, np.nan, inplace=True)
        pmi_df = np.log2(ratio_df)
        return pmi_df

    @classmethod
    def from_CORM(cls, coorm):
        return cls.from_COPRM(COPRM.from_CORM(coorm))

    @classmethod
    def from_COPRM(cls, probrm):
        pmi_df = cls.probability_df2pmi_df(probrm.df)
        return cls(pmi_df, num_datapoints=probrm.num_datapoints,
                   name=probrm.name)

    def tightness(self, probrm, *args, raw=False):
        """Takes a COPRM as an argument.
        Returns the 'tightness' of the provided COPRM according to this PMIRM variant.
        If raw=False (default), the provided probability matrix is applied in a formatted
        way, where only events/pattern with probability over 0 count, the rest goes completely ignored:
        If raw=True, there is no provided probability matrix formatting at all:
        tightness = mean(COPRM * NormalizedPMIRM)"""
        if raw:
            scores = probrm.df * self.df
        else:
            scores = probrm.df.clip(lower=0, upper=1).replace(0, np.NaN) * self.df
        meanscore = scores.stack().mean(skipna=True)
        return meanscore

    def fp_tightness(self, fingerprint, *args, fpformat=None):
        """Returns the 'tightness' of a provided fingerprint relative to this PMIRM variant
        by converting it to PROBRM and using the tightness function."""
        probrm = COPRM.from_CORM(
                    CORM.from_fingerprints([fingerprint], fpformat=fpformat))
        tightness = self.tightness(probrm)
        return tightness


class ZPMIRM(PMIRM):
    # ZPMI = Zscores(PMI), dataframe-wide, i.e. not by columns
    # ZPMI = (PMI - mean(PMI)) / std(PMI)

    @classmethod
    def from_CORM(cls, coorm):
        return cls.from_PMIRM(PMIRM.from_CORM(coorm))

    @classmethod
    def from_COPRM(cls, probrm):
        return cls.from_PMIRM(PMIRM.from_COPRM(probrm))

    @classmethod
    def from_PMIRM(cls, pmirm):
        offdiag_values = tuple(pmirm.offdiag_values(skipna=True))
        mean = float(sum(offdiag_values)) / len(offdiag_values)
        std = np.std(offdiag_values)
        zdf = (pmirm.df - mean) / std
        np.fill_diagonal(zdf.values, 0)
        return cls(zdf, num_datapoints=pmirm.num_datapoints,
                name=pmirm.name)
