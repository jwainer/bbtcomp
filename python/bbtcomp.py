""" A Python implementation of the BBT model com compare multiple
    algorithms on multiple data sets.

    The program requires cmdstanpy as an interface to Stan.
    Find it how to install it at  https://mc-stan.org/cmdstanpy/

    VERSION = 0.3.0

"""

import os
import numpy as np
import pandas as pd
import cmdstanpy


def allpairs(n):
    k = 0
    for i in range(n-1):
        for j in range(i + 1, n):
            yield (i, j, k)
            k += 1


class Wintable():
    def print(self):
        """
        Generate a DataFrame with the names of two algorithms, the number of wins for each algorithm, and return it.
        """
        name1 = np.array(self.alg_names)[self.table[:, 0]]
        name2 = np.array(self.alg_names)[self.table[:, 1]]
        return pd.DataFrame({"alg1": name1, "alg2": name2,
                             "win1": self.table[:, 2],
                             "win2": self.table[:, 3]})


class BBTmodel():
    def diagnose(self):
        print(self.model.diagnose())


#
# bbtcomp
#

def bbtcomp(dat,
            datasd=None,
            dbcol=0,
            lrope=True,
            paired=True,
            lrope_value=0.4,
            deal_with_ties="s",
            **kwargs):
    """bbtcomp: generate a bbt model

    This is the main function in the package. It generates a BBTmodel
    which can be further queried to obtain the probabilities distributions and
    other queries.

    Parameters
    ----------

    dat: numpy array or pandas DataFrame
        The input data. Each column is (in general) an algorithm
        Each line is a data set.
        One particular optional column (dbcol) is an indicator of the db name or id.
        The point of the dbcol is that multiple lines with the same dbcol name
        are folds of the same data set (for the computations of the local rope).
        If dbcol is -1 or None, it is assumed that there is no dbcol,
        and thus each entry is the average of the fold values for the
        algorithm and for the data set.

    datasd: numpy array or pandas DataFrame
        If dat has no dbcol, then to compute the local ROPE, one needs
        the standard deviation for each
        fold and algorithm. This is the purpose of the datasd

    dbcol: int or None
        The index of the column that contains the db name or id. If -1 or None,
        then there is no such column.

    lrope: bol
        Whether to use the local Rope (if possible) to compute wins and losses.
        Local ROPE can only be calculated
        if either there is fold information (and thus a dbcol column) or the
        datasd information. If local ROPE cannot
        be conputed, this parameter will be silently ignored.

    paired: bol
        Whether to use the paired version of the local ROPE. Only possible is
        dbcol is present. Silently ingored otherwise.

    deal_with_ties: str in ["s", "a", "f", "d"]
        How to deal with ties.
        -  a: add
        -  s: spread (default)
        -  f: forget
        -  d: use Davidson model


    Returns
    -------

    a BBTmodel

    """
   
    if deal_with_ties == "d":
        use_davidson = True
    else:
        use_davidson = False
    
    w = makewintable(dat, datasd, dbcol, lrope, paired, lrope_value,
                     deal_with_ties)
    mod = mcmcbbt(w, use_davidson, **kwargs)                 )
    return mod


#
# makewintable
#


def makewintable(dat, datsd=None, dbcol=0, lrope=True,
                 paired=True, lrope_value=0.4,
                 deal_with_ties="s"):

    names = list(dat.columns)
    if dbcol is None or dbcol == -1:
        dbcol = -1
    elif type(dbcol) is str:
        dbname = dbcol
        names.remove(dbname)
    else:
        dbname = names[dbcol]
        names.remove(dbname)

    nalgs = len(names)

    if dbcol == -1 or len(np.unique(dat[dbname])) == len(dat[dbname]):
        out = compute_dif_no_paired(dat, datsd, lrope_value)
    else:
        if not lrope:
            datmean = dat.groupby(dbname).mean()
            out = compute_dif_no_paired(datmean, None, lrope_value)

        else:
            out = np.zeros((nalgs * (nalgs - 1) // 2, 5), dtype=int)
            for db in np.unique(dat[dbname]):
                inx = dat[dbname] == db
                x = dat[inx].drop(dbname, axis=1)
                for i, j, k in allpairs(nalgs):
                    dif = x.iloc[:, i] - x.iloc[:, j]
                    mean = dif.mean()
                    sd = dif.std()
                    win1 = int(mean > lrope_value * sd)
                    win2 = int(mean < -lrope_value * sd)
                    ties = int(mean >= -lrope_value * sd and
                               mean <= lrope_value * sd)
                    out[k, 0] = i
                    out[k, 1] = j
                    out[k, 2] += win1
                    out[k, 3] += win2
                    out[k, 4] += ties

    obj = Wintable()
    obj.table = proc_ties(out, deal_with_ties)
    obj.table_pre = out
    obj.alg_names = names
    obj.lrope = lrope
    obj.paired = paired
    obj.lrope_value = lrope_value
    obj.ties_proc = deal_with_ties
    return obj


def compute_dif_no_paired(dmean, dstd, lrope_value):
    n = dmean.shape[1]
    out = np.empty((n*(n-1)//2, 5), dtype=int)
    for i, j, k in allpairs(n):
        delta = dmean.iloc[:, i] - dmean.iloc[:, j]
        if dstd is None:
            sd = 0.0
        else:
            sd = (dstd.iloc[:, i]**2 + dstd.iloc[:, j]**2)/2
        win1 = np.sum(delta > lrope_value * sd)
        win2 = np.sum(delta < -lrope_value * sd)
        ties = np.sum(np.logical_and(delta >= -lrope_value * sd,
                                     delta <= lrope_value * sd))
        out[k, :] = (i, j, win1, win2, ties)

    return out


def proc_ties(out, deal_with_ties):
    if deal_with_ties == "d":
        return out
    if deal_with_ties == "a":
        newout = out[:, 0:4]
        newout[:, 2] += out[:, 4]
        newout[:, 3] += out[:, 4]
        return newout
    if deal_with_ties == "s":
        newout = out[:, 0:4]
        newout[:, 2] += np.ceil(out[:, 4]/2).astype(int)
        newout[:, 3] += np.ceil(out[:, 4]/2).astype(int)
        return newout

#
#   mcmcbbt
#


def mcmcbbt(win, use_davidson=False, **kwargs):
      
    tab = win.table
    N = tab.shape[0]
    K = len(win.alg_names)
    data = dict(K=K, N=N,
                player1=tab[:, 0] + 1,
                player2=tab[:, 1] + 1,
                win1=tab[:, 2],
                win2=tab[:, 3])

    mod_file = "./bbt-full.stan"
    data.update(dict(use_davidson=int(use_davidson),
                         ties=tab[:, 4]))

    mod = cmdstanpy.CmdStanModel(stan_file=mod_file)
    fit = mod.sample(data=data, **kwargs)

    obj = BBTmodel()
    obj.model = fit
    obj.davidson = use_davidson
    obj.wintable = win
    return obj

#
# table_pwin
#


def table_pwin(modout,
                 control=None,
                 selected=None,
                 short=True,
                 rope=(0.45, 0.55),
                 columns=("median", "mean", "low", "high", "delta",
                          "above.50", "in.rope"),
                 hdi=0.89,
                 ndigits=2):
    if short and len(columns) == 7:
        columns = ("mean", "delta", "above.50", "in.rope")

    outcolumns = ["pairs"] + [c for c in columns if c in
                              ("median", "mean", "low", "high",
                               "delta", "above.50", "in.rope")]
    zz, outnames = get_pwin(modout, control, selected)
    out = dict()
    out["pairs"] = outnames
    if "median" in outcolumns:
        out["median"] = np.around(np.median(zz, axis=0), decimals=ndigits)
    if "mean" in outcolumns:
        out["mean"] = np.around(np.mean(zz, axis=0), decimals=ndigits)
    if "above.50" in outcolumns:
        out["above.50"] = np.around(np.sum(zz > 0.50, axis=0)/zz.shape[0],
                                    decimals=ndigits)
    if "in.rope" in outcolumns:
        out["in.rope"] = np.around(np.sum(np.logical_and(zz >= rope[0],
                                                         zz <= rope[1]),
                                          axis=0)/zz.shape[0],
                                   decimals=ndigits)
    if "delta" in outcolumns or "low" in outcolumns or "high" in outcolumns:
        hdix = np.apply_along_axis(lambda x: newhdi(x, hdi), 0, zz)
        if "low" in outcolumns:
            out["low"] = round(hdix[0, :], ndigits)
        if "high" in outcolumns:
            out["high"] = round(hdix[1, :], ndigits)
        if "delta" in outcolumns:
            out["delta"] = np.around(np.apply_along_axis(lambda x: x[1] - x[0],
                                                         0, hdix),
                                     decimals=ndigits)

    return pd.DataFrame(out)[outcolumns]


def get_pwin(mod, control=None, selected=None):

    def auxprob(a, b):
        return a/(a+b)

    zz = mod.model.stan_variable("beta")
    nalg = zz.shape[1]
    names = mod.wintable.alg_names
    ord = np.argsort(-np.mean(zz, axis=0))
    outnames = []
    ordnames = np.array(names)[ord]
    zz = np.exp(zz[:, ord])
    if selected is not None:
        oo = list(ordnames)
        ordnames = np.array([x for x in ordnames if x in selected])
        pp = [oo.index(c) for c in selected]
        zz = zz[:, pp]
        nalg = len(selected)
    if control is None or control not in ordnames:
        out = np.empty((zz.shape[0], int(nalg * (nalg - 1) // 2)))
        for i, j, k in allpairs(nalg):
            out[:, k] = auxprob(zz[:, i], zz[:, j])
            outnames.append(ordnames[i] + " > " + ordnames[j])
        return (out, outnames)
    else:
        j = np.where(ordnames == control)[0][0]
        out = np.empty((zz.shape[0], nalg - 1))
        k = 0
        for i in range(j):
            out[:, k] = auxprob(zz[:, i], zz[:, j])
            k += 1
            outnames.append(ordnames[i] + " > " + ordnames[j])
        for i in range(j + 1, nalg):
            out[:, k] = auxprob(zz[:, j], zz[:, i])
            k += 1
            outnames.append(ordnames[j] + " > " + ordnames[i])
        return (out, outnames)
#
#  ppc_summary
#


def table_ppc(mod):

    tab = mod.wintable.table
    y = tab[:, 2]
    yrep = mod.model.stan_variable("win1_rep")
    h50 = np.apply_along_axis(lambda x: newhdi(x, 0.5), 0, yrep)
    h90 = np.apply_along_axis(lambda x: newhdi(x, 0.9), 0, yrep)
    h95 = np.apply_along_axis(lambda x: newhdi(x, 0.95), 0, yrep)
    h100 = np.apply_along_axis(lambda x: (np.min(x), np.max(x)), 0, yrep)
    x50 = round(np.mean(np.logical_and(y >= h50[0, :], y <= h50[1, :])), 2)
    x90 = round(np.mean(np.logical_and(y >= h90[0, :], y <= h90[1, :])), 2)
    x95 = round(np.mean(np.logical_and(y >= h95[0, :], y <= h95[1, :])), 2)
    x100 = round(np.mean(np.logical_and(y >= h100[0, :], y <= h100[1, :])), 2)
    return pd.DataFrame({"hdi": [0.5, 0.9, 0.95, 1.0],
                         "proportion": [x50, x90, x95, x100]})  


def newhdi(arr, credmass):
    # based on R's HDInterval hdivector
    import math
    x = np.sort(arr)
    n = len(x)
    exclude = n - math.floor(n * credmass) - 1
    low_poss = x[0:exclude]
    upp_poss = x[(n - exclude):n]
    best = np.argmin(upp_poss - low_poss)
    result = (low_poss[best], upp_poss[best])
    return result
