r"""
Cache of computations of top xi evaluations and various evaluation of tautological classes.
"""

from ast import literal_eval
import os
import os.path

from sage.env import DOT_SAGE
import sage.misc.persist

from admcycles.diffstrata.adm_eval_cache import ADM_EVAL_CACHE
from admcycles.diffstrata.xi_cache import XI_TOP_CACHE


class Cache(dict):
    r"""
    Dictionary with synchronization in a ``.sobj`` file.

    TESTS:

    The synchronization files stores only diffs::

        sage: from admcycles.diffstrata.cache import Cache
        sage: filename = tmp_filename(ext='.sobj')
        sage: C = Cache({2: 'two'}, filename=filename)
        sage: C[3] = 'three'
        sage: D = Cache(filename=filename)
        sage: D
        {3: 'three'}

    An example without synchronization file::

        sage: C = Cache({2: 'two'})
        sage: C[3] = 'three'
        sage: C
        {2: 'two', 3: 'three'}
    """

    def __init__(self, values=None, filename=None):
        self._filename = filename
        if values is None:
            self._default = {}
            dict.__init__(self)
        else:
            self._default = values.copy()
            dict.__init__(self, values)
        self.update_from_sobj_file(self._filename)

    def reset_default(self):
        r"""
        Reset the content of the dictionary to the set of default values.

        TESTS::

            sage: from admcycles.diffstrata.cache import Cache
            sage: C = Cache({1: 'ein', 2: 'zwei'})
            sage: C[3] = 'drei'
            sage: sorted(C.keys())
            [1, 2, 3]
            sage: C.reset_default()
            sage: sorted(C.keys())
            [1, 2]
        """
        self.clear()
        self.update(self._default)

    def update_from_sobj_file(self, filename=None):
        r"""
        Update the values of this dictionary from a ``.sobj`` file.

        Missing or empty files are ignored. The cache file is synchronized.

        EXAMPLES::

            sage: from admcycles.diffstrata.cache import Cache
            sage: filename1 = tmp_filename(ext='.sobj')
            sage: filename2 = tmp_filename(ext='.sobj')
            sage: C1 = Cache({2: 'two'}, filename=filename1)
            sage: C1[3] = 'three'
            sage: C2 = Cache(filename=filename2)
            sage: C2
            {}
            sage: C2.update_from_sobj_file(filename1)
            sage: Cache(filename=filename2)
            {3: 'three'}
        """
        if filename is None:
            filename = self._filename
        if filename is not None and os.path.isfile(
                filename) and os.stat(filename).st_size:
            self.update(sage.misc.persist.load(filename))
            if filename != self._filename:
                self.save_to_sobj_file()

    def diff(self):
        cache = self.copy()
        for key in self._default:
            del cache[key]
        return cache

    def save_to_sobj_file(self, filename=None, only_diff=True):
        r"""
        Save the values of this dictionary to a ``.sobj`` file.
        """
        if filename is None:
            filename = self._filename
        if filename is not None:
            cache = self.diff() if only_diff else self.copy()
            if cache:
                sage.misc.persist.save(cache, filename)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self.save_to_sobj_file()


class FakeCache(Cache):
    r"""
    A cache where ``__setitem__`` does not do anything.

    This class is intended for benchmarking without caching.

    TESTS::

        sage: from admcycles.diffstrata.cache import FakeCache
        sage: C = FakeCache()
        sage: C[3] = 'trois'
        sage: C
        {}
    """

    def __setitem__(self, key, value):
        pass


ADM_EVALS_FILENAME = os.path.join(DOT_SAGE, 'adm_evals.sobj')
TOP_XIS_FILENAME = os.path.join(DOT_SAGE, 'top_xis.sobj')

ADM_EVALS = Cache(ADM_EVAL_CACHE, ADM_EVALS_FILENAME)
TOP_XIS = Cache(XI_TOP_CACHE, TOP_XIS_FILENAME)


def list_top_xis(xi_dict=TOP_XIS):
    """
    A generator for decyphering LevelStratum.dict_key

    By default, the xi cache is used, alternatively a different
    cache dictionary may be specified.

    Args:
        xi_dict (dict, optional): cache dictionary. Defaults to ``TOP_XIS``.

    Yields:
        tuple: signature list, residue conditions, top xi evaluation
    """
    for key in sorted(xi_dict):
        sig_list, rcs = key
        yield (sig_list, rcs, xi_dict[key])


def print_top_xis(xi_dict=TOP_XIS):
    """
    The top xi powers in the cache are printed in human-readable form.

    Alternatively, any valid xi cache dictionary may be specified.

    The dimensions of the table might have to be adapted to console size and
    complexity of the involved strata.

    Args:
        xi_dict (dict, optional): xi cache dictionary. Defaults to ``TOP_XIS``.

    EXAMPLES::

        sage: from admcycles.diffstrata import print_top_xis
        sage: from admcycles.diffstrata.cache import TOP_XIS
        sage: TOP_XIS.reset_default()
        sage: print_top_xis()
        Stratum                                 | Residue Conditions                          | xi^dim
        --------------------------------------------------------------------------------------------------
        (-8, -2, 8)                             | [(0, 0), (0, 1)]                            | 1
        (-8, 0, 6)                              | [(0, 0)]                                    | 1
        (-8, 1, 1, 2, 2)                        | [(0, 0)]                                    | 49
        ...
        (3, 3)                                  | ()                                          | 0
        (4,)                                    | ()                                          | 305/580608
        (6,)                                    | ()                                          | -87983/199065600
        (8,)                                    | ()                                          | 339849/504627200
    """
    LEFT = 40
    MIDDLE = 45
    RIGHT = 13
    print('{:<{width}}'.format('Stratum', width=LEFT) + '|' +
          '{:<{width}}'.format(' Residue Conditions', width=MIDDLE) + '|' +
          '{:<{width}}'.format(' xi^dim', width=RIGHT))
    print('-' * (LEFT + MIDDLE + RIGHT))
    for sig_list, res_conds, xi in list_top_xis(xi_dict):
        if len(sig_list) == 1:
            sig = sig_list[0]
        else:
            sig = list(sig_list)
        if not res_conds:
            rcs = tuple()
        else:
            if len(res_conds) == 1:
                rcs = list(res_conds[0])
            else:
                rcs = [list(c) for c in res_conds]
        print('{:<{width}}'.format(str(sig), width=LEFT) + '| ' +
              '{:<{width}}'.format(str(rcs), width=MIDDLE - 1) + '| ' +
              '{:<{width}}'.format(str(xi), width=RIGHT - 1))


def print_adm_evals(adm_dict=ADM_EVALS):
    """
    The cached evaluations are printed in human-readable form.

    Alternatively, any valid adm cache dictionary may be specified (see load_adm_evals
    for details on the format).

    The dimensions of the table might have to be adapted to console size and
    complexity of the involved strata.

    Args:
        adm_dict (dict, optional): evaluations cache dictionary. Defaults to None.

    EXAMPLES::

        sage: from admcycles.diffstrata import print_adm_evals
        sage: from admcycles.diffstrata.cache import ADM_EVALS
        sage: ADM_EVALS.reset_default()
        sage: print_adm_evals()
        Stratum                                 | Psis                                        | eval
        --------------------------------------------------------------------------------------------------
        ...
        (2,)                                    | {1: 3}                                      | 1/1920
        (2, 2)                                  | {1: 1, 2: 5}                                | 71/322560
        (2, 2)                                  | {1: 2, 2: 4}                                | 13/53760
        (2, 2)                                  | {1: 3, 2: 3}                                | 13/53760
        (2, 2)                                  | {1: 6}                                      | 43/322560
        (4,)                                    | {1: 5}                                      | 13/580608
        (6,)                                    | {1: 7}                                      | 281/199065600
    """
    LEFT = 40
    MIDDLE = 45
    RIGHT = 13
    print('{:<{width}}'.format('Stratum', width=LEFT) + '|' +
          '{:<{width}}'.format(' Psis', width=MIDDLE) + '|' +
          '{:<{width}}'.format(' eval', width=RIGHT))
    print('-' * (LEFT + MIDDLE + RIGHT))
    for key in sorted(adm_dict):  # sort by stratum
        sig, psis = key
        value = adm_dict[key]
        dpsis = dict(psis)
        print('{:<{width}}'.format(str(sig), width=LEFT) + '| ' +
              '{:<{width}}'.format(str(dpsis), width=MIDDLE - 1) + '| ' +
              '{:<{width}}'.format(str(value), width=RIGHT - 1))


def jsonify_dict(d):
    json_dict = {}
    for k, v in d.items():
        json_dict[str(k)] = v
    return json_dict


def unjsonify_dict(json_dict):
    restored_dict = {}
    for k, v in json_dict.items():
        restored_dict[literal_eval(k)] = v
    return restored_dict


def import_adm_evals(filename):
    """
    Import a dictionary of the form adm_key : value into the cache.

    Args:
        filename (String): sobj filename
    """
    ADM_EVALS.update_from_sobj_file(filename)


def import_top_xis(filename):
    """
    Import a dictionary of the form LevelStratum.dict_key : value into the cache.

    Args:
        filename (String): sobj filename
    """
    TOP_XIS.update_from_sobj_file(filename)


def load_adm_evals():
    """
    The cache dictionary for the admcycles evaluations (via GeneralisedStratum.adm_evaluate)

    The cache dictionary is of the form:

    adm_key -> GeneralisedStratum.adm_evaluate(adm_key)  (a rational number)

    Returns:
        dict: cache dictionary
    """
    return ADM_EVALS


def load_xis():
    """
    The cache dictionary for the top xi powers (evaluated)

    The cache dictionary is of the form:

    L.dict_key -> (L.xi)^L.dim().evaluate()  (a rational number)

    for a LevelStratum L.

    Returns:
        dict: cache dictionary
    """
    return TOP_XIS
