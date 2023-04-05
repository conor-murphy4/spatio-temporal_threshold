# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import os
import abc
import copy
import json


class _ParamUpdate(abc.ABCMeta):
    """Metaclass for automating the updating of default parameters when inheriting parameter reader classes."""
    def __new__(mcs, name, bases, attrs):
        cls = super().__new__(mcs, name, bases, attrs)
        def_params = [base._def_params for base in bases if hasattr(base, '_def_params')] + [cls._def_params]
        cls._def_params = dict()
        for def_param in def_params:
            cls._def_params.update(def_param)
        return cls


class Parameters(metaclass=_ParamUpdate):
    _def_params = dict()

    def __init__(self, path=None, **kwargs):
        """Read and store parameters from JSON file.

        Args:
            path (str, None): Path to JSON with input parameters
            **kwargs: Keyword arguments can be used to override parameters from input file
        """
        self.__params = copy.deepcopy(self._def_params)
        if path:
            with open(path, 'r') as fh:
                params = json.load(fh)
                self._update(params, kwargs)
        else:
            params = kwargs
        self._update(self.__params, params)
        if params:
            raise ValueError('Unrecognized input fields: {}'.format(', '.join(params)))
        self._dirname = os.path.dirname(os.path.realpath(path)) if path else '.'

    def _access(self, key):
        """Reads the parameter value for the given key.

        For nested dictionaries, the key should be given with each level separated by `/`.

        Args:
            key (str): Name/location of parameter to read

        Returns:
            any: Value corresponding to the key
        """
        keys = key.split('/')
        out = self.__params
        for k in keys:
            out = out[k]
        return out

    @classmethod
    def _update(cls, original, update):
        """Updates a given nested dictionary based on values from a second dictionary.

        Only entries that already exist in the original dictionary will be updated. The only exception is if the key
        corresponds to ``dof``, in which case it will be fully overwritten.

        Updated key-value pairs are removed from the second dictionary, the original dictionary is modified in-place.

        Args:
            original (dict): Dictionary to be updated
            update (dict): Dictionary with key-value pairs to update
        """
        for key in original:
            if key not in update:
                continue
            if isinstance(original[key], dict) and key != 'dof':
                cls._update(original[key], update.pop(key))
                continue
            original[key] = update.pop(key)
