import importlib
import os
import sys


def _bootstrap_dependencies():
    '''
    Ensure external bindings expose the legacy module names expected by the
    generated Cython code (e.g. ``TACS`` instead of ``tacs.TACS``).
    '''
    dependency_aliases = {
        'TACS': 'tacs.TACS',
        'constitutive': 'tacs.constitutive',
        'elements': 'tacs.elements',
        'functions': 'tacs.functions',
    }

    for alias, target in dependency_aliases.items():
        if alias in sys.modules:
            continue

        try:
            module = importlib.import_module(target)
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                f"stacs requires '{target}' so it can register the '{alias}' "
                "module name expected by the compiled extension."
            ) from exc

        sys.modules.setdefault(alias, module)


_bootstrap_dependencies()


def _load_extension():
    '''
    Import the compiled ``stacs.STACS`` module and register the legacy
    top-level name ``STACS`` for generated Cython code that expects it.
    '''
    module = importlib.import_module(f'{__name__}.STACS')
    sys.modules.setdefault('STACS', module)
    return module


STACS = _load_extension()


def get_cython_include():
    '''
    Get the include directory for the Cython .pxd files in PSPACE
    '''
    return [os.path.abspath(os.path.dirname(__file__))]

def get_include():
    '''
    Get the include directory for the Cython .pxd files in PSPACE
    '''
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_inc_dirs = ['src/include']

    inc_dirs = []
    for path in rel_inc_dirs:
    	inc_dirs.append(os.path.join(root_path, path))

    return inc_dirs

def get_libraries():
    '''
    Get the library directories
    '''
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_lib_dirs = ['lib']
    libs = ['stacs']
    lib_dirs = []
    for path in rel_lib_dirs:
    	lib_dirs.append(os.path.join(root_path, path))

    return lib_dirs, libs
