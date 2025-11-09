PYTHON = python3

# Optional: support building with a pinned legacy Cython version
# Override via environment if needed
LEGACY_CYTHON_VERSION ?= 0.29.36
LEGACY_PYTHON_SITE ?= $(CURDIR)/python_legacy

default:
	./build.sh

interface:
	${PYTHON} setup.py build_ext --inplace

complex_interface:
	${PYTHON} setup.py build_ext --inplace --define PSPACE_USE_COMPLEX

.PHONY: legacy_cython
legacy_cython:
	@${PYTHON} tools/ensure_legacy_cython.py "$(LEGACY_PYTHON_SITE)" "$(LEGACY_CYTHON_VERSION)"

.PHONY: interface_legacy
interface_legacy: legacy_cython
	PYTHONPATH=$(LEGACY_PYTHON_SITE)$${PYTHONPATH:+:$$PYTHONPATH} ${PYTHON} setup.py build_ext --inplace

.PHONY: complex_interface_legacy
complex_interface_legacy: legacy_cython
	PYTHONPATH=$(LEGACY_PYTHON_SITE)$${PYTHONPATH:+:$$PYTHONPATH} ${PYTHON} setup.py build_ext --inplace --define PSPACE_USE_COMPLEX
