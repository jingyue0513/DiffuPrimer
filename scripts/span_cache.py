"""Byte-preserving cache for deterministic sequence-only span calculations."""
import dataclasses
import fcntl
import hashlib
import inspect
import json
import os
from pathlib import Path

import numpy as np
import primer3

from diffuprimer import constraints, guidance, sequence


ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "results/span_cache"


def install():
    CACHE.mkdir(parents=True, exist_ok=True)
    dependencies = [inspect.getsource(constraints.Conditions), inspect.getsource(constraints.oligo_properties),
                    inspect.getsource(guidance.interval_table), inspect.getsource(sequence.reverse_complement), primer3.__version__]

    def wrap(function, temperature=False):
        identity = hashlib.sha256("\n".join(dependencies + [inspect.getsource(function)]).encode()).hexdigest()

        def cached(template, conditions, *extra):
            key = hashlib.sha256(json.dumps({"function": identity, "template": hashlib.sha256(template.encode()).hexdigest(),
                "conditions": dataclasses.asdict(conditions)}, sort_keys=True).encode())
            if extra:
                intervals, allowed = extra[0]
                for array in [intervals, *allowed]:
                    key.update(array.dtype.str.encode())
                    key.update(str(array.shape).encode())
                    key.update(array.tobytes())
            path = CACHE / (key.hexdigest() + ".npz")
            with path.with_suffix(".lock").open("a") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                if path.exists():
                    with np.load(path, allow_pickle=False) as data:
                        if str(data["identity"]) != identity:
                            raise ValueError("Span cache identity changed")
                        return data["temperatures"].copy() if temperature else (data["intervals"].copy(), [data["F"].copy(), data["R"].copy()])
                result = function(template, conditions, *extra)
                temp = path.with_suffix(".tmp." + str(os.getpid()))
                with temp.open("wb") as handle:
                    if temperature:
                        np.savez_compressed(handle, identity=identity, temperatures=result)
                    else:
                        np.savez_compressed(handle, identity=identity, intervals=result[0], F=result[1][0], R=result[1][1])
                temp.replace(path)
                return result
        return cached

    guidance.feasible_oligo_spans = wrap(guidance.feasible_oligo_spans)
    if hasattr(guidance, "compatible_span_temperatures"):
        guidance.compatible_span_temperatures = wrap(guidance.compatible_span_temperatures, temperature=True)
