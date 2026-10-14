"""Reproducible random streams.

Each grid point / task gets its own stream derived from (base seed, integer key path) with
numpy SeedSequence, so results do not depend on execution order or the number of workers.
"""
from __future__ import annotations

import hashlib

import numpy as np


def key_to_ints(key) -> list[int]:
    """Map a key (int, str, or tuple of those) to a list of uint32 words, stably across runs."""
    if isinstance(key, (tuple, list)):
        out = []
        for k in key:
            out.extend(key_to_ints(k))
        return out
    if isinstance(key, (int, np.integer)):
        return [int(key) & 0xFFFFFFFF, (int(key) >> 32) & 0xFFFFFFFF]
    h = hashlib.sha256(str(key).encode()).digest()
    return [int.from_bytes(h[i:i + 4], "little") for i in range(0, 16, 4)]


def stream(seed: int, *key) -> np.random.Generator:
    """Independent Generator for (seed, key...)."""
    ss = np.random.SeedSequence(entropy=int(seed), spawn_key=tuple(key_to_ints(key)))
    return np.random.Generator(np.random.PCG64(ss))


def numba_seed(seed: int, *key) -> int:
    """A 63-bit integer seed for numba kernels (np.random.seed inside njit), derived the same way."""
    ss = np.random.SeedSequence(entropy=int(seed), spawn_key=tuple(key_to_ints(key)))
    return int(ss.generate_state(2, dtype=np.uint32).astype(np.uint64) @ np.array([1, 2**31], dtype=np.uint64)) & (2**62 - 1)
