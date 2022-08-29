#!/usr/bin/env python3


import hashlib
import uuid

import sys
import os

# Generate reproducible UUIDs using hashed file name(s) as seed(s)
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print(f"Usage: {sys.argv[0]} file1.txt ...")
        sys.exit(1)

    for f in sys.argv[1:]:
        hasher = hashlib.md5()
        hasher.update(os.path.basename(f).encode("utf-8"))

        print(uuid.UUID(hasher.hexdigest()))
