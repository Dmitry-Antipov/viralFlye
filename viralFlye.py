#!/usr/bin/env python3


"""
This script sets up environment paths
and invokes viralFlye without installation.
"""

import os
import sys

if sys.version_info < (3,6):
    raise SystemExit("Requires Python 3.6 or higher")

#Setting environment for local run without installation
viralflye_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, viralflye_root)

VIRALCOMPLETE = "viralComplete"
vc_absolute = os.path.join(viralflye_root, VIRALCOMPLETE)
os.environ["PATH"] = vc_absolute + os.pathsep + os.environ["PATH"]

#entry point
from viralflye.main import main
sys.exit(main())
