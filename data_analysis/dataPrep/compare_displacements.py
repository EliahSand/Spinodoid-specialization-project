from __future__ import annotations

from .common import DISPLACEMENTS
from .group_runner import run_component_group


if __name__ == "__main__":
    run_component_group("displacements", DISPLACEMENTS)

