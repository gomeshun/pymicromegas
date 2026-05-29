from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymicromegas.pymicromegas import MicrOmegas


class MicrOmegasProcess:
    def __init__(self, model: MicrOmegas):
        self.model = model

    @property
    def lib(self):
        return self.model.lib