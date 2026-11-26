"""pKa request/response provenance and representation sensitivity; no predictor is bundled."""

import hashlib
import json
import math
from typing import Literal, Protocol

from pydantic import Field

from .ast import AnalysisRecord, AtomID, Identifier, Node


class PKAQuery(Node):
    query_id: Identifier
    acid_state_id: Identifier
    acid_view_id: Identifier
    conjugate_base_state_id: Identifier
    conjugate_base_view_id: Identifier
    site_atom_id: AtomID
    equilibrium_type: Literal["microscopic", "macroscopic"]
    solvent: str = Field(min_length=1)
    temperature_kelvin: float = Field(gt=0, allow_inf_nan=False)
    model_id: str = Field(min_length=1)
    model_version: str = Field(min_length=1)
    raw_input: dict

    def equilibrium_key(self) -> str:
        # Only the resonance view IDs are excluded. Sites, tautomers (state IDs),
        # protonation equilibria, solvent, temperature and models remain distinct.
        data = self.model_dump(exclude={"query_id", "acid_view_id", "conjugate_base_view_id", "raw_input"})
        return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


class PKAModel(Protocol):
    model_id: str
    model_version: str

    def predict(self, raw_input: dict) -> dict: ...


class UnknownPKAModel:
    model_id, model_version = "unconfigured", "none"

    def predict(self, raw_input: dict) -> dict:
        return {"status": "unknown", "value": None, "reason": "No pKa predictor or source is configured"}


def query_pka(query: PKAQuery, model: PKAModel) -> AnalysisRecord:
    if (model.model_id, model.model_version) != (query.model_id, query.model_version):
        raise ValueError("pKa query model identity/version mismatch")
    # Isolate the submitted input from mutating adapters.
    raw_output = model.predict(json.loads(json.dumps(query.raw_input, allow_nan=False)))
    value = raw_output.get("value")
    if value is not None and (type(value) not in {int, float} or not math.isfinite(value)):
        raise ValueError("pKa response must be a finite numeric value or null")
    return AnalysisRecord(analysis_id=query.query_id, operation="PKA_QUERY", state_id=query.acid_state_id,
                          view_id=query.acid_view_id, status="unknown" if value is None else "complete",
                          input=query.model_dump(mode="json"),
                          output={"equilibrium_key": query.equilibrium_key(), "raw_output": raw_output,
                                  "pka": value, "is_experimental_measurement": False},
                          uncertainties=["Prediction applicability and source reliability require separate review."])


def representation_sensitivity(records: list[AnalysisRecord]) -> list[AnalysisRecord]:
    groups = {}
    for record in records:
        if record.operation != "PKA_QUERY":
            continue
        query = PKAQuery.model_validate(record.input)
        key = query.equilibrium_key()
        groups.setdefault(key, []).append(record)
    result = []
    for index, (key, values) in enumerate(sorted(groups.items())):
        predictions = [v.output.get("pka") for v in values]
        numbers = [x for x in predictions if x is not None]
        representations = {(v.input["acid_view_id"], v.input["conjugate_base_view_id"]) for v in values}
        differing = len(set(numbers)) > 1 and len(representations) > 1
        result.append(AnalysisRecord(
            analysis_id=f"sensitivity.{index}", operation="REPRESENTATION_SENSITIVITY",
            status="complete" if len(numbers) >= 2 and len(representations) >= 2 else "unknown",
            input={"equilibrium_key": key, "query_ids": [v.analysis_id for v in values]},
            output={"representation_sensitive": differing, "all_predictions": predictions,
                    "range": max(numbers) - min(numbers) if numbers else None,
                    "selected_prediction": None, "representations": len(representations)},
            uncertainties=["Differences are retained; no preferred representation or prediction is selected."],
        ))
    return result
