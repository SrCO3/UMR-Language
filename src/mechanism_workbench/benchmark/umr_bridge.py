"""Version dispatch; all chemistry stays in the existing UMR implementations."""
import hashlib
import json

from ..umr.ast import MechanismDocument as Document01
from ..umr.v02.ast import MechanismDocument as Document02
from ..umr.v03.ast import MechanismDocument as Document03
from ..umr.v04.ast import MechanismDocument as Document04

from ..umr.v05.ast import MechanismDocument as Document05
from ..umr.v06.ast import MechanismDocument as Document06
from ..umr.v07.ast import MechanismDocument as Document07
DOCUMENTS={"umr/0.1":Document01,"umr/0.2":Document02,"umr/0.3":Document03,"umr/0.4":Document04,"umr/0.5":Document05,"umr/0.6":Document06}
DOCUMENTS['umr/0.7']=Document07
from ..umr.v08.ast import MechanismDocument as Document08
DOCUMENTS['umr/0.8']=Document08
from ..umr.v09.ast import MechanismDocument as Document09
DOCUMENTS['umr/0.9']=Document09
from ..umr.v10.ast import MechanismDocument as Document10
DOCUMENTS['umr/1.0']=Document10


def canonical_json(value):return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(",",":"),allow_nan=False)
def digest(value):return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def parse_document(data):
    version=data.get("schema_version")
    if version not in DOCUMENTS:raise ValueError(f"Unsupported UMR schema version: {version}")
    return DOCUMENTS[version].model_validate(data)
